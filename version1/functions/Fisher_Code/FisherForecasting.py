import os
import numpy as np
from scipy.integrate import quad
from classy import Class
from getdist import MCSamples, plots
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import linregress

# Define the cosmological parameters
params = {
    'output': 'mPk',           
    'omega_b': 0.0223828,      
    'omega_cdm': 0.1201075,    
    'h': 0.6736,               
    'A_s': 2.1e-9,             
    'n_s': 0.9665,             
    'tau_reio': 0.0543,        
    'Omega_Lambda': 0.6825,    
    'Omega_k': 0.0,            
    'z_max_pk': 0.75,          
}

# Define the survey parameters
survey_params = [
    {
        "name": "MeerKLASS L-band",
        "V_survey": 1.3e9,
        "f_sky": 0.10,
        "nu_min": 900.0,
        "nu_max": 1185.0,
        "t_obs": 4000.0,
        "N_dish": 64,
        "k_min": 3.3
    },
    {
        "name": "MeerKLASS UHF-band",
        "V_survey": 10.8e9,
        "f_sky": 0.10,
        "nu_min": 580.0,
        "nu_max": 1000.0,
        "t_obs": 4000.0,
        "N_dish": 64,
        "k_min": 1.6
    },
    {
        "name": "SKA-MID Band 1",
        "V_survey": 221.6e9,
        "f_sky": 0.48,
        "nu_min": 350.0,
        "nu_max": 1050.0,
        "t_obs": 20000.0,
        "N_dish": 197,
        "k_min": 0.53
    },
]

# Define constants
eta = 1.0
nu21 = 1420
N_pol = 2.0
T_spl = 3.0
T_CMB = 2.73
T_rx_coeff = 7.5
T_rx_freq_coeff = 10.0


def hubble_parameter(z):
    H0 = params['h'] /2997.9 
    matter = params['omega_b'] * (1 + z)**3 + params['omega_cdm'] * (1 + z)**3
    dark_energy = params['Omega_Lambda']
    curvature = params['Omega_k'] * (1 + z)**2
    return H0 * np.sqrt(matter + dark_energy + curvature)

def calculate_noise_power(survey, k_values, z):
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()

    T_gal = 25.0 * (408.0 / nu21) ** 2.75


    T_rx=T_rx_coeff + T_rx_freq_coeff * ((nu21 * 1e-3) - 0.75) ** 2

    T_sys = (T_rx + T_gal + T_spl + T_CMB)

    def chi(z):
        chi, _ = quad(lambda x: 1 / hubble_parameter(x), 0, z)
        return chi

    lambda_z = (0.21 / 3.083e24) * (1 + z)
    noise_power = (T_sys ** 2 * chi(z)**2 * ((z+1)/hubble_parameter(z)) * (lambda_z * 4 * np.pi * survey["f_sky"]) /
                   (eta * N_pol * survey["N_dish"] * survey["t_obs"] * 3600))
    return noise_power

def calculate_fisher_matrix(k, z, survey):
    P0 = 37300
    kT0 = 0.0165
    m = 0.40
    n = 0.30
    x = (np.log(k) - np.log(kT0))/np.log(kT0)

    def P_fit(k):
        return np.where(k < kT0, P0 ** (1-m*x), P0 ** (1-n*x))

    def dP_fit_dm(k):
        return np.where(k < kT0, -P_fit(k) * x, 0)

    def dP_fit_dn(k):
        return np.where(k >= kT0, -P_fit(k), 0)

    def dP_fit_dP0(k):
        return P_fit(k) / P0

    def dP_fit_dkT0(k):
        return np.where(k < kT0, -m * P_fit(k) / kT0, -n * P_fit(k) / kT0)

    errors = (P_fit(k_values) + calculate_noise_power(survey, k_values, z)) / np.sqrt(
        survey["V_survey"] * (4 * np.pi * k_values**2 * survey["k_min"] * 1e-3) / (2 * np.pi)**3)

    Fisher_matrix = np.zeros((4, 4))
    dP_dm = dP_fit_dm(k_values)
    dP_dn = dP_fit_dn(k_values)
    dP_dP0 = dP_fit_dP0(k_values)
    dP_dkT0 = dP_fit_dkT0(k_values)

    for i in range(len(k_values)):
        Fisher_matrix[0, 0] += dP_dm[i]**2 / errors[i]**2
        Fisher_matrix[0, 1] += dP_dm[i] * dP_dn[i] / errors[i]**2
        Fisher_matrix[0, 2] += dP_dm[i] * dP_dP0[i] / errors[i]**2
        Fisher_matrix[0, 3] += dP_dm[i] * dP_dkT0[i] / errors[i]**2
        Fisher_matrix[1, 1] += dP_dn[i]**2 / errors[i]**2
        Fisher_matrix[1, 2] += dP_dn[i] * dP_dP0[i] / errors[i]**2
        Fisher_matrix[1, 3] += dP_dn[i] * dP_dkT0[i] / errors[i]**2
        Fisher_matrix[2, 2] += dP_dP0[i]**2 / errors[i]**2
        Fisher_matrix[2, 3] += dP_dP0[i] * dP_dkT0[i] / errors[i]**2
        Fisher_matrix[3, 3] += dP_dkT0[i]**2 / errors[i]**2

    Fisher_matrix = np.triu(Fisher_matrix) + np.triu(Fisher_matrix, 1).T
    return Fisher_matrix

# Define k_values and z_value
k_values = np.linspace(0.005, 0.025, 1000)
z_value = 0.5

# Calculate Fisher matrices and covariance matrices
Fisher_matrices = [calculate_fisher_matrix(k_values, z_value, survey) for survey in survey_params]
covariance_matrices = [np.linalg.inv(F) for F in Fisher_matrices]

# Define the fiducial values for the parameters
fiducial_values = [0.40, 0.3, 37300, 0.0165]

# Convert covariance matrices to MCSamples objects
samples_list = [MCSamples(samples=np.random.multivariate_normal(mean=fiducial_values, cov=cov_mat, size=100000),
                          names=[r'$\alpha$', r'$\beta$', r'$10^4P_0$', r'$k_{0}$'])
                for cov_mat in covariance_matrices]

# Compute the relative uncertainties for each parameter
relative_errors_list = [np.sqrt(np.diag(cov_mat)) / fiducial_values for cov_mat in covariance_matrices]

# Create LaTeX output table for parameter constraints with relative errors
latex_table = r"\begin{table}[h!]\n"
latex_table += r"\begin{center}\n"
latex_table += r"\begin{tabular}{|c|c|c|c|c|}\n"
latex_table += r"\hline\n"
latex_table += r"Survey & $\sigma(\alpha)/\alpha$ & $\sigma(\beta)/\beta$ & $\sigma(10^5 P_0)/(10^5 P_0)$ & $\sigma(K_{0})/K_{0}$ \\\n"
latex_table += r"\hline\n"

for i, (samples, rel_errors) in enumerate(zip(samples_list, relative_errors_list)):
    survey_name = survey_params[i]["name"]
    latex_table += f"{survey_name} & {rel_errors[0]:.4f} & {rel_errors[1]:.4f} & {rel_errors[2]:.6f} & {rel_errors[3]:.4f} \\\n"

latex_table += r"\hline\n"
latex_table += r"\end{tabular}\n"
latex_table += r"\end{center}\n"
latex_table += r"\caption{Relative parameter uncertainties for each survey.}\n"
latex_table += r"\label{tab:relative_parameter_uncertainties}\n"
latex_table += r"\end{table}\n"

# Print LaTeX table
print(latex_table)

# Set larger plot size
plt.figure(figsize=(15, 12))

# Plot parameter constraints using plot_triangle
g = plots.get_subplot_plotter(subplot_size=4, width_inch=10)

# Adjust font size and quality settings
g.settings.axes_fontsize = 16
g.settings.legend_fontsize = 16
g.settings.lab_fontsize = 18
g.settings.alpha_filled_add = 0.80
g.settings.figure_legend_loc = 'upper right'
g.settings.linewidth = 2.5
g.settings.tight_layout = True


# Plot triangle plot with specified colors
colors = ['blue', 'green', 'red']
g.triangle_plot(
    samples_list, 
    filled=True, 
    legend_labels=[survey["name"] for survey in survey_params], 
    use_math_text=True, 
    contour_colors=colors
)

# Define a formatter to scale down the tick labels of P_0 by 1e4
def scale_down_by_1e4(x, pos):
    return f"{x / 1e4:.1f}"

# Apply the formatter to the x-axis of the specific subplot for P_0
ax_x = g.subplots[-1, 2]
ax_x.xaxis.set_major_formatter(ticker.FuncFormatter(scale_down_by_1e4))

# Apply the formatter to the y-axis of the specific subplot for P_0
ax_y = g.subplots[2, 1]
if ax_y is not None:
    ax_y.yaxis.set_major_formatter(ticker.FuncFormatter(scale_down_by_1e4))

# Ensure ticks are bold and thick
for ax_row in g.subplots:
    for ax in ax_row:
        if ax is not None:
            ax.tick_params(axis='both', which='major', labelsize=14, width=2)
            for spine in ax.spines.values():
                spine.set_linewidth(2.5)

# Save the figure to the home directory with high resolution
output_path = '/home/yolanda/Downloads/Parameter_constraints_high_quality.pdf'
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

# Get the subplot for m and n parameters
ax = g.subplots[1, 0]

# Get the data points
x_data = samples_list[0].samples[:, 0]  # m values
y_data = samples_list[0].samples[:, 1]  # n values

# Fit a linear regression line
slope, intercept, r_value, p_value, std_err = linregress(x_data, y_data)

# Print the equation for the line
print(f"n = {slope:.4f} * m + {intercept:.4f}")
