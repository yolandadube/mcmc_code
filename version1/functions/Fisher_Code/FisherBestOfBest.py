# Relevant python modules
import os
import sys
import emcee
import numpy as np
from scipy.integrate import quad
from classy import Class
import getdist.plots as gdplots
import matplotlib.pyplot as plt
import matplotlib.font_manager
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
from getdist import plots, MCSamples

# Define the cosmological parameters
cosmo = Class()
params = {
    'output': 'mPk , mTk',  # Request matter power spectrum calculation
    'omega_b': 0.0223828,  # Baryon density
    'omega_cdm': 0.1201075,  # Cold Dark Matter density
    'h': 0.6736,  # Hubble parameter
    'A_s': 2.1e-9,  # Amplitude of primordial power spectrum
    'n_s': 0.9665,  # Spectral index of primordial power spectrum
    'tau_reio': 0.0543,  # Optical depth to reionization
    'Omega_Lambda': 0.6825,  # Present-day density parameter for dark energy
    'Omega_k': 0.0,  # Present-day density parameter for curvature
    'z_max_pk': 0.75,  # Maximum redshift for storing power spectrum
}

cosmo.set(params)

# Survey parameters

# Add paths to your custom modules
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/surveySpecifications/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/surveyRedshifts/')

# Import functions from the custom modules
from specifications import surveySpecs
from surveyRedshiftRange import redshiftRange

# Define functions to calculate cosmological quantities
def E(z, cosmo):
    Omega_m0 = params['omega_b'] + params['omega_cdm']
    Omega_Lambda0 = params['Omega_Lambda']
    Omega_k0 = params['Omega_k']
    return np.sqrt(Omega_m0 * (1 + z) ** 3 + Omega_Lambda0 + Omega_k0 * (1 + z) ** 2)

def H(z, cosmo):
    H0 = 100 * params['h']  # H0 in km/s/Mpc
    E_z = E(z, cosmo)
    return H0 * E_z

c = 299792.458

def comoving_distance(z, cosmo):
    integrand = lambda z_prime: 1.0 / H(z_prime, cosmo)
    r, _ = quad(integrand, 0, z)
    return (c / (100 * params['h'])) * r

# Function to calculate the survey volume
def survey_volume(survey, specsSurvey, cosmo):
    z_min, z_max = redshiftRange(survey, specsSurvey)
    N_dish, V_survey, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey)

    r_min = comoving_distance(z_min, cosmo)
    r_max = comoving_distance(z_max, cosmo)

    if survey_area_sky is not None:
        solid_angle = 4 * np.pi * survey_area_sky
    else:
        solid_angle = 4 * np.pi

    volume = (solid_angle / 3) * (r_max ** 3 - r_min ** 3)

    return volume

# Defining k_min
def l_x(z, specsSurvey, cosmo):
    N_dish, V_survey, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey)
    lx = comoving_distance(z, cosmo) * np.sqrt(survey_area_sky)
    return lx


def l_z(survey, specsSurvey, cosmo):
    z_min, z_max = redshiftRange(survey, specsSurvey)
    r_min = comoving_distance(z_min, cosmo)
    r_max = comoving_distance(z_max, cosmo)
    l_z = r_max - r_min
    return l_z

def kmin(z, specsSurvey, cosmo):
    l_x_val = l_x(z, specsSurvey, cosmo)  # Corrected the function call
    l_z_val = l_z(survey, specsSurvey, cosmo)  # Assuming survey is a global variable
    k_min = 2 * np.pi / np.sqrt(2 * l_x_val ** 2 + l_z_val ** 2)
    return k_min



def Nmodes(k, z, survey, specsSurvey, cosmo):
    delta = 2 * kmin(z, specsSurvey, cosmo)
    volume = survey_volume(survey, specsSurvey, cosmo)
    N_modes = (volume * k ** 2 * delta) / (2 * np.pi) ** 3
    return N_modes

def Tsys(z):
    T_spl = 2
    nu21 = 1420
    T_CMB = 2.73
    T_gal = 25 * (408 / nu21) ** 2.75
    T_rx = 7.5 + 10 * (nu21 / 1e3 - 0.75) ** 2
    T_sys = T_spl + T_CMB + T_gal + T_rx
    return T_sys

def P_N(z, specsSurvey):
    eta = 2
    N_pol = 1
    N_dish, V_survey, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey)
    lambda21 = 0.21 * (1 + z)
    P_N = Tsys(z) * comoving_distance(z, cosmo) ** 2 * lambda21 * ((1 + z) / H(z, cosmo)) * (
                4 * np.pi * survey_area_sky / eta * N_pol * N_dish * t_total)
    return P_N


# Function for the power spectrum turnover region
def P_fit(k, kT0, m, n, P0):
    x = (np.log(k) - np.log(kT0)) / np.log(kT0)
    return np.where(k < kT0, P0 ** (1 - m * x), P0 ** (1 - n * x))

def error(k, z, survey, specsSurvey, cosmo, kT0, m, n, P0):
    err = (P_fit(k, kT0, m, n, P0) + P_N(z, specsSurvey)) / np.sqrt(Nmodes(k, z, survey, specsSurvey, cosmo))
    return err

m = 0.3  # Replace with your fiducial value for alpha
n = 0.45  # Replace with your fiducial value for beta
P0 = 3.73e4  # Replace with your fiducial value for P0
kT0 = 0.0167  # Replace with your fiducial value for k0

fiducial_values = (kT0, m, n, P0)

def calculate_fisher_matrix(k_values, z, survey, specsSurvey, cosmo, kT0, m, n, P0):
    Fisher_matrix = np.zeros((4, 4))

    for i, k in enumerate(k_values):
        dP_dm = np.where(k < kT0, -P0 ** (1 - m * ((np.log(k) - np.log(kT0)) / np.log(kT0))) * ((np.log(k) - np.log(
            kT0)) / np.log(kT0)), 0)
        dP_dn = np.where(k >= kT0, -P0 ** (1 - n * ((np.log(k) - np.log(kT0)) / np.log(kT0))), 0)
        dP_dP0 = P_fit(k, kT0, m, n, P0) / P0
        dP_dkT0 = np.where(k < kT0, -m * P_fit(k, kT0, m, n, P0) / kT0, -n * P_fit(k, kT0, m, n, P0))

        for j, k_prime in enumerate(k_values):
            err = error(k_prime, z, survey, specsSurvey, cosmo, kT0, m, n, P0)
            Fisher_matrix[0, 0] += dP_dm * dP_dm / err ** 2
            Fisher_matrix[0, 1] += dP_dm * dP_dn / err ** 2
            Fisher_matrix[0, 2] += dP_dm * dP_dP0 / err ** 2
            Fisher_matrix[0, 3] += dP_dm * dP_dkT0 / err ** 2

            Fisher_matrix[1, 1] += dP_dn * dP_dn / err ** 2
            Fisher_matrix[1, 2] += dP_dn * dP_dP0 / err ** 2
            Fisher_matrix[1, 3] += dP_dn * dP_dkT0 / err ** 2

            Fisher_matrix[2, 2] += dP_dP0 * dP_dP0 / err ** 2
            Fisher_matrix[2, 3] += dP_dP0 * dP_dkT0 / err ** 2

            Fisher_matrix[3, 3] += dP_dkT0 * dP_dkT0 / err ** 2

    Fisher_matrix = np.triu(Fisher_matrix) + np.triu(Fisher_matrix, 1).T
    return Fisher_matrix


# Define k_values and z_value
k_values = np.linspace(0.005, 0.025, 1000)
z_value = 0.5

# Define the list of surveys you want to include
surveys = ["HI IM"]

# Calculate Fisher matrices for each survey
Fisher_matrices = []
for survey in surveys:
    specsSurvey = "skaOIMBand1"
    Fisher_matrix = calculate_fisher_matrix(k_values, z_value, survey, specsSurvey, cosmo, *fiducial_values)
    Fisher_matrices.append(Fisher_matrix)

# Convert covariance matrices to MCSamples objects
samples_list = []
for Fisher_matrix in Fisher_matrices:
    covariance_matrix = np.linalg.inv(Fisher_matrix)
    samples = MCSamples(samples=np.random.multivariate_normal(mean=fiducial_values, cov=covariance_matrix, size=100000),
                        names=[r'$\alpha$', r'$\beta$', r'$P0$', r'$k_{0}$'])
    samples_list.append(samples)

# Set up plot settings
plt.figure(figsize=(15, 12))
g = plots.get_subplot_plotter(subplot_size=4, width_inch=10)
g.settings.axes_fontsize = 16
g.settings.legend_fontsize = 16
g.settings.lab_fontsize = 18
g.settings.alpha_filled_add = 0.80
g.settings.figure_legend_loc = 'upper right'
g.settings.figure_legend_ncol = 2

# Plot results
for i, survey in enumerate(surveys):
    g.triangle_plot(samples_list[i], filled=True, contour_colors=['blue'], title_limit=1,
                    legend_labels=[survey], plot_3d_with_param='P_0')

plt.tight_layout()
plt.savefig('fisher_forecast.pdf')
plt.show()
