import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Enable LaTeX for all text in the plot
plt.rcParams.update({
    'text.usetex': True,    # Use LaTeX for text rendering
    'font.size': 28,            # Default text size
    'axes.titlesize': 30,    # Axes title size
    'axes.labelsize': 30,   # Axes labels size
    'xtick.labelsize': 30,   # X-axis tick labels size
    'ytick.labelsize': 30,   # Y-axis tick labels size
    'legend.fontsize': 28,  # Legend font size
    'text.latex.preamble': r'\usepackage{amsmath}' # Use amsmath package for LaTeX formatting
})

# Define the file paths (update them according to your local file paths)
beam_file = "/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z0.5.dat"
nobeam_file = "/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_NoBeam_z0.5.dat"

# Read the data from each file
beam_data = np.loadtxt(beam_file)
nobeam_data = np.loadtxt(nobeam_file)

# Extract the k0 values (second column) from each dataset
beam_k0 = beam_data[:, 1]  # Assuming k0 is the second column
nobeam_k0 = nobeam_data[:, 1]  # Assuming k0 is the second column

# 2. Overlayed Gaussian Fit
plt.figure(figsize=(15, 12))

# Fit Gaussian to Beam data
mu_beam, std_beam = norm.fit(beam_k0)
x_beam = np.linspace(min(beam_k0), max(beam_k0), 100)
p_beam = norm.pdf(x_beam, mu_beam, std_beam)

# Fit Gaussian to No Beam data
mu_nobeam, std_nobeam = norm.fit(nobeam_k0)
x_nobeam = np.linspace(min(nobeam_k0), max(nobeam_k0), 100)
p_nobeam = norm.pdf(x_nobeam, mu_nobeam, std_nobeam)

# Plot histograms
plt.hist(beam_k0, bins=50, density=True, alpha=0.6, color='green', edgecolor='black', label='Considering Beam Effect')
plt.hist(nobeam_k0, bins=50, density=True, alpha=0.6, color='orange', edgecolor='black', label='Ignoring Beam Effect')

# Overlay Gaussian curves
plt.plot(x_beam, p_beam, 'k--', linewidth=2, label=f'$\\mu={mu_beam:.4f}$, $\\sigma={std_beam:.4f}$')
plt.plot(x_nobeam, p_nobeam, 'r--', linewidth=2, label=f'$\\mu={mu_nobeam:.4f}$, $\\sigma={std_nobeam:.4f}$')

# Add labels and title
plt.xlabel(r'$k_0$')
plt.ylabel('Probability Density')
plt.legend()

# Customize spines for a bold look
ax = plt.gca()
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)

# Customize tick parameters
ax.tick_params(axis='both', which='major', length=10, width=3)  # Increase length and width of major ticks
ax.tick_params(axis='both', which='minor', length=10, width=3)   # Customize minor ticks if necessary

# Save the plot to a file with high resolution
plt.savefig("/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/BeamGaussianFit.pdf", bbox_inches='tight', dpi=300)

# Show the plot
#plt.title('Gaussian Fit Over Histogram')
#plt.show()
