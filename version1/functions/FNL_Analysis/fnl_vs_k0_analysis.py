import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Set global plot parameters
def set_plot_params():
    plt.rcParams.update({
        'text.usetex': True,  # Use LaTeX for text rendering
        'font.size': 18,      # Default text size
        'axes.titlesize': 22, # Axes title size
        'axes.labelsize': 20, # Axes labels size
        'xtick.labelsize': 18,# X-axis tick labels size
        'ytick.labelsize': 18,# Y-axis tick labels size
        'legend.fontsize': 18,# Legend font size
        'text.latex.preamble': r'\usepackage{amsmath}', # Use amsmath package for additional LaTeX features
        'figure.dpi': 300,    # High resolution for printing
        'lines.solid_joinstyle': 'miter',
        'lines.antialiased': True,
        'ps.useafm': True,
        'pdf.use14corefonts': True,
        'axes.linewidth': 3,  # Spine thickness
        'xtick.major.size': 5,  # Major tick size
        'xtick.major.width': 2,  # Major tick width
        'ytick.major.size': 5,  # Major tick size
        'ytick.major.width': 2,  # Major tick width
        'xtick.minor.size': 5,  # Minor tick size
        'xtick.minor.width': 2,  # Minor tick width
        'ytick.minor.size': 5,  # Minor tick size
        'ytick.minor.width': 2   # Minor tick width
    })

# Set plot parameters
set_plot_params()

# Define the file path pattern to get all relevant files
file_path_pattern = "/home/yolanda/Documents/mcmc_code/version1/result/k0_vs_fnl/MCMC_skaOIMBand1skaOIMBand1_0.001_fnl*.dat"

# Find all files matching the pattern
file_list = glob.glob(file_path_pattern)

# Sort the file list to ensure files are processed in order
file_list.sort()

# Initialize lists to store f_NL values and corresponding k_0 means
fnl_values = []
k0_means = []

# Iterate over the files
for file_path in file_list:
    # Extract the f_NL value from the file name
    file_name = os.path.basename(file_path)
    # Detect negative or positive f_NL value from the filename
    if "m" in file_name.split('fnl')[-1].split('.')[0]:
        fnl_value = -int(file_name.split('fnl')[-1].split('.')[0].replace("m", ""))
    else:
        fnl_value = int(file_name.split('fnl')[-1].split('.')[0])
    
    fnl_values.append(fnl_value)
    
    # Load the data from the file
    data = np.loadtxt(file_path)
    
    # Extract the k0 values (assuming k0 is the second column)
    k0_values = data[:, 1]  # Adjust if k0 is not the second column
    
    # Compute the mean k0 for this file
    k0_mean = np.mean(k0_values)
    k0_means.append(k0_mean)

# Plot the fnl vs. k0 mean
plt.figure(figsize=(12, 8))

# Separate positive and negative f_NL values for different colors
positive_indices = [i for i, x in enumerate(fnl_values) if x >= 0]
negative_indices = [i for i, x in enumerate(fnl_values) if x < 0]

positive_fnl_values = [fnl_values[i] for i in positive_indices]
positive_k0_means = [k0_means[i] for i in positive_indices]

negative_fnl_values = [fnl_values[i] for i in negative_indices]
negative_k0_means = [k0_means[i] for i in negative_indices]

# Plot positive f_NL values
plt.plot(positive_fnl_values, positive_k0_means, marker='o', linestyle='solid', color='blue', linewidth=6, markersize=10, label=r'$\mathrm{Positive\ f_{NL}}$')

# Plot negative f_NL values
plt.plot(negative_fnl_values, negative_k0_means, marker='o', linestyle='solid', color='red', linewidth=6, markersize=10, label=r'$\mathrm{Negative\ f_{NL}}$')

# Add labels and title
plt.xlabel(r'$f_{NL}$', fontsize=22)
plt.ylabel(r'$k_0$', fontsize=22)
#plt.title(r'Effect of $f_{NL}$ on $k_0$', fontsize=24)
plt.legend()

# Customize ticks and spines
ax = plt.gca()

# Customize tick parameters
ax.tick_params(axis='both', which='major', length=5, width=2)
ax.tick_params(axis='both', which='minor', length=5, width=2)

plt.tight_layout()  # Adjust layout for better spacing

# Save the plot
plt.savefig("fnl_vs_k0_latex_plot.pdf", bbox_inches='tight')

# Show the plot
plt.show()

# Print basic statistics for each f_NL
for fnl, k0_mean in zip(fnl_values, k0_means):
    print(f"$f_{{NL}}$: {fnl}, Mean $k_0$: {k0_mean:.4f}")
