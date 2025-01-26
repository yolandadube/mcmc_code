import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Enable LaTeX for all text in the plot
plt.rcParams.update({
    'text.usetex': True,        # Use LaTeX for text rendering
    'font.size': 18,            # Default text size
    'axes.titlesize': 22,       # Axes title size
    'axes.labelsize': 20,       # Axes labels size
    'xtick.labelsize': 18,      # X-axis tick labels size
    'ytick.labelsize': 18,      # Y-axis tick labels size
    'legend.fontsize': 18,      # Legend font size
    'text.latex.preamble': r'\usepackage{amsmath}' # Use amsmath package for additional LaTeX features
})

# Define the file path pattern to get all relevant files
file_path_pattern = "/home/yolanda/Documents/mcmc_code/version1/result/Fnl_vs_Redshift/*.dat"

# Find all files matching the pattern
file_list = glob.glob(file_path_pattern)

# Sort the file list to ensure files are processed in order
file_list.sort()

# Initialize a dictionary to store data for each fnl value
data = {1: [], 2: [], 5: [], 10: []}

# Regular expression to extract fnl and z from the filename
regex = r'_fnl(\d+)_z(\d+\.\d+)\.dat'

# Iterate over the files and extract data
for file_path in file_list:
    file_name = os.path.basename(file_path)
    match = re.search(regex, file_name)
    
    if match:
        fnl_value = int(match.group(1))
        z_value = float(match.group(2))
        
        # Load the data from the file
        file_data = np.loadtxt(file_path)
        
        # Assuming k0 is the second column
        k0_values = file_data[:, 1]  # Adjust if k0 is not the second column
        k0_mean = np.mean(k0_values)
        
        # Store the mean k0 value for this fnl and z
        if fnl_value in data:
            data[fnl_value].append((z_value, k0_mean))

# Plotting
plt.figure(figsize=(12, 8), dpi=300)

# Color palette for different fnl values
colors = {1: 'blue', 2: 'green', 5: 'red', 10: 'purple'}

# Plot each fnl series
for fnl_value, values in data.items():
    values.sort(key=lambda x: x[0])  # Sort by redshift
    z_values, k0_means = zip(*values)
    plt.plot(z_values, k0_means, marker='o', linestyle='-', color=colors[fnl_value], linewidth=2, markersize=8, label=f'$f_{{NL}}={fnl_value}$')

# Add labels and title
plt.xlabel(r'$z$', fontsize=20)
plt.ylabel(r'$k_0$', fontsize=20)
#plt.title(r'Variation of $k_0$ with Redshift for Different $f_{NL}$ Values', fontsize=24)
plt.legend()

# Customize ticks and spines
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.grid(True, which='both', linestyle='--', linewidth=1)

# Make spines bold
ax = plt.gca()
ax.spines['top'].set_linewidth(4)
ax.spines['right'].set_linewidth(4)
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)

# Customize tick parameters
ax.tick_params(axis='both', which='major', length=8, width=2.5)
ax.tick_params(axis='both', which='minor', length=8, width=2.5)

# Save the plot
plt.savefig("fnl_vs_redshift_latex_plot.pdf", bbox_inches='tight', dpi=300)

# Show the plot
plt.show()

# Print data for verification
for fnl_value, values in data.items():
    print(f"\n$f_{{NL}}={fnl_value}$:")
    for z, k0_mean in values:
        print(f"Redshift $z={z:.2f}$, Mean $k_0={k0_mean:.4f}$")
