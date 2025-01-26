import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Enable LaTeX for all text in the plot
plt.rcParams.update({
    'text.usetex': True,   # Use LaTeX for text rendering
    'font.size': 18,       # Default text size
    'axes.titlesize': 22,  # Axes title size
    'axes.labelsize': 20,  # Axes labels size
    'xtick.labelsize': 18, # X-axis tick labels size
    'ytick.labelsize': 18, # Y-axis tick labels size
    'legend.fontsize': 18, # Legend font size
    'text.latex.preamble': r'\usepackage{amsmath}'  # Use amsmath package for additional LaTeX features
})

# Define the file path pattern to get all relevant files
file_path_pattern = "/home/yolanda/Documents/mcmc_code/version1/result/Fnl_vs_Redshift/*.dat"

# Find all files matching the pattern
file_list = glob.glob(file_path_pattern)

# Sort the file list to ensure files are processed in order
file_list.sort()

# Initialize a dictionary to store data for the heatmap
data = {'fNL': [], 'Redshift': [], 'k0': []}

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
        k0_mean = np.mean(file_data[:, 1])  # Adjust if k0 is not the second column
        
        # Store the k0 mean for this fnl and z
        data['fNL'].append(fnl_value)
        data['Redshift'].append(z_value)
        data['k0'].append(k0_mean)

# Convert to DataFrame
df = pd.DataFrame(data)

# Pivot table for heatmap
pivot_df = df.pivot_table(values='k0', index='fNL', columns='Redshift', aggfunc='mean')

# Plotting
plt.figure(figsize=(12, 8), dpi=300)

# Create a heatmap
sns.heatmap(pivot_df, annot=True, cmap='coolwarm', cbar_kws={'label': r'$k_0$'})

# Add labels and title
plt.xlabel(r'Redshift $z$', fontsize=22)
plt.ylabel(r'$f_{NL}$', fontsize=22)
#plt.title(r'Heatmap of $k_0$ for Different $f_{NL}$ and Redshift', fontsize=24)

# Customize ticks and spines
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)

# Make spines bold
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)

# Customize tick parameters
ax.tick_params(axis='both', which='major', length=5, width=2)
ax.tick_params(axis='both', which='minor', length=3, width=1)

# Save the plot
plt.savefig("fnl_vs_redshift_heatmap.pdf", bbox_inches='tight', dpi=300)

# Show the plot
plt.show()
