import sys, os, glob, re
import numpy as np
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.font_manager
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
from getdist import plots, MCSamples
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from matplotlib.ticker import LogLocator, FuncFormatter
from scipy.interpolate import CubicSpline
from matplotlib.ticker import LogFormatter 
from matplotlib.pyplot import rc, rcParams
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
import pandas as pd  # Ensure pandas is imported

#--------------------------------------------------------------------------------------------------------------------------

# Set global plot parameters for A4 document
plt.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'font.size': 20,
    'axes.titlesize': 35,
    'axes.labelsize': 30,
    'xtick.labelsize': 30,
    'ytick.labelsize': 30,
    'legend.fontsize': 28,
    'xtick.major.size': 5,
    'xtick.major.width': 2.5,
    'ytick.major.size': 5,
    'ytick.major.width': 2.5,
    'xtick.major.pad': 5,
    'ytick.major.pad': 5,
    'xtick.minor.size': 3,
    'xtick.minor.width': 2,
    'ytick.minor.size': 3,
    'ytick.minor.width': 2,
    'figure.dpi': 500,
    'lines.solid_joinstyle': 'miter',
    'lines.antialiased': True,
    'ps.useafm': True,
    'pdf.use14corefonts': True,
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsmath, amssymb}',
    'axes.linewidth': 2.5,
})

fig = plt.figure(figsize=(10, 8))

plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2.5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2.5
plt.rcParams['xtick.major.pad'] = 5
plt.rcParams['ytick.major.pad'] = 5
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.width'] = 2

plt.rcParams['figure.dpi'] = 300
plt.rcParams['lines.solid_joinstyle'] = 'miter'
plt.rcParams['lines.antialiased'] = True
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.tick_params(axis='both', which='both', direction='in')
plt.rcParams['axes.linewidth'] = 2.5  # set the value globally

# Define the file paths
beam_files_pattern = "/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z*.dat"
nobeam_files_pattern = "/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_NoBeam_z*.dat"
beam_files = glob.glob(beam_files_pattern)
nobeam_files = glob.glob(nobeam_files_pattern)

# Extract k0 values and corresponding redshifts
beam_k0_data = []
nobeam_k0_data = []

# Function to extract k0 and z from a file
def extract_k0_and_z(file_path):
    data = np.loadtxt(file_path)
    k0_values = data[:, 1]  # Assuming k0 is the second column
    
    # Extract redshift z from the filename using regex
    match = re.search(r'_z(\d+\.\d+)', file_path)
    z_value = float(match.group(1)) if match else None
    return k0_values, z_value

# Process beam files
for file in beam_files:
    k0_values, z = extract_k0_and_z(file)
    beam_k0_data.append((k0_values, z))

# Process no beam files
for file in nobeam_files:
    k0_values, z = extract_k0_and_z(file)
    nobeam_k0_data.append((k0_values, z))

# Combine data into a DataFrame for easy plotting
data = []
for k0_values, z in beam_k0_data:
    for k0 in k0_values:
        data.append({'k0': k0, 'z': z, 'Beam': 'Considering Beam Effect'})
for k0_values, z in nobeam_k0_data:
    for k0 in k0_values:
        data.append({'k0': k0, 'z': z, 'Beam': 'Ignoring Beam Effect'})

df = pd.DataFrame(data)

# Define the specific z values and their corresponding labels
z_values = [0.40, 0.50, 0.55, 0.70, 1.00, 2.50]
z_labels = [r'$0.40$', r'$0.50$', r'$0.55$', r'$0.70$', r'$1.00$', r'$2.50$']

# Plotting using Seaborn
#plt.figure(figsize=(8, 6), dpi=100)

# Choose a bright color palette for better visibility
palette = sns.color_palette("bright", 2)

# Boxplot with appropriate line widths and vibrant colors
sns.boxplot(x='z', y='k0', hue='Beam', data=df, palette=palette, linewidth=1.5)

# Add a dotted line for the fiducial value of k0 = 0.0163
plt.axhline(0.0163, color='grey', linestyle='--', linewidth=1.5, label=r'Fiducial $k_0=0.0163$')

# Use LaTeX formatting for x and y labels
plt.xlabel(r'$z$', fontsize=30)
plt.ylabel(r'$k_0$', fontsize=30)

# Format the tick labels with LaTeX and fixed decimal places
def latex_ticks_x(x, pos):
    return r'${:.2f}$'.format(z_values[int(x)])

def latex_ticks_y(y, pos):
    return r'${:.4f}$'.format(y)

ax = plt.gca()
ax.xaxis.set_major_formatter(FuncFormatter(latex_ticks_x))
ax.yaxis.set_major_formatter(FuncFormatter(latex_ticks_y))

# Customize spines and ticks for a more prominent look
ax.tick_params(axis='both', which='major', length=3, width=2.5)
ax.tick_params(axis='both', which='minor', length=3, width=2.5)

# Apply spine settings directly to the axes object
for spine in ax.spines.values():
    spine.set_linewidth(2.5)

# Set the xticks to the specific redshift values
ax.set_xticks(range(len(z_values)))
ax.set_xticklabels(z_labels)

plt.tight_layout()

# Save the plot to a file with high resolution
plt.savefig("Boxplot_k0_BeamAnalysis_A4.pdf", bbox_inches='tight')

# Show the plot
plt.show()

# Print basic statistics for comparison
for k0_values, z in beam_k0_data:
    k0_mean = np.mean(k0_values)
    k0_std = np.std(k0_values)
    k0_median = np.median(k0_values)
    print(f"Beam Analysis z={z:.2f} - Mean $k_0$: {k0_mean:.4f}, Median $k_0$: {k0_median:.4f}, Std Dev: {k0_std:.4f}")

for k0_values, z in nobeam_k0_data:
    k0_mean = np.mean(k0_values)
    k0_std = np.std(k0_values)
    k0_median = np.median(k0_values)
    print(f"No Beam Analysis z={z:.2f} - Mean $k_0$: {k0_mean:.4f}, Median $k_0$: {k0_median:.4f}, Std Dev: {k0_std:.4f}")
