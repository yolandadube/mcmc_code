import numpy as np
import matplotlib.pyplot as plt
import os

# Set global plot parameters for A4 document
def set_plot_params():
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
        'font.size': 14,
        'axes.titlesize': 16,
        'axes.labelsize': 20,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 22,
        'xtick.major.size': 5,
        'xtick.major.width': 2.5,
        'ytick.major.size': 5,
        'ytick.major.width': 2.5,
        'xtick.major.pad': 5,
        'ytick.major.pad': 5,
        'xtick.minor.size': 3,
        'xtick.minor.width': 1,
        'ytick.minor.size': 3,
        'ytick.minor.width': 1,
        'figure.dpi': 300,
        'lines.solid_joinstyle': 'miter',
        'lines.antialiased': True,
        'ps.useafm': True,
        'pdf.use14corefonts': True,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath, amssymb}',
        'axes.linewidth': 1.5,
    })

# Set plot parameters
set_plot_params()

def plot_k0_vs_kfg_for_surveys(base_path, k_fg_values, destination_folder, surveys, statistic='median', log_scale='both'):
    """
    :param log_scale: Can be 'x', 'y', or 'both' to set the respective axes to logarithmic scale.
    """
    # Ensure the destination folder exists
    os.makedirs(destination_folder, exist_ok=True)

    # Initialize the plot with a smaller size suitable for A4 page
    plt.figure(figsize=(10, 8), dpi=300)

    # Define line styles, markers, and colors for each survey
    line_styles = ['-', '--', '-.']
    markers = ['o', 's', '^']
    colors = ['red', 'green', 'blue']

    for index, survey in enumerate(surveys):
        k0_statistics = []
        k0_lower = []
        k0_upper = []

        # Loop through each k_fg value and extract k0 statistics
        for k_fg in k_fg_values:
            file_path = os.path.join(base_path, survey['file_template'].format(k_fg=k_fg))
            data = np.loadtxt(file_path)
            k0_samples = data[:, 1]

            # Compute the median for k0
            k0_stat = np.median(k0_samples)
            k0_statistics.append(k0_stat)

            # Compute the 1-sigma error bounds
            k0_16, k0_84 = np.percentile(k0_samples, [16, 84])
            k0_lower.append(k0_stat - k0_16)
            k0_upper.append(k0_84 - k0_stat)

        # Customize the line style, marker, and color for each survey
        plt.errorbar(k_fg_values, k0_statistics, yerr=[k0_lower, k0_upper], fmt=markers[index], linestyle=line_styles[index], color=colors[index], capsize=5, elinewidth=1.5, linewidth=2.5, label=rf'{survey["name"]}')

    # Set axis labels and legend
    plt.xlabel(r'$k_{fg} \, [h/\text{Mpc}]$', fontsize=20)
    plt.ylabel(r'$k_0 \, [h/\text{Mpc}]$', fontsize=20)
    plt.legend(fontsize=16)

    # Optional: Set axes to logarithmic scale
    if 'x' in log_scale:
        plt.xscale('log')
    if 'y' in log_scale:
        plt.yscale('log')

    # Other plot customizations
    plt.axhline(y=0.0163, color='grey', linestyle='--', linewidth=2, label='k0 = 0.0163')

    # Apply spine settings directly to the axes object
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)  # Standard frame thickness for A4

    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(destination_folder, 'combined_k0_vs_kfg_plot_log_scale.pdf'), dpi=300)
    plt.close()

# Example usage
base_path = '/home/yolanda/Documents/mcmc_code/version1/functions/ResidualData/'
k_fg_values = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]  # Ensure these are > 0 for log scale
destination_folder = os.path.join(base_path, 'k0_vs_kfg_values_log_scale')

surveys = [
    {"name": "SKA-MID Band 1", "file_template": "MCMC_skaOIMBand1skaOIMBand1_{k_fg:.3f}.dat"},
    {"name": "MeerKAT L-Band", "file_template": "MCMC_meerKATLBandmeerKATLBand_{k_fg:.3f}.dat"},
    {"name": "MeerKAT UHF-Band", "file_template": "MCMC_meerKATUHFBandmeerKATUHFBand_{k_fg:.3f}.dat"}
]

plot_k0_vs_kfg_for_surveys(base_path, k_fg_values, destination_folder, surveys, statistic='median', log_scale='both')
