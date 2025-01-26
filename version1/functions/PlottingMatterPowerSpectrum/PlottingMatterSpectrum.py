import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from scipy import interpolate

# Set global plot parameters
def set_plot_params():
    plt.rcParams.update({
        'figure.facecolor': (186/255, 215/255, 230/255),
        'axes.facecolor': (186/255, 215/255, 230/255),
        'savefig.facecolor': (186/255, 215/255, 230/255),
        'font.size': 200,
        'axes.titlesize': 220,
        'axes.labelsize': 220,
        'xtick.labelsize': 200,
        'ytick.labelsize': 200,
        'legend.fontsize': 220,
        'xtick.major.size': 30,
        'xtick.major.width': 25,
        'ytick.major.size': 30,
        'ytick.major.width': 25,
        'xtick.major.pad': 30,
        'ytick.major.pad': 30,
        'xtick.minor.size': 20,
        'xtick.minor.width': 8,
        'ytick.minor.size': 25,
        'ytick.minor.width': 8,
        'figure.dpi': 300,
        'lines.solid_joinstyle': 'miter',
        'lines.antialiased': True,
        'ps.useafm': True,
        'pdf.use14corefonts': True,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{amsmath, amssymb}',
        'axes.linewidth': 6,
    })

# Function to check if the files exist
def check_files_exist(file_paths):
    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            return False
    return True

def plot_power_spectrum(file_paths, output_path, dpi=300, k_values=None, custom_k_values=None, smooth_curve=True):
    """
    Plot and save the power spectrum as a high-quality image.

    Parameters:
    - file_paths: List of paths to the files containing the power spectrum data.
    - output_path: Path to save the high-quality image of the plot.
    - dpi: Dots per inch (resolution of the saved plot). Default is 300.
    - k_values: Array of wave numbers from the data. Default is None.
    - custom_k_values: Array of custom wave numbers. Default is None.
    - smooth_curve: Boolean indicating whether to plot a smooth curve. Default is True.
    """
    
    # Append the directory containing the module to the Python path
    module_directory = '/home/yolanda/Documents/mcmc_code/version1/functions/matterPowerSpectrum/'
    sys.path.append(module_directory)

    # Import the pk function directly from matterPk module
    from matterPk import pk

    # Check if the files exist
    if not check_files_exist(file_paths):
        return

    # Set plot parameters
    set_plot_params()

    # Initialize plot
    fig, ax = plt.subplots(figsize=(60, 50))  # Increased size for A0 poster

    # Plot power spectrum for each file
    for i, file_path in enumerate(file_paths):
        # Load the data
        if k_values is None:
            data = np.loadtxt(file_path)
            k_values = data[:, 0]  # Wave numbers

        # If custom_k_values are provided, use them
        if custom_k_values is not None:
            k_values = custom_k_values

        # Get power spectrum values using the provided function
        Pk_values = pk(file_path, k_values)

        # Perform interpolation to create a smooth curve
        if smooth_curve:
            interp_func = interpolate.interp1d(k_values, Pk_values, kind='cubic')
            smooth_k_values = np.linspace(k_values.min(), k_values.max(), 1000)
            smooth_Pk_values = interp_func(smooth_k_values)
        else:
            smooth_k_values, smooth_Pk_values = k_values, Pk_values

        # Find the peak point
        peak_index = np.argmax(Pk_values)
        peak_k = k_values[peak_index]
        peak_Pk = Pk_values[peak_index]

        # Print the peak values
        print(f"Peak at k = {peak_k:.5f}, P(k) = {peak_Pk:.5e}")

        # Plotting
        label = 'Linear Power Spectrum' if i == 0 else 'Non-Linear Power Spectrum'
        linestyle = '-' if i == 0 else '--'
        
        ax.plot(smooth_k_values, smooth_Pk_values, label=label, linestyle=linestyle, linewidth=30.0)  # Increased line width
        if i == 1:
            ax.scatter(peak_k, peak_Pk, color='black', s=500, label='Peak')  # Plot the peak point with increased size
            # Adjust annotation position to ensure it's within the plot area
            ax.annotate('Turnover', xy=(peak_k, peak_Pk), xytext=(peak_k * 1.5, peak_Pk * 1.5),
                         arrowprops=dict(facecolor='red', arrowstyle='->'), fontsize=100, color='red')  # Red arrow and text

    # Configure plot
    ax.set_xlabel('$k$ [h/Mpc]')
    ax.set_ylabel('$P(k)$ [$h^{-3}$ Mpc$^3$]')
    ax.legend()
    ax.set_xscale('log')  # Set x-axis to logarithmic scale
    ax.set_yscale('log')  # Set y-axis to logarithmic scale

    # Apply spine settings directly to the axes object
    for spine in ax.spines.values():
        spine.set_linewidth(40)  # Increased frame thickness

    plt.tight_layout()  # Adjust layout for better spacing

    # Save the plot as a high-quality image
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()  # Close the plot to free up memory

# Define custom k_values array
k_values = np.logspace(-3, 0)

# Example usage
file_paths = [
    '/home/yolanda/Documents/mcmc_code/version1/result/PowerSpectrumDatFiles/CREATED00_pk_fnl10.dat',
    #'/home/yolanda/Documents/mcmc_code/version1/result/PowerSpectrumDatFiles/CREATED00_pk_nl.dat'
]
output_path = '/home/yolanda/Documents/mcmc_code/version1/result/PowerSpectrumDatFiles/power_spectrum_plot.pdf'
plot_power_spectrum(file_paths, output_path, k_values=k_values, smooth_curve=True)
