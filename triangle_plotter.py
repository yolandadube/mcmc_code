import os
import json
import shutil
import getdist
import numpy as np
import pandas as pd
from scipy.linalg import inv
import matplotlib.pyplot as plt
from matplotlib import rcParams
from getdist import plots, MCSamples


# Set plotting parameters for high-quality output
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family']     = 'serif'
rcParams['font.serif']      = ['Times New Roman']
rcParams['font.size']       = 36
rcParams['axes.labelsize']  = 38
rcParams['axes.linewidth']  = 3.5
rcParams['xtick.labelsize'] = 34
rcParams['ytick.labelsize'] = 34
rcParams['legend.fontsize'] = 34
rcParams['figure.figsize']  = (15, 12)
rcParams['text.usetex']     = True  # Enable LaTeX formatting

# Define standard colors for each survey
SURVEY_COLORS = {
    "SKA-MID Band 1"  : "red",
    "MeerKAT L-band"  : "blue",
    "MeerKAT UHF-band": "green"
}

# Define survey redshift ranges - adding this missing definition
SURVEY_REDSHIFT_RANGES = {
    "SKA-MID Band 1"  : (0.35, 3.00),
    "MeerKAT L-band"  : (0.20, 0.58),
    "MeerKAT UHF-band": (0.40, 1.45)
}

# Define the k_fg values to process - now includes lower values
k_fg_values = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01]

# Base directories
input_dir = "MCMC_Fisher_FixedWidth_Results"  # Updated to match new directory structure
output_dir = "Fisher_Triangle_Analysis"
bin_width = 0.1  # Fixed bin width used

# Create output directory
os.makedirs(output_dir, exist_ok=True)

def correlation_coefficient(cov, i, j):
    """Calculate correlation coefficient between parameters i and j"""
    return cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])

def read_direct_from_cov_file(file_path):
    """
    Read covariance matrix directly from text file format.
    Format example:
        alpha: 7.625170e-02  -4.138278e-03  -3.453290e-05
        beta: -4.138278e-03  1.687200e-03  1.328934e-05
          k0: -3.453290e-05  1.328934e-05  1.108962e-07
    """
    try:
        # Read all lines from file
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find the lines with the covariance matrix
        data_lines = []
        start_reading = False
        for line in lines:
            line = line.strip()
            if start_reading and line.startswith('#'):
                break  # Stop at the next comment after matrix
            
            # Start reading at the alpha line
            if 'alpha:' in line:
                start_reading = True
                data_lines.append(line)
            elif start_reading:
                data_lines.append(line)
        
        # Create empty 3x3 matrix
        cov_matrix = np.zeros((3, 3))
        
        # Fill the matrix from data lines
        for i, line in enumerate(data_lines):
            # Get entries after the colon
            entries = line.split(':')[1].strip().split()
            for j, entry in enumerate(entries):
                cov_matrix[i, j] = float(entry)
        
        return cov_matrix
    
    except Exception as e:
        print(f"Error reading cov file {file_path}: {e}")
        return None

# Add this function to combine Fisher matrices with volume weighting
def combine_fisher_matrices(covariance_matrices, epsilon=1e-10, bin_volumes=None):
    """
    Combine results by first converting covariances to Fisher matrices,
    then applying volume weighting before summing.

    Parameters:
    -----------
    covariance_matrices : list of 2D arrays
        List of covariance matrices to combine
    epsilon : float
        Small value to add to diagonal for numerical stability
    bin_volumes : list of float, optional
        Volumes of each redshift bin for proper weighting

    Returns:
    --------
    combined_cov : 2D array
        Combined covariance matrix
    """
    fisher_matrices = []
    good_indices = []  # Track which indices contain valid matrices

    # Print diagnostic information about input matrices
    print(f"  Combining {len(covariance_matrices)} matrices")
    
    # If no volumes provided, use uniform weighting
    if (bin_volumes is None):
        bin_volumes = [1.0] * len(covariance_matrices)
    
    # Compute total volume for normalization
    total_volume = sum(bin_volumes)
    print(f"  Total survey volume: {total_volume:.2e}")
    
    # Check the parameter correlations in input matrices
    for i, cov_matrix in enumerate(covariance_matrices):
        cov = np.array(cov_matrix)  # FIXED: Changed 'cov' to 'cov_matrix'

        # Skip invalid matrices
        if np.any(np.isnan(cov)) or np.linalg.det(cov) < 1e-15:
            print(f"  Skipping invalid matrix {i}")
            continue

        # Print correlation between alpha and beta (parameters 0 and 1)
        alpha_beta_corr = correlation_coefficient(cov, 0, 1)
        print(f"  Matrix {i}: alpha-beta correlation = {alpha_beta_corr:.4f}, volume = {bin_volumes[i]:.2e}")

        # Condition the matrix for numerical stability
        diag_vals = np.diag(cov)
        regularization = np.diag(diag_vals * epsilon)
        cov_reg = cov + regularization

        # Check eigenvalues for degeneracy diagnosis
        eigenvalues = np.linalg.eigvals(cov_reg)
        condition_number = np.max(np.abs(eigenvalues)) / np.min(np.abs(eigenvalues))
        print(f"  Matrix {i}: condition number = {condition_number:.2e}")

        # Convert to Fisher matrix
        try:
            fisher = inv(cov_reg)
            fisher_matrices.append(fisher)
            good_indices.append(i)  # Track which indices were used
        except np.linalg.LinAlgError:
            print(f"  Matrix inversion failed for matrix {i}")
            continue

    if not fisher_matrices:
        return None

    # Get the volumes for the matrices we're actually using
    used_volumes = [bin_volumes[i] for i in good_indices]
    
    # Combine Fisher matrices using volume weighting
    combined_fisher = np.zeros_like(fisher_matrices[0])
    
    for i, fisher in enumerate(fisher_matrices):
        # Weight by volume fraction of the good matrices
        volume_weight = used_volumes[i] / sum(used_volumes)
        print(f"  Matrix {good_indices[i]}: volume weight = {volume_weight:.4f}")
        combined_fisher += fisher * volume_weight
    
    # Convert back to covariance matrix
    try:
        combined_cov = inv(combined_fisher)

        # Print diagnostic information about the combined matrix
        alpha_beta_corr = correlation_coefficient(combined_cov, 0, 1)
        print(f"  Combined matrix: alpha-beta correlation = {alpha_beta_corr:.4f}")

        eigenvalues = np.linalg.eigvals(combined_cov)
        condition_number = np.max(np.abs(eigenvalues)) / np.min(np.abs(eigenvalues))
        print(f"  Combined matrix: condition number = {condition_number:.2e}")

        # Print the uncertainties from the combined matrix
        uncertainties = np.sqrt(np.diag(combined_cov))
        print(f"  Combined uncertainties: {uncertainties}")

        return combined_cov

    except np.linalg.LinAlgError:
        print("  Failed to invert combined Fisher matrix")
        return None

# In the read_covariance_files function, add ability to return bin volumes
def read_covariance_files(bin_dir, survey_file_name):
    """
    Read and parse all covariance files for a survey in a directory.
    Tries different file formats.
    
    Returns:
    --------
    tuple: (covariance_matrices, bin_volumes)
        The covariance matrices and corresponding bin volumes
    """
    cov_matrices = []
    bin_volumes = []
    
    # 1. First try the JSON file
    json_file = os.path.join(bin_dir, f"{survey_file_name}_results.json")
    if (os.path.exists(json_file)):
        try:
            with open(json_file, 'r') as f:
                results = json.load(f)
            
            # Extract both matrices and volumes if available
            matrices = results.get("covariance_matrices", [])
            volumes = results.get("bin_volumes", [1.0] * len(matrices))
            
            # If volumes list is shorter than matrices, pad with last value
            if (len(volumes) < len(matrices)):
                volumes.extend([volumes[-1]] * (len(matrices) - len(volumes)))
                
            return matrices, volumes
        except Exception as e:
            print(f"Error reading JSON file {json_file}: {e}")
    
    # 2. Try NPZ file as backup
    npz_file = os.path.join(bin_dir, f"{survey_file_name}_results.npz")
    if (os.path.exists(npz_file)):
        try:
            results_npz = np.load(npz_file, allow_pickle=True)
            if ("covariance_matrices" in results_npz):
                matrices = results_npz["covariance_matrices"].tolist()
                
                # Try to extract volumes, default to uniform if not found
                if ("bin_volumes" in results_npz):
                    volumes = results_npz["bin_volumes"].tolist()
                else:
                    volumes = [1.0] * len(matrices)
                    
                return matrices, volumes
        except Exception as e:
            print(f"Error reading NPZ file {npz_file}: {e}")
    
    # 3. Look for individual covariance text files
    print(f"  Looking for individual covariance files...")
    for file in os.listdir(bin_dir):
        if (file.startswith(survey_file_name) and file.endswith("_cov.txt")):
            file_path = os.path.join(bin_dir, file)
            print(f"  Reading covariance from {file_path}")
            cov_matrix = read_direct_from_cov_file(file_path)
            if (cov_matrix is not None):
                cov_matrices.append(cov_matrix)
                
                # Try to extract bin volume from file name or content
                volume = 1.0  # Default volume
                
                # Try to extract volume from the file content
                try:
                    with open(file_path, 'r') as f:
                        for line in f:
                            if ("Bin volume =" in line):
                                vol_str = line.split("=")[1].strip().split()[0]
                                volume = float(vol_str)
                                break
                except:
                    pass
                
                bin_volumes.append(volume)
    
    if (cov_matrices):
        print(f"  Found {len(cov_matrices)} individual covariance files")
        return cov_matrices, bin_volumes
    
    # If all methods fail
    print(f"No covariance data found for {survey_file_name}")
    return [], []

# In the process_k_fg function, update the covariance combination
def process_k_fg(k_fg):
    """Process results for a specific k_fg value with fixed bin width"""
    # Create output directory for this k_fg value
    k_fg_output_dir = os.path.join(output_dir, f"k_fg_{k_fg:.4f}")
    os.makedirs(k_fg_output_dir, exist_ok=True)

    # Path to the bin directory for this k_fg value (using fixed bin width)
    bin_dir = os.path.join(input_dir, f"k_fg_{k_fg:.4f}", f"dz{bin_width:.2f}")

    if not os.path.exists(bin_dir):
        print(f"Directory not found: {bin_dir}")
        return None

    # Dictionary to store combined results for each survey
    survey_results = {}

    # Process each survey
    for survey_name in SURVEY_COLORS.keys():
        # Clean survey name for filenames
        survey_file_name = survey_name.replace(' ', '_')

        print(f"\nProcessing {survey_name} in k_fg={k_fg:.4f}, bin width={bin_width:.2f}")
        
        # Try to read all covariance matrices using various formats
        cov_matrices, bin_volumes = read_covariance_files(bin_dir, survey_file_name)
        
        if not cov_matrices:
            print(f"No covariance matrices found for {survey_name}")
            continue
        
        print(f"Found {len(cov_matrices)} covariance matrices for {survey_name}")
        
        # Print a sample of the first covariance matrix for verification
        if len(cov_matrices) > 0:
            print("First covariance matrix:")
            print(np.array(cov_matrices[0]))
            print(f"Determinant: {np.linalg.det(np.array(cov_matrices[0]))}")
            
            # Extract and print uncertainties from the first matrix
            cov0 = np.array(cov_matrices[0])
            unc0 = np.sqrt(np.diag(cov0))
            print(f"Uncertainties from first matrix: {unc0}")
            print(f"Volume of first bin: {bin_volumes[0] if bin_volumes else 'N/A'}")

        # Combine the covariance matrices with volume weighting
        combined_cov = combine_fisher_matrices(cov_matrices, bin_volumes=bin_volumes)

        if combined_cov is None:
            print(f"Failed to combine matrices for {survey_name}")
            continue

        # Get uncertainties
        uncertainties = np.sqrt(np.diag(combined_cov))

        # Calculate detection significance (Cunnington method: alpha/sigma(alpha))
        # Using fiducial alpha value of 1.0 from the "fiducial" array below
        detection_significance = 1.0 / uncertainties[0]  # alpha=1.0 / sigma(alpha)

        # Store results
        survey_results[survey_name] = {
            "covariance": combined_cov,
            "uncertainties": uncertainties,
            "parameters": ["alpha", "beta", "k0"],
            "fiducial": [1.0, 0.8, 0.0164],
            "detection_significance": detection_significance  # Add this field
        }

        print(f"  {survey_name}: Combined {len(cov_matrices)} matrices")
        print(f"    σ(α) = {uncertainties[0]:.6f}")
        print(f"    σ(β) = {uncertainties[1]:.6f}")
        print(f"    σ(k₀) = {uncertainties[2]:.6f}")
        print(f"    Detection significance: {detection_significance:.2f}σ")

        # Save the combined covariance matrix for later analysis
        output_cov_file = os.path.join(k_fg_output_dir, f"{survey_file_name}_combined_cov.npy")
        np.save(output_cov_file, combined_cov)

    # Create triangle plots using getdist
    create_triangle_plots(survey_results, k_fg_output_dir, k_fg)

    return survey_results

def create_triangle_plots(survey_results, output_dir, k_fg):
    """Create triangle plots for all surveys at a given k_fg value using GetDist"""
    if not survey_results:
        print("No survey results to plot!")
        return
        
    # Define parameter labels
    params = ["alpha", "beta", "k0"]
    latex_labels = [r"\alpha", r"\beta", r"k_0"]

    # Generate GetDist sample objects for each survey
    samples_dict = {}

    for survey_name, results in survey_results.items():
        # Get covariance matrix and fiducial values
        cov = results["covariance"]
        fiducial = results["fiducial"]

        # Generate mock samples from the covariance matrix (for better density estimates)
        # Use more samples to better capture the degeneracy
        np.random.seed(42)  # For reproducibility
        sample_points = np.random.multivariate_normal(fiducial, cov, size=500000)
        
        # Create a direct scatter plot to verify correlations
        fig, ax = plt.subplots(figsize=(8, 8))
        plt.scatter(sample_points[:5000, 0], sample_points[:5000, 1], alpha=0.3, s=3)
        plt.xlabel(r'$\alpha$', fontsize=14)
        plt.ylabel(r'$\beta$', fontsize=14)
        corr = correlation_coefficient(cov, 0, 1)
        plt.title(f"Correlation: {corr:.4f}")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{survey_name.replace(' ', '_')}_direct_correlation.pdf"))
        plt.close()

        # Create MCSamples object
        samples = MCSamples(
            samples=sample_points,
            names=params,
            labels=[r'$%s$' % l for l in latex_labels],
            label=survey_name
        )

        samples_dict[survey_name] = samples

    # Create individual plots for each survey
    for name, samples in samples_dict.items():
        # Create the triangle plot
        g = plots.get_subplot_plotter(width_inch=10)
        
        # Use only the most basic settings that should be available in all GetDist versions
        g.settings.figure_legend_frame = True
        g.settings.axis_marker_lw = 1.5
        
        # Use only num_plot_contours which is a core setting
        g.settings.num_plot_contours = 3
        
        # Reduce contour line thickness
        g.settings.linewidth_contour = 0.8  # Added this line to make contours thinner
        
        # Create the triangle plot with minimal settings
        g.triangle_plot(
            [samples], 
            params, 
            filled=True,
            contour_colors=[SURVEY_COLORS[name]],
            markers={'alpha': 1.2, 'beta': 0.8, 'k0': 0.0164}
        )

        # After plot is created, manually adjust each axis
        for ax in plt.gcf().get_axes():
            ax.grid(False)
            
            # Manually add axis labels if needed (since we're not using no_title_or_label)
            if ax.get_xlabel() == '' and ax.get_ylabel() != '':
                # This is an x-axis without a label but with a y-label
                y_label = ax.get_ylabel().strip('$').replace('\\', '')
                # Special handling for k_0 vs k0
                if y_label == 'k_0':
                    param_idx = params.index('k0')
                else:
                    param_idx = params.index(y_label)
                ax.set_xlabel(f'${latex_labels[param_idx]}$')
            
            if ax.get_ylabel() == '' and ax.get_xlabel() != '':
                # This is a y-axis without a label but with an x-label
                x_label = ax.get_xlabel().strip('$').replace('\\', '')
                # Special handling for k_0 vs k0
                if x_label == 'k_0':
                    param_idx = params.index('k0')
                else:
                    param_idx = params.index(x_label)
                ax.set_ylabel(f'${latex_labels[param_idx]}$')

        # Save the plot with high resolution
        survey_file_name = name.replace(' ', '_')
        plt.savefig(os.path.join(output_dir, f"{survey_file_name}_k_fg_{k_fg:.4f}.pdf"), 
                    bbox_inches='tight', dpi=600)
        plt.close()

    # Create a combined plot with all surveys - using minimal settings
    g = plots.get_subplot_plotter(width_inch=12)
    g.settings.figure_legend_frame = True
    g.settings.axis_marker_lw = 1.5
    g.settings.legend_fontsize = 24
    g.settings.num_plot_contours = 3

    # Ensure SKA-MID is last in the list so it's on top
    ordered_names = []
    ordered_samples = []

    # First add non-SKA surveys
    for name in samples_dict.keys():
        if "SKA" not in name:
            ordered_names.append(name)
            ordered_samples.append(samples_dict[name])

    # Then add SKA (so it's on top)
    for name in samples_dict.keys():
        if "SKA" in name:
            ordered_names.append(name)
            ordered_samples.append(samples_dict[name])

    # Customize colors for each survey in the correct order
    colors = [SURVEY_COLORS[name] for name in ordered_names]

    # Create the plot with all surveys - with minimal settings
    g.triangle_plot(
        ordered_samples, 
        params, 
        filled=True,
        contour_colors=colors,
        legend_labels=ordered_names,
        legend_loc='upper right'
    )

    # After plot is created, turn off all grid lines with matplotlib directly
    for ax in plt.gcf().get_axes():
        ax.grid(False)

    # Save the plot
    plt.savefig(os.path.join(output_dir, f"triangle_plot_k_fg_{k_fg:.4f}.pdf"), 
                bbox_inches='tight', dpi=600)
    plt.close()

    print(f"Triangle plots created for k_fg = {k_fg:.4f}")

def create_comparison_plots(all_results):
    """Create comparison plots showing how results change with k_fg"""
    comparison_dir = os.path.join(output_dir, "Comparisons")
    os.makedirs(comparison_dir, exist_ok=True)

    # Extract parameter uncertainties for each survey and k_fg
    data = []
    for k_fg, results in all_results.items():
        for survey_name, survey_results in results.items():
            uncertainties = survey_results["uncertainties"]
            # Calculate Fisher Figure of Merit (FoM) as inverse of uncertainty volume
            fom = 1.0 / (uncertainties[0] * uncertainties[1] * uncertainties[2])

            
            # Calculate detection significance
            detection_significance = 1.0 / uncertainties[0]  # alpha=1.0 / sigma(alpha)

            # Also extract parameter correlations
            cov = survey_results["covariance"]
            alpha_beta_corr = correlation_coefficient(cov, 0, 1)
            alpha_k0_corr = correlation_coefficient(cov, 0, 2)
            beta_k0_corr = correlation_coefficient(cov, 1, 2)

            data.append({
                "k_fg": k_fg,
                "Survey": survey_name,
                "sigma(alpha)": uncertainties[0],
                "sigma(beta)": uncertainties[1],
                "sigma(k0)": uncertainties[2],
                "corr(alpha,beta)": alpha_beta_corr,
                "corr(alpha,k0)": alpha_k0_corr,
                "corr(beta,k0)": beta_k0_corr,
                "FoM": fom,
                "Detection_Significance": detection_significance
            })

    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(data)

    # Save the complete results with correlations
    df.to_csv(os.path.join(comparison_dir, "all_results_with_correlations.csv"), index=False)

    # Plot 1: Parameter uncertainties vs k_fg for each survey
    param_names = [r"\alpha", r"\beta", r"k_0"]
    param_columns = ["sigma(alpha)", "sigma(beta)", "sigma(k0)"]

    plt.figure(figsize=(15, 10))

    for i, (param, col) in enumerate(zip(param_names, param_columns)):
        plt.subplot(2, 2, i+1)

        for survey in df["Survey"].unique():
            survey_data = df[df["Survey"] == survey]
            # Sort by k_fg to ensure lines connect properly
            survey_data = survey_data.sort_values(by="k_fg")

            # Scatter plot
            plt.scatter(
                survey_data["k_fg"], 
                survey_data[col],
                color=SURVEY_COLORS[survey],
                label=survey, 
                s=120,
                edgecolors='black',
                linewidths=1.0,
                alpha=0.9
            )

            # Connect points with lines
            plt.plot(
                survey_data["k_fg"], 
                survey_data[col],
                color=SURVEY_COLORS[survey],
                alpha=0.5,
                linestyle='--'
            )

        plt.xlabel(r"$k_{fg}$ [h/Mpc]", fontsize=16)
        plt.ylabel(r"$\sigma(" + param + r")$", fontsize=16)
        plt.grid(True, linestyle='--', alpha=0.7)
        if i == 0:  # Only add legend to the first subplot
            plt.legend(loc='best', frameon=True, framealpha=0.9, edgecolor='grey')

    # Plot 4: Parameter correlations vs k_fg
    plt.subplot(2, 2, 4)
    for survey in df["Survey"].unique():
        survey_data = df[df["Survey"] == survey]
        # Sort by k_fg to ensure lines connect properly
        survey_data = survey_data.sort_values(by="k_fg")
        
        plt.scatter(
            survey_data["k_fg"], 
            survey_data["corr(alpha,beta)"],
            color=SURVEY_COLORS[survey],
            label=survey, 
            s=120,
            edgecolors='black',
            linewidths=1.0,
            alpha=0.9
        )
        plt.plot(
            survey_data["k_fg"], 
            survey_data["corr(alpha,beta)"],
            color=SURVEY_COLORS[survey],
            alpha=0.5,
            linestyle='--'
        )
        plt.xlabel(r"$k_{fg}$ [h/Mpc]", fontsize=16)
    plt.ylabel(r"Correlation($\alpha$, $\beta$)", fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.7)

    # Use a simpler tight_layout approach to prevent errors
    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, "parameter_uncertainties_vs_kfg.pdf"), 
                dpi=600, bbox_inches='tight')
    plt.close()

    # New plot to visualize correlations
    plt.figure(figsize=(15, 5))

    # Plot correlations for alpha-beta, alpha-k0, beta-k0
    corr_pairs = [
        ("corr(alpha,beta)", r"Correlation($\alpha$, $\beta$)"),
        ("corr(alpha,k0)", r"Correlation($\alpha$, $k_0$)"),
        ("corr(beta,k0)", r"Correlation($\beta$, $k_0$)")
    ]

    for i, (col, label) in enumerate(corr_pairs):
        plt.subplot(1, 3, i+1)

        for survey in df["Survey"].unique():
            survey_data = df[df["Survey"] == survey]
            # Sort by k_fg
            survey_data = survey_data.sort_values(by="k_fg")

            plt.scatter(
                survey_data["k_fg"], 
                survey_data[col],
                color=SURVEY_COLORS[survey],
                label=survey, 
                s=100,
                edgecolors='black',
                linewidths=1.0,
                alpha=0.9
            )

            plt.plot(
                survey_data["k_fg"], 
                survey_data[col],
                color=SURVEY_COLORS[survey],
                alpha=0.5,
                linestyle='--'
            )

        plt.xlabel(r"$k_{fg}$ [h/Mpc]", fontsize=16)
        plt.ylabel(label, fontsize=16)
        plt.grid(True, linestyle='--', alpha=0.7)
        if i == 0:
            plt.legend(loc='best', frameon=True, framealpha=0.9, edgecolor='grey')
        plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, "parameter_correlations.pdf"), 
                dpi=600, bbox_inches='tight')
    plt.close()

    # Create a summary plot of best k_fg for each parameter
    best_k_fg = {}

    for survey in df["Survey"].unique():
        survey_data = df[df["Survey"] == survey]

        best_k_fg[survey] = {
            "alpha": survey_data.loc[survey_data["sigma(alpha)"].idxmin()]["k_fg"],
            "beta": survey_data.loc[survey_data["sigma(beta)"].idxmin()]["k_fg"],
            "k0": survey_data.loc[survey_data["sigma(k0)"].idxmin()]["k_fg"],
            "FoM": survey_data.loc[survey_data["FoM"].idxmax()]["k_fg"]
        }

    # Create a visualization of the optimal k_fg value for each parameter and survey
    fig, ax = plt.subplots(figsize=(10, 6))

    surveys = list(best_k_fg.keys())
    params = ["alpha", "beta", "k0", "FoM"]
    param_labels = [r"$\alpha$", r"$\beta$", r"$k_0$", "FoM"]

    x = np.arange(len(surveys))
    width = 0.2

    for i, param in enumerate(params):
        values = [best_k_fg[survey][param] for survey in surveys]
        ax.bar(x + (i-1.5)*width, values, width, label=param_labels[i], alpha=0.8)

    ax.set_ylabel('Optimal $k_{fg}$ [h/Mpc]', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(surveys, fontsize=14)
    ax.legend(loc='best', frameon=True, framealpha=0.9, edgecolor='grey')
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.savefig(os.path.join(comparison_dir, "optimal_kfg_by_parameter.pdf"), 
                dpi=600, bbox_inches='tight')
    plt.close()

    # Generate a summary table
    summary_data = []
    for survey in df["Survey"].unique():
        survey_data = df[df["Survey"] == survey]
        best_row = survey_data.loc[survey_data["FoM"].idxmax()]

        summary_data.append({
            "Survey": survey,
            "Best k_fg": best_row["k_fg"],
            "sigma(alpha)": best_row["sigma(alpha)"],
            "sigma(beta)": best_row["sigma(beta)"],
            "sigma(k0)": best_row["sigma(k0)"],
            "corr(alpha,beta)": best_row["corr(alpha,beta)"],
            "FoM": best_row["FoM"]
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(comparison_dir, "best_constraints_summary.csv"), index=False)

    # Also create a pretty LaTeX table
    latex_table = summary_df.to_latex(
        index=False,
        formatters={
            "Best k_fg": lambda x: f"{x:.4f}",
            "sigma(alpha)": lambda x: f"{x:.6f}",
            "sigma(beta)": lambda x: f"{x:.6f}",
            "sigma(k0)": lambda x: f"{x:.6f}",
            "corr(alpha,beta)": lambda x: f"{x:.4f}",
            "FoM": lambda x: f"{x:.2f}"
        },
        caption="Summary of Best Constraints for Each Survey",
        label="tab:best_constraints"
    )

    with open(os.path.join(comparison_dir, "best_constraints_summary.tex"), 'w') as f:
        f.write(latex_table)

    print("Comparison plots and summary tables created")
    return df

def plot_k0_vs_redshift():
    """
    Plot σ(k_0) vs redshift for each survey in each k_fg value folder.
    Uses the correct redshift range for each survey.
    """
    # Create output directory for these plots
    k0_plots_dir = os.path.join(output_dir, "K0_Analysis")
    os.makedirs(k0_plots_dir, exist_ok=True)

    # For each k_fg value
    for k_fg in k_fg_values:
        print(f"Processing k_fg = {k_fg:.4f} for redshift analysis")

        # Path to the bin directory for this k_fg value (using fixed bin width)
        bin_dir = os.path.join(input_dir, f"k_fg_{k_fg:.4f}", f"dz{bin_width:.2f}")

        if not os.path.exists(bin_dir):
            print(f"Directory not found: {bin_dir}")
            continue

        # Create plot
        plt.figure(figsize=(12, 10))

        # Process each survey
        for survey_name in SURVEY_COLORS.keys():
            # Clean survey name for filenames
            survey_file_name = survey_name.replace(' ', '_')

            # Get the correct redshift range for this survey
            z_min, z_max = SURVEY_REDSHIFT_RANGES[survey_name]

            # Try to load the results JSON file
            json_file = os.path.join(bin_dir, f"{survey_file_name}_results.json")
            if not os.path.exists(json_file):
                print(f"Results file not found: {json_file}")
                # Try to find individual text files
                continue

            with open(json_file, 'r') as f:
                results = json.load(f)

            # Extract redshifts and uncertainties
            redshifts = np.array(results["redshifts"])
            
            # Check if uncertainties are directly available
            if "uncertainties" in results:
                uncertainties = results["uncertainties"]
                # Extract k_0 uncertainty for each redshift bin (k_0 is the 3rd parameter, index 2)
                k0_errors = [unc[2] for unc in uncertainties]
            else:
                # Need to compute uncertainties from covariance matrices
                cov_matrices = results["covariance_matrices"]
                k0_errors = []
                for cov in cov_matrices:
                    cov_array = np.array(cov)
                    # Check if the matrix is valid
                    if np.any(np.isnan(cov_array)) or np.linalg.det(cov_array) < 1e-15:
                        k0_errors.append(np.nan)  # Use NaN for invalid matrices
                    else:
                        k0_error = np.sqrt(cov_array[2, 2])
                        k0_errors.append(k0_error)

            # Filter redshifts to only include the valid range for this survey
            valid_range_mask = (redshifts >= z_min) & (redshifts <= z_max)
            filtered_redshifts = redshifts[valid_range_mask]
            filtered_k0_errors = np.array(k0_errors)[valid_range_mask]

            # Filter out NaN values
            valid_mask = ~np.isnan(filtered_k0_errors)
            valid_redshifts = filtered_redshifts[valid_mask]
            valid_k0_errors = filtered_k0_errors[valid_mask]

            if len(valid_redshifts) == 0:
                print(f"No valid data points for {survey_name} in redshift range {z_min}-{z_max}")
                continue

            # Plot this survey with dots only (no connecting lines)
            plt.scatter(
                valid_redshifts,
                valid_k0_errors,
                color=SURVEY_COLORS[survey_name],
                label=f"{survey_name}",
                s=80,
                edgecolors='grey',
                linewidths=0.8,
                alpha=0.9
            )

        # Add labels and legend
        plt.xlabel("Redshift", fontsize=28)
        plt.ylabel(r"$\sigma(k_0)$ [h/Mpc]", fontsize=28)
        plt.legend(loc='best', frameon=True, framealpha=0.9, edgecolor='grey')
        plt.grid(True, linestyle='--', alpha=0.7)

        # Use scientific notation for y-axis
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

        # Save the plot
        plt.tight_layout()
        plt.savefig(os.path.join(k0_plots_dir, f"k0_error_vs_redshift_k_fg_{k_fg:.4f}.pdf"), 
                    dpi=600, bbox_inches='tight')
        plt.close()

    print("K0 vs redshift plots created")

def plot_k0_vs_kfg():
    """
    Plot σ(k_0) vs k_fg for each survey.
    """
    # Create output directory for these plots
    k0_plots_dir = os.path.join(output_dir, "K0_Analysis")
    os.makedirs(k0_plots_dir, exist_ok=True)

    # Dictionary to store results
    results_by_survey = {survey: [] for survey in SURVEY_COLORS.keys()}

    # For each k_fg value
    for k_fg in k_fg_values:
        # Look for the combined covariance matrix saved by the process_k_fg function
        for survey_name in SURVEY_COLORS.keys():
            # Clean survey name for filenames
            survey_file_name = survey_name.replace(' ', '_')

            # Look for the combined covariance matrix saved by the main script
            cov_file = os.path.join(output_dir, f"k_fg_{k_fg:.4f}", f"{survey_file_name}_combined_cov.npy")

            if os.path.exists(cov_file):
                # Load the combined covariance matrix
                combined_cov = np.load(cov_file)

                # Extract k_0 uncertainty (k_0 is the 3rd parameter, index 2)
                k0_error = np.sqrt(combined_cov[2, 2])

                # Store result
                results_by_survey[survey_name].append({
                    "k_fg": k_fg,
                    "k0_error": k0_error
                })

    # Create plot
    plt.figure(figsize=(12, 10))

    # Plot for each survey
    for survey_name, results in results_by_survey.items():
        if not results:
            continue

        # Sort by k_fg to ensure proper ordering
        results.sort(key=lambda x: x["k_fg"])

        # Extract k_fg values and k0 errors
        k_fg_list = [r["k_fg"] for r in results]
        k0_error_list = [r["k0_error"] for r in results]

        # Plot this survey with dots and connecting lines
        plt.scatter(
            k_fg_list,
            k0_error_list,
            color=SURVEY_COLORS[survey_name],
            label=f"{survey_name}",
            s=80,
            edgecolors='grey',
            linewidths=0.8,
            alpha=0.9
        )
        
        # Connect points with a line
        plt.plot(
            k_fg_list,
            k0_error_list,
            color=SURVEY_COLORS[survey_name],
            linestyle='--',
            alpha=0.6
        )

    # Add labels and legend
    plt.xlabel(r"$k_{fg}$ [h/Mpc]", fontsize=28)
    plt.ylabel(r"$\sigma(k_0)$ [h/Mpc]", fontsize=28)
    plt.legend(loc='best', frameon=True, framealpha=0.9, edgecolor='grey')
    plt.grid(True, linestyle='--', alpha=0.7)

    # Use scientific notation for y-axis
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    # Save the plot
    plt.tight_layout()
    plt.savefig(os.path.join(k0_plots_dir, f"k0_error_vs_kfg.pdf"), 
                dpi=600, bbox_inches='tight')
    plt.close()

    print("K0 vs k_fg plot created")

def plot_combined_constraints():
    """
    Create a combined plot showing the parameter constraints for all surveys
    """
    # Create output directory
    k0_plots_dir = os.path.join(output_dir, "Combined_Analysis")
    os.makedirs(k0_plots_dir, exist_ok=True)

    # Track best results
    best_results = {}

    # For each survey
    for survey_name in SURVEY_COLORS.keys():
        survey_file_name = survey_name.replace(' ', '_')
        best_k0_error = float('inf')
        best_k_fg = None
        best_cov = None

        # Find the best k_fg for this survey (the one with lowest k0 error)
        for k_fg in k_fg_values:
            cov_file = os.path.join(output_dir, f"k_fg_{k_fg:.4f}", f"{survey_file_name}_combined_cov.npy")

            if os.path.exists(cov_file):
                combined_cov = np.load(cov_file)
                k0_error = np.sqrt(combined_cov[2, 2])

                if k0_error < best_k0_error:
                    best_k0_error = k0_error
                    best_k_fg = k_fg
                    best_cov = combined_cov

        if best_cov is not None:
            best_results[survey_name] = {
                "k_fg": best_k_fg,
                "cov": best_cov,
                "uncertainties": np.sqrt(np.diag(best_cov)),
                "fiducial": [1.2, 0.8, 0.0164]
            }

    # Create two-panel plot: constraints vs survey and summary table
    plt.figure(figsize=(16, 8))

    # Panel 1: Bar chart of constraints
    plt.subplot(1, 2, 1)

    surveys = list(best_results.keys())
    x = np.arange(len(surveys))
    width = 0.3

    # Extract uncertainties
    alpha_errors = [best_results[s]["uncertainties"][0] for s in surveys]
    beta_errors = [best_results[s]["uncertainties"][1] for s in surveys]
    k0_errors = [best_results[s]["uncertainties"][2] for s in surveys]

    # Scale to make visible on same plot
    scale_factor = 10

    # Plot bars
    plt.bar(x - width, alpha_errors, width, label=r'$\sigma(\alpha)$', alpha=0.7, color='blue')
    plt.bar(x, beta_errors, width, label=r'$\sigma(\beta)$', alpha=0.7, color='green') 
    plt.bar(x + width, [err * scale_factor for err in k0_errors], width, 
           label=r'$\sigma(k_0) \times ' + str(scale_factor) + '$', alpha=0.7, color='red')

    plt.xlabel('Survey', fontsize=16)
    plt.ylabel('Parameter Uncertainty', fontsize=16)
    plt.xticks(x, surveys, fontsize=12)
    plt.legend(fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Panel 2: Text summary - Use plain text column labels without Greek characters
    plt.subplot(1, 2, 2)
    plt.axis('off')

    table_data = []
    for i, survey in enumerate(surveys):
        result = best_results[survey]
        table_data.append([
            survey, 
            f"{result['k_fg']:.4f}",
            f"{result['uncertainties'][0]:.6f}",
            f"{result['uncertainties'][1]:.6f}",
            f"{result['uncertainties'][2]:.6f}"
        ])

    # Use plain ASCII text for column labels (no Greek letters)
    column_labels = ["Survey", "Best k_fg", "sigma(alpha)", "sigma(beta)", "sigma(k0)"]

    table = plt.table(
        cellText=table_data,
        colLabels=column_labels,
        loc='center',
        cellLoc='center',
        colColours=['lightgrey']*5
    )

    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.5)

    plt.tight_layout()
    plt.savefig(os.path.join(k0_plots_dir, "best_survey_constraints.pdf"), 
                dpi=600, bbox_inches='tight')
    plt.close()

    # Also create a standalone table in LaTeX format for direct use in papers
    table_file = os.path.join(k0_plots_dir, "best_survey_constraints.tex")
    with open(table_file, 'w') as f:
        f.write("\\begin{table}\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{lccccc}\n")  # Added one more column
        f.write("\\hline\n")
        f.write("Survey & Best $k_{fg}$ & $\\sigma(\\alpha)$ & $\\sigma(\\beta)$ & $\\sigma(k_0)$ & Detection Significance \\\\\n")
        f.write("\\hline\n")

        for i, survey in enumerate(surveys):
            result = best_results[survey]
            detection = 1.0 / result['uncertainties'][0]  # alpha/sigma(alpha)
            f.write(f"{survey} & {result['k_fg']:.4f} & {result['uncertainties'][0]:.6f} & ")
            f.write(f"{result['uncertainties'][1]:.6f} & {result['uncertainties'][2]:.6f} & {detection:.1f}$\\sigma$ \\\\\n")

        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\caption{Summary of best constraints and detection significance for each survey}\n")
        f.write("\\label{tab:best_constraints}\n")
        f.write("\\end{table}\n")

    print("Combined constraints plot and LaTeX table created")

def create_animation_script():
    """Create a script to combine triangle PDFs into an animation"""
    script_path = os.path.join(output_dir, "create_animation.sh")

    with open(script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# This script creates an animation of triangle plots as k_fg changes\n\n")

        # First ensure the Animation directory exists
        f.write("mkdir -p Animation\n\n")

        # Copy all the combined triangle plots to numbered frames for the animation
        f.write("# Copy triangle plots to numbered frames\n")
        for i, k_fg in enumerate(k_fg_values):
            f.write(f"cp k_fg_{k_fg:.4f}/triangle_plot_k_fg_{k_fg:.4f}.pdf Animation/frame_{i:02d}.pdf\n")

        f.write("\n# Convert PDFs to PNG (requires ImageMagick)\n")
        f.write("cd Animation\n")
        f.write("for pdf in frame_*.pdf; do\n")
        f.write("    png=${pdf%.pdf}.png\n")
        f.write("    convert -density 150 $pdf -quality 100 $png\n")
        f.write("done\n\n")

        f.write("# Create animated GIF\n")
        f.write("convert -delay 100 -loop 0 frame_*.png triangle_animation.gif\n\n")

        f.write("echo 'Animation created at Animation/triangle_animation.gif'\n")

    # Make the script executable
    os.chmod(script_path, 0o755)
    print(f"Animation script created at {script_path}")

def main():
    """Main function to process all data and create plots"""
    print(f"Processing Fisher results with fixed bin width of {bin_width}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")

    # Process each k_fg value
    all_results = {}
    processed_k_fg_values = []  # Keep track of which k_fg values actually work

    for k_fg in k_fg_values:
        print(f"\nProcessing k_fg = {k_fg:.4f}")
        bin_dir = os.path.join(input_dir, f"k_fg_{k_fg:.4f}", f"dz{bin_width:.2f}")

        if not os.path.exists(bin_dir):
            print(f"⚠️ Directory not found: {bin_dir}")
            print(f"⚠️ Did you run FisherVaryingBinVolume.py for k_fg = {k_fg:.4f}?")
            continue

        results = process_k_fg(k_fg)
        if results:
            all_results[k_fg] = results
            processed_k_fg_values.append(k_fg)

    # Create comparison plots with the k_fg values that actually worked
    if all_results:
        print(f"\nSuccessfully processed {len(processed_k_fg_values)} k_fg values: {processed_k_fg_values}")
        df = create_comparison_plots(all_results)

        # Create k0-focused analysis plots with available data
        plot_k0_vs_redshift()
        plot_k0_vs_kfg()
        plot_combined_constraints()

        create_animation_script()
        print("\nAll processing complete!")
        print(f"Results saved to {output_dir}/")
    else:
        print("\n❌ No valid results to process.")
        print("Please ensure you've run FisherVaryingBinVolume.py for at least one of these k_fg values:")
        for k_fg in k_fg_values:
            print(f"  - {k_fg:.4f}")

if __name__ == "__main__":
    main()
