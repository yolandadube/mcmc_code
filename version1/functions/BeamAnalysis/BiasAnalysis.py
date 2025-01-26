import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Enable LaTeX style for matplotlib
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Define the file paths for the different cases and surveys
file_paths = {
    "SKA Mid Band 1": {
        "Ignoring Effect of Beam and Foregrounds": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_skaOIMBand1skaOIMBand1_NoBeam_NoFg_0.001.dat",
        "Ignoring Foregrounds Only": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_skaOIMBand1skaOIMBand1_Beam_NoFg_0.001.dat",
        "Ignoring Effect of Beam Only": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_skaOIMBand1skaOIMBand1_NoBeam_Fg_0.001.dat"
    },
    "MeerKAT L Band": {
        "Ignoring Effect of Beam and Foregrounds": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_meerKATLBandmeerKATLBand_NoBeam_NoFg_0.001.dat",
        "Ignoring Foregrounds Only": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_meerKATLBandmeerKATLBand_Beam_NoFg_0.001.dat",
        "Ignoring Effect of Beam Only": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_meerKATLBandmeerKATLBand_NoBeam_Fg_0.001.dat"
    },
    "MeerKAT UHF Band": {
        "Ignoring Effect of Beam and Foregrounds": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_meerKATUHFBandmeerKATUHFBand_NoBeam_NoFg_0.001.dat",
        "Ignoring Foregrounds Only": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_meerKATUHFBandmeerKATUHFBand_Beam_NoFg_0.001.dat",
        "Ignoring Effect of Beam Only": "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias/MCMC_ParabolicModel_meerKATUHFBandmeerKATUHFBand_NoBeam_Fg_0.001.dat"
    }
}

# Create the output directory if it doesn't exist
output_dir = "/home/yolanda/Documents/mcmc_code/version1/result/Cause_Of_Bias"
output_file_path = os.path.join(output_dir, "detected_k0_values_plot.pdf")
os.makedirs(output_dir, exist_ok=True)

# Load the data for each case and survey, using the second column (index 1) for k_0
data_frames_corrected = []
for survey, cases in file_paths.items():
    for case, file_path in cases.items():
        df = pd.read_csv(file_path, delim_whitespace=True, header=None, usecols=[1], names=["k_0"])  # Second column (index 1)
        df["Survey"] = survey
        df["Case"] = case
        data_frames_corrected.append(df)

# Concatenate all dataframes
combined_df_corrected = pd.concat(data_frames_corrected, ignore_index=True)

# Set the fiducial value
fiducial_value = 0.0167

# Plot the corrected boxplot with the fiducial value and LaTeX math formatting
plt.figure(figsize=(15, 12))

# Create the boxplot with specific red, green, and blue colors
sns.boxplot(x="Survey", y="k_0", hue="Case", data=combined_df_corrected, 
            palette={"Ignoring Effect of Beam and Foregrounds": "red", "Ignoring Foregrounds Only": "green", "Ignoring Effect of Beam Only": "blue"})

# Add the fiducial value line as grey
plt.axhline(y=fiducial_value, color='grey', linestyle='--', label=r"Fiducial Value: $k_0 = 0.0168$")

# Adjusting font sizes
plt.xlabel(r"\textbf{Surveys}", fontsize=28)
plt.ylabel(r"$k_0$", fontsize=32)
plt.xticks(rotation=0, fontsize=28)
plt.yticks(fontsize=30)

# Add legend for the fiducial value
plt.legend(fontsize=20, title_fontsize=20)

# Save the figure as a high-quality PDF to the specified path
plt.tight_layout()
plt.savefig(output_file_path, format='pdf', dpi=300)

# Show the plot
#plt.show()
