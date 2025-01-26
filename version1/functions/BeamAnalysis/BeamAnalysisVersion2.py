import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# File paths for Beam and NoBeam data (adjusted with the correct file path)
beam_files = [
    '/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z0.4.dat',
    '/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z0.5.dat',
    '/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z0.7.dat',
    '/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z0.53.dat',
    '/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z1.0.dat',
    '/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_Beam_z2.5.dat'
]

nobeam_file = '/home/yolanda/Documents/mcmc_code/version1/result/Effect_Of_Beam_Analysis/MCMC_skaOIMBand1skaOIMBand1_0.001_NoBeam_z0.5.dat'

# Load data
beam_data = [np.loadtxt(f) for f in beam_files]
nobeam_data = np.loadtxt(nobeam_file)

# Assuming k_0 is in the second column (index 1), calculate the mean
beam_k0_means = [np.mean(data[:, 1]) for data in beam_data]
nobeam_k0_mean = np.mean(nobeam_data[:, 1])

# Calculate the bias (difference between Beam and NoBeam values)
biases = [beam_k0 - nobeam_k0_mean for beam_k0 in beam_k0_means]

# Redshift labels corresponding to the files
redshifts = ['z0.4', 'z0.5', 'z0.7', 'z0.53', 'z1.0', 'z2.5']

# Create a DataFrame for the results
bias_df = pd.DataFrame({
    'Redshift': redshifts,
    'Beam k_0': beam_k0_means,
    'NoBeam k_0': nobeam_k0_mean,
    'Bias (Beam - NoBeam)': biases
})

# Display the results in a table
print(bias_df)

# Optionally, plot the bias as a function of redshift
plt.figure(figsize=(8, 6))
plt.plot(redshifts, biases, marker='o', linestyle='-', color='b')
plt.title('Beam Bias in Detected Turnover (k_0)')
plt.xlabel('Redshift')
plt.ylabel('Bias (Beam - NoBeam)')
plt.grid(True)
plt.show()
