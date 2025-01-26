import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# main path directory
mainPath = "/home/yolanda/Documents/mcmc_code/version1/"

# Import necessary functions from other parts of your project
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/nModes/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/surveySpecifications/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/surveyRedshifts/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/surveyVolume/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/Theory/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/noise/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/backgroundCosmology/')

from theory import model
from numberOfModes import Nmodes
from specifications import surveySpecs
from thermalNoise import HIThermalNoise
from surveyRedshiftRange import redshiftRange
from comovingDistance import comoving_distance
from shotNoise import shot_noise  # Import the shot noise function

def fisherForecasting(mainPath, dictionary, cosmologY, z, delta_z, kMin, kMax, surveys, specsSurveys, modeSurveys):

    k = np.linspace(kMin, kMax, 1000, endpoint=True)

    parameters = {'alpha': 0.45, 'beta': 0.30, 'P0': 37300, 'k0': 0.0163}  # Fiducial parameters

    # Directory to save the inverse Fisher matrices (covariance matrices)
    save_dir = os.path.join(mainPath, 'result', 'contour', 'FisherforecastingDATA')

    # Ensure the directory exists
    os.makedirs(save_dir, exist_ok=True)

    def compute_fisher_matrix():
        
        for idx, specsSurvey in enumerate(specsSurveys):

            print(f"Processing survey: {specsSurvey}")

            # Retrieve survey specifications
            N_dish, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey)
            print(f"Survey specifications retrieved for: {specsSurvey}")

            # Calculate redshift range
            z_min, z_max = redshiftRange("galaxy", specsSurvey)
            print(f"Redshift range for {specsSurvey}: z_min={z_min}, z_max={z_max}")

            # Calculate the number of modes
            N_modes = Nmodes(mainPath, k, z, delta_z, cosmologY, surveys, specsSurveys, modeSurveys)
            print(f"Number of modes calculated for {specsSurvey}")

            # Initialize the Fisher matrix assuming 4 parameters (alpha, beta, P0, k0)
            F_matrix = np.zeros((4, 4))

            for i, ki in enumerate(k):

                x = (np.log(ki) / np.log(parameters['k0'])) - 1

                for j, zj in enumerate(z):
                    # Compute derivative terms ∂P/∂θ
                    dP_dm = -x**2 * parameters['P0']**(1 - parameters['alpha'] * x**2) * np.log(parameters['P0'])
                    dP_dn = -x**2 * parameters['P0']**(1 - parameters['beta'] * x**2) * np.log(parameters['P0'])
                    dP_dP0 = (1 - parameters['alpha'] * x**2) * parameters['P0']**(-parameters['alpha'] * x**2)
                    dP_dkto = parameters['alpha'] * x * parameters['P0']**(1 - parameters['beta'] * x**2) * (np.log(ki) + 1) /\
                    (np.log(parameters['k0'])**2) * np.log(parameters['P0'])

                    # Calculate the power spectrum and noise
                    P = model(ki, parameters['P0'], parameters['k0'], parameters['alpha'], parameters['beta'])
                    
                    # Use HIThermalNoise for "HI IM" and shot_noise for "galaxy"
                    if surveys[idx] == "HI IM":
                        PN = HIThermalNoise(mainPath, 0, ki, zj, dictionary, cosmologY, surveys[idx], specsSurvey, modeSurveys[idx], \
                                            survey_area_sky, 1.0, N_dish, D_dish, t_total, D_res)  # Set missing parameters as default
                    elif surveys[idx] == "galaxy":
                        PN = shot_noise(mainPath, zj, cosmologY, specsSurvey)

                    # Calculate the comoving distance using cosmologY
                    chi = comoving_distance(zj, cosmologY)

                    Cov_P = ((P + PN) / np.sqrt(N_modes[i]))**2

                    # Compute Fisher matrix components
                    F_matrix[0, 0] += (dP_dm**2) / Cov_P
                    F_matrix[1, 1] += (dP_dn**2) / Cov_P
                    F_matrix[2, 2] += (dP_dP0**2) / Cov_P
                    F_matrix[3, 3] += (dP_dkto**2) / Cov_P

                    # Off-diagonal terms
                    F_matrix[0, 1] += (dP_dm * dP_dn) / Cov_P
                    F_matrix[0, 2] += (dP_dm * dP_dP0) / Cov_P
                    F_matrix[0, 3] += (dP_dm * dP_dkto) / Cov_P
                    F_matrix[1, 2] += (dP_dn * dP_dP0) / Cov_P
                    F_matrix[1, 3] += (dP_dn * dP_dkto) / Cov_P
                    F_matrix[2, 3] += (dP_dP0 * dP_dkto) / Cov_P

                    # Symmetrize the Fisher matrix
                    F_matrix[1, 0] = F_matrix[0, 1]
                    F_matrix[2, 0] = F_matrix[0, 2]
                    F_matrix[3, 0] = F_matrix[0, 3]
                    F_matrix[2, 1] = F_matrix[1, 2]
                    F_matrix[3, 1] = F_matrix[1, 3]
                    F_matrix[3, 2] = F_matrix[2, 3]

            # Invert the Fisher matrix to get the covariance matrix
            cov_matrix = np.linalg.inv(F_matrix)

            # Save the covariance matrix to a .dat file
            output_filename = os.path.join(save_dir, f"cov_matrix_{specsSurvey}.dat")
            np.savetxt(output_filename, cov_matrix)
            print(f"Covariance matrix saved to: {output_filename}")

    # Run the Fisher matrix computation
    compute_fisher_matrix()
