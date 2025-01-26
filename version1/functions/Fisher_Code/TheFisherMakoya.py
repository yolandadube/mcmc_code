# Relevant python modules
import os
import sys
import emcee
import numpy as np
from scipy.integrate import quad
from classy import Class
import getdist.plots as gdplots
import matplotlib.pyplot as plt
import matplotlib.font_manager
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
from getdist import plots, MCSamples
from colossus.cosmology import cosmology

# Define the cosmological parameters
cosmo = Class()
params = {
    'output'      : 'mPk , mTk',         # Request matter power spectrum calculation
    'omega_b'     : 0.0223828,           # Baryon density
    'omega_cdm'   : 0.1201075,           # Cold Dark Matter density
    'h'           : 0.6736,              # Hubble parameter
    'A_s'         : 2.1e-9,              # Amplitude of primordial power spectrum
    'n_s'         : 0.9665,              # Spectral index of primordial power spectrum
    'tau_reio'    : 0.0543,              # Optical depth to reionization
    'Omega_Lambda': 0.6825,              # Present-day density parameter for dark energy
    'Omega_k'     : 0.0,                 # Present-day density parameter for curvature
    'z_max_pk'    : 0.75,                # Maximum redshift for storing power spectrum
}

cosmo.set(params)

# Add paths to your custom modules
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/surveySpecifications/')
sys.path.append('/home/yolanda/Documents/mcmc_code/version1/functions/surveyRedshifts/')

from surveyRedshiftRange import redshiftRange
from specifications import surveySpecs

def fisherCode(params, k, alpha, beta, P0, k0, z, survey, specsSurvey):

    # Define functions to calculate cosmological quantities
    def E(z, cosmo):
        Omega_k0 = params['Omega_k']
        Omega_Lambda0 = params['Omega_Lambda']
        Omega_m0 = params['omega_b'] + params['omega_cdm']
        return np.sqrt(Omega_m0 * (1 + z) ** 3 + Omega_Lambda0 + Omega_k0 * (1 + z) ** 2)

    def H(z, cosmo):
        H0 = params['h']/2997.92458 # H0 in h/Mpc
        E_z = E(z, cosmo)
        return H0 * E_z

    def comoving_distance(z, cosmo):
        integrand = lambda z_prime: 1.0 / H(z_prime, cosmo)
        r, _ = quad(integrand, 0, z)
        c = 2997.92458
        return (c / (100 * params['h'])) * r

    def survey_volume(survey, specsSurvey, cosmo):
        # Retrieve redshift range
        z_min, z_max = redshiftRange(survey, specsSurvey)

        # Retrieve survey specifications
        N_dish, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey)

        # Calculate comoving distances
        r_min = comoving_distance(z_min, cosmo)  # Comoving distance in Mpc/h
        r_max = comoving_distance(z_max, cosmo)  # Comoving distance in Mpc/h

        # Calculate the fraction of sky covered by the survey
        total_area_of_sky = 4.0 * np.pi * (180.0 / np.pi) ** 2.0     # Total area of the sky in square degrees
        f_sky = survey_area_sky / total_area_of_sky

        # Calculate the survey volume in (Mpc/h)^3
        volume_mpc_h3 = 221.6

        print(f"Survey Volume (Mpc/h)^3: {volume_mpc_h3}")

        return volume_mpc_h3


    # Defining k_min
    def l_x(z, specsSurvey, cosmo):
        N_dish, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey)
        lx = comoving_distance(zeff, cosmo) * np.sqrt(survey_area_sky)
        return lx

    def l_z(survey, specsSurvey, cosmo):
        z_min, z_max = redshiftRange(survey, specsSurvey)
        r_min = comoving_distance(z_min, cosmo)
        r_max = comoving_distance(z_max, cosmo)
        l_z = r_max - r_min
        return l_z

    def kmin(z, specsSurvey, cosmo):
        l_x_val = l_x(z, specsSurvey, cosmo)
        l_z_val = l_z(survey, specsSurvey, cosmo)
        k_min = 2 * np.pi / np.sqrt(2 * l_x_val ** 2 + l_z_val ** 2)
        return k_min

    def Nmodes(k, z, survey, specsSurvey, cosmo):
        delta = 2 * kmin(z, specsSurvey, cosmo)
        volume = survey_volume(survey, specsSurvey, cosmo)
        N_modes = (4*np.pi*volume * k ** 2 * delta) / (2 * np.pi) ** 3
        return N_modes

    def Tsys(z):
        T_spl = 2
        nu21 = 1420
        T_CMB = 2.73
        T_gal = 25 * (408 / nu21) ** 2.75
        T_rx = 7.5 + 10 * (nu21 / 1e3 - 0.75) ** 2
        T_sys = T_spl + T_CMB + T_gal + T_rx
        return T_sys

    def P_N(z, specsSurvey):
        eta = 2
        N_pol = 1
        lambda21 = 0.21 * (1 + z)
        N_dish, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey)
        P_N = Tsys(z)**2 * comoving_distance(z, cosmo) ** 2 * lambda21 * ((1 + z) / H(z, cosmo)) * (4 * np.pi * survey_area_sky / eta * N_pol * N_dish * t_total)
        return P_N

    # Function for the power spectrum turnover region
    def P_fit(k, k0, alpha, beta, P0):
        x = (np.log(k) - np.log(k0)) / np.log(k0)
        return np.where(k < k0, P0 ** (1 - alpha * x**2), P0 ** (1 - beta * x**2))

    def error(k, z, survey, specsSurvey, cosmo, k0, alpha, beta, P0):
        err = (P_fit(k, k0, alpha, beta, P0) + P_N(z, specsSurvey)) / 2*np.sqrt(Nmodes(k, z, survey, specsSurvey, cosmo))
        return err

    def calculate_fisher_matrix(k_vals, z, survey, specsSurvey, cosmo, k0, alpha, beta, P0):

        Fisher_matrix = np.zeros((4, 4))

        for i, k in enumerate(k_vals):

            x = (np.log(k) - np.log(k0)) / np.log(k0)

            dP_dm = -x**2 * (P0)**(1 - alpha * x**2) * np.log(P0)
            dP_dn = -x**2 * (P0)**(1 - beta * x**2) * np.log(P0)
            dP_dP0 = np.where(k < k0, (1 - alpha * x**2) * (P0)**(-alpha * x**2), (1 - beta * x**2) * (P0)**(-beta * x**2))
            dP_dk0 = np.where(k < k0, -alpha * x * (P0)**(1 - alpha * x**2) * -((np.log(k) + 1) / ((np.log(k0))**2)) * np.log(P0), -beta * x * (P0)**(1 - beta * x**2) * -((np.log(k) + 1) / ((np.log(k0))**2)) * np.log(P0))

            err = error(k, z, survey, specsSurvey, cosmo, k0, alpha, beta, P0)

            Fisher_matrix[0, 0] += np.sum(dP_dm * dP_dm / err ** 2)
            Fisher_matrix[0, 1] += np.sum(dP_dm * dP_dn / err ** 2)
            Fisher_matrix[0, 2] += np.sum(dP_dm * dP_dP0 / err ** 2)
            Fisher_matrix[0, 3] += np.sum(dP_dm * dP_dk0 / err ** 2)

            Fisher_matrix[1, 1] += np.sum(dP_dn * dP_dn / err ** 2)
            Fisher_matrix[1, 2] += np.sum(dP_dn * dP_dP0 / err ** 2)
            Fisher_matrix[1, 3] += np.sum(dP_dn * dP_dk0 / err ** 2)

            Fisher_matrix[2, 2] += np.sum(dP_dP0 * dP_dP0 / err ** 2)
            Fisher_matrix[2, 3] += np.sum(dP_dP0 * dP_dk0 / err ** 2)

            Fisher_matrix[3, 3] += np.sum(dP_dk0 * dP_dk0 / err ** 2)

        Fisher_matrix = np.triu(Fisher_matrix) + np.triu(Fisher_matrix, 1).T
        return Fisher_matrix

    def invert_fisher_matrix(Fisher_matrix):
        # Check for eigenvalues before inverting
        eigenvalues, _ = np.linalg.eigh(Fisher_matrix)
        if np.any(eigenvalues <= 0):
            print("Warning: Fisher matrix has non-positive eigenvalues. Regularization may be required.")
            # Add a small regularization term if needed
            Fisher_matrix += np.eye(Fisher_matrix.shape[0]) * 1e-6

        # Invert the Fisher matrix to obtain the covariance matrix
        covariance_matrix = np.linalg.inv(Fisher_matrix)
        return covariance_matrix

    def generate_triangle_plot(covariance_matrix, labels):
        # Check if covariance matrix is positive-semidefinite
        if np.all(np.linalg.eigvals(covariance_matrix) > 0):
            # Generate samples from the covariance matrix
            samples = np.random.multivariate_normal(mean=[0] * covariance_matrix.shape[0], cov=covariance_matrix, size=10000)

            # Create a GetDist MCSamples object
            gdist_samples = MCSamples(samples=samples, names=labels, labels=labels)

            # Create a triangle plot
            g = plots.get_subplot_plotter()
            g.triangle_plot(gdist_samples, filled=True)
            g.export('triangle_plot.png')  # Save the plot as a PNG file
        else:
            print("Covariance matrix is not positive-semidefinite; cannot generate triangle plot.")

    # Calculate Fisher matrix
    Fisher_matrix = calculate_fisher_matrix(k, z, survey, specsSurvey, cosmo, k0, alpha, beta, P0)

    # Invert Fisher matrix to get covariance matrix
    covariance_matrix = invert_fisher_matrix(Fisher_matrix)

    # Generate triangle plot
    generate_triangle_plot(covariance_matrix, labels=[r'm', r'n', r'P_0', r'k_{T0}'])

# Example usage:
k = np.linspace(0.005, 0.025, 1000)  # Example k values
z = 0.5  # Example redshift
survey = "HI IM"
specsSurvey = "skaOIMBand1"
cosmo = cosmo  # Use the cosmology object defined earlier
k0 = 0.0163  # Example k0
alpha = 0.75  # Example alpha
beta = 0.85  # Example beta
P0 = 37300  # Example P0

# Example usage of the function
fisherCode(params, k, alpha, beta, P0, k0, z, survey, specsSurvey)
