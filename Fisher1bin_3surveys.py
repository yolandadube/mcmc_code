import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import classy  # Import the classy module
from scipy.optimize import minimize_scalar
import getdist
from getdist import plots, MCSamples  # Import getdist modules
import matplotlib.pyplot as plt
import os  # To handle file paths

# -----------------------------------------------------------------------------
#  Define Cosmological Parameters
# -----------------------------------------------------------------------------

cosmo_params = {
    'H0'         : 67.4,          # Hubble constant (km/s/Mpc)
    'Om0'        : 0.315,         # Total matter density
    'Ob0'        : 0.0486,        # Baryon density
    'Ocdm0'      : 0.315 - 0.0486,  # CDM Density
    'A_s'        : 2.100549e-9,   # Primordial power spectrum amplitude
    'n_s'        : 0.9665,        # Primordial power spectrum index
    'N_ur'       : 3.046,         # Number of ultra-relativistic species
    'k_pivot'    : 0.05,          # Pivot wavenumber (Mpc^-1)
}

# Unpack for global use, but access via the dictionary for clarity.
H0 = cosmo_params['H0']
Om0 = cosmo_params['Om0']
c  = 299792.458  # km/s

# Redshift bin parameters
z_c     = 0.5
Delta_z = 0.1

# Parameters for foreground and beam factors (adjust as needed)
k_fg = 0.001  # Foreground scale

# Define the SAVE PATH
save_path = "/home/yolanda/Documents/mcmc_code/version1/result/contour/FisherforecastingDATA/"

# -----------------------------------------------------------------------------
#  Functions
# -----------------------------------------------------------------------------

def E(z):
    """Dimensionless Hubble parameter."""
    return np.sqrt(cosmo_params['Om0'] * (1 + z)**3 + (1 - cosmo_params['Om0']))

def H(z):
    """Hubble parameter in km/s/Mpc."""
    return cosmo_params['H0'] * E(z)

def comoving_radial_distance(z):
    """Comoving radial distance in Mpc."""
    integral, err = quad(lambda z_prime: 299792.458 / H(z_prime), 0, z)
    return integral

def b_HI(z):
    """HI bias as a function of redshift. Cunnington 2022"""
    return 0.842 + 0.693*z - 0.046*z**2

def T_HI(z):
    """Mean HI temperature as a function of redshift (in K). From the PDF paper. Equation (3)"""
    return (56.4e-6 * (1+z)**0.62)  # Kelvin

def foreground_factor(k, mu, k_fg):
    """Foreground factor D_fg(k, mu). From the PDF paper. Equation (3.4)"""
    return 1 - np.exp(-(mu * k / k_fg)**2)

def beam_factor(k, mu, z, D_dish):
    """Beam factor D_b(k, mu, z). From the PDF Paper. Equation (3.5)"""
    lambda_21 = 0.21106*3.24078e-23  # Wavelength of 21 cm line in meters
    theta_b   = 1.22 * lambda_21 * (1 + z) / D_dish
    r_z       = comoving_radial_distance(z)
    return np.exp(-(k**2 * (1 - mu**2) * r_z**2 * theta_b**2) / (16 * np.log(2)))

def P_HI_original(k, z, P_matter_interp, mu, f_z, k_fg_val, D_dish_val):
    """HI power spectrum with RSD, foregrounds, and beam."""
    b_HI_z  = b_HI(z)
    T_HI_z  = T_HI(z)
    D_fg    = foreground_factor(k, mu, k_fg_val)
    D_b     = beam_factor(k, mu, z, D_dish_val)
    return (T_HI_z * (b_HI_z + f_z * mu**2))**2 * P_matter_interp(k) * D_fg * D_b * (H(z)/100)**3

def P_HI_monopole(k, z, P_matter_interp, f_z, k_fg_val, D_dish_val):
    """HI power spectrum monopole (integral over mu)."""
    def integrand(mu):
        return P_HI_original(k, z, P_matter_interp, mu, f_z, k_fg_val, D_dish_val)
    integral, err = quad(integrand, -1, 1)
    return 0.5 * integral

def sigma_P(k, z, P_HI_func, P_N_func, N_modes):
    """
    Variance in the power spectrum measurement.

    P_HI_func: Function that returns the HI power spectrum (P(k)).
    P_N_func: Function that returns the noise power spectrum.
    """
    P_HI_val = P_HI_func(k,z)
    P_N_val  = P_N_func(k,z)
    return (P_HI_val + P_N_val)**2 / N_modes

def dP_dalpha(P0, x, alpha):
    return -x**2 * P0**(1 - alpha * x**2) * np.log(P0)

def dP_dbeta(P0, x, beta):
    return -x**2 * P0**(1 - beta * x**2) * np.log(P0)

def dP_dk0(P0, x, k, k0, alpha, beta):
    log_ratio = np.log(k + 1e-10) / (np.log(k0+ 1e-10) ** 2) #Added a small number to prevent log(0) errors.
    if k < k0:
        return alpha * x * P0**(1 - alpha * x**2) * log_ratio * np.log(P0)
    else:
        return beta * x * P0**(1 - beta * x**2) * log_ratio * np.log(P0)

def fisher_matrix_element(k, z, P_HI_func, P_N_func, N_modes, sigma_P_func,
                           param1_index, param2_index, k0, alpha, beta, P0):
    """
    Calculates a single element of the Fisher matrix using analytical derivatives.

    P_HI_func: Function that returns the HI power spectrum (P(k)).
    P_N_func: Function that returns the noise power spectrum.
    N_modes: Number of Fourier modes.
    sigma_P_func: Function that returns the variance in the power spectrum.
    param1_index: Index of the first parameter.
    param2_index: Index of the second parameter.
    k0       : Value of k0 at which to evaluate derivatives.
    alpha    : Value of alpha at which to evaluate derivatives.
    beta     : Value of beta at which to evaluate derivatives.
    P0       : The peak value of the power spectrum.
    """
    x = np.log(k)/np.log(k0) - 1  #Dimensionless wavenumber

    # Calculate derivatives analytically
    if param1_index == 0:  # k0
        dP_dparam1 = dP_dk0(P0, x, k, k0, alpha, beta)
    elif param1_index == 1:  # alpha
        dP_dparam1 = dP_dalpha(P0, x, alpha)
    elif param1_index == 2:  # beta
        dP_dparam1 = dP_dbeta(P0, x, beta)
    else:
        raise ValueError("Invalid param1_index.")

    if param2_index == 0:  # k0
        dP_dparam2 = dP_dk0(P0, x, k, k0, alpha, beta)
    elif param2_index == 1:  # alpha
        dP_dparam2 = dP_dalpha(P0, x, alpha)
    elif param2_index == 2:  # beta
        dP_dparam2 = dP_dbeta(P0, x, beta)
    else:
        raise ValueError("Invalid param2_index.")

    return (1 / sigma_P_func(k, z, P_HI_func, P_N_func, N_modes)) * dP_dparam1 * dP_dparam2

def survey_volume(zc, Delta_z, survey_area):
    """Comoving survey volume in Mpc^3."""
    r_max = comoving_radial_distance(zc + Delta_z/2)
    r_min = comoving_radial_distance(zc - Delta_z/2)
    return survey_area * (r_max**3 - r_min**3) / 3

def delta_k(zc, z_eff, Delta_z, survey_area):
    """
    Fundamental mode in the radial direction based on Cunnington's method.
    
    Parameters:
    zc           : float - Central redshift of the survey
    Delta_z      : float - Redshift bin width
    survey_area  : float - Survey area in steradians
    
    Returns:
    float - Fundamental mode Δk
    """
    # Comoving distances
    r_min = comoving_radial_distance(zc - Delta_z / 2)
    r_max = comoving_radial_distance(zc + Delta_z / 2)
    l_z = r_max - r_min  # Comoving distance between min and max redshift

    # Effective comoving distance to median redshift
    r_eff = comoving_radial_distance(zc)

    # Transverse dimensions based on angular sky coverage
    l_x = r_eff * np.sqrt(survey_area)
    l_y = r_eff * np.sqrt(survey_area)

    # Fundamental mode Δk
    return 2 * np.pi / np.sqrt(l_x**2 + l_y**2 + l_z**2)


def calculate_Nmodes(k, z, delta_k, survey_volume):
    """Number of independent Fourier modes."""
    return survey_volume * 4 * np.pi * k**2 * delta_k / (2 * np.pi)**3

def T_sys(nu):
    """Total system temperature as a function of frequency (nu in GHz)."""
    nu_GHz  = nu  # Input is already in GHz in this function
    T_RX    = 7.5 + 10 * (nu_GHz)**-2.75
    T_spl   = 3
    T_gal   = 25 * (408e-3 / nu_GHz)**2.75  # Converts MHz to GHz
    return T_RX + T_spl + 2.73 + T_gal  # T_CMB = 2.73 K

def P_N_SD(k, z, T_sys_func, survey_area, N_pol, N_dish, t_survey):
    """Thermal noise power spectrum for single-dish."""
    nu_21cm = 1.42040575177  # GHz (rest-frame)
    nu      = nu_21cm / (1 + z)  # Observed frequency
    T_sys   = T_sys_func(nu)
    D_C     = comoving_radial_distance(z)
    lambda_val = (299792458 / (nu * 1e9))*3.24078e-23  # Wavelength in Mpc
    V_pix   = (D_C * lambda_val)**2 * (299792.458 / H(z))
    S_area  = survey_area  # in steradians
    return T_sys**2 * V_pix / (N_pol * N_dish * t_survey * 3600 * eta) * (H(z)/100)**3  # (K)^2 (Mpc/h)^3

def get_survey_params(survey_name):
    """
    Returns survey parameters based on the survey name.
    """
    if survey_name == "MeerKAT_UHF":
        # MeerKAT UHF Band
        sky_area  = 4000 * (np.pi / 180)**2   # steradians (4000 deg^2)
        Ndish     = 64                       # Number of dishes
        t_survey  = 1000                     # hours (Total survey time)
        nu_min    = 580e6                    # Hz (UHF band)
        nu_max    = 1015e6                   # Hz
        eta       = 1.0                      # observing efficiency
        Npol      = 2                        # Number of polarizations
        D_dish    = 13.5                      # Dish diameter
        z_eff     = 0.925
    elif survey_name == "MeerKAT_L":
        # MeerKAT L Band
        sky_area  = 4000 * (np.pi / 180)**2   # steradians (4000 deg^2)
        Ndish     = 64                       # Number of dishes
        t_survey  = 1000                     # hours (Total survey time)
        nu_min    = 900e6                    # Hz (L Band)
        nu_max    = 1670e6                   # Hz
        eta       = 1.0                      # observing efficiency
        Npol      = 2                        # Number of polarizations
        D_dish    = 13.5                      # Dish diameter
        z_eff     = 0.39
    elif survey_name == "SKAO":
        # SKAO
        sky_area  = 20000 * (np.pi / 180)**2  # steradians (20000 deg^2)
        Ndish     = 197                      # Number of dishes
        t_survey  = 10000                    # hours (Total survey time)
        nu_min    = 50e6                     # Hz
        nu_max    = 350e6                    # Hz
        eta       = 1.0                      # observing efficiency
        Npol      = 2                        # Number of polarizations
        D_dish    = 15                        # Dish diameter
        z_eff     = 1.675
    else:
        raise ValueError("Invalid survey_name.")

    return sky_area, Ndish, t_survey, nu_min, nu_max, eta, Npol, D_dish, z_eff

 # -----------------------------------------------------------------------------
#  Main part of the code
# -----------------------------------------------------------------------------

def fisher_forecast(params, k_min, k_max, n_k_bins, P_matter_interp, survey_name, f_z):
    """
    Performs the Fisher forecast using analytical derivatives.

    params:  List of parameters to forecast (k0, alpha, beta).
    k_min: Minimum k value (Mpc^-1).
    k_max: Maximum k value (Mpc^-1).
    n_k_bins: Number of k bins to integrate over.
    P_matter_interp: Interpolated matter power spectrum function.
    survey_name: Name of the survey to use.
    f_z: Growth rate at redshift z.
    """

    k0, alpha, beta = params

    # 1.  Calculate survey parameters
    sky_area, Ndish, t_survey, nu_min, nu_max, eta, Npol, D_dish, z_eff = get_survey_params(survey_name)
    V_survey = survey_volume(z_c, Delta_z, sky_area)
    dk = delta_k(z_c, z_eff, Delta_z, sky_area)

    # 2.  Set up k values
    k_values = np.logspace(np.log10(k_min), np.log10(k_max), n_k_bins) #Log-spaced k values
    #Fixed values,
    k_fg_val = 0.001 #Set fixed k_fg value.
    k0 = 0.016 #h/Mpc
    # 3.  Find the peak of the HI power spectrum (monopole)
    def P_HI_neg(k): #Negative of the power spectrum for minimization
        return -P_HI_monopole(k, z_c, P_matter_interp, f_z, k_fg, D_dish)

    #The peak *location* is approximately k0.  Search around that value.
    result = minimize_scalar(P_HI_neg, bounds=(k0*0.5, k0*1.5), method='bounded')
    k_peak = result.x
    P0 = P_HI_monopole(k_peak, z_c, P_matter_interp, f_z, k_fg, D_dish)  #Peak of the monopole

    print(f"Peak of HI power spectrum found at k = {k_peak:.4f} Mpc^-1, P0 = {P0:.4e} K^2 (Mpc/h)^3")

    # 4. Define functions that depend on k and z
    # Now, these functions are *fixed* at the peak.
    def P_HI_k_z(k, z):
      return P_HI_monopole(k_peak, z, P_matter_interp, f_z, k_fg, D_dish) #EVALUATE AT THE PEAK, not at k.

    def P_N_k_z(k, z):
      return P_N_SD(k, z, T_sys, sky_area, Npol, Ndish, t_survey) #Pass T_sys as a function

    def sigma_P_k_z(k, z, P_HI_func, P_N_func, N_modes):
      P_HI_val = P_HI_func(k,z)
      P_N_val  = P_N_func(k,z)
      return (P_HI_val + P_N_val)**2 / N_modes

    # 5.  Calculate the Fisher matrix
    n_params = len(params)
    Fisher_matrix = np.zeros((n_params, n_params))

    for i in range(n_params):
        for j in range(n_params):
            for k in k_values:
                N_modes = calculate_Nmodes(k, z_c, dk, V_survey)
                Fisher_matrix[i, j] += fisher_matrix_element(k, z_c, P_HI_k_z, P_N_k_z, N_modes, sigma_P_k_z, i, j, k0, alpha, beta, P0)

    # 6.  Invert the Fisher matrix
    try:
        covariance_matrix = np.linalg.inv(Fisher_matrix)
        inv_check = Fisher_matrix @ covariance_matrix
    except np.linalg.LinAlgError:
        print("Fisher matrix is singular (non-invertible).  Check your parameters and survey setup.")
        return None, None

    # 7.  Calculate the forecasted uncertainties
    uncertainties = np.sqrt(np.diag(covariance_matrix))

    return covariance_matrix, uncertainties

 # -----------------------------------------------------------------------------
#  Example Usage (YOU NEED TO MODIFY THIS)
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # 1. Define parameters to forecast (k0, alpha, beta) - REPLACE WITH YOUR PARAMETERS
    params_to_forecast = [0.016, 1.0, 2.0]  # Example values: [k0, alpha, beta]
    param_names        = ["k0", "alpha", "beta"] #For clearer output

    # 2. Define Surveys
    survey_names = ["MeerKAT_UHF", "MeerKAT_L", "SKAO"]  # List of surveys to iterate over
    
    # 3.  Set k range and number of bins
    k_min_val   = 0.005  # Mpc^-1
    k_max_val   = 0.025    # Mpc^-1
    n_k_bins_val     = 500

    # 4. Get P_matter from CLASS
    # Initialize CLASS parameters
    params_CLASS = {
        'output'         : 'mPk',
        'z_max_pk'       : 3.0, #Adjust this as needed
        'P_k_max_1/Mpc'  : 5.0, #Adjust this as needed
        'H0'             : cosmo_params['H0'],
        'Omega_b'        : cosmo_params['Ob0'], #Planck 2018 values
        'Omega_cdm'      : cosmo_params['Ocdm0'],  #Om0 defined above
        'A_s'            : cosmo_params['A_s'], #Planck 2018
        'n_s'            : cosmo_params['n_s'], #Planck 2018
        'N_ur'           : cosmo_params['N_ur'],
        'k_pivot'        : cosmo_params['k_pivot'],
        'lensing'        : 'no',
        'non linear'     : 'halofit',
    }

    # Create a CLASS instance
    cosmo = classy.Class()
    cosmo.set(params_CLASS)
    cosmo.compute()

    # Create the interpolated matter power spectrum *function*
    def P_matter_CLASS(k, z):
        """
        Returns the matter power spectrum at redshift z obtained from CLASS.
        k is in Mpc^-1, P(k) is in (Mpc/h)^3
        """
        return cosmo.pk(k*cosmo_params['H0']/100, z) * (cosmo_params['H0']/100)**3 # k in h/Mpc, P(k) in (Mpc/h)^3

   #Manually compute the growth rate f = Omega_m(z)^gamma, gamma ~ 0.55
    z = z_c #Redshift of interest
    Omega_m_z = cosmo_params['Om0'] * (1 + z)**3 / (np.sqrt(cosmo_params['Om0'] * (1 + z)**3 + (1 - cosmo_params['Om0'])))**2
    f_z = Omega_m_z**0.55
    print(f"Growth rate f(z={z_c}) = {f_z:.4f}")

     #Vectorize to speed things up
    k_CLASS = np.logspace(-3, -1, 200)  # Focus on turnover region
    P_matter_values = np.array([P_matter_CLASS(k, z_c) for k in k_CLASS]) #Evaluate at z_c

    P_matter_interp = interp1d(k_CLASS, P_matter_values, kind='cubic', fill_value="extrapolate") #Linear is fine

    # 4. Iterate over surveys
    for survey_name in survey_names:
        print(f"\n-------------------- Survey: {survey_name} --------------------")

        # Get survey parameters (including D_dish)
        sky_area, Ndish, t_survey, nu_min, nu_max, eta, Npol, D_dish, z_eff = get_survey_params(survey_name)

        # Perform Fisher forecasting
        covariance_matrix, uncertainties = fisher_forecast(params_to_forecast, k_min_val, k_max_val, n_k_bins_val, P_matter_interp, survey_name, f_z)

        # Print Results
        print("Covariance Matrix:\n", covariance_matrix)
        print("\nForecasted Uncertainties:")
        if uncertainties is not None:
            for i, unc in enumerate(uncertainties):
                print(f"  sigma({param_names[i]}): {unc:.4f}")

        # Construct file name for covariance matrix
        filename = os.path.join(save_path, f"covariance_{survey_name}.dat")

        # Check if directory exists
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        # Saving data
        if covariance_matrix is not None:
            np.savetxt(filename, covariance_matrix)
            print(f"Covariance matrix saved to {filename}")
        else:
            print(f"No covariance matrix to save for {survey_name}")
