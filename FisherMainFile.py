import glob
import time
import os, sys
import numpy as np
from tqdm import tqdm
from natsort import natsorted
from concurrent.futures import ThreadPoolExecutor

# main path directory
mainPath = "/home/yolanda/Documents/mcmc_code/version1/"

# import functions
# Add paths to different function directories
sys.path.append(mainPath + "functions/daTa/")
sys.path.append(mainPath + "functions/noise/")
sys.path.append(mainPath + "functions/nModes/")
sys.path.append(mainPath + "functions/Theory/")
sys.path.append(mainPath + "functions/FisherCodes/")
sys.path.append(mainPath + "functions/planckBestFit/")
sys.path.append(mainPath + "functions/pyCLASSWrapper/")
sys.path.append(mainPath + "functions/surveyRedshifts/")
sys.path.append(mainPath + "functions/surveyVolume/")
sys.path.append(mainPath + "functions/backgroundCosmology/")
sys.path.append(mainPath + "functions/surveySpecifications/")
sys.path.append(mainPath + "functions/markovChainMonteCarlo/")
sys.path.append(mainPath + "functions/tracerPowerSpectrum/monopole/")

# Import specific functions
from MCMC import mcmcAlgorithm
from compileCLASS import pyCompileCLASS
from planckCosmologicalParameters import planckValues
from TurnoverFisher import fisherForecasting

# compile CLASS
compilation = pyCompileCLASS(mainPath)

# dictionary fiducial Planck cosmological parameters
dictionary    = {'H0'         :67.66,     # hubble parameter today
                 'omb'        :0.02242,   # physical omega baryonic today, \omega_b = \Omega_b0 * h^{2} 
                 'omcdm'      :0.11933,   # physical omega cold dark matter today, \omega_cdm = \Omega_cdm0 * h^{2}
                 'omega_k0'   :0.0,       # omega curvature today
                 'omega_n0'   :0.0,       # omega neutrino today
                 'n_s'        :0.9665,    # spectral index of primordial power spectrum
                 'A_s'        :2.105e-09, # A_s 
                 'sigma80'    :0.8102,    # sigma_8 today
                 'gamma'      :0.55,      # linear growth index
                 'w'          :-1.0,      # equation of state parameter for dark energy
                 'fnl'        :0.0        # primordial non-Gaussianity parameter
                }

PlanckData = planckValues(dictionary)

# gauge of perturbations: 'newtonian' or 'synchronous'
gauge      = 'synchronous'

# 1st order contributions in the source number count fluctuations
NewtonianO1 = True

# MCMC parameters list
mcmcParameters  = ['p0', 'k0', 'alpha', 'beta']

# priors for MCMC parameters
priors          = {'p0'    :[1.0e+03, 1.0e+06],
                   'k0'    :[5.0e-03, 5.0e-02], 
                   'alpha' :[-5, 5],
                   'beta'  :[-5, 5]    
                  }

# fingers-of-God effect
fingerOfGod = False

# foreground effect
foreGround  = True

# beam effect
beam        = True

# non-linear cut-off [h Mpc^{-1}]
k_NL0       = 0.08 

# wedge parameter (Nwedge = 0, 1, 3)
Nwedge      = 1.0

# do non-linear matter power spectrum
nonLinear   = False

#--------------------------------------------------------------------------------------------------------------------------

# number of Markov chains
chains      = 55

# maximum number of walks per Markov chain
walks       = 10000

# redshift bin size
delta_z     = 0.1

# redshift bin centre
z             = np.array([0.25])

#LBAND
#z               = np.array([0.25, 0.35, 0.45, 0.55])

#UHF
#z = np.array([0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35])

#SKAO
#z = np.array([
    #0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00,
    #1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80,
    #1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60,
    #2.70, 2.80, 2.90, 3.00])


#ELG
#z = np.array([0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65])

#BGS
#z = np.array([0.05, 0.15, 0.25, 0.35, 0.45])

#HAlpha
#z = np.array([0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75])



# generate the array of k
kMin, kMax       = [6.0e-03, 5.0e-02]
k                = np.linspace(kMin, kMax, 700, endpoint=True)

# Define a list of survey configurations

# specifications of survey: 
# GALAXY:                 : "HAlpha", "desiBGS", "desiELG",
#                           "skaOGalBand2", "megaMapLBG"
# HI IM                   : "skaOIMBand1", "skaOIMBand2", 
#                           "meerKATLBand", "meerKATUHFBand", "dsa2000",  
#                           "hirax256", "hirax1024", "puma5k", 'puma32k'

survey_configurations = [
    {
        "surveys": ["HI IM", "HI IM"],
        "modeSurveys": ["single dish", "single dish"],
        "specsSurveys": ["meerKATLBand", "meerKATLBand"],
        "k_fg_values": [0.005]
    }
]

# Keep track of processed (survey, z) combinations
processed_bins = set()

def run_fisher_forecasting(z_value, k_fg_value, config):
    print(f"Running Fisher forecasting for z = {z_value:.2f}, k_fg = {k_fg_value:.3f}")
    fisherForecasting(mainPath, dictionary, PlanckData, np.array([z_value]), delta_z, k, 
                      config["surveys"], config["specsSurveys"], config["modeSurveys"],
                      gauge, nonLinear, k_NL0, k_fg_value, Nwedge, fingerOfGod, 
                      foreGround, NewtonianO1, beam)

for config in survey_configurations:
    for k_fg_value in config["k_fg_values"]:
        start = time.time()

        # Commented out Fisher forecasting
        """
        print(f"Running Fisher forecasting for k_fg = {k_fg_value:.3f}")
        fisherForecasting(mainPath, dictionary, PlanckData, z, delta_z, k, 
                          config["surveys"], config["specsSurveys"], config["modeSurveys"],
                          gauge, nonLinear, k_NL0, k_fg_value, Nwedge, fingerOfGod, 
                          foreGround, NewtonianO1, beam)
        """

        # Run MCMC algorithm for all redshift values in z array
        print(f"Running MCMC for k_fg = {k_fg_value:.3f}")
        flatSample = mcmcAlgorithm(mainPath, mcmcParameters, priors, chains, walks, PlanckData, gauge, nonLinear, k, 
                                   k_NL0, k_fg_value, Nwedge, z, delta_z, dictionary, 
                                   config["surveys"], config["specsSurveys"], config["modeSurveys"], 
                                   fingerOfGod, foreGround, beam, NewtonianO1, file_suffix=f"kfg{k_fg_value:.3f}")

        end = time.time()

        timeTaken = end - start

        if timeTaken > 60.0:
            print("Survey: {}, k_fg: {:.3f}, Total run time = {:.1f} min".format(config["surveys"], k_fg_value, (end - start) / 60.0))
        else:
            print("Survey: {}, k_fg: {:.3f}, Total run time = {:.1f} s".format(config["surveys"], k_fg_value, end - start))
