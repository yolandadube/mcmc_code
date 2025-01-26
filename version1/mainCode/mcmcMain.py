import glob
import time
import os, sys
import numpy as np
from classy import Class
from tqdm import tqdm
from natsort import natsorted

# main path directory
mainPath = "/home/yolanda/Documents/mcmc_code/version1/"

# import functions
sys.path.append(mainPath+"functions/planckBestFit/")
sys.path.append(mainPath+"functions/pyCLASSWrapper/")
sys.path.append(mainPath+"functions/surveyRedshifts/")
sys.path.append(mainPath+"functions/markovChainMonteCarlo/")
sys.path.append(mainPath+"functions/Fisher_Code/")

#---#
from MCMC import mcmcAlgorithm
from compileCLASS import pyCompileCLASS
from FisherForecastingNew import fisherForecasting
#from surveyRedshiftRange import redshiftRange
from planckCosmologicalParameters import planckValues

# compile CLASS
compilation = pyCompileCLASS(mainPath)


# dictionary fiducial Planck cosmological parameters
dictionary    = {'H0'         :67.66,     # hubble parameter today
                 'omb0'             :0.02242,   # is the baryonic matter density parameter
                 'h'                     :0.6766 ,   # 
                 'omm0'            :0.3130,    # is the total matter density parameter
                 'omcdm0'        :0.2640,    # is the cold dark matter density parameter
                 'omega_k0'     :0.0,       # omega curvature today
                 'omega_n0'     :0.0,       # omega neutrino today
                 'n_s'                 :0.9665,    # spectral index of primordial power spectrum
                 'A_s'                 :2.105e-09, # A_s 
                 'sigma80'         :0.8102,    # sigma_8 today
                 'gamma'           :0.55,      # linear growth index
                 'w'                     :-1.0,      # equation of state parameter for dark energy
                 'fnl'                    :0.0        # primordial non-Gaussianity parameter
                }


PlanckData = planckValues(dictionary)

# gauge of perturbations: 'newtonian' or 'synchronous'
gauge      = 'synchronous'

# 1st order contributions in the source number count fluctuations
NewtonianO1 = True

# MCMC parameters list
mcmcParameters  = ['p0', 'k0', 'alpha', 'beta']
mcmcParameters1 = ['A', 'k0', 'n', 'm', 'beta']    # Parameters for the empirical model
mcmcParameters2 = ['A', 'k0', 'n', 'm', 'Delta']   # Parameters for the smoothly broken power law

# priors for MCMC parameters
priors  = {    'p0'    :[1.0e+03, 1.0e+06],
               'k0'         :[5.0e-03, 5.0e-02], 
               'alpha'    :[-10, 10],
               'beta'      :[-10, 10],  
            }

# Priors for the empirical model
priors1 = {
                'A'      : [6.0e04, 5.0e+06],     # Amplitude range
                'k0'     : [1.0e-03, 3.5e-02],    # Turnover point range
                'n'      : [0.30, 2.0],              # Exponent before turnover
                'm'      : [1.0, 6.0],           # Controls sharpness of transition
                'beta'   : [0, 10.0],                # Exponent after turnover
            }

# Priors for the smoothly broken power law
priors2 = {
                'A'      : [2e04, 8e+4],        # Amplitude range
                'k0'     : [0.010, 0.025],          # Turnover point range
                'n'      : [0.0, 5.0],            # Exponent before turnover
                'm'      : [-0.8, 0.0],            # Exponent after turnover
                'Delta'  : [0.0, 3.0],             # Controls smoothness of transition
            }


# fingers-of-God effect
fingerOfGod = False

# foreground effect
foreGround  = True

# beam effect
beam        = True

# non-linear cut-off [h Mpc^{-1}]
k_NL0       = 0.08 

# foreground scale cut-off [h Mpc^{-1}]
#k_fg       = 0.005

# wedge parameter (Nwedge = 0, 1, 3)
Nwedge      = 1.0

# do non-linear matter power spectrum
nonLinear   = False

#--------------------------------------------------------------------------------------------------------------------------

# number of Markov chains
chains           = 75

# maximum number of walks per Markov chain
walks            = 10000

# redshift bin size
delta_z          = 1.325

# redshift bin centre
z                = np.array([1.675])

# generate the array of k
kMin, kMax       = [5.0e-03, 3.5e-02]
k                = np.linspace(kMin, kMax, 5000, endpoint=True)

# Define a list of survey configurations

# specifications of survey: 
# GALAXY:                 : "HAlpha", "desiBGS", "desiELG",
#                           "skaOGalBand2", "megaMapLBG"
# HI IM                   : "skaOIMBand1", "skaOIMBand2", 
#                           "meerKATLBand", "meerKATUHFBand", "dsa2000",  
#                          "hirax256", "hirax1024", "puma5k", 'puma32k'
'''

    {
        "surveys"     : ["HI IM", "HI IM"],
        "modeSurveys" : ["single dish", "single dish"],
        "specsSurveys": ["skaOIMBand1", "skaOIMBand1"],
        "k_fg_values" : [0.001]
    },
    {
        "surveys": ["HI IM", "HI IM"],
        "modeSurveys": ["single dish", "single dish"],
        "specsSurveys": ["meerKATLBand", "meerKATLBand"],
        "k_fg_values": [0.001, 0.005]
        # "k_fg_values": np.linspace(0.01, 0.05, 5).tolist()  # Alternative way to generate k_fg_values
    },
    {
        "surveys": ["HI IM", "HI IM"],
        "modeSurveys": ["single dish", "single dish"],
        "specsSurveys": ["meerKATUHFBand", "meerKATUHFBand"],
        "k_fg_values": [0.001, 0.005]
        # "k_fg_values": np.linspace(0.01, 0.05, 5).tolist()  # Alternative way to generate k_fg_values
    },
    {
        "surveys"     : ["galaxy", "galaxy"],
        "modeSurveys" : ["single dish", "single dish"],
        "specsSurveys": ["HAlpha", "HAlpha"],
        "k_fg_values" : [0.001]
    },
    {
        "surveys": ["galaxy", "galaxy"],
        "modeSurveys": ["single dish", "single dish"],
        "specsSurveys": ["desiBGS", "desiBGS"],
        "k_fg_values": [0.001]
        # "k_fg_values": np.linspace(0.01, 0.05, 5).tolist()  # Alternative way to generate k_fg_values
    },
    {
        "surveys": ["galaxy", "galaxy"],
        "modeSurveys": ["single dish", "single dish"],
        "specsSurveys": ["desiELG", "desiELG"],
        "k_fg_values": [0.001]
        # "k_fg_values": np.linspace(0.01, 0.05, 5).tolist()  # Alternative way to generate k_fg_values
    }
'''

survey_configurations = [
    {
        "surveys"     : ["HI IM", "HI IM"],
        "modeSurveys" : ["single dish", "single dish"],
        "specsSurveys": ["skaOIMBand1", "skaOIMBand1"],
        "k_fg_values" : [0.001]
    },   
]

for config in survey_configurations:

    fisherForecasting(mainPath, dictionary, PlanckData, z, delta_z, kMin, kMax, config["surveys"], config["specsSurveys"], config["modeSurveys"])

    # Iterate over each k_fg value for the current survey configuration

    
    for k_fg_value in config["k_fg_values"]:
        
        start = time.time()

        # Run MCMC algorithm with the current k_fg value
        flatSample = mcmcAlgorithm(mainPath, mcmcParameters, priors, chains, walks, PlanckData, gauge, nonLinear, k, \
                                   k_NL0, k_fg_value, Nwedge, z, delta_z, dictionary, config["surveys"], config["specsSurveys"], config["modeSurveys"], \
                                   fingerOfGod, foreGround, beam, NewtonianO1, file_suffix=k_fg_value)

        end = time.time()

        timeTaken = end - start

        if timeTaken > 60.0:
            print("Survey: {}, k_fg: {:.3f}, run time = {:.1f} min".format(config["surveys"], k_fg_value, (end - start) / 60.0))
        else:
            print("Survey: {}, k_fg: {:.3f}, run time = {:.1f} s".format(config["surveys"], k_fg_value, end - start))
    
