import glob
import emcee
import os, sys
import numpy as np
import arviz as az
import seaborn as sns
import matplotlib.cm as cm
import multiprocessing as mp
from natsort import natsorted
import matplotlib.font_manager
import matplotlib.pyplot as pl
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator
from matplotlib.ticker import LogFormatter 
from matplotlib.pyplot import rc, rcParams
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
#--------------------------------------------------------------------------------------------------------------------------

fig = pl.figure(figsize = (15.0, 12.0))

pl.rcParams['xtick.major.size']  = 3
pl.rcParams['xtick.major.width'] = 2
pl.rcParams['ytick.major.size']  = 3
pl.rcParams['ytick.major.width'] = 2
pl.rcParams['xtick.major.pad']   = 6
pl.rcParams['ytick.major.pad']   = 6
# ...
pl.rcParams['xtick.minor.size']  = 3
pl.rcParams['xtick.minor.width'] = .3
pl.rcParams['ytick.minor.size']  = 3
pl.rcParams['ytick.minor.width'] = .3

# ...
pl.rcParams['figure.dpi']            = 80
pl.rcParams['lines.solid_joinstyle'] = 'miter'  
pl.rcParams['lines.antialiased']     = True
pl.rc('xtick', labelsize = 20)
pl.rc('ytick', labelsize = 20)

pl.rcParams['ps.useafm']          = True
pl.rcParams['pdf.use14corefonts'] = True
pl.rcParams['text.usetex']        = True

#pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pl.tick_params(axis = 'both', which = 'both', direction = 'in')
pl.rcParams['axes.linewidth'] = 1.4 #set the value globally
#--------------------------------------------------------------------------------------------------------------------------

def mcmcAlgorithm2(mainPath, mcmcParameters2, priors2, chains, walks, cosmologY, gauge, nonLinear, k, k_NL0, k_fg, \
                  Nwedge, z, delta_z, dictionary, surveys, specsSurveys, modeSurveys, fingerOfGod, \
                  foreGround, beam, NewtonianO1, file_suffix=""):

    # chains     ---> number of Markov chains
    # walks      ---> number of walks per Markov chain
    # priors     ---> dictionary that contains prior 
    #               information on each fitting parameter
    # mcmcParameters ---> list of parameters for MCMC

    sys.path.append(mainPath+"functions/daTa/")
    sys.path.append(mainPath+"functions/guess/")
    sys.path.append(mainPath+"functions/latex/")
    sys.path.append(mainPath+"functions/probability/")
    sys.path.append(mainPath+"functions/pyCLASSWrapper/")
    #---#
    from runCLASS import pyCLASS
    from Data import dataPTracer
    from guess2 import initialGuess2
    from latexFonts2 import latexDictionary2
    from probability2 import logProbability2
    #---#

    # surveys
    survey1, survey2           = surveys
    modeSurvey1, modeSurvey2   = modeSurveys
    specsSurvey1, specsSurvey2 = specsSurveys

    # compute matter pk and tk from CLASS
    pyCLASS(mainPath, cosmologY, gauge, nonLinear, z)

    # get all the data files for pk and tk
    dataFiles      = natsorted(glob.glob(mainPath+'CLASSDataFile/outputs/*'))
    pkFile, tkFile = dataFiles

    # simulate data + err for power spectrum
    PData          = dataPTracer(mainPath, cosmologY, gauge, nonLinear, pkFile, tkFile, \
                                 k, k_NL0, k_fg, Nwedge, z, delta_z, dictionary, \
                                 surveys, specsSurveys, modeSurveys, \
                                 fingerOfGod, foreGround, beam, \
                                 NewtonianO1)

    # initial guess for model fitting parameters
    mcmcGuess      = initialGuess2(mainPath, cosmologY, k, PData[0])

    # write MCMC initial 'guess' to '.dat' file
    filePath       = f"{mainPath}result/"
    dataFile1      = open(filePath+f"MCMC_guess.dat", 'w')
    for i in range(len(mcmcGuess)):
        dataFile1.write("%e  %s"%(mcmcGuess[i],''))
    dataFile1.close()

    # model fit parameters in parameter space
    #models   = np.array(mcmcGuess) + (1.0e4*np.random.randn(chains, len(mcmcParameters2)))
    step_sizes = np.array([1e4, 2e-3, 0.1, 1.5, 0.5])  # Step sizes for [A, k0, n, m, beta]
    models = np.array(mcmcGuess) + step_sizes * np.random.randn(chains, len(mcmcParameters2))

    
    with mp.Pool() as pool:
        sampler  = emcee.EnsembleSampler(chains, len(mcmcParameters2), logProbability2, pool=pool, \
                                         args=(mainPath, k, PData, priors2))
        sampler.run_mcmc(models, walks, progress=True)

    # Monitor the acceptance rate
    acceptance_fraction = sampler.acceptance_fraction
    print(f"Mean acceptance fraction: {np.mean(acceptance_fraction):.3f}")

    # Use convergence point to filter outliers
    chain_data = sampler.chain  # Shape: (nwalkers, nsteps, ndim)

    # Use the last 10% of the samples to identify convergence
    convergence_data = chain_data[:, int(0.9 * walks):, :]  # Take the last 10% of steps
    median_convergence = np.median(convergence_data, axis=1)  # Median of the final positions
    std_convergence = np.std(convergence_data, axis=1)  # Std deviation of the final positions

    # Define a threshold based on standard deviation
    threshold = 2 * std_convergence  # You can adjust the multiplier as needed

    # Identify good walkers (those that converge near the median)
    good_walkers = np.all(np.abs(chain_data[:, -1, :] - median_convergence) < threshold, axis=1)

    # Ensure some walkers are retained
    if np.sum(good_walkers) < chains / 2:
        print("Warning: Too many walkers filtered out. Retaining all walkers.")
        filtered_chain_data = chain_data  # Keep original data if too many are removed
    else:
        filtered_chain_data = chain_data[good_walkers, :, :]

    # burnIn walks = 40% x walks
    burnIn = int(0.7 * walks)  # Adjusted burn-in

    # discard the burnIn walks and flatten the filtered chain data
    flatSample_filtered = filtered_chain_data[:, burnIn:, :].reshape(-1, filtered_chain_data.shape[2])

    # trace plot for each parameter to check for convergence
    fig, axes = pl.subplots(len(mcmcParameters2), sharex=True)
    axesLabels = [latexDictionary2(mcmcParameters2[i]) for i in range(len(mcmcParameters2))]
    for i in range(len(mcmcParameters2)):
        axes[i].plot(filtered_chain_data[:,:,i].T, color=sns.xkcd_rgb["black"], linestyle='-', linewidth=1.0, alpha=0.3)
        axes[i].set_ylabel(axesLabels[i], fontsize=10)
        axes[i].tick_params(axis='both', labelsize=10)
    axes[-1].set_xlabel(r'$\mathrm{number\;of\;walks}$', fontsize=10)
    pl.tight_layout()
    pl.savefig(f"{mainPath}result/convergence/chains_{specsSurvey1}{specsSurvey2}.pdf", format='pdf', dpi=300)
    pl.show(block=False)
    pl.close('all')

    # Reshape the chain data: (n_walkers, n_steps, n_params) -> (n_chains, n_draws, n_params)
    chain_data = sampler.get_chain(discard=0, flat=False)  # Get the chain without burn-in
    chain_data = np.swapaxes(chain_data, 0, 1)  # Swap the first two dimensions
    
    # Convert to InferenceData for ArviZ
    data = az.convert_to_inference_data(chain_data)
    
    # Run Gelman-Rubin diagnostic test
    r_hat_values = az.rhat(data)
    print(f"Gelman-Rubin R-hat values:\n{r_hat_values}")

    # Check if any R-hat values are significantly greater than 1.1
    if np.any(r_hat_values.to_array().values > 1.1):  # Convert to array to check values
        print("Warning: Some parameters have not yet converged (R-hat > 1.1).")

    # Write the cleaned MCMC 'flatSamples' to the '.dat' file
    dataFile2 = open(filePath + f"MCMC_Model2_{specsSurvey1}{specsSurvey2}_{file_suffix}.dat", 'w')
    for i in range(len(flatSample_filtered)):
        for j in range(len(flatSample_filtered[i])):
            dataFile2.write("%e  %s" % (flatSample_filtered[i][j], ''))
        dataFile2.write("\n")
    dataFile2.close()

    return flatSample_filtered

