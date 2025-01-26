import glob
import emcee
import os, sys
import numpy as np
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

def mcmcAlgorithm(mainPath, mcmcParameters, priors, chains, walks, cosmologY, gauge, nonLinear, k, k_NL0, k_fg, \
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
    from guess import initialGuess
    from latexFonts import latexDictionary
    from probability import logProbability
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
    mcmcGuess      = initialGuess(mainPath, cosmologY, k, PData[0])

    # write MCMC initial 'guess' to '.dat' file
    filePath       = f"{mainPath}result/"
    dataFile1      = open(filePath+f"MCMC_guess.dat", 'w')
    for i in range(len(mcmcGuess)):
        dataFile1.write("%e  %s"%(mcmcGuess[i],''))
    dataFile1.close()

    # model fit parameters in parameter space
    models   = np.array(mcmcGuess) + (1.0e-04*np.random.randn(chains, len(mcmcParameters)))
    
    with mp.Pool() as pool:
        sampler  = emcee.EnsembleSampler(chains, len(mcmcParameters), logProbability, pool=pool, \
                                         args=(mainPath, k, PData, priors))
        sampler.run_mcmc(models, walks, progress=True)

    # Monitor the acceptance rate
    acceptance_fraction = sampler.acceptance_fraction
    print(f"Mean acceptance fraction: {np.mean(acceptance_fraction):.3f}")

    # trace plot for each parameter to check for convergence
    fig, axes      = pl.subplots(len(mcmcParameters), sharex=True)
    axesLabels     = [ latexDictionary(mcmcParameters[i]) for i in range(len(mcmcParameters)) ]
    for i in range(len(mcmcParameters)):
        axes[i].plot(sampler.chain[:,:,i].T, color=sns.xkcd_rgb["black"], linestyle='-', linewidth=1.0, alpha=0.3)
        axes[i].set_ylabel(axesLabels[i], fontsize=10)
        axes[i].tick_params(axis='both', labelsize=10)
    axes[-1].set_xlabel(r'$\mathrm{number\;of\;walks}$', fontsize=10)

    pl.tight_layout()
    pl.savefig(f"{mainPath}result/convergence/chains_{specsSurvey1}{specsSurvey2}.pdf", format='pdf', dpi=300)
    pl.show(block=False)
    pl.close('all')

    # burnIn walks = 30% x walks
    burnIn         = int(0.3*walks)

    # discard the burnIn walks
    flatSample     = sampler.get_chain(discard=burnIn, thin=15, flat=True)
    
    # write MCMC 'flatSamples' to '.dat' file with k_fg value appended to the file name
    dataFile2      = open(filePath+f"MCMC_ParabolicModel_{specsSurvey1}{specsSurvey2}_{file_suffix}.dat", 'w')
    
    for i in range(len(flatSample)):
        for j in range(len(flatSample[i])):
            dataFile2.write("%e  %s"%(flatSample[i][j],''))
        dataFile2.write("\n")
    dataFile2.close()

    return flatSample
