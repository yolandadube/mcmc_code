import sys, os
import numpy as np
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.font_manager
import matplotlib.pyplot as pl
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
from getdist import plots, MCSamples
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from matplotlib.ticker import LogLocator
from scipy.interpolate import CubicSpline
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
pl.rc('xtick', labelsize = 15)
pl.rc('ytick', labelsize = 15)

pl.rcParams['ps.useafm']          = True
pl.rcParams['pdf.use14corefonts'] = True
pl.rcParams['text.usetex']        = True

#pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pl.tick_params(axis = 'both', which = 'both', direction = 'in')
pl.rcParams['axes.linewidth'] = 1.4 #set the value globally
pl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
#--------------------------------------------------------------------------------------------------------------------------

# main path directory
mainPath = "/home/yolanda/Documents/mcmc_code/version1/"
#--------------------------------------------------------------------------------------------------------------------------

sys.path.append(mainPath+"functions/latex/")
from latexFonts import latexDictionary
#--------------------------------------------------------------------------------------------------------------------------

# MCMC parameters list
mcmcParameters  = ['p0', 'k0', 'alpha', 'beta', 'A', 'f_NL']
#--------------------------------------------------------------------------------------------------------------------------

filePath    = f"{mainPath}result/contour/run1/"
filePath1    = f"{mainPath}result/DatafilesToBePlotted/"
#---#

# load initial guess data file for MCMC parameters
guessFile   = np.loadtxt(filePath1+'MCMC_guess.dat')

# load MCMC data file
file1 = np.loadtxt(filePath1+'MCMC_skaOIMBand1skaOIMBand1_fnl2.dat')
file2 = np.loadtxt(filePath1+'MCMC_skaOIMBand1skaOIMBand1_fnl5.dat')
file3 = np.loadtxt(filePath1+'MCMC_skaOIMBand1skaOIMBand1_fnl10.dat')
#file4 = np.loadtxt(filePath1+'MCMC_skaOGalBand2skaoGalBand2_0.001.dat')
#file5 = np.loadtxt(filePath1+'MCMC_megaMapLBGmegaMapLBG_0.001.dat')
#--------------------------------------------------------------

# axis labels for triangular plot
axesNames   = ["x%s"%i for i in range(len(mcmcParameters))]
axesLabels  = [ latexDictionary(mcmcParameters[i]) for i in range(len(mcmcParameters)) ]

# convert MCMC data files to GetDist MCSamples objects
samplesMCMC1 = MCSamples(samples=file1, names=axesNames, labels=axesLabels)
samplesMCMC2 = MCSamples(samples=file2, names=axesNames, labels=axesLabels)
samplesMCMC3 = MCSamples(samples=file3, names=axesNames, labels=axesLabels)
#samplesMCMC4 = MCSamples(samples=file4, names=axesNames, labels=axesLabels)
#samplesMCMC5 = MCSamples(samples=file4, names=axesNames, labels=axesLabels)
#--------------------------------------------------------------

# plotting
#g = plots.get_subplot_plotter(rc_sizes=True, width_inch=None)
g = plots.get_single_plotter(rc_sizes=True, width_inch=None)

# customizing subplots
g.settings.scaling                    = True
g.settings.linewidth                  = 3.0
g.settings.line_styles                = ['-', '-',]  # ['-', ...]
#g.settings.line_labels                = False # remove legend
g.settings.prob_y_ticks               = True
g.settings.tight_layout               = True
g.settings.axes_fontsize              = 15
g.settings.axes_labelsize             = 15
g.settings.legend_fontsize            = 15
g.settings.norm_prob_label            = 15
g.settings.norm_1d_density            = True
g.settings.num_plot_contours          = 2
g.settings.figure_legend_loc          = "upper right"
g.settings.subplot_size_ratio         = 1
g.settings.figure_legend_frame        = True
g.settings.axis_tick_max_labels       = 3
g.settings.alpha_factor_contour_lines = 1

# colours, line styles, labels
# cases      = [samplesMCMC1]
# colors     = [sns.xkcd_rgb["blue"]]
# linestyles = ['-']
# labelS     = [r'$\mathrm{SKAO\;Band}\;1$']

cases      = [samplesMCMC1, samplesMCMC2, samplesMCMC3]
colors = [sns.xkcd_rgb["orange"], sns.xkcd_rgb["purple"], sns.xkcd_rgb["green"]]
linestyles = ['-', '-', '-']
#labelS = ["MeerKLASS L-band", "MeerKLASS UHF-Band", "SKA Mid-Band 1"]
#labelS = [r'$\text{desiBGS}$', r'$\text{desiELG}$']
#labelS     = [r'$\mathrm{SKAO\;Band}\;1$', r'$\mathrm{HIRAX\;256}$']
#labelS = [r'$\mathrm{MeerKLASS\ L\text{-}band}$', 
          #r'$\mathrm{MeerKLASS\ UHF\text{-}Band}$', 
          #r'$\mathrm{SKA\ Mid\text{-}Band\ 1}$']
labelS = [
    "MeerKLASS L-band, $f_{\mathrm{NL}} = 10$",
    "MeerKLASS UHF-Band, $f_{\mathrm{NL}} = 10$",
    "SKA Mid-Band 1, $f_{\mathrm{NL}} = 10$"
]



'''
g.triangle_plot(cases, ['x0', 'x1', 'x2', 'x3', 'x4', 'x5'], filled_compare=True, \
				filled=False, contour_colors=colors, \
				line_args=[{'ls':linestyles[i], 'color':colors[i]} \
							for i in range(len(cases))], \
				legend_labels=labelS, alpha=1.0)
'''
g.triangle_plot(cases, filled_compare=True, filled=False, contour_colors=colors, \
				line_args=[{'ls':linestyles[i], 'color':colors[i]} \
							for i in range(len(cases))], \
				legend_labels=labelS, alpha=1.0)

# k0_fiducial value
k0_fiducial = 0.0163                                                             # Taken from Cunnington

# Loop over potential subplot locations where k0 might be on the x-axis
for i in range(len(mcmcParameters)):                                             # Assuming k0 is in the 'k0' position in mcmcParameters
    for j in range(len(mcmcParameters)):
        # Check if the subplot at this location exists
        if g.subplots[j, i] is not None:
            # Check if this subplot has k0 on the x-axis
            if i == mcmcParameters.index('k0'):
                # Add the vertical line for k0 fiducial value
                g.subplots[j, i].axvline(x=k0_fiducial, color='grey', linestyle='--', linewidth=1)
            # If k0 is on the y-axis of the subplot, and you also want a horizontal line, you can add:
            # if j == mcmcParameters.index('k0'):
            #     g.subplots[j, i].axhline(y=k0_fiducial, color='grey', linestyle='--', linewidth=1)


'''
# customize subplots
for i in range(len(mcmcParameters)):
	for j in range(i,len(mcmcParameters)):
		# axes 
		ax = g.subplots[j,i]
		if (i==j):
			# insert subplot label 
			#ax.set_title(f'{latexDictionary(mcmcParameters[i])}', fontsize = 15)
			# insert vertical line
			ax.vlines(guessFile[i], -1.0e+06, 1.0e+06, color=sns.xkcd_rgb["grey"], \
				 	  linestyle='--', linewidth=1.0, alpha=1.0)
		else:
			# insert vertical line
			ax.vlines(guessFile[i], -1.0e+06, 1.0e+06, color=sns.xkcd_rgb["grey"], \
			 	  	  linestyle='--', linewidth=1.0, alpha=1.0)
			# insert horizontal line
			ax.hlines(guessFile[j], -1.0e+06, 1.0e+06, color=sns.xkcd_rgb["grey"], \
			 	  	  linestyle='--', linewidth=1.0, alpha=1.0)
                         
'''

# saving as pdf
g.export(f'{filePath}MCMCcontour10.pdf')
#--------------------------------------------------------------------------------------------------------------------------
