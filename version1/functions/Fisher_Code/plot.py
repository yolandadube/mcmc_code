import warnings
warnings.filterwarnings('ignore')
#-------------------#
import matplotlib
matplotlib.use('Agg')
#-------------------#

import os, sys
import numpy as np
import matplotlib.cm as cm
import matplotlib.font_manager
import matplotlib.pyplot as pl
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
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
#---------------------------------------------------------------------------------------------------------------------

fig = pl.figure(figsize = (15.0, 12.0))

pl.rcParams['xtick.major.size']  = 3
pl.rcParams['xtick.major.width'] = 2
pl.rcParams['ytick.major.size']  = 3
pl.rcParams['ytick.major.width'] = 2
pl.rcParams['xtick.major.pad']= 6
pl.rcParams['ytick.major.pad']= 6
# ...
pl.rcParams['xtick.minor.size']  = 3
pl.rcParams['xtick.minor.width'] = .3
pl.rcParams['ytick.minor.size']  = 3
pl.rcParams['ytick.minor.width'] = .3

# ...
pl.rcParams['figure.dpi'] = 80
pl.rcParams['lines.solid_joinstyle'] = 'miter'  
pl.rcParams['lines.antialiased'] = True
pl.rc('xtick', labelsize = 40.0)
pl.rc('ytick', labelsize = 40.0)

pl.rcParams['ps.useafm'] = True
pl.rcParams['pdf.use14corefonts'] = True
pl.rcParams['text.usetex'] = True

#pl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pl.tick_params(axis = 'both', which = 'both', direction = 'in')
pl.rcParams['axes.linewidth'] = 1.4 #set the value globally
#--------------------------------------------------------------------------------------------------------------------------

# code main path
mainPath = "/home/yolanda/Documents/mcmc_code/version1/"
#--------------------------------------------------------------------------------------------------------------------------

# dictionary fiducial Planck cosmological parameters
dictionary    = {'H0'       :67.66,     # hubble parameter today
                 'omb0'          :0.02242,   # is the baryonic matter density parameter
                 'h'                  :0.6766 ,   # 
                 'omm0'         :0.3130,    # is the total matter density parameter
                 'omcdm0'     :0.2640,    # is the cold dark matter density parameter
                 'omega_k0'   :0.0,       # omega curvature today
                 'omega_n0'   :0.0,       # omega neutrino today
                 'n_s'               :0.9665,    # spectral index of primordial power spectrum
                 'A_s'               :2.105e-09, # A_s 
                 'sigma80'      :0.8102,    # sigma_8 today
                 'gamma'        :0.55,      # linear growth index
                 'w'                  :-1.0,      # equation of state parameter for dark energy
                 'fnl'                 :0.0        # primordial non-Gaussianity parameter
                }
#--------------------------------------------------------------------------------------------------------------------------

# do contour plots

sys.path.append(mainPath+"functions/Fisher_Code/")
from FisherPlotContour import plot_contours
from FisherPlotContourGalaxy import plot_contours_galaxy
#---#
conTours = plot_contours_galaxy(mainPath, dictionary)
#--------------------------------------------------------------------------------------------------------------------------