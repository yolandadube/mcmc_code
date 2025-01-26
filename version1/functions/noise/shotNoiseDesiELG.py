# importing packages
import os, sys
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline 
#----------------------------------------------------------------------

# arXiv: 1611.00036
def shot_noiseDesiELG(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from comovingDistance import comoving_distance
	from conformalHubbleParameter import curlyH, curlyHprime

	# number density in deg^{-2}
	redshifts    = np.array([0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, \
							1.55, 1.65])
	dNg_dz       = np.array([309.0, 2269.0, 1923.0, 2094.0, 1441.0, 1353.0, 1337.0, \
							523.0, 466.0, 329.0, 126.0]) 
	dNg_dzInterp = interpolate.interp1d(redshifts, dNg_dz)
	# conversion from per sq degrees to per steradians
	convert     = (180.0/np.pi)**2.0
	if ((z<np.min(redshifts)) or (z>np.max(redshifts))):
		Ng_bar  = convert*dNg_dzSpline(z)
	else:	
		Ng_bar  = convert*dNg_dzInterp(z)

	# comoving distance in h^{-1} Mpc
	chi = comoving_distance(z, cosmologY)

	# conformal hubble parameter in h Mpc^{-1}
	curlH = curlyH(mainPath, z, cosmologY)

	# arXiv: 2107.14057
	# ng in h^{3} Mpc^{-3}
	ng = ((1.0+z)*curlH*Ng_bar) / chi**2.0
	
	# shot noise
	shot_noisE = 1.0/ng

	return shot_noisE
#----------------------------------------------------------------------




