# importing packages
import numpy as np
from scipy import interpolate
#----------------------------------------------------------------------

def shot_noiseMegaMapLBG(z):

	# MegaMapper survey: arXiv: 1903.09208 (Table 1), 2106.09713
	redshifts = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
	# comoving galaxy number density in h^{3} Mpc^{-3}
	n = np.array([25.0, 12.0, 6.0, 3.0, 1.5, 0.8, 0.4]) * 10.0**(-4.0) 
	nInterp = interpolate.interp1d(redshifts, n)
	# shot noise in h^{-3} Mpc^{3}
	shot_noisE = 1.0 / nInterp(z)	

	return shot_noisE
#----------------------------------------------------------------------