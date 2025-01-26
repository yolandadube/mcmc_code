# importing packages
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline
#----------------------------------------------------------------------

def shot_noiseHAlpha(z):

	# H-alpha stage IV spectroscopic survey: arXiv: 1807.07076
	redshifts = np.array([0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])
	# comoving galaxy number density in h^{3} Mpc^{-3}
	n       = np.array([1.60, 1.29, 1.04, 0.84, 0.68, 0.55, 0.45, 0.37, 0.30, 0.24]) * 10.0**(-3.0) 
	nInterp = interpolate.interp1d(redshifts, n)
	nSpline = CubicSpline(redshifts, n)
	# shot noise in h^{-3} Mpc^{3}
	if ((z<np.min(redshifts)) or (z>np.max(redshifts))):
		shot_noisE = 1.0 / nSpline(z)
	else:
		shot_noisE = 1.0 / nInterp(z)	

	return shot_noisE
#----------------------------------------------------------------------