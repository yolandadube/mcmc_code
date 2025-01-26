# importing packages
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline
#----------------------------------------------------------------------

def shot_noiseSkaOGalBand2(z):

	# arXiv: 1811.02743 (Table 8)
	redshifts = np.array([0.05, 0.15, 0.25, 0.35, 0.45])
	# comoving galaxy number density in h^{3} Mpc^{-3}
	n = np.array([27.3, 4.93, 0.949, 0.223, 0.0644]) * 10.0**(-3.0) 
	nInterp = interpolate.interp1d(redshifts, n)
	nSpline = CubicSpline(redshifts, n)
	# shot noise in h^{-3} Mpc^{3}
	if ((z<np.min(redshifts)) or (z>np.max(redshifts))):
		shot_noisE = 1.0 / nSpline(z)	
	else:
		shot_noisE = 1.0 / nInterp(z)

	return shot_noisE
#----------------------------------------------------------------------