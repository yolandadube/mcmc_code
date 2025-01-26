# importing packages
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline
#----------------------------------------------------------------------

def bMegaMapLBG(z):

	# MegaMapper survey: arXiv: 1903.09208 (Table 1), 2106.09713 (eq 2.3)
	redshifts = np.array([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])

	# b1
	b1       = np.array([2.5, 3.3, 4.1, 4.9, 5.8, 6.6, 7.4])
	b1Interp = interpolate.interp1d(redshifts, b1)
	b1Spline = CubicSpline(redshifts, b1)
	
	# b2
	b2       = 0.0
	b2Interp = 0.0
	b2Spline = 0.0
	
	# bs
	bs       = 0.0
	bsInterp = 0.0
	bsSpline = 0.0
	
	# sigma in h^{-1} Mpc
	sigma       = 0.0
	sigmaInterp = 0.0
	sigmaSpline = 0.0
	
	# evolution bias
	be       = 0.0
	beInterp = 0.0
	beSpline = 0.0
	
	# maginification bias
	Q       = 0.0
	QInterp = 0.0
	QSpline = 0.0
	
	if ((z<np.min(redshifts)) or (z>np.max(redshifts))):
		return b1Spline(z), b2Spline, bsSpline, sigmaSpline, beSpline, QSpline
	else:
		return b1Interp(z), b2Interp, bsInterp, sigmaInterp, beInterp, QInterp
#----------------------------------------------------------------------