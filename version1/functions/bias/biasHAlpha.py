# importing packages
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline
#----------------------------------------------------------------------

def bHAlpha(z):

	# H-alpha stage IV spectroscopic survey: arXiv: 1807.07076
	redshifts = np.array([0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])

	# b1
	b1       = np.array([1.33, 1.40, 1.47, 1.54, 1.61, 1.68, 1.75, 1.82, 1.89, 1.96])
	b1Interp = interpolate.interp1d(redshifts, b1)
	b1Spline = CubicSpline(redshifts, b1)
	
	# b2
	b2       = np.array([-0.749, -0.737, -0.721, -0.703, -0.682, -0.658, -0.631, -0.600, -0.567, -0.531])
	b2Interp = interpolate.interp1d(redshifts, b2)
	b2Spline = CubicSpline(redshifts, b2)
	
	# bs
	bs       = np.array([-0.189, -0.229, -0.269, -0.309, -0.349, -0.389, -0.429, -0.469, -0.509, -0.549])
	bsInterp = interpolate.interp1d(redshifts, bs)
	bsSpline = CubicSpline(redshifts, bs)
	
	# sigma in h^{-1} Mpc
	sigma       = np.array([4.62, 4.51, 4.39, 4.27, 4.15, 4.03, 3.92, 3.81, 3.70, 3.61]) 
	sigmaInterp = interpolate.interp1d(redshifts, sigma)
	sigmaSpline = CubicSpline(redshifts, sigma)
	
	# evolution bias
	be       = np.array([-6.13, -5.94, -5.75, -5.54, -5.34, -5.14, -4.94, -4.74, -4.54, -4.35])
	beInterp = interpolate.interp1d(redshifts, be)
	beSpline = CubicSpline(redshifts, be) 
	
	# maginification bias
	Q       = np.array([1.97, 2.07, 2.17, 2.26, 2.33, 2.41, 2.47, 2.52, 2.57, 2.61])
	QInterp = interpolate.interp1d(redshifts, Q)
	QSpline = CubicSpline(redshifts, Q)
	
	if ((z<np.min(redshifts)) or (z>np.max(redshifts))):
		return b1Spline(z), b2Spline(z), bsSpline(z), sigmaSpline(z), beSpline(z), QSpline(z)
	else:
		return b1Interp(z), b2Interp(z), bsInterp(z), sigmaInterp(z), beInterp(z), QInterp(z)
#----------------------------------------------------------------------