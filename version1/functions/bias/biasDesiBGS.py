# importing packages
import os, sys
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline
#----------------------------------------------------------------------

# arXiv: 2004.12981
def bDesiBGS(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	from growthFactor import D
	#---#

	# b1
	b1 = 1.34 / D(z, cosmologY)
	#b1 = 0.99 + (0.73*z) - (1.29*z**2.0) + (10.21*z**3.0)

	# b2
	b2 = 0.0

	# bs
	bs = (4.0/7.0)*(1.0-b1)

	# FoG parameter
	sigma = 0.0

	# arXiv: 2107.14057
	# be
	be = -2.25 -(4.02*z) + (0.318*z**2.0) - (14.6*z**3.0)

	# Q
	Q = 0.282 + (2.36*z) + (2.27*z**2.0) + (11.1*z**3.0) 

	return b1, b2, bs, sigma, be, Q
#----------------------------------------------------------------------