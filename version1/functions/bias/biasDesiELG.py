# importing packages
import sys, os
import numpy as np
#----------------------------------------------------------------------

# arXiv: 2004.12981
def bDesiELG(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	from growthFactor import D
	#---#

	# b1
	b1 = 0.84 / D(z, cosmologY)

	# b2
	b2 = 0.0

	# bs
	bs = (4.0/7.0)*(1.0-b1)

	# FoG parameter
	sigma = 0.0

	# be
	be = 0.0

	# Q
	Q = 0.0 

	return b1, b2, bs, sigma, be, Q
#----------------------------------------------------------------------