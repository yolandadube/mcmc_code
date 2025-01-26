# importing packages
import numpy as np
#----------------------------------------------------------------------

def bSKAOGalBand2(z):

	# arXiv: 1811.02743 (Table 8)
	
	# b1
	b1 = 0.625 + (0.550*z) + (0.462*z**2.0)
	
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