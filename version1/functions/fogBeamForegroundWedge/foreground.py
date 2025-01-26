# importing packages
import numpy as np
#----------------------------------------------------------------------

def foreGroundRemoval(mu, k, k_fg, survey):

	# foreground in case of "HI IM"
	if (survey=="HI IM"):
		# k_parallel
		k_parallel = mu*k
		'''
		# heaviside function
		if (np.abs(k_parallel)>k_fg):
			checkFg = 1.0
		else:
			checkFg = 0.0
		'''
		# exponential function
		checkFg = 1.0 - np.exp( -(k_parallel/k_fg)**2.0 )
	# no foreground in case of "galaxy"
	else:
		checkFg = 1.0

	return checkFg
#----------------------------------------------------------------------