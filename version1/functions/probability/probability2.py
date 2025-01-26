import os, sys
import numpy as np
#----------------------------------------------------------------------

def logProbability2(PModel, mainPath, k, PData, priors):


	# PData      ---> simulated data + error for the power spectrum

	# PModel     ---> fitting parameters in our assumed model

	# priors     ---> dictionary that contains prior 
	#			      information on each parameter

	sys.path.append(mainPath+"functions/prior/")
	sys.path.append(mainPath+"functions/likelihood/")
	#---#
	from prior import logPrior
	from likelihood2 import logLikelihood2

	# log prior
	lp = logPrior(PModel, priors)

	# check if prior is not infinite
	if not np.isfinite(lp):
		return -np.inf

	# log likelihood
	lh = logLikelihood2(mainPath, k, PData, PModel)

	return lp+lh
#----------------------------------------------------------------------