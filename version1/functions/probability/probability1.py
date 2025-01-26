import os, sys
import numpy as np
#----------------------------------------------------------------------

def logProbability1(PModel, mainPath, k, PData, priors):


	# PData      ---> simulated data + error for the power spectrum

	# PModel     ---> fitting parameters in our assumed model

	# priors     ---> dictionary that contains prior 
	#			      information on each parameter

	sys.path.append(mainPath+"functions/prior/")
	sys.path.append(mainPath+"functions/likelihood/")
	#---#
	from prior import logPrior
	from likelihood1 import logLikelihood1

	# log prior
	lp = logPrior(PModel, priors)

	# check if prior is not infinite
	if not np.isfinite(lp):
		return -np.inf

	# log likelihood
	lh = logLikelihood1(mainPath, k, PData, PModel)

	return lp+lh
#----------------------------------------------------------------------