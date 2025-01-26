import functools 
import numpy as np
#----------------------------------------------------------------------

def logPrior(PModel, priors):

	# we assume flat prior probabilities

	# PModel     ---> fitting parameters in our assumed model

	# priors     ---> dictionary that contains prior 
	#			      information on each parameter

	Priors       = np.zeros(len(PModel), dtype=np.float64)
	for i in range(len(PModel)):
		# parameter name
		parameter   = list(priors.keys())[i]
		# accepted minimum and maximum values of parameter 
		miN, maX    = priors[parameter]

		if (PModel[i]>miN and PModel[i]<maX):
			Priors[i] = 0.0
		else:
			Priors[i] = -np.inf

	# log of total priors
	logPriors = np.sum(Priors)

	return logPriors
#----------------------------------------------------------------------