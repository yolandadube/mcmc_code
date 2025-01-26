import os, sys
import numpy as np
#----------------------------------------------------------------------

def logLikelihood1(mainPath, k, PData, PModel):

	# PData      ---> simulated data + error for the power spectrum

	# PModel     ---> fitting parameters in our assumed model

	sys.path.append(mainPath+"functions/theory/")
	#---#
	from EmpiricalFittingModel import empirical_fitting_model
	#---#

	# unpack tracerPk + error ---> data
	P, errP  = PData
	# compute the sum of chi square
	A, k0, n, m, beta  = PModel
	# theoretical model for power spectrum
	Pth      = empirical_fitting_model(k, A, k0, n, m, beta)
	# chi square 
	chiSq    = np.sum( ((P-Pth)/errP)**2.0 )

	# Gaussian log likelihood
	gaussLogLikelihood = -0.5*chiSq

	return gaussLogLikelihood
#----------------------------------------------------------------------