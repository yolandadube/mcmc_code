import os, sys
import numpy as np
#----------------------------------------------------------------------

def logLikelihood(mainPath, k, PData, PModel):

	# PData      ---> simulated data + error for the power spectrum

	# PModel     ---> fitting parameters in our assumed model

	sys.path.append(mainPath+"functions/theory/")
	#---#
	from theory import model
	#---#

	# unpack tracerPk + error ---> data
	P, errP  = PData
	# compute the sum of chi square
	P0, k0, alpha, beta  = PModel
	# theoretical model for power spectrum
	Pth      = model(k, P0, k0, alpha, beta)
	# chi square 
	chiSq    = np.sum( ((P-Pth)/errP)**2.0 )

	# Gaussian log likelihood
	gaussLogLikelihood = -0.5*chiSq

	return gaussLogLikelihood
#----------------------------------------------------------------------