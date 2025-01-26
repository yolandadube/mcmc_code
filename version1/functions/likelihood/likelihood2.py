import os, sys
import numpy as np
#----------------------------------------------------------------------

def logLikelihood2(mainPath, k, PData, PModel):

	# PData      ---> simulated data + error for the power spectrum

	# PModel     ---> fitting parameters in our assumed model

	sys.path.append(mainPath+"functions/theory/")
	#---#
	from SmoothlyBrokenPowerModel import smoothly_broken_power_law
	#---#

	# unpack tracerPk + error ---> data
	P, errP  = PData
	# compute the sum of chi square
	A, k0, n, m, Delta  = PModel
	# theoretical model for power spectrum
	Pth      = smoothly_broken_power_law(k, A, k0, n, m, Delta)
	# chi square 
	chiSq    = np.sum( ((P-Pth)/errP)**2.0 )

	# Gaussian log likelihood
	gaussLogLikelihood = -0.5*chiSq

	return gaussLogLikelihood
#----------------------------------------------------------------------