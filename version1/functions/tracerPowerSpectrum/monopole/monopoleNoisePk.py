# importing packages
import sys, os
import numpy as np
from scipy import integrate
#----------------------------------------------------------------------

def noiseP0(k, k_NL0, k_fg, pk, tk, mainPath, z, Nwedge, fingerOfGod, foreGround, beam, dictionary, cosmologY, \
			surveys, specsSurveys, modeSurveys, NewtonianO1):
	
	sys.path.append(mainPath+"functions/tracerPowerSpectrum/")
	#---#
	from noisePk import noisePTracer
	#--#

	# precision of integration 
	precision = 1.49e-08

	# mu-range [-1, 1] for integration 
	mu        = np.linspace(-0.99999, 0.99999, 1000, endpoint=True)

	monopole  = np.zeros(len(k), dtype=np.float64)
	# compute the monopole
	for i in range(len(k)):
		
		# trapezoidal rule
		pNoise  = noisePTracer(mu, k_NL0, k[i], k_fg, pk[i], tk[i], mainPath, z, Nwedge, fingerOfGod, \
							   foreGround, beam, dictionary, cosmologY, surveys, specsSurveys, \
							   modeSurveys, NewtonianO1)

		p0      = np.trapz(pNoise, x=mu)
		

		'''
		# quadrature integration
		p0 = integrate.quad(noisePTracer, -1.0, 1.0, \
							args=(k[i], k_NL0, k_fg, pk[i], tk[i], mainPath, z, Nwedge, fingerOfGod, \
							   	  foreGround, beam, dictionary, cosmologY, surveys, specsSurveys, \
							   	  modeSurveys, NewtonianO1), \
							epsabs=precision, epsrel=precision, limit=55)[0]
		'''
		
		'''
		# romberg integration
		p0  = integrate.romberg(noisePTracer, -1.0, 1.0, \
								args=(k[i], k_NL0, k_fg, pk[i], tk[i], mainPath, z, Nwedge, fingerOfGod, \
							   		  foreGround, beam, dictionary, cosmologY, surveys, specsSurveys, \
							   		  modeSurveys, NewtonianO1), \
								tol=precision, rtol=precision, divmax=15)
		'''
		monopole[i] = 0.5*p0

	return monopole
#--------------------------------------------------------------------------------------------------------