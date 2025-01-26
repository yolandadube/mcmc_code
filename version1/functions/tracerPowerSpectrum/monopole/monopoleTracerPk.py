# importing packages
import sys, os
import numpy as np
from scipy import integrate

#----------------------------------------------------------------------

def P0(k, pk, tk, k_NL0, k_fg, Nwedge, z, mainPath, dictionary, cosmologY, surveys, specsSurveys, modeSurveys, \
	   fingerOfGod, foreGround, NewtonianO1, beam):
	
	sys.path.append(mainPath+"functions/tracerPowerSpectrum/")
	#---#
	from pkTracer import Ptracer
	#--#

	# precision of integration 
	precision = 1.49e-08

	# mu-range [-1, 1] for integration 
	mu        = np.linspace(-0.99999, 0.99999, 10000, endpoint=True)

	monopole  = np.zeros(len(k), dtype=np.float64)
	# compute the monopole
	for i in range(len(k)):
		
		# trapezoidal rule
		pTracer = Ptracer(mu, k[i], k_NL0, k_fg, pk[i], tk[i], mainPath, z, Nwedge, fingerOfGod, \
						  foreGround, dictionary, cosmologY, surveys, specsSurveys, modeSurveys, \
						  NewtonianO1, beam)
		
		p0      = np.trapz(pTracer, x=mu)
		

		'''
		# quadrature integration
		p0 = integrate.quad(Ptracer, -1.0, 1.0, \
							args=(k[i], k_NL0, k_fg, pk[i], tk[i], mainPath, z, Nwedge, fingerOfGod, \
								  foreGround, dictionary, cosmologY, surveys, specsSurveys, \
								  modeSurveys, NewtonianO1), \
							epsabs=precision, epsrel=precision, limit=55)[0]
		'''

		'''
		# romberg integration
		p0  = integrate.romberg(Ptracer, -1.0, 1.0, \
								args=(k[i], k_NL0, k_fg, pk[i], tk[i], mainPath, z, Nwedge, fingerOfGod, \
								  	  foreGround, dictionary, cosmologY, surveys, specsSurveys, \
								  	  modeSurveys, NewtonianO1), \
								tol=precision, rtol=precision, divmax=15)
		'''
		monopole[i] = 0.5*p0
	
	return monopole
#--------------------------------------------------------------------------------------------------------