# importing packages
import sys, os
import numpy as np
#----------------------------------------------------------------------

def deltaC(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from omegaMatter import omega_m
	#---#

	omega_m_ = omega_m(mainPath, z, cosmologY)

	return ( (3.0*(12.0*np.pi)**(2.0/3.0)) / 20.0 ) * ( 1.0 + (0.0123*np.log10(omega_m_)) )
#----------------------------------------------------------------------

def deltaCprime(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from omegaMatter import omega_m
	from omegaMatterConformalDerivative import omega_mPrime
	#---#

	omega_m_      = omega_m(mainPath, z, cosmologY)
	omega_m_prime = omega_mPrime(mainPath, z, cosmologY)
	deltaCPrime_  = ( (3.0*(12.0*np.pi)**(2.0/3.0)) / 20.0 ) * 0.0123 * ( omega_m_prime/(omega_m_*np.log(10.0)) )

	return deltaCPrime_ 
#----------------------------------------------------------------------