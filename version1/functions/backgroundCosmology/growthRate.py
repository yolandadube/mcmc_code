# importing packages
import sys, os
#----------------------------------------------------------------------

# linear growth rate
def f(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from omegaMatter import omega_m
	#---#

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY

	omega_m_ = omega_m(mainPath, z, cosmologY)

	return omega_m_**gamma
#----------------------------------------------------------------------

# f' = df / f_eta
def fprime(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from omegaMatter import omega_m
	from omegaMatterConformalDerivative import omega_mPrime
	#---#

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY

	omega_m_      = omega_m(mainPath, z, cosmologY)
	omega_m_prime = omega_mPrime(mainPath, z, cosmologY)
	
	return gamma*omega_m_**(gamma-1.0)*omega_m_prime
#----------------------------------------------------------------------