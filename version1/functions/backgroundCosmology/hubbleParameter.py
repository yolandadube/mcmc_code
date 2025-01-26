# importing packages
import numpy as np
#----------------------------------------------------------------------

def H(z, cosmologY):

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY
	
	# omega dark energy
	omega_lambda_0 = 1.0 - omm0
	# Hubble constant defined at time t = today, i.e. redshift z = 0 in h Mpc^{-1}
	H0             = h / 2997.9 
	# Hubble parameter at any redshift z in h Mpc^{-1}
	Hubble         = (H0/h) * np.sqrt( (omm0*(1.0+z)**(3.0)) + omega_lambda_0)
	
	return Hubble
#----------------------------------------------------------------------