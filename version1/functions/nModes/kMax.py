import numpy as np
#--------------------------------------------------------------------------------------------------------------------------

def k_Max(z, k_NL0, cosmologY):

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY

	# non-linear cut-off, k_max [h Mpc^{-1}]
	k_max = k_NL0 * (1.0+z)**(2.0/(2.0+n_s))

	return k_max
#--------------------------------------------------------------------------------------------------------------------------