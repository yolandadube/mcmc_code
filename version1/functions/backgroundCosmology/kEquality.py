# importing packages
import sys, os
import numpy as np
#----------------------------------------------------------------------

def k_eq(mainPath, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from scaleFactor import a
	#---#

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY
	
	# Hubble constant defined at time t = today, i.e. redshift z = 0
	H0   = h / 2997.9 

	# redshift of matter-radiation equality
	z_eq = 3512.881472

	# a_equality
	a_eq = a(z_eq)

	return np.sqrt( (2.0*omm0*(H0/h)**2.0)/a_eq )
#----------------------------------------------------------------------