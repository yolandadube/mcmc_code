# importing packages
import sys, os 
#----------------------------------------------------------------------

def omega_m(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from hubbleParameter import H
	#---#

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY
	
	Hubble = H(z, cosmologY)
	# Hubble constant defined at time t = today, i.e. redshift z = 0
	H0     = h / 2997.9 

	# omegaM
	omegaM = omm0 * (1.0+z)**(3.0) * ( (H0/h) / Hubble)**2.0

	return omegaM
#----------------------------------------------------------------------