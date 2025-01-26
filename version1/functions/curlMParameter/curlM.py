# importing packages
import sys, os
#----------------------------------------------------------------------

def curlyM(mainPath, k, Tk, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from growthRate import f
	from growthFactor import D
	from omegaMatter import omega_m
	#---#

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY

	# hubble constant defined at time t = today, i.e. redshfit z = 0
	H0       = h / 2997.9
	# growth factor
	D_       = D(z, cosmologY)
	# linear growth rate
	fgrowth  = f(mainPath, z, cosmologY)
	# omega_m
	omega_m_ = omega_m(mainPath, z, cosmologY)
	# g_in
	g_in     = (3.0/5.0) * (1.0+z) * D_ * ( 1.0 + ((2.0*fgrowth)/(3.0*omega_m_)) )

	return (2.0 / 3.0) * (1.0/g_in) * ( (D_* Tk) / omm0 ) * ( k / (H0/h) )**2.0 
#----------------------------------------------------------------------