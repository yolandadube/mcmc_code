# importing packages
from colossus.cosmology import cosmology
#----------------------------------------------------------------------

# growth factor 
def D(z, cosmologY):

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY

	# parameters dictionary for cosmology
	parameters = {'flat': True, 'H0': H0, 'Om0': omm0, 'Ob0': omb0, 'sigma8': sigma80, 'ns': n_s}
	cosmology.addCosmology('myCosmology', **parameters)
	cosmo      = cosmology.setCosmology('myCosmology')

	return cosmo.growthFactor(z)
#----------------------------------------------------------------------