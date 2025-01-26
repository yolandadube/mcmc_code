def planckValues(dictionary):

	# Planck 2018 (1807.06209) Table 2 (cosmological parameters)
	# Planck 2018 (1905.05697) (primordial non-Gaussianity)

	# dimensionless hubble parameter today
	H0          = dictionary['H0']
	h           = dictionary['h'] 
	# omega baryonic today
	omb0        = dictionary['omb0']
	# omega cold dark matter today
	omcdm0      = dictionary['omcdm0']
	# Omega matter today
	omm0        = omb0 + omcdm0
	# omega curvature today
	omega_k0    = dictionary['omega_k0']
	# omega neutrino today
	omega_n0    = dictionary['omega_n0']
	# spectral index of primordial power spectrum
	n_s         = dictionary['n_s']
	# A_s
	A_s         = dictionary['A_s']
	# sigma_8
	sigma80     = dictionary['sigma80']
	# linear growth index
	gamma       = dictionary['gamma']
	# equation of state parameter for dark energy
	w           = dictionary['w']
	# primordial non-Gaussianity
	fnl         = dictionary['fnl']

	return A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, omega_n0, n_s, gamma, w, fnl
#-----------------------------------------------------------------------------------------