import numpy as np
from scipy import interpolate
#--------------------------------------------------------------------------------------------------------------------------

# matter power spectrum from CLASS
def pk(pkFile, k):

	# open pk data file
	dataFile = np.loadtxt(pkFile)
	# wavenumber K [h Mpc^{-1}]
	K        = dataFile[:,0]
	# matter power spectrum Pk [h^{-3} Mpc^{3}]
	PK       = dataFile[:,1]
	# interpolate PK
	Pk       = interpolate.CubicSpline(K, PK)

	return Pk(k)
#--------------------------------------------------------------------------------------------------------------------------