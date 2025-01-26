import numpy as np
from scipy import interpolate
#--------------------------------------------------------------------------------------------------------------------------

# transfer function at z=0 from CLASS
def tk(tkFile, k):

	# open pk data file
	dataFile = np.loadtxt(tkFile)
	# wavenumber K [h Mpc^{-1}]
	K        = dataFile[:,0]
	# Phi transfer function 
	TKPhi    = dataFile[:,6]
	# interpolate TK
	TkPhi    = interpolate.CubicSpline(K, TKPhi)
	# normalize Tk to unity on large scale
	Tk       = TkPhi(k)/TKPhi[0]

	return Tk
#--------------------------------------------------------------------------------------------------------------------------