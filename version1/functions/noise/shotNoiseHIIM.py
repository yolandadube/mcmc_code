# importing packages
import numpy as np
from scipy import interpolate
#----------------------------------------------------------------------

def shot_noiseHIIM(z):

	# shot noise for HI galaxies
	
	# shot noise in h^{-3} Mpc^{3} from Table 5 of arXiv: 1804.09180
	redshifts   = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
	Pshot       = np.array([104.0, 124.0, 65.0, 39.0, 14.0, 7.0])
	
	PshotSpline = interpolate.CubicSpline(redshifts, Pshot)
	# shot noise in h^{-3} Mpc^{3}
	shot_noisE  = PshotSpline(z)	

	return shot_noisE
#----------------------------------------------------------------------