# importing packages
import sys, os
import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as pl
#----------------------------------------------------------------------

# arXiv: 1911.03964
def nbL(mainPath, mu, k, z, cosmologY, specsSurvey):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	sys.path.append(mainPath+"functions/surveySpecifications/")
	# importing functions
	from specifications import surveySpecs
	from comovingDistance import comoving_distance

	# comoving distance in h^{-1} Mpc
	chi   = comoving_distance(z, cosmologY)
	# survey specification for dish diameter in m
	D_dish = surveySpecs(specsSurvey)[1]
	
	# frequency of 21 cm in Hz
	frequency_21 = 1420.0 * 10.0**(6.0)
	# speed of light in m s^{-1}
	c            = 2997.9 * 10.0**(5.0)
	# wavelength of 21 cm in m
	lambda_21    = c / frequency_21
	# redshift 21 cm lambda
	lambda_      = lambda_21*(1.0+z)
	# k_perp
	k_perp       = k*np.sqrt(1.0-mu**2.0)
	# l = u*lambda, u = (k_perp*chi)/2*pi
	l            = (chi*lambda_*k_perp)/(2.0*np.pi)

	if (specsSurvey == 'hirax256'):
		a  = 0.4847
		b  = -0.3300
		c  = 1.3157
		d  = 1.5974
		e  = 6.8390
		Ns = 16.0
		# antenna distribution 
		L = Ns*D_dish
		nb_term1 = (Ns/D_dish)**2.0
		nb_term2 = a + (b*(l/L))
		nb_term3 = 1.0 + (c*(l/L)**d)
		nb_term4 = np.exp( -(l/L)**e)
		nb_l     = nb_term1*(nb_term2/nb_term3)*nb_term4
	#---#
	if (specsSurvey == 'hirax1024'):
		a  = 0.4847
		b  = -0.3300
		c  = 1.3157
		d  = 1.5974
		e  = 6.8390
		Ns = 32.0
		# antenna distribution 
		L = Ns*D_dish
		nb_term1 = (Ns/D_dish)**2.0
		nb_term2 = a + (b*(l/L))
		nb_term3 = 1.0 + (c*(l/L)**d)
		nb_term4 = np.exp( -(l/L)**e)
		nb_l     = nb_term1*(nb_term2/nb_term3)*nb_term4
	#---#
	# if(specsSurvey=="chime"):
	# 	# antena distribution parameters
	# 	A = 48.511
	# 	B = 60.693
	# 	C = 2.4797
	# 	# antenna distribution
	# 	nb_l = A * np.exp( -(l/B)**C )
	#---#
	if (specsSurvey == 'puma5k'):
		a  = 0.5698
		b  = -0.5274
		c  = 0.8358
		d  = 1.6635
		e  = 7.3177
		Ns = 100.0
		# antenna distribution 
		L = Ns*D_dish
		nb_term1 = (Ns/D_dish)**2.0
		nb_term2 = a + (b*(l/L))
		nb_term3 = 1.0 + (c*(l/L)**d)
		nb_term4 = np.exp( -(l/L)**e)
		nb_l     = nb_term1*(nb_term2/nb_term3)*nb_term4
	#---#
	if (specsSurvey == 'puma32k'):
		a  = 0.5698
		b  = -0.5274
		c  = 0.8358
		d  = 1.6635
		e  = 7.3177
		Ns = 253.0
		# antenna distribution 
		L = Ns*D_dish
		nb_term1 = (Ns/D_dish)**2.0
		nb_term2 = a + (b*(l/L))
		nb_term3 = 1.0 + (c*(l/L)**d)
		nb_term4 = np.exp( -(l/L)**e)
		nb_l     = nb_term1*(nb_term2/nb_term3)*nb_term4
	#---#
	if (specsSurvey == 'dsa2000'):
		# use data file arXiv: 2311.00896
		nbData   = np.loadtxt(mainPath+"functions/noise/DSA_2000baselinedistribuion.dat")
		# baseline length
		bL       = nbData[:,0]
		# baseline count
		nb       = nbData[:,1]
		# spline interpolation
		nbSpline = CubicSpline(bL, nb)
		nb_l     = nbSpline(l) / lambda_**2.0
		
	return nb_l
#----------------------------------------------------------------------