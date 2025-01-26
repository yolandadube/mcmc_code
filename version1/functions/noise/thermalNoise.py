# importing packages
import sys, os
import numpy as np
#--------------------------------------------------------------------------------------------------

def HIThermalNoise(mainPath, mu, k, z, dictionary, cosmologY, survey, specsSurvey, modeSurvey, \
					surveyArea, systemTemperature, N_dish, D_dish, t_total, beam):

	sys.path.append(mainPath+"functions/noise/")
	sys.path.append(mainPath+"functions/planckBestFit/")
	sys.path.append(mainPath+"functions/HITemperature/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	sys.path.append(mainPath+"functions/fogBeamForegroundWedge/")
	# importing functions
	from hubbleParameter import H
	from beamFactor import DP_beam
	from HIbrightnessTemperature import THI
	from antennaBaselineDistribution import nbL
	from conformalHubbleParameter import curlyH
	from comovingDistance import comoving_distance
	from planckCosmologicalParameters import planckValues

	# comoving distance in h^{-1} Mpc
	chi               = comoving_distance(z, cosmologY)
	# conformal hubble parameter h Mpc^{-1}
	curlH             = curlyH(mainPath, z, cosmologY)
	# total area of sky in square degrees
	total_area_of_sky = 4.0*np.pi*(180.0/np.pi)**2.0
	# fraction of sky covered
	f_sky             = surveyArea/total_area_of_sky 
	# background brightness temperature in K
	Tb                = THI(z) * 10.0**(-3.0)
	# number of dishes
	N_Dish            = N_dish
	# total time of observation in seconds
	t_Total           = t_total * 3600.0
	# frequency of 21 cm in MHz
	frequency_21_MHz  = 1420.0
	# frequency of 21 cm in Hz
	frequency_21      = 1420.0 * 10.0**(6.0)
	# speed of light in m s^{-1}
	c                 = 2997.9 * 10.0**(5.0)
	# wavelength of 21 cm in m
	lambda_21         = c / frequency_21
	# wavelength at z
	lambdA            = lambda_21*(1.0+z)
	# effective beam area
	Ad                = (np.pi/4.0)*D_dish**2.0
	Ae                = 0.7*Ad

	# alpha term
	if (modeSurvey=="single dish"):
		alpha    = 1.0 / N_Dish
	elif (modeSurvey=="interferometer"):
		# nb
		nbTerm   = lambdA**2.0*nbL(mainPath, mu, k, z, cosmologY, specsSurvey)
		alpha    = (lambdA**2.0/Ae)**2.0 * (1.0/nbTerm)

	# beam term
	if (beam==True):
		beta     = DP_beam(mainPath, k, mu, z, cosmologY, survey, specsSurvey, modeSurvey)
	else:
		beta     = 1.0

	# Pnoise
	Pnoise_term1 = (2.0*np.pi*f_sky)/(N_Dish*t_Total)
	Pnoise_term2 = (chi**2.0)
	Pnoise_term3 = (systemTemperature)**2.0
	Pnoise_term4 = lambdA
	Pnoise_term5 = (1+z)/H(z, cosmologY)

	PnOise       = Pnoise_term1*Pnoise_term2*Pnoise_term3*Pnoise_term4*Pnoise_term5
	
	return PnOise
#--------------------------------------------------------------------------------------------------