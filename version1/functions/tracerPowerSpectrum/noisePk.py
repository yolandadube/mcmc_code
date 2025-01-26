# importing packages
import sys, os
import numpy as np
#----------------------------------------------------------------------

def noisePTracer(mu, k, k_NL0, k_fg, pk, tk, mainPath, z, Nwedge, fingerOfGod, foreGround, beam, dictionary, \
				 cosmologY, surveys, specsSurveys, modeSurveys, NewtonianO1):

	sys.path.append(mainPath+"functions/noise/")
	sys.path.append(mainPath+"functions/systemTemperature/")
	sys.path.append(mainPath+"functions/tracerPowerSpectrum/")
	sys.path.append(mainPath+"functions/surveySpecifications/")
	#---#
	from pkTracer import Ptracer
	from shotNoise import shot_noise
	from systemTemPerature import Tsys
	from specifications import surveySpecs
	from thermalNoise import HIThermalNoise
	
	# surveys
	survey1, survey2           = surveys
	modeSurvey1, modeSurvey2   = modeSurveys
	specsSurvey1, specsSurvey2 = specsSurveys

	# tracer power spectrum
	pTracer = Ptracer(mu, k, k_NL0, k_fg, pk, tk, mainPath, z, Nwedge, fingerOfGod, foreGround, dictionary, \
					  cosmologY, surveys, specsSurveys, modeSurveys, NewtonianO1, beam)
	# single-tracer
	if (specsSurvey1==specsSurvey2):

		# system temperature in Kelvin
		tSys = Tsys(mainPath, z, survey1, specsSurvey1)

		# survey specifications
		N_dish, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff = surveySpecs(specsSurvey1)

		# noise term
		if ( survey1=="HI IM" ):
			# shot noise
			pShotNoise    = 0.0
			# thermal noise
			pThermalNoise = HIThermalNoise(mainPath, mu, k, z, dictionary, cosmologY, survey1, \
											specsSurvey1, modeSurvey1, survey_area_sky, tSys, N_dish, \
											D_dish, t_total, beam)
		else:
			# shot noise
			pShotNoise    = shot_noise(mainPath, z, cosmologY, specsSurvey1)
			# thermal noise
			pThermalNoise = 0.0
	# multi-tracer ---> cross noise is 0
	else:
		# shot noise
		pShotNoise        = 0.0
		# thermal noise
		pThermalNoise     = 0.0

	return pTracer + pShotNoise + pThermalNoise
#----------------------------------------------------------------------