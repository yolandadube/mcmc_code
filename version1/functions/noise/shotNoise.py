# importing packages
import os, sys
#----------------------------------------------------------------------

def shot_noise(mainPath, z, cosmologY, specsSurvey):

	sys.path.append(mainPath+"functions/noise")
	# importing functions
	from shotNoiseHIIM import shot_noiseHIIM
	from shotNoiseHAlpha import shot_noiseHAlpha
	from shotNoiseDesiBGS import shot_noiseDesiBGS
	from shotNoiseDesiELG import shot_noiseDesiELG
	from shotNoiseMegaMapLBG import shot_noiseMegaMapLBG
	from shotNoiseSKAOGalBand2 import shot_noiseSkaOGalBand2
	#---#

	# H\alpha euclid-like
	if (specsSurvey=="HAlpha"):
		shot_noisE = shot_noiseHAlpha(z)
	# DESI BGS
	elif (specsSurvey=="desiBGS"):
		shot_noisE = shot_noiseDesiBGS(mainPath, z, cosmologY)
	# DESI ELG
	elif (specsSurvey=="desiELG"):
		shot_noisE = shot_noiseDesiELG(mainPath, z, cosmologY)
	# skaOGalBand2
	elif (specsSurvey=="skaOGalBand2"):
		shot_noisE = shot_noiseSkaOGalBand2(z)
	# MegaMapper LBG
	elif (specsSurvey=="megaMapLBG"):
		shot_noisE = shot_noiseMegaMapLBG(z)
	else:
		# HI IM surveys
		shot_noisE = shot_noiseHIIM(z)
		
	return shot_noisE
#----------------------------------------------------------------------