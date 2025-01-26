# importing packages
import os, sys 
#----------------------------------------------------------------------

def thetaFOV(mainPath, z, specsSurvey):

	sys.path.append(mainPath+"functions/surveySpecifications/")
	# importing functions
	from specifications import surveySpecs
	#---#

	# frequency of 21 cm in Hz
	frequency_21 = 1420.0 * 10.0**(6.0)
	# speed of light in m s^{-1}
	c            = 2997.9 * 10.0**(5.0)
	# wavelength of 21 cm in m
	lambda_21    = c / frequency_21
	# wavelength at z
	lambdA = lambda_21*(1.0+z)
	# survey specifications for dish diameter
	D_dish = surveySpecs(specsSurvey)[2]

	# field-of-view
	theta_FOV = 1.22*(lambdA/D_dish)

	return theta_FOV
#----------------------------------------------------------------------