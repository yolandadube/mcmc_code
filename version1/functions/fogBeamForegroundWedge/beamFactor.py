# importing packages
import sys, os
import numpy as np
#----------------------------------------------------------------------------------

def DP_beam(mainPath, k, mu, z, cosmologY, survey, specsSurvey, modeSurvey):

	sys.path.append(mainPath+"functions/fieldOfView/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")

	# importing functions
	from fov import thetaFOV
	from comovingDistance import comoving_distance

	# single dish mode
	if (modeSurvey=="single dish"):

		if (survey=="HI IM"):
			# comoving distance in h^{-1} Mpc
			chi     = comoving_distance(z, cosmologY)
			# field-of-view
			theta_b = thetaFOV(mainPath, z, specsSurvey)
			# beam
			k_perp  = k*np.sqrt(1.0-mu**2.0)
			beAm    = np.exp(-(k_perp**2.0*chi**2.0*theta_b**2.0)/(16.0*np.log(2.0)))
		else:
			beAm    = 1.0
	# interferometer mode	
	else:
		beAm = thetaFOV(mainPath, z, specsSurvey)

	return beAm
#------------------------------------------------------------------------------------