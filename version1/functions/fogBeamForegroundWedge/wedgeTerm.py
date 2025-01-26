# importing packages
import sys, os
import numpy as np
#----------------------------------------------------------------------

def wEdge(mainPath, cosmologY, mu, k, z, surveys, modeSurveys, specsSurveys, Nwedge):

	sys.path.append(mainPath+"functions/fieldOfView/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from fov import thetaFOV
	from hubbleParameter import H
	from comovingDistance import comoving_distance
	#---#

	# surveys
	survey1, survey2           = surveys
	modeSurvey1, modeSurvey2   = modeSurveys
	specsSurvey1, specsSurvey2 = specsSurveys

	# comoving distance in h^{-1} Mpc
	chi        = comoving_distance(z, cosmologY)
	# Hubble parameter in h Mpc^{-1}
	H_         = H(z, cosmologY)
	# magnitude of k_parallel component
	# k_parallel
	k_parallel = mu*k
	# magnitude of k_perp component
	# k perp
	k_perp     = k*np.sqrt(1.0-mu**2.0)

	if (survey1=="galaxy" and survey2=="galaxy"):
		checkWedge = 1.0
	elif ((survey1=="HI IM" and survey2=="HI IM") and
			(modeSurvey1=="single dish" and modeSurvey2=="single dish")):
		checkWedge = 1.0
	elif ((survey1=="galaxy" and survey2=="HI IM") and
			(modeSurvey1=="single dish" and modeSurvey2=="single dish")):
		checkWedge = 1.0
	elif ((survey1=="HI IM" and survey2=="galaxy") and
		(modeSurvey1=="single dish" and modeSurvey2=="single dish")):
		checkWedge = 1.0
	else:
		if ((survey1=="galaxy" and survey2=="HI IM") and
			(modeSurvey1=="single dish" and modeSurvey2=="interferometer")):
			# field-or-view
			FOV   = thetaFOV(mainPath, z, specsSurvey2)
		elif ((survey1=="HI IM" and survey2=="galaxy") and
				(modeSurvey1=="interferometer" and modeSurvey2=="single dish")):
			# field-or-view
			FOV   = thetaFOV(mainPath, z, specsSurvey1)
		elif ((survey1=="HI IM" and survey2=="HI IM") and 
				(modeSurvey1=="interferometer" and modeSurvey2=="interferometer")):
			# field-or-view
			FOV   = thetaFOV(mainPath, z, specsSurvey1)
			
		# k_wedge
		kWedge     = ((chi*H_)/(1.0+z))*np.sin(1.22*Nwedge*(FOV/2.0))
		checkWedge = np.heaviside((np.abs(k_parallel)-(kWedge*k_perp)), 0.0)

	return checkWedge
#----------------------------------------------------------------------