# python package
import os, sys
import numpy as np
#-------------------------------------------------------------------------------------------

def Ptracer(mu, k, k_NL0, k_fg, pk, tk, mainPath, z, Nwedge, fingerOfGod, foreGround, dictionary, cosmologY, \
			surveys, specsSurveys, modeSurveys, NewtonianO1, beam):

	sys.path.append(mainPath+"functions/kernel/")
	sys.path.append(mainPath+"functions/HITemperature/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	sys.path.append(mainPath+"functions/fogBeamForegroundWedge/")

	from wedgeTerm import wEdge
	from beamFactor import DP_beam
	from fingersOfGod import DP_FOG
	from foreground import foreGroundRemoval
	from comovingDistance import comoving_distance
	from O1NewtonianKernel import O1NewtonianFourierKernel

	# surveys
	survey1, survey2           = surveys
	modeSurvey1, modeSurvey2   = modeSurveys
	specsSurvey1, specsSurvey2 = specsSurveys

	# Fourier kernel
	O1NSurvey1  = O1NewtonianFourierKernel(mainPath, k, mu, tk, z, dictionary['fnl'], cosmologY, \
											survey1, specsSurvey1, NewtonianO1)
	O1NSurvey2  = O1NewtonianFourierKernel(mainPath, k, mu, tk, z, dictionary['fnl'], cosmologY, \
											survey2, specsSurvey2, NewtonianO1)
	kernelSq    = O1NSurvey1*np.conjugate(O1NSurvey2)

	# finger-of-god effect
	if (fingerOfGod==True):
		DP_fog  = DP_FOG(mainPath, k, mu, z, cosmologY, surveys, specsSurveys)
	else:
		DP_fog  = 1.0

	if (foreGround==True):
		# apply foreground check
		fgCheck1    = foreGroundRemoval(mu, k, k_fg, survey1)
		fgCheck2    = foreGroundRemoval(mu, k, k_fg, survey2)
		fgCheck     = fgCheck1*fgCheck2
	else:
		fgCheck     = 1.0	

	if (beam==True):
		Db1 = DP_beam(mainPath, k, mu, z, cosmologY, survey1, specsSurvey1, modeSurvey1)
		Db2 = DP_beam(mainPath, k, mu, z, cosmologY, survey2, specsSurvey2, modeSurvey2)
		Db  = (Db1) * (Db2)  # Fixed line
	else:
		Db  = 1.0



	# case for interferometer - apply wedge check
	wedgeCheck  = wEdge(mainPath, cosmologY, mu, k, z, surveys, modeSurveys, \
						specsSurveys, Nwedge)
	
	# tracer power spectrum
	PT          = DP_fog*Db*fgCheck*wedgeCheck*kernelSq*pk

	return PT
#-------------------------------------------------------------------------------------------