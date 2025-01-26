# importing packages
import sys, os
import numpy as np
#----------------------------------------------------------------------

# arXiv: 9311057
def DP_FOG(mainPath, k, mu, z, cosmologY, surveys, specsSurveys):

	sys.path.append(mainPath+"functions/bias/")
	# importing functions
	from biasParameters import clusteringBias
	#---#
	
	# surveys
	survey1, survey2           = surveys
	specsSurvey1, specsSurvey2 = specsSurveys

	# single-tracer
	if(specsSurvey1==specsSurvey2):
		# damping parameter sigma_P
		sigma          = clusteringBias(mainPath, z, cosmologY, survey1, specsSurvey1)[3]
		dampingTermFoG = np.exp( -0.5*k**2.0*mu**2.0*sigma**2.0 )
	# multi-tracer
	else:
		# damping parameter sigma_P
		sigma1         = clusteringBias(mainPath, z, cosmologY, survey1, specsSurvey1)[3]
		sigma2         = clusteringBias(mainPath, z, cosmologY, survey2, specsSurvey2)[3] 
		dampingTermFoG = np.exp( -0.5*k**2.0*mu**2.0*(sigma1**2.0+sigma2**2.0) )

	return dampingTermFoG
#----------------------------------------------------------------------