# importing packages
import sys, os
#----------------------------------------------------------------------

def bNonGaussian(mainPath, z, p, cosmologY, survey, specsSurvey):

	sys.path.append(mainPath+"functions/bias/")
	sys.path.append(mainPath+"functions/criticalDensityParameter/")
	# importing functions
	from criticalDensity import deltaC
	from biasParameters import clusteringBias
	#---#

	# delta_c
	delta_c = deltaC(mainPath, z, cosmologY)

	# gaussian bias 
	b10, b20, bs, sigma, be, Q = clusteringBias(mainPath, z, cosmologY, survey, specsSurvey)

	# non-gaussian bias
	b01 = 2.0*delta_c*(b10-p)

	b11 = 4.0*((delta_c*b20)+((((13.0/21.0)*delta_c)-1.0)*(b10-p))+1.0)

	bn  = 4.0*((delta_c*(p-b10))+1.0)

	b02 = 4.0*delta_c*((delta_c*b20)-((2.0*(((4.0/21.0)*delta_c)+1.0))*(b10-p)))

	return b01, b11, bn, b02
#----------------------------------------------------------------------