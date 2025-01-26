# importing packages
import sys, os
#----------------------------------------------------------------------

def O1NewtonianFourierKernel(mainPath, k, mu, Tk, z, fnl, cosmologY, survey, specsSurvey, NewtonianO1):

	sys.path.append(mainPath+"functions/bias/")
	sys.path.append(mainPath+"functions/pParameter/")
	sys.path.append(mainPath+"functions/curlMParameter/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from P import p
	from curlM import curlyM
	from growthRate import f
	from nonGaussianBias import bNonGaussian
	from biasParameters import clusteringBias
	#---#

	if (NewtonianO1==True):
		# Gaussian
		b1      = clusteringBias(mainPath, z, cosmologY, survey, specsSurvey)[0]
		# linear growth rate
		fgrowth = f(mainPath, z, cosmologY)
		# bias
		term1   = b1
		# rsd
		term2   = fgrowth*mu**2.0
		# nuisance p-parameter
		pTerm   = p(survey)
		# non-Gaussian
		curlM_k = curlyM(mainPath, k, Tk, z, cosmologY)
		b01     = bNonGaussian(mainPath, z, pTerm, cosmologY, survey, specsSurvey)[0]
		term3   = b01*(fnl/curlM_k)
		
	else:
		term1 = 0.0
		term2 = 0.0
		term3 = 0.0

	return term1+term2+term3
#----------------------------------------------------------------------