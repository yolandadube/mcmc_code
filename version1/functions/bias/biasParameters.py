# importing packages
import os, sys
import numpy as np
#----------------------------------------------------------------------

def clusteringBias(mainPath, z, cosmologY, survey, specsSurvey):

	b1=b2=bs=sigma=be=Q=None

	sys.path.append(mainPath+"functions/bias/")
	sys.path.append(mainPath+"functions/HITemperature/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from biasHIIM import bHIIM
	from biasHAlpha import bHAlpha
	from biasDesiBGS import bDesiBGS
	from biasDesiELG import bDesiELG
	from biasMegaMapLBG import bMegaMapLBG
	from HIbrightnessTemperature import THI
	from biasSKAOGalBand2 import bSKAOGalBand2
	from conformalHubbleParameter import curlyH, curlyHprime
	#---#

	# fiducial Planck cosmological parameters
	A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, \
	omega_n0, n_s, gamma, w, fnl = cosmologY

	if (survey=="galaxy"):
		
		# H\alpha euclid-like
		if (specsSurvey=="HAlpha"):
			b1, b2, bs, sigma, be, Q = bHAlpha(z)
		# DESI BGS
		elif (specsSurvey=="desiBGS"):
			b1, b2, bs, sigma, be, Q = bDesiBGS(mainPath, z, cosmologY)
		# DESI ELG
		elif (specsSurvey=="desiELG"):
			b1, b2, bs, sigma, be, Q = bDesiELG(mainPath, z, cosmologY)
		# skaOGalBand2
		elif (specsSurvey=="skaOGalBand2"):
			b1, b2, bs, sigma, be, Q = bSKAOGalBand2(z)
		# MegaMapper LBG
		elif (specsSurvey=="megaMapLBG"):
			b1, b2, bs, sigma, be, Q = bMegaMapLBG(z)
			
	if (survey=="HI IM"):
		
		# b1, b2 
		b1, b2 = bHIIM(z)
		
		# bs
		bs = (4.0/7.0)*(1.0-b1)

		# sigma, arXiv: 1906.07032
		sigma0 = 11.0*h
		sigma  =  sigma0 * (1.0+z)**(-1.9) * np.exp(-(z/11.0)**2.0)
		#sigma  = (1.596*h*np.cos((0.833*z) - 7.68) + 1.44) * np.exp(-0.1689*z)
		 
		# evolution bias
		Tb_HI     = THI(z)
		dTb_HI_dz = 2.3242*10.0**(-1.0) - (4.8272*10.0**(-2.0)*z)
		curlH     = curlyH(mainPath, z, cosmologY)
		curlH1    = curlyHprime(mainPath, z, cosmologY)
		term1     = 1.0 + (curlH1/curlH**2.0)
		term2     = -((1.0+z)/Tb_HI)*dTb_HI_dz
		be        = term1+term2

		# magnification bias
		Q = 1.0

	return b1, b2, bs, sigma, be, Q
#----------------------------------------------------------------------
