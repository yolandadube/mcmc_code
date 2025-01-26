# importing packages
import sys, os 
#----------------------------------------------------------------------

# omega_m' = domega_m / d_eta
def omega_mPrime(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from omegaMatter import omega_m
	from conformalHubbleParameter import curlyH, curlyHprime
	#---#

	omegaM = omega_m(mainPath, z, cosmologY)
	curlH  = curlyH(mainPath, z, cosmologY) 
	curlH1 = curlyHprime(mainPath, z, cosmologY)

	return -omegaM * ( curlH + ((2.0*curlH1)/curlH) )
#----------------------------------------------------------------------