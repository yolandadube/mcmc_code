# importing packages
import sys, os 
#----------------------------------------------------------------------

def curlyH(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from scaleFactor import a
	from hubbleParameter import H
	#---#

	return a(z) * H(z, cosmologY)

def curlyHprime(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from omegaMatter import omega_m
	#---#

	curlH = curlyH(mainPath, z, cosmologY) 
	omega_m_ = omega_m(mainPath, z, cosmologY)

	return curlH**2.0 * ( 1.0 - (1.5* omega_m_) ) 

def curlyHdoublePrime(mainPath, z, cosmologY):

	sys.path.append(mainPath+"functions/backgroundCosmology/")
	# importing functions
	from omegaMatter import omega_m
	#---#

	curlH    = curlyH(mainPath, z, cosmologY) 
	omega_m_ = omega_m(mainPath, z, cosmologY)
	
	return  2.0 * curlH**3.0 * ( 1.0 - ((3.0*omega_m_)/4.0) )
#----------------------------------------------------------------------