import os, sys
#--------------------------------------------------------------------------------------------------------------------------

# request power spectra and transfer functions
def spectra(mainPath, cosmologY, gauge, nonLinear, z):

	sys.path.append(mainPath+"functions/pyCLASSWrapper/")
	#---#
	from runCLASS import pyCLASS
	#---#

	# compute pk and tk from CLASS
	outputs = pyCLASS(mainPath, cosmologY, gauge, nonLinear, z)
	
	return None
#--------------------------------------------------------------------------------------------------------------------------