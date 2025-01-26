import os, sys 
import numpy as np
#----------------------------------------------------------------------

# system temperature in K
def Tsys(mainPath, z, survey, specsSurvey):

	sys.path.append(mainPath+"functions/skyTemperature/")
	sys.path.append(mainPath+"functions/systemTemperature/")

	#---#
	from TSky import T_sky
	from instrumentalTemperature import Tinst
	
	#---#
	if (survey=='HI IM'):
		systemTemperaturE = Tinst(specsSurvey) + T_sky(z)
	else:
		systemTemperaturE = 0.0

	return systemTemperaturE
#----------------------------------------------------------------------