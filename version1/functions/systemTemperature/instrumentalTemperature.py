import os, sys 
import numpy as np
#----------------------------------------------------------------------

# intrumental temperature in K
def Tinst(specsSurvey):

	# arXiv: 1704.01941
	if (specsSurvey=="skaOIMBand1"):
		TInstrument = 25.0

	if (specsSurvey=="skaOIMBand2"):
		TInstrument = 25.0

	if (specsSurvey=="meerKATLBand"):
		TInstrument = 25.0

	if (specsSurvey=="meerKATUHFBand"):
		TInstrument = 25.0

	if (specsSurvey=="hirax256"):
		TInstrument = 50.0

	if (specsSurvey=="hirax1024"):
		TInstrument = 50.0

	if (specsSurvey=="puma5k"):
		# arXiv: 1810.09572
		eta_c       = 0.9
		eta_s       = 0.9
		Tampl       = 50.0  # K
		Tground     = 300.0 # K
		TInstrument = (Tampl/(eta_c*eta_s)) + ((Tground*(1.0-eta_s))/eta_s)
		#TInstrument = 50.0

	if (specsSurvey=="puma32k"):
		# arXiv: 1810.09572
		eta_c       = 0.9
		eta_s       = 0.9
		Tampl       = 50.0  # K
		Tground     = 300.0 # K
		TInstrument = (Tampl/(eta_c*eta_s)) + ((Tground*(1.0-eta_s))/eta_s)
		#TInstrument = 50.0

	if (specsSurvey=="dsa2000"):
		TInstrument = 25.0

	return TInstrument
#----------------------------------------------------------------------