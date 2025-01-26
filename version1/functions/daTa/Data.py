import os, sys
import numpy as np
#----------------------------------------------------------------------

# arXiv: 2202.13828 
def dataPTracer(mainPath, cosmologY, gauge, nonLinear, pkFile, tkFile, k, k_NL0, k_fg, Nwedge, \
				z, delta_z, dictionary, surveys, specsSurveys, modeSurveys, fingerOfGod, \
				foreGround, beam, NewtonianO1):

	sys.path.append(mainPath+"functions/nModes/")
	sys.path.append(mainPath+"functions/theory/")
	sys.path.append(mainPath+"functions/transferFunction/")
	sys.path.append(mainPath+"functions/matterPowerSpectrum/")
	sys.path.append(mainPath+"functions/tracerPowerSpectrum/monopole/")
	#---#
	from matterPk import pk
	from transfer import tk
	from monopoleTracerPk import P0
	from numberOfModes import Nmodes
	from monopoleNoisePk import noiseP0

	# compute matter transfer function at z=0
	Tk             = tk(tkFile, k)

	# get matter pk at z from CLASS
	Pk             = pk(pkFile, k)

	# compute the monopole of the tracer power spectrum
	pTracer0       =  P0(k, Pk, Tk, k_NL0, k_fg, Nwedge, z, mainPath, dictionary, cosmologY, \
						 surveys, specsSurveys, modeSurveys, fingerOfGod, foreGround, NewtonianO1, beam)

	# compute the monopole of the tracer power spectrum with noise
	pNoiseTracer0  =  noiseP0(k, k_NL0, k_fg, Pk, Tk, mainPath, z, Nwedge, fingerOfGod, \
							  foreGround, beam, dictionary, cosmologY, surveys, \
							  specsSurveys, modeSurveys, NewtonianO1)

	# number of modes in each k-bin
	N_modes        = Nmodes(mainPath, k, z, delta_z, cosmologY, surveys, specsSurveys, \
							modeSurveys)
	# compute error in pTracer
	errPTracer0    = pNoiseTracer0/np.sqrt(N_modes)
		
	return pTracer0, errPTracer0
#------------------------------------------------------------------------------------------------------------

