import os, sys
import numpy as np
from scipy.optimize import curve_fit
#----------------------------------------------------------------------

def initialGuess2(mainPath, cosmologY, k, Pdata):
	'''
	sys.path.append(mainPath+"functions/theory/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	#---#
	from SmoothlyBrokenPowerModel import smoothly_broken_power_law
	from kEquality import k_eq
	#---#

	# equality scale (reference for the turnover scale)
	kEq        = k_eq(mainPath, cosmologY)

	# fit model to P (data) given the inital guess 'p0=[A, k0, n, m, Delta]'
	popt, pcov = curve_fit(smoothly_broken_power_law, k, Pdata, p0=[3.0e+04, kEq, 1.75, 1.85, 1.5], maxfev=20000)
	
	# fitting parameters for model
	A, k0, n, m, Delta= popt
	Pmodel              = smoothly_broken_power_law(k, A, k0, n, m, Delta)
	'''
	A  = 4.7e4
	k0 = 0.019
	n  = 0.65
	m  = -0.65
	Delta = 2.2

	return [A, k0, n, m, Delta]
	
#----------------------------------------------------------------------