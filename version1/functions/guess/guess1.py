import os, sys
import numpy as np
from scipy.optimize import curve_fit
#----------------------------------------------------------------------

def initialGuess1(mainPath, cosmologY, k, Pdata):

	'''

	sys.path.append(mainPath+"functions/theory/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	#---#
	from EmpiricalFittingModel import empirical_fitting_model
	from kEquality import k_eq
	#---#

	# equality scale (reference for the turnover scale)
	kEq        = k_eq(mainPath, cosmologY)

	# fit model to P (data) given the inital guess 'p0=[A, k0, n, m, beta]'
	popt, pcov = curve_fit(empirical_fitting_model, k, Pdata, p0=[8.9e+04, kEq, 0.1, 6.5, 0.2], maxfev=20000)
	
	# fitting parameters for model
	A, k0, n, m, beta= popt
	Pmodel              = empirical_fitting_model(k, A, k0, n, m, beta)

	'''

	A = 3e5
	k0 = 0.018
	n = 0.65
	m = 4.0
	beta = 0.52

	return [A, k0, n, m, beta]
#------------------------------------------------------------------------------------------------------