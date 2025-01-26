import os, sys
import numpy as np
from scipy.optimize import curve_fit
#----------------------------------------------------------------------

def initialGuess(mainPath, cosmologY, k, Pdata):

	sys.path.append(mainPath+"functions/theory/")
	sys.path.append(mainPath+"functions/backgroundCosmology/")
	#---#
	from theory import model
	from kEquality import k_eq
	#---#

	# equality scale (reference for the turnover scale)
	kEq        = k_eq(mainPath, cosmologY)
	print('k_Equality=', kEq)

	# fit model to P (data) given the inital guess 'p0=[P0, k0, alpha, beta]'
	popt, pcov = curve_fit(model, k, Pdata, p0=[3.0e+04, kEq, 1.75, 1.85], maxfev=20000)
	
	# fitting parameters for model
	P0, k0, alpha, beta= popt
	Pmodel              = model(k, P0, k0, alpha, beta)

	return [P0, k0, alpha, beta]
#----------------------------------------------------------------------