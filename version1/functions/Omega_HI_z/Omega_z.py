import numpy as mp 

def Omega_HI(z):
	Omega = 0.00048 + 0.00039*z - 0.000065*z**2
	return Omega