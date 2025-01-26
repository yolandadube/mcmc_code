# importing packages
from colossus.cosmology import cosmology
import numpy as np

def comoving_distance(z, cosmologY):
    '''
    # Check if z is a single-element array and convert to scalar
    if isinstance(z, np.ndarray) and z.size == 1:
        z = z.item()
    elif z is None or not isinstance(z, (int, float)):
        raise TypeError("Invalid input: 'z' must be a numerical value.")
        
    '''

    # fiducial Planck cosmological parameters
    A_s, sigma80, H0, h, omm0, omb0, omcdm0, omega_k0, omega_n0, n_s, gamma, w, fnl = cosmologY

    # parameters dictionary for cosmology
    parameters = {'flat': True, 'H0': H0, 'Om0': omm0, 'Ob0': omb0, 'sigma8': sigma80, 'ns': n_s}
    cosmology.addCosmology('myCosmology', params=parameters)
    cosmo = cosmology.setCosmology('myCosmology')

    # comoving distance in h^{-1} Mpc at redshift z
    chi = cosmo.comovingDistance(z_min=0.0, z_max=z, transverse=False)

    return chi
#------------------------------------------------------------------------------------------------
