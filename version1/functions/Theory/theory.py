import numpy as np

def model(k, P0, k0, alpha, beta):
    """
    Fit for the monopole of the power spectrum around the turnover scale k0.
    
    Parameters:
    k (array-like): Wavenumbers.
    P0 (float): Monopole amplitude around turnover for the fit.
    k0 (float): Turnover scale.
    alpha (float): Parameter that regulates the shape of the parabola for k < k0.
    beta (float): Parameter that regulates the shape of the parabola for k >= k0.
    
    Returns:
    fit (array-like): The fitted monopole power spectrum.
    """
    # Logarithmic scaling
    x = np.log(k) / np.log(k0) - 1.0
    
    # Fitting formula
    fit = np.where(
        k < k0,
        P0 ** (1.0 - (alpha * x ** 2.0)),
        P0 ** (1.0 - (beta * x ** 2.0))
    )
    
    return fit
