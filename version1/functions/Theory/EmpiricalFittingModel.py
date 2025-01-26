def empirical_fitting_model(k, A, k0, n, m, beta):
    """
    Empirical Fitting Model for detecting turnover in data.

    Parameters:
    k (float or ndarray): The independent variable, typically a wavenumber or scale.
    A (float): Amplitude of the function.
    k_turn (float): The turnover point where the behavior of the function changes.
    n (float): Exponent for the power-law behavior before the turnover.
    m (float): Controls the sharpness of the transition at the turnover point.
    beta (float): Exponent controlling the slope after the turnover.

    Returns:
    float or ndarray: The value of the empirical fitting model at k.
    """
    return A * k**n / (1 + (k / k0)**m)**beta
