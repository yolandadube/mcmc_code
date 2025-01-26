def smoothly_broken_power_law(k, A, k0, n, m, Delta):
    """
    Smoothly Broken Power Law Model for detecting turnover in data.

    Parameters:
    k (float or ndarray): The independent variable, typically a wavenumber or scale.
    A (float): Amplitude of the function.
    k_turn (float): The turnover point where the behavior of the function changes.
    n1 (float): Exponent for the power-law behavior before the turnover.
    n2 (float): Exponent for the power-law behavior after the turnover.
    Delta (float): Controls the smoothness of the transition at the turnover point.

    Returns:
    float or ndarray: The value of the smoothly broken power law model at k.
    """
    return A * (k / k0)**n * (1 + (k / k0)**Delta)**((m - n) / Delta)
