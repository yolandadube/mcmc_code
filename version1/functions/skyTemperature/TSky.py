# sky temperature in K, arXiv: 1810.09572
def T_sky(z):

    # frequency of 21 cm emission in MHz
    nu21 = 1420.0

    T_spl = 3                                                         # spill over
    T_CMB = 2.73                                                      # CMB contributions
    T_gal = 25*((408.0/nu21)*(1+z))**2.75
    T_rx  = 7.5 + 10.0*((nu21/(1e3*(1+z))) - 0.75)**2

    return T_spl + T_CMB + T_gal + T_rx
