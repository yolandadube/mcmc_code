import numpy as np
import matplotlib.pyplot as plt
from classy import Class

# Define cosmological parameters (Planck 2018 values for Lambda-CDM model)
cosmo_params = {
    'h': 0.674,
    'Omega_b': 0.02237 / 0.674**2,
    'Omega_cdm': 0.1200 / 0.674**2,
    'n_s': 0.965,
    'A_s': 2.1e-9,
    'tau_reio': 0.0544,
    'output': 'mPk',
    'P_k_max_1/Mpc': 10.0,  # Maximum k in units of h/Mpc
    'z_max_pk': 0.0         # Maximum redshift for matter power spectrum (z=0 for present time)
}

# Initialize Class and compute
cosmo = Class()
cosmo.set(cosmo_params)
cosmo.compute()

# Define a range of wavenumbers k in h/Mpc
k = np.logspace(-3, 0, 1000)  # From 10^-3 to 10^0 h/Mpc

# Calculate the matter power spectrum P(k) at redshift z = 0
Pk = np.array([cosmo.pk(kk, 0.0) for kk in k])

# Find the peak (turnover) by identifying the maximum of P(k)
peak_index = np.argmax(Pk)
k_turnover = k[peak_index]
P_turnover = Pk[peak_index]

# Plotting the matter power spectrum
plt.figure(figsize=(15, 12))
plt.loglog(k, Pk, label=r'Matter Power Spectrum $P(k)$', color='b', lw=3.5)

# Shaded region for the galaxy survey probing range (k ~ 0.01 to 0.2 h/Mpc)
plt.fill_between(k, Pk.min(), Pk.max(), where=((k >= 0.01) & (k <= 0.2)), 
                 color='pink', alpha=0.5, label=r'Galaxy Survey Probed Region')

# Shaded region for SKA probing range (k ~ 0.001 to 0.05 h/Mpc)
plt.fill_between(k, Pk.min(), Pk.max(), where=((k >= 0.001) & (k <= 0.05)), 
                 color='lightblue', alpha=0.6, label=r'SKA Probed Region')

# Highlighting the power spectrum turnover at the peak
plt.axvline(x=k_turnover, color='red', linestyle='--', lw=3)

# Labeling the turnover
plt.text(k_turnover * 1.1, 5e3, r'Turnover Scale', color='red', fontsize=32, ha='left')

# LaTeX-formatted labels
plt.xlabel(r'$k \, [h/\mathrm{Mpc}]$', fontsize=35)
plt.ylabel(r'$P(k) \, [(\mathrm{Mpc}/h)^3]$', fontsize=35)

# Increase the font size of tick labels
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

# Grid and legend settings
plt.legend(fontsize=30)

# Save the figure as a high-quality PDF
plt.savefig("matter_power_spectrum.pdf", format="pdf", dpi=300)

# Show the plot
#plt.show()

# Clean up
cosmo.struct_cleanup()
cosmo.empty()
