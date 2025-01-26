import numpy as np
import matplotlib.pyplot as plt

# Define the function P(k)
def P_k(k, k0, P0, alpha, beta):
    x = (np.log10(k) - np.log10(k0)) / np.log10(k0)
    P = np.where(k < k0, P0**(1 - alpha * x**2), P0**(1 - beta * x**2))
    return P

# Parameters
k0 = 0.0165    # Reference scale
P0 = 37300    # Amplitude of power spectrum
k = np.logspace(-3, 0, 1000)  # k values

# Different pairs of alpha and beta, with increasing values
alpha_beta_pairs = [(0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]

# Plotting
plt.figure(figsize=(15, 12))

for alpha, beta in alpha_beta_pairs:
    P = P_k(k, k0, P0, alpha, beta)
    plt.plot(k, P, label=fr'$\alpha={{{alpha}}}, \ \beta={{{beta}}}$', linewidth=3.5)  # Increased line width

# Add LaTeX equation in the middle of the plot
plt.text(0.01, 5e3, r'$P(k) = P_0^{1-\alpha x^2} \ \text{for} \ k<k_0$' + "\n" +
         r'$P(k) = P_0^{1-\beta x^2} \ \text{for} \ k \geq k_0$', 
         fontsize=32, ha='center', va='center')

# Use LaTeX for all text in the plot with proper math typesetting
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$k$', fontsize=35)
plt.ylabel(r'$P(k)$', fontsize=35)

# Increase font size for the legend and set title size
plt.legend(fontsize=32, title_fontsize=22)

# Optionally re-enable grid for better readability
#plt.grid(True, which="both", ls="--")

# Increase the size of tick labels
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

# Save the figure as a high-quality PDF
plt.savefig("Fitting_Model_plot.pdf", format="pdf", dpi=300)

plt.show()
