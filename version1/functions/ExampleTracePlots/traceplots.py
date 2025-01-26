import numpy as np
import matplotlib.pyplot as plt

# Define the target distribution (1D Gaussian)
def target_distribution(x):
    return np.exp(-0.5 * x**2)

# Metropolis-Hastings MCMC sampler
def metropolis_hastings_1d(num_samples, initial_position, proposal_width, target_dist):
    samples = np.zeros(num_samples)
    samples[0] = initial_position

    for i in range(1, num_samples):
        current_position = samples[i - 1]
        proposed_position = current_position + proposal_width * np.random.randn()
        
        # Calculate acceptance probability
        acceptance_prob = min(1, target_dist(proposed_position) / target_dist(current_position))
        
        if np.random.rand() < acceptance_prob:
            samples[i] = proposed_position
        else:
            samples[i] = current_position

    return samples

# Parameters for the MCMC
num_samples = 20000
initial_position = 0

# Case 1: Converging Trace Plot (Appropriate Proposal Width)
proposal_width_converging = 0.5
samples_converging = metropolis_hastings_1d(num_samples, initial_position, proposal_width_converging, target_distribution)

# Case 2: Poorly Mixing / Diverging Trace Plot (Too Small Proposal Width)
proposal_width_poor_mixing = 0.05  # Too small, causing the sampler to move slowly
samples_poor_mixing = metropolis_hastings_1d(num_samples, initial_position, proposal_width_poor_mixing, target_distribution)

# Case 3: Chaotic Diverging Trace Plot (Very Large Proposal Width)
proposal_width_diverging = 75.0  # Extremely large, causing chaotic behavior
samples_diverging = metropolis_hastings_1d(num_samples, initial_position, proposal_width_diverging, target_distribution)

# Create the figure and subplots
fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=True)  # Share y-axis

# Set font sizes and line widths for better visibility
plt.rc('font', family='serif', size=20)
plt.rc('axes', linewidth=2)

# Adjust spine, tick thickness, and length
spine_thickness = 3.5  # Set the desired spine thickness
tick_size = 20
tick_width = 3.5
tick_length = 5  # Adjust the length of the ticks

for ax in axs:
    # Adjust the spines
    for spine in ax.spines.values():
        spine.set_linewidth(spine_thickness)
    
    # Adjust the tick parameters
    ax.tick_params(axis='both', which='major', labelsize=tick_size, width=tick_width, length=tick_length)

# Plot the converging trace plot
axs[0].plot(samples_converging, color='blue')
axs[0].set_title('Converging Trace Plot', fontsize=22)
axs[0].set_xlabel('Iteration', fontsize=20)
axs[0].set_ylabel('Sample Value', fontsize=20)

# Plot the poorly mixing trace plot
axs[1].plot(samples_poor_mixing, color='orange')
axs[1].set_title('Poorly Mixing Trace Plot', fontsize=22)
axs[1].set_xlabel('Iteration', fontsize=20)

# Plot the chaotic diverging trace plot
axs[2].plot(samples_diverging, color='red')
axs[2].set_title('Chaotic Diverging Trace Plot', fontsize=22)
axs[2].set_xlabel('Iteration', fontsize=20)

# Tight layout for better spacing
plt.tight_layout()

# Save the figure as a high-quality PDF
plt.savefig('/home/yolanda/Documents/mcmc_code/version1/functions/ExampleTracePlots/mcmc_trace_plots_high_quality.pdf', format='pdf', bbox_inches='tight', dpi=300)

# Show the plot
plt.show()
