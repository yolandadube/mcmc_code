import numpy as np
import matplotlib.pyplot as plt

def plot_model_from_file(file_path):
    # Import the model function from the module
    from theory import model

    file_path = "/home/yolanda/Documents/mcmc_code/version1/result/MCMC_guess.dat"
    
    # Load parameters from file
    P0, k0, alpha, beta, A, f_NL = np.loadtxt(file_path)
    
    # Define k values
    k = np.logspace(-3, -1, 100) 
    
    # Call the model function
    P_model = model(k, P0, k0, alpha, beta, A, f_NL)
    
    # Plot the model
    plt.figure(figsize=(8, 6))
    plt.plot(k, P_model, label='Model')
    plt.xlabel('k')
    plt.ylabel('P(k)')
    plt.title('Model Power Spectrum')
    plt.legend()
    plt.show()

# Example usage
file_path = 'Guess_values.dat'
plot_model_from_file(file_path)
