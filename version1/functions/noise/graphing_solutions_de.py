import numpy as np
import matplotlib.pyplot as plt

# Define the equation as a function
def plot_equation(A_values, x_range):
    x = np.linspace(x_range[0], x_range[1], 400)

    # Create a plot for each value of A
    for A in A_values:
        y = np.tan(np.log(A * (x - 1) / x))

        # Plot the function
        plt.plot(x, y, label=f'A = {A}')

    # Adding labels and legend
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Plot of y = tan(ln(A(x-1)/x))')
    plt.legend()

    # Show the plot
    plt.show()

# Values of A to plot
A_values = [0.5, 1, 1.5, 2]

# Range of x values to consider
# Avoid x = 0 to prevent division by zero and x = 1 to avoid log(0)
x_range = [1.1, 3]

# Call the function to plot the equation for different A values
plot_equation(A_values, x_range)
