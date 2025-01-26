from vpython import*
import numpy as np
from scipy.linalg import eigh

# Define the system parameters
m1 = 1.0  # Mass of the first mass (kg)
m2 = 1.0  # Mass of the second mass (kg)
k1 = 2.0  # Spring constant for the first spring (N/m)
k2 = 1.0  # Spring constant for the second spring (N/m)

# Mass matrix
M = np.array([[m1, 0],
              [0, m2]])

# Stiffness matrix
K = np.array([[k1 + k2, -k2],
              [-k2, k2]])

# Solve the eigenvalue problem
eigenvalues, eigenvectors = eigh(K, M)

# Natural frequencies (rad/s)
omega = np.sqrt(eigenvalues)

# Time span for simulation
t = np.linspace(0, 20, 1000)

# Simulate the normal modes
def simulate_mode(v, omega, t):
    return np.outer(v, np.cos(omega * t))

# Initialize VPython objects
mass1 = sphere(pos=vector(0, 0, 0), radius=0.1, color=color.red)
mass2 = sphere(pos=vector(0, -1, 0), radius=0.1, color=color.blue)
spring1 = cylinder(pos=vector(0, 1, 0), axis=mass1.pos - vector(0, 1, 0), radius=0.05)
spring2 = cylinder(pos=mass1.pos, axis=mass2.pos - mass1.pos, radius=0.05)

# Choose a mode to visualize (0 or 1)
mode_to_visualize = 0
mode = simulate_mode(eigenvectors[:, mode_to_visualize], omega[mode_to_visualize], t)

# Animation loop
for i in range(len(t)):
    #rate(50)  # Adjust the rate to control the speed of the animation
    mass1.pos.y = mode[0, i]
    mass2.pos.y = mode[1, i]
    spring1.axis = mass1.pos - spring1.pos
    spring2.pos = mass1.pos
    spring2.axis = mass2.pos - mass1.pos
