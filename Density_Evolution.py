import numpy as np
from scipy.special import iv
import Parameters as pm

# Constants
M_disk_initial = pm.M_p    # Initial mass of the disk (mass of black hole)
R_disk_initial = pm.R_disk # Initial radius of the disk 
viscosity = pm.viscosity_disk   # Viscosity 

def sigma(t, R):
    """Compute the surface density of the disk at radius R and time t."""
    x = R / R_disk_initial
    tau = t * viscosity / R_disk_initial**2

    # Avoid division by zero with a small epsilon
    epsilon = 1e-10
    tau = max(tau, epsilon)
    
    return (M_disk_initial / (np.pi * R_disk_initial**2)) * tau**-1 * x**-0.25 * np.exp(-(1 + x**2) / tau) * iv(0.25, 2*x/tau)

# Time array
t = np.linspace(0, 1, 100)   # Time (arbitrary units)

# Radial array
R = np.linspace(0, 2*R_disk_initial, 100)   # Radial distance from the center (from 0 to twice the initial disk radius)

# Initial condition for surface density
sigma_initial = np.zeros_like(R)
for i, r in enumerate(R):
    sigma_initial[i] = sigma(0, r)

# Time step for finite difference method
dt = t[1] - t[0]
dr = R[1] - R[0]

# Create a matrix to hold the surface density at each time and radius
sigma_matrix = np.zeros((len(t), len(R)))
sigma_matrix[0, :] = sigma_initial

# Apply the explicit finite difference method
for n in range(len(t)-1):
    # Handle the interior points
    sigma_matrix[n+1, 1:-1] = sigma_matrix[n, 1:-1] + 3 * dt * (sigma_matrix[n, 2:] - 2*sigma_matrix[n, 1:-1] + sigma_matrix[n, :-2]) / dr**2 / R[1:-1]
    
    # Handle the boundaries (using first-order forward and backward differences)
    sigma_matrix[n+1, 0] = sigma_matrix[n, 0] + 3 * dt * (-2*sigma_matrix[n, 0] + 2*sigma_matrix[n, 1]) / dr**2 / R[0]
    sigma_matrix[n+1, -1] = sigma_matrix[n, -1] + 3 * dt * (-2*sigma_matrix[n, -1] + 2*sigma_matrix[n, -2]) / dr**2 / R[-1]
