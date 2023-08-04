import numpy as np
from Parameters import G, c, sigma_T, M_p, epsilon, R_disk, viscosity_disk, omega_disk, L_dot

# Eddington Luminosity
def eddington_luminosity(M_p, G, c, sigma_T):
    return (4 * np.pi * G * M_p * c) / sigma_T

# Eddington Accretion Rate
def eddington_accretion_rate(L_edd, epsilon, c):
    return L_edd / (epsilon * c**2)

# Accretion Rate in a 3-Body System (using conservation of mass)
def accretion_rate_mass(R, viscosity_disk):
    return np.pi * R_disk**2 * viscosity_disk

# Accretion rate in a 3-Body System (using conservation of mass and angular momentum)
def accretion_rate_md_dot(L_dot, R, omega):
    """
    Calculate the accretion rate of the dust population.

    Parameters:
        L_dot (float): Rate of change of the angular momentum of the accretion disk.
        R_disk (float): Radius of the accretion disk in m.
        omega_disk (float): Angular velocity of the accretion disk in rad/s.

    Returns:
        float: Accretion rate of the dust population in m/s.
    """
    return L_dot / (2 * np.pi * R_disk**2 * omega_disk)

L_Dot = accretion_rate_md_dot(L_dot, R_disk, omega_disk)
