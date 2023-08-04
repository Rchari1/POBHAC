import sympy as sp
import numpy as np
import Parameters as prm

# Define the symbols
R = sp.symbols('R')

# Gravitational potential for a mass M_p at distance R
Phi = -prm.G * prm.M_p / R

# Assume the angular velocity is the one of the disk
Omega = prm.omega_disk

def circular_velocity(R, Phi):
    """Calculate the circular velocity of a particle in the disk.
    
    Args:
        R: The radial distance from the center of the disk.
        Phi: The gravitational potential function.
    
    Returns:
        v_phi: The circular velocity of the particle.
    """
    dPhi_dR = sp.diff(Phi, R)
    v_phi = sp.sqrt(-R * dPhi_dR)
    return v_phi

def rate_of_shearing(R, Omega):
    """Calculate the rate of shearing in the disk.
    
    Args:
        R: The radial distance from the center of the disk.
        Omega: The angular velocity function.
    
    Returns:
        A: The rate of shearing in the disk.
    """
    dOmega_dR = sp.diff(Omega, R)
    A = R * dOmega_dR
    return A
