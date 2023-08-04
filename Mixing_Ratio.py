import numpy as np
from scipy import constants as const
from Accretion_Rates import eddington_luminosity, eddington_accretion_rate, accretion_rate_mass, accretion_rate_md_dot

from Parameters import *

def calculate_mixing_ratio(f_pl):
    """
    Calculate the mixing ratio below the meteoric source.

    Parameters:
        f_pl (float): Abundance by number of the element in the dust particles.

    Returns:
        float: Mixing ratio (n/nb).
    """
    F = -6.5e9 * f_pl * (epsilon * accretion_rate_md_dot(L_dot, R_disk, omega_disk) / 1e8) * (R_disk / R_exomoon)**(-2)  # Particle flux

    return 4.1e-3 * f_pl * (epsilon * accretion_rate_md_dot(L_dot, R_disk, omega_disk) / 1e8) * (M_exoplanet / 1.898e30)**(-1) * (T / 1000)**(0.5) * ((mb * mu) / (mu * (mu + mb)))**(0.5)


