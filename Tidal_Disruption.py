import Parameters as p

def tidal_radius():
    """Calculate the tidal radius of an exomoon around a black hole.
    
    Parameters from parameters.py:
        p.M_p (float): The mass of the black hole in solar masses.
        p.M_exomoon (float): The mass of the exomoon in lunar masses.
        p.R_exomoon (float): The radius of the exomoon in lunar radii.

    Returns:
        float: The tidal radius in lunar radii.
    """
    r_T = p.R_exomoon / 2 * (p.M_p / p.M_exomoon)**(1/3)
    return r_T

def capture_rate(N_moon, sigma_moon):
    """Calculate the capture rate of exomoons by a black hole.
    
    Parameters:
        N_moon (float): The number density of exomoons in the vicinity of the black hole in exomoons per cubic lunar radius.
        sigma_moon (float): The velocity dispersion of the exomoons in kilometers per second.

    Uses the tidal radius calculated in the tidal_radius function.

    Returns:
        float: The capture rate in exomoons per year.
    """
    from math import pi
    r_T = tidal_radius()
    Gamma = N_moon * sigma_moon / (2 * pi * r_T**2)
    return Gamma
