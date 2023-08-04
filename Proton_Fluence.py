# main.py
from Parameters import n_be, n_o, sigma

def calculate_proton_fluence(n_be, n_o, sigma):
    """
    Calculate the proton fluence (F_p) required to obtain a given Be/O ratio.

    Parameters:
        n_be (float): Number density of beryllium (n_{Be}).
        n_o (float): Number density of oxygen (n_O).
        sigma (float): Cross-section for the spallation reaction.

    Returns:
        float: Proton fluence (F_p) in units of cm^-2.
    """
    return n_be / (sigma * n_o)

proton_fluence = calculate_proton_fluence(n_be, n_o, sigma)
#"cm^-2"
