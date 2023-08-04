import scipy.constants as const

# Constants
G = const.G             # gravitational constant in m^3 kg^-1 s^-2
c = const.c             # speed of light in m/s
sigma_T = 6.65e-29      # Thomson scattering cross-section in m^2
epsilon = 0.1           # radiative efficiency for a Schwarzschild black hole

M_p = 1.989e30          # mass of the black hole in kg, assuming it's similar to the mass of the sun
M_p *= 10**8            # scaling up to represent the mass of the supermassive black hole


# Computationl Parameters

M_exomoon = 7.35e22     # mass of the exomoon in kg, assuming it's similar to the moon's mass
R_exomoon = 1.737e6     # radius of the exomoon in m, assuming it's similar to the moon's radius

# Exoplanet Parameters (assumptions)
M_exoplanet = 5.972e24  # mass of the exoplanet in kg, assuming it's similar to the Earth's mass
R_exoplanet = 6.371e6   # radius of the exoplanet in m, assuming it's similar to the Earth's radius

# Accretion Disk Parameters (assumptions)
R_disk = 1e12           # radius of the accretion disk in m, assuming it's 0.01 parsec
R_disk_cm = R_disk/100 			# radius of the accretion disk in cm
viscosity_disk = 1e-3   # viscosity of the disk in m^2/s, assuming it's similar to water
omega_disk = 1e-7       # angular velocity of the disk in rad/s, arbitrary
L_dot = 1.5e30 # Rate of change of the angular momentum of the accretion disk in g cm^2/s^2


# Proton Fluence (assumptions)
n_be = 1e15  # Example number density of beryllium in cm^-3
n_o = 1e16   # Example number density of oxygen in cm^-3
sigma = 2e-27  # Example cross-section for spallation reaction in cm^2

# Mixing Ratio's (assumptions)
T = 1000  # Temperature in K (assumed)
mb = 16  # Mass of the background element (oxygen) in AMU (assumed)
mu = 9.012  # Mean mass of the elements making up the dust particles in AMU (assumed)
f_pl = 0.01 
