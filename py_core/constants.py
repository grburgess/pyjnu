"""
constants.py
Author:
    -Stephan Meighen-Berger
Collection of physical constants. The values
here will probably not be changed.
"""
phys_const = {
    'sigmaT': 6.6524e-25,
    'c': 2.99792458e10,
    'c2': 8.9875517873681764e20,
    'qe': 4.8032068e-10,
    'me': 9.1093897e-28,
    'mec2': 8.1871111680068e-7,
    'mp': 1.6749286e-24,
    'Msun': 1.989e33,
    'h': 6.6260755e-27,
    'hbar': 1.05457266e-27,
    'kb': 1.380657e-16,
    'sigmaSB': 5.67051e-5,
    'ath': 7.5657e-15,
    'Ggrav': 6.67259e-8,
    'Bcrit': 4.414e13,
    'prefactor nuB': 2.7992483604657531e6,

    # Cosmological Constants
    'Hubble': 69.6,  # km s^-1 Mpc^-1

    # Mass in units
    'me_u': 1.0,
    'mmu_u': 206.768289933,

    # Particle masses
    'mass_11': 0.,
    'mass_22': 0.,
    'mass_22_local': 0.,
    'mass_22_D': 0.
}

unit_conv = {
    'Mpc_to_cm': 3.08568e24,  # km s^-1 Mpc^-1
    's_to_year': 3.17098e-8,
    'Hz_to_eV': 4.13566553855e-15,  # Convert Hz to eV
}

kinetic_const = {
    'alphaf': 0.0072973525664,
    'lambdac': 2.4263102367e-10,
    'cyclosynchrotron': 2.7992483604657531e6
}

pdg_id_lib = {
    'electron': '11',
    'photon': '22'
}
