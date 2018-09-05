"""
config.py
Authors:
    -Stephan Meighen-Berger
Config file for the pyjnu runs.
Changes should for now be made here for
differing runs. Once typical variables
are established they should become class
variables of the PyRun class object.
"""

config = {
    # Run parameters
    # Currently each particle has a unique
    # energy grid which is idiotic.
    # This config file provides standards
    'grid_22': 150,
    'grid_11': 100,
    'grid_22_local': 150,
    'emin_22': 1e-2,
    'emax_22': 1e30,
    'emin_22_local': 0.01,
    'emax_22_local': 1e30,
    'emin_11': 10.,
    'emax_11': 1e8,
    'p1': 2.7,
    'prec': 1e-7,
    'loga': 0,
    'node_number': 15,

    # Geometry parameters
    'Bfield': 0.088,
    'delta': 230,
    'R': 0.19e15,
    'gc': 3.1e4,
    'gmin': 1e3,
    'gmax': 1.3e5,
    'gb': 3.1e4,
    'Ke': 2e40,

    # Position parameters
    'd': 540.0,  # Distance in Mpc
    'z': 0.116
}
