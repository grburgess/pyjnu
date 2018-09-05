"""
geometry.py
Authors:
    - Stephan Meighen-Berger
    - Martina Karl
Geometry file of the package. Here the treatment of
the geometry of the object will be dealt with.
"""
# External imports
from math import pi
# Local imports
from config import config
from constants import phys_const
from logger import Logger


class geometry(Logger):
    """
    class: geometry
    Geometry class of the package. This will create the
    object.
    Parameters:
        -None
    Returns:
        -None
    """

    def __init__(self, Bfield, delta, R, d, z):
        """
        function: __init__
        Initialized the class and sets parameters
        Parameters:
            -float Bfield:
                The magnetic field of the object
            -float delta:
                -
            -float R:
                The radius of the object
            -float d:
                The distance to the object (in MPc)
            -float z:
                The redshift of the object
        Returns:
            -None
        """
        self.logger.info('Creating the geometry object...')
        self.B = Bfield
        self.delta = delta
        self.R = R
        self.gc = config['gc']
        self.gmax = config['gmax']
        self.Ke = config['Ke']
        self.V = 4. / 3. * pi * self.R**3
        self.nuB = phys_const['prefactor nuB'] * self.B
        self.A = (4.0 * pi * 3.**(1. / 3.) * phys_const['qe']**2. /
                  phys_const['c'] / phys_const['h'])
        self.ASSA = (4.0 * pi * 3.**(1. / 3.) * phys_const['qe']**2. *
                     self.nuB / phys_const['c'])
        self.tauCst = self.R / (8.0 * pi * phys_const['me'] * self.nuB**2.)
        self.d = d
        self.z = z
        self.K = self.Ke * self.delta**3 / self.V
        self.logger.info('Finished creating the geometry object...')
