"""
geometry.py
Authors:
    - Stephan Meighen-Berger
    - Martina Karl
Geometry file of the package. Here the treatment of
the geometry of the object will be dealt with.
"""
# External imports
import numpy as np

from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import  z_at_value
import astropy.constants as constants


# Local imports
from pyjnu.config import pyjnu_config
from constants import phys_const


_43pi = 4./3. * np.pi 

class geometry(object):
    """
    class: geometry
    Geometry class of the package. This will create the
    object.
    Parameters:
        -None
    Returns:
        -None
    """

    def __init__(self, Bfield, delta, radius, luminosity_distance=None, z=None):
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

        assert (z is not None) or (luminosity_distance is not None), 'distance or redshift must be specified'

        if z is not None:

            assert luminosity_distance is None, 'specify a luminosity_distance or z not both'

            luminosity_distance = cosmo.luminosity_distance(z)
        else:

            z = z_at_value(cosmo.luminosity_distance, luminosity_distance)
            
        
        self._B = Bfield
        self._delta = delta
        self._radius = radius
        self._luminosity_distance = luminosity_distance
        self._z = z

        self._compute_geometry()


    @staticmethod
    def _compute_volume(radius):
        return _43pi * radius**3

    def _compute_geometry(self):
        """
        COMMENTS!!
        """
        self._gc = pyjnu_config['geometry']['gc']

        self._gmax = pyjnu_config['geometry']['gmax']

        self._Ke = pyjnu_config['geometry']['Ke']

        self._volume = self._compute_volume(self._radius)

        self._nuB = phys_const['prefactor nuB'] * self._B

        # self._A = (4.0 * pi * 3.**(1. / 3.) * phys_const['qe']**2. /
        #           phys_const['c'] / phys_const['h'])

        self._A = (4.0 * pi * 3.**(1. / 3.) * constants.q_e**2. /
                  constants.c / constants.h)

        # self._ASSA = (4.0 * pi * 3.**(1. / 3.) * phys_const['qe']**2. *
        #              self._nuB / phys_const['c'])

         self._ASSA = (4.0 * pi * 3.**(1. / 3.) * constants.q_e**2. *
                       self._nuB / constants.c)
        
        self._tauCst = self._radius / (8.0 * pi * constants.m_e * self._nuB**2.)
        
        self._K = self._Ke * self._delta**3 / self._V

    @property
    def volume(self):
        return self._volume

    @property
    def luminosity_distance(self):
        return self._luminosity_distance
    
        
