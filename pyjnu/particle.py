"""
particle.py
Authors:
    -Martina Karl
    -Stephan Meighen-Berger
Deals with the different interactions and particles
of the models.
"""

import numpy as np
from config import config
from constants import phys_const

from pyjnu.utils.interval import IntervalSet



class Particle(object):
    """
    class: particle
    Class to create particles (the fluxes)
    Parameters:
        -None
    Returns:
        -None
    """

    def __init__(self, emin, emax, n_grid_points, mass):
        """
        function: __init__
        Function to initialize the instance.
        Parameters:
            -str PDG_ID:
                The PDG_ID of the particle
        Returns:
            -None
        """
        #self.logger.info('Creating particle ' + PDG_ID)



        self._mass = mass
        self._grid = IntervalSet.from_logrithmic(emin, emax, size)
        
        # self._step = np.exp(np.log(self._emax / self._emin) / self._size)
        # self._e_grid = np.logspace(np.log10(self._emin), np.log10(self._emax),
        #                           self._size, base=np.e, endpoint=False)
        # self._e_borders = self._e_grid * sqrt(self._step)
        # # First position in the borders
        # self._e_borders = np.insert(self._e_borders, 0,
        #                            self._emin / sqrt(self._step))
        # self._e_diff = np.diff(self._e_borders)

        self._flux = {}
        self._dflux = {}
        #self._logger.info('Finished particle ' + PDG_ID)

    @property
    def energy_grid(self):
        return self._grid


class Electron(Particle):

    def __init__(self):

        # self._mass = phys_const['mass_' + PDG_ID]
        # self._emin = config['emin_' + PDG_ID]
        # self._emax = config['emax_' + PDG_ID]
        # self._size = config['grid_' + PDG_ID]


        super(Electron, self).__init__(emin, emax, n_grid_points, mass)

class Photon(Particle):

    def __init__(self):

        # self._mass = phys_const['mass_' + PDG_ID]
        # self._emin = config['emin_' + PDG_ID]
        # self._emax = config['emax_' + PDG_ID]
        # self._size = config['grid_' + PDG_ID]


        super(Photon, self).__init__(emin, emax, n_grid_points, mass)
