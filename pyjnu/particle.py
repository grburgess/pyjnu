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
from math import sqrt



class particle(object):
    """
    class: particle
    Class to create particles (the fluxes)
    Parameters:
        -None
    Returns:
        -None
    """

    def __init__(self, PDG_ID):
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
        self._mass = phys_const['mass_' + PDG_ID]
        self._emin = config['emin_' + PDG_ID]
        self._emax = config['emax_' + PDG_ID]
        self._size = config['grid_' + PDG_ID]

        self._step = np.exp(np.log(self.emax / self.emin) / self.size)
        self._e_grid = np.logspace(np.log10(self.emin), np.log10(self.emax),
                                  self.size, base=np.e, endpoint=False)
        self._e_borders = self.e_grid * sqrt(self.step)
        # First position in the borders
        self._e_borders = np.insert(self.e_borders, 0,
                                   self.emin / sqrt(self.step))
        self._e_diff = np.diff(self.e_borders)

        self._flux = {}
        self._dflux = {}
        #self.logger.info('Finished particle ' + PDG_ID)
