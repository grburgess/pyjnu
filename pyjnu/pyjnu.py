"""
pyjnu.py
Authors:
    -Stephan Meighen-Berger
Main run file for the python implementation of
pyjnu.
"""
# External
import numpy as np

import cPickle
from math import pi
# Internal
from particle import particle
from geometry import geometry
from config import config
from constants import pdg_id_lib, unit_conv, phys_const
from nodes import nodes
from scipy.integrate import quad

import h5py

from pyjnu.io.package_data import get_path_of_data_file
from pyjnu.config import pyjnu_config

class PyRun(object):
    """
    class: PyRun
    PyRun is the interface to the package for the user.
    Parameters:
        -None
    Returns:
        -None
    """

    def __init__(self):
        """
        function: __init__
        Function to initialize the class.
        This will prepare the instance for solving
        Parameters:
            -kwargs:
                Settings that deviate from the config file
        Returns:
            -None
        """


        # Loading rates
        
        self._read_rates()

        
        # User variables

    
        
        # Creating particles
        
        self.particles = {}
        for key in pdg_id_lib:
            self.particles[pdg_id_lib[key]] = particle(pdg_id_lib[key])
        self.particles['22_local'] = particle('22_local')
        
        


        def _read_rates(self):
            """
            """
            with h5py.File(get_path_of_data_file('rates.h5')) as f:

                self._selfynch_rate = f['synchrotron_rate']
                self._opacity = f['opacitiy']
                self._ic_rate = f['ic_rate']

    def _init_spectrum(self, PDG_ID):
        """
        function: _init_spectrum
        Function to create the initial spectrum of particles.
        This is currently a version for the stable state case
        Parameters:
            -str PDG_ID:
                PDG_ID of the desired particle
        Returns:
            -None
        """
        # Shorthand
        part = self.particles[PDG_ID]
        f = lambda x: self._int_bpl(x, [config['p1'], config['gb']])
        tmp_flux = np.zeros(self.particles[PDG_ID].size)
        for i in range(0, self.particles[PDG_ID].size):
            if (part.e_borders[i + 1] > config['gmin'] and part.e_borders[i] < config['gmax']):
                Imin = np.max([part.e_borders[i], config['gmin']])
                Imax = np.min([part.e_borders[i + 1], config['gmax']])
                # Need to compare results with and without refinement...
                # or just use the scipy version below?
                tmp_flux[i] = (self.geom.K * quad(f, Imin, Imax)[0])
        # Normalizing
        total = np.sum(tmp_flux)

        self.particles[PDG_ID].flux['0'] = tmp_flux * self.geom.K / total

    def _int_bpl(self, x, param):
        """
        function: __int_pbl
        Returns exponential depending on the input
        Parameters:
            -float x:
                Number
            -np.array param:
                Parameter list with 2 elements
        Returns:
            -float res:
                The exponential
        """
        if (x < param[1]):
            return ((x / param[1])**(-param[0]))
        else:
            return ((x / param[1])**(-param[0] - 1))

   

    def rescale(self, ID1, ID2):
        """
        function: rescale
        Rescales the energy grid from one particle
        to the other.
        Parameters:
            -str ID1:
                PDG_ID of the particle to use as a baseline
            -str ID2:
                PDG_ID of particle to overwrite
        Returns:
            -None
        """

        nubar = self.particles[ID1].e_borders[0] * self.geom.nuB
        part1 = self.particles[ID1]
        part2 = self.particles[ID2]
        # Case A : nubar < ph2.Eb[0],
        # that is to say we have to forget some of the points
        if nubar < part2.e_borders[0]:

            #  TODO: still need to be rewritten using slicing
            i = 0
            
            while nubar < self.particles[ID2].e_borders[0]:
                i += 1
                nubar = (self.particles[ID1].e_borders[i] * self.geom.nuB)
            i = i - 1
            # Now  ph1.Eb[i] < ph2.Eb[0] <= ph1.Eb[i+1]
            for j in range(i, config['grid_' + ID1] - 1):
                self.particles[ID2].flux['0'][j - 1] = (
                    part1.flux['0'][j] * (part1.e_borders[j + 1] * self.geom.nuB - part2.e_borders[j - i]) /
                    (part1.e_diff[j] * self.geom.nuB) +
                    part1.flux['0'][j + 1] * (part2.e_borders[j - i + 1] - part1.e_borders[j + 1] * self.geom.nuB) /
                    (part1.e_diff[j + 1] * self.geom.nuB))
        # Case B : nubar > ph2.Eb[0]
        else:
            
            i = np.where(nubar < part2.e_borders)[0][0] - 1
            # Now ph2.Eb[i] < ph1.Eb[0] * nuB < ph2.Eb[i+1]
            self.particles[ID2].flux['0'][i] = (
                part1.flux['0'][0] * (part2.e_borders[i + 1] - part1.e_borders[0] * self.geom.nuB) / (
                    (part1.e_diff[0]) * self.geom.nuB))
            low = i + 1
            up = config['grid_' + ID1] - 1
            self.particles[ID2].flux['0'][low:up] = (
                part1.flux['0'][low - i:up - i] *
                (part2.e_borders[low + 1:up + 1] - part1.e_borders[low - i:up - i] * self.geom.nuB) /
                (part1.e_diff[low - i:up - i] * self.geom.nuB) + part1.flux['0'][low - i - 1:up - i - 1] *
                (-part2.e_borders[low:up] + part1.e_borders[low - i:up - i] * self.geom.nuB) / (
                    (part1.e_diff[low - i - 1:up - i - 1]) * self.geom.nuB))

    def solve_steady(self):
        """
        function: solve_steady
        Solves for the steady state solution
        Parameters:
            -None
        Returns:
            -None
        """
        
        self.geom = geometry(pyjnu_config['Bfield'], pyjnu_config['delta'], pyjnu_config['R'], pyjnu_config['d'],
                             pyjnu_config['z'])
        
        
        self._init_spectrum('11')
        
        # Synchrotron
        
        A = self.geom.A
        # Move below part to particl
        self.particles['22_local'].flux['0'] = np.zeros(config['grid_22_local'])

        # transition from electron_j to gamma_i

        self.particles['22_local'].flux['0'] = \
            A * np.dot(self._synch_rate,
                       self.particles['11'].flux['0'] /
                       self.particles['11'].e_diff) / self.particles['22'].e_diff

        # There is a slight difference in the particle spectra here
        # debug_out = [(self.particles['22_local'].e_grid[i] * self.geom.delta / (1. + self.geom.z),
        #               self.particles['22_local'].flux['0'][i])
        #              for i in range(0, config['grid 22_local'])]
        # print(debug_out)
        # This is not implemented in Damien's code
        # Also not sure what type of function onemexpdexp is,
        # it is not implemented in Damien's code
        # Opacity
        # ASSA = self.geom.ASSA
        # tauCst = self.geom.tauCst
        # for i in range(0, config['grid 22']):
        #     summtau = 0.
        #     for j in range(0, config['grid 11']):
        #         summtau += (self.opacity_rate[i][j] *
        #                     (self.particles['11'].flux['0'][j] /
        #                      (self.particles['11'].e_borders[j+1] -
        #                       self.particles['11'].e_borders[j])))
        #     summtau = ASSA*tauCst*summtau
        #     self.particles['22'].flux['0'][i] = (
        #         self.particles['22'].flux['0'][i] * onemexpdexp(summtau)
        #     )

        # IC
        
        self.particles['22_local'].flux['2'] = np.zeros(config['grid_22_local'])
        self.particles['22'].flux['0'] = np.zeros(config['grid_22'])
        self.particles['22'].flux['2'] = np.zeros(config['grid_22'])
        # Rescaling the photon grid according to the local version
        # print('The local photon before rescaling')
        # print([(self.particles['22_local'].e_grid[i] * self.geom.delta / (1. + self.geom.z),
        #         self.particles['22_local'].flux['0'][i],
        #         self.particles['22_local'].flux['2'][i])
        #         for i in range(0, len(self.particles['22_local'].e_grid))])
        self.rescale('22_local', '22')
        # print('The resulting photon after rescaling')
        # print([(self.particles['22'].e_grid[i] * self.geom.delta / (1. + self.geom.z),
        #         self.particles['22'].flux['0'][i],
        #         self.particles['22'].flux['2'][i])
        #         for i in range(0, len(self.particles['22'].e_grid))])
        # Results differ slightly compared to the original but not too bad
        # debug_out = [(self.particles['22'].e_grid[i] * self.geom.delta / (1. + self.geom.z),
        #               self.particles['22'].flux['0'][i])
        #              for i in range(0, config['grid 22'])]
        # print(debug_out)
        # Need to rework this part

        summIC = np.zeros(config['grid_22'])
        for i in range(0, config['grid_22']):

            # TODO: Do we need the mask? or is it always one?
            mask = self._ic_rate[i, :, :, 2] == 1.

            # This is now a matrix vector multiplication
            term = (np.dot(self._ic_rate[i, :, :, 0] * mask,
                           self.particles['11'].flux['0'] / self.particles['11'].e_diff))
            summIC[i] = np.sum(self.particles['22'].flux['0'] * term)

        summICfinal = (0.75 * phys_const['c'] * phys_const['sigmaT'] * summIC / self.particles['22'].e_diff)

        # Multiply the emissivity by the volume.
        self.particles['22'].flux['2'] = (
            (4. * pi * self.geom.R**3 * self.geom.delta**3 /
             (self.geom.d**2 * unit_conv['Mpc_to_cm']**2 * 4. * pi * 3.)) * (
                 (self.geom.R / phys_const['c']) * summICfinal) * phys_const['h'] * self.particles['22'].e_grid)
        self.particles['22'].flux['0'] = (
            (4. * pi * self.geom.R**3 * self.geom.delta**3 /
             (self.geom.d**2 * unit_conv['Mpc_to_cm']**2 * 4. * pi * 3.)) * self.particles['22'].flux['0'] *
            phys_const['h'] * self.particles['22'].e_grid)

        # debug_out = [(self.particles['22'].e_grid[i] * self.geom.delta / (1. + self.geom.z),
        #               self.particles['22'].flux['2'][i])
        #              for i in range(0, config['grid 22'])]
        # print(debug_out)

        # Pair
        # Is not implemented
        # print('The resulting photon after IC')
        # print([(self.particles['22'].e_grid[i] * self.geom.delta / (1. + self.geom.z),
        #         self.particles['22'].flux['0'][i],
        #         self.particles['22'].flux['2'][i])
        #        for i in range(0, len(self.particles['22'].e_grid))])
        # Normalizing pre-storage

        # Todo: Why? Think this is causes a problem when
        # changing parameters and refitting

        self.particles['22'].e_grid = (self.particles['22'].e_grid * self.geom.delta / (1. + self.geom.z))
        # Storing results



        # what is this?
        saveString = '../data/results.pkl'
        with open(saveString, 'wb') as f:
            cPickle.dump(self.particles, f, protocol=cPickle.HIGHEST_PROTOCOL)
        # Load with
        # with open(saveString, 'rb') as f:
        #             self.particles = cPickle.load(f)
        
        
