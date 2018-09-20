"""
pyjnu.py
Authors:
    -Stephan Meighen-Berger
    -Theo Glauch
Main run file for the python implementation of
pyjnu.
"""
# External
import numpy as np
import csv
import cPickle
from math import pi
# Internal
from particle import particle
from geometry import geometry
from config import config
from constants import pdg_id_lib, unit_conv, phys_const
from nodes import nodes
from logger import Logger
from scipy.integrate import quad


class PyRun(Logger):
    """
    class: PyRun
    PyRun is the interface to the package for the user.
    Parameters:
        -None
    Returns:
        -None
    """

    def __init__(self, **kwargs):
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
        self.logger.info('Initializing the class...')
        # User variables

        self.config = config
        self.logger.info('Setting user defined variables..')
        for key in kwargs.keys():
            self.config[key] = kwargs[key]

        self.logger.info('Finished setting the variables.')
        # Creating particles
        self.logger.info('Creating the particle instances...')
        self.particles = {}
        for key in pdg_id_lib:
            self.particles[pdg_id_lib[key]] = particle(pdg_id_lib[key])
        self.particles['22_local'] = particle('22_local')
        self.logger.info('Finished particle creation')
        # Loading rates
        self.logger.info('Loading the rates...')
        self.emissivity()
        self.logger.info('Finished loading')
        self.logger.info('Finished initialization')
        return

    def __init_spectrum(self, PDG_ID):
        """
        function: __init_spectrum
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
        f = lambda x: self.__int__bpl(x, [config['p1'], config['gb']])
        tmp_flux = np.zeros(self.particles[PDG_ID].size)
        for i in range(0, self.particles[PDG_ID].size):
            if (part.e_borders[i + 1] > config['gmin'] and
                    part.e_borders[i] < config['gmax']):
                Imin = np.max([part.e_borders[i], config['gmin']])
                Imax = np.min([part.e_borders[i + 1], config['gmax']])
                # Need to compare results with and without refinement...
                # or just use the scipy version below?
                tmp_flux[i] = (self.geom.K * quad(f, Imin, Imax)[0])
        # Normalizing
        total = np.sum(tmp_flux)

        self.particles[PDG_ID].flux['0'] = tmp_flux * self.geom.K / total

    def __int__bpl(self, x, param):
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

    def emissivity(self):
        """
        function: emissivity
        Calculates the emissions of the object
        Parameters:
            -None
        Returns:
            -None
        """
        # Synchrotron

        # Added pickle storage to maybe decrease loading time
        saveString = '../data/synchrotron_rate.pkl'
        try:
            with open(saveString, 'rb') as f:
                    self.synch_rate = cPickle.load(f)
        except:
            tmp = []
            with open('../data/synchrotron_rate.txt', 'r') as f:
                reader = csv.reader(f, delimiter='\t',
                                    quoting=csv.QUOTE_NONNUMERIC)
                for row in reader:
                    tmp.append(row)
            tmp = np.array(tmp)
            self.synch_rate = np.split(tmp[:, 2], len(tmp) / config['grid_11'])
            with open(saveString, 'wb') as f:
                        cPickle.dump(self.synch_rate,
                                     f, protocol=cPickle.HIGHEST_PROTOCOL)
            tmp = None
        # Load with
        # with open(saveString, 'rb') as f:
        #             self.particles = cPickle.load(f)
        # tmp = []
        # with open('../data/synchrotron_rate.txt', 'r') as f:
        #     reader = csv.reader(f, delimiter='\t',
        #                         quoting=csv.QUOTE_NONNUMERIC)
        #     for row in reader:
        #         tmp.append(row)
        # tmp = np.array(tmp)
        # self.synch_rate = np.split(tmp[:, 2], len(tmp) / config['grid_11'])

        # Opacity

        # Added pickle storage to maybe decrease loading time
        saveString = '../data/tau.pkl'
        try:
            with open(saveString, 'rb') as f:
                    self.opacity_rate = cPickle.load(f)
        except:
            tmp = []
            with open('../data/tau.txt', 'r') as f:
                reader = csv.reader(f, delimiter='\t',
                                    quoting=csv.QUOTE_NONNUMERIC)
                for row in reader:
                    tmp.append(row)
            tmp = np.array(tmp)
            self.opacity_rate = np.split(tmp[:, 2], len(tmp) / config['grid_11'])
            with open(saveString, 'wb') as f:
                        cPickle.dump(self.opacity_rate,
                                     f, protocol=cPickle.HIGHEST_PROTOCOL)
            tmp = None
        # tmp = []
        # with open('../data/tau.txt', 'r') as f:
        #     reader = csv.reader(f, delimiter='\t',
        #                         quoting=csv.QUOTE_NONNUMERIC)
        #     for row in reader:
        #         tmp.append(row)
        # tmp = np.array(tmp)
        # self.opacity_rate = np.split(tmp[:, 2], len(tmp) / config['grid_11'])

        # Inverse Compton

        # Added pickle storage to maybe decrease loading time
        saveString = '../data/IC_rate.pkl'
        try:
            with open(saveString, 'rb') as f:
                    self.IC_rate = cPickle.load(f)
        except:
            self.IC_rate = np.zeros((config['grid_22'],
                                    config['grid_22'],
                                    config['grid_11'],
                                    2))
            with open('../data/IC_rate.txt', 'r') as f:
                reader = csv.reader(f, delimiter='\t',
                                    quoting=csv.QUOTE_NONNUMERIC)
                for row in reader:
                    self.IC_rate[int(row[0])][int(row[1])][int(row[2])] = (
                        [row[3], row[4]]
                    )  
            with open(saveString, 'wb') as f:
                        cPickle.dump(self.IC_rate,
                                     f, protocol=cPickle.HIGHEST_PROTOCOL)

        # self.IC_rate = np.zeros((config['grid_22'],
        #                          config['grid_22'],
        #                          config['grid_11'],
        #                          2))
        # with open('../data/IC_rate.txt', 'r') as f:
        #     reader = csv.reader(f, delimiter='\t',
        #                         quoting=csv.QUOTE_NONNUMERIC)
        #     for row in reader:
        #         self.IC_rate[int(row[0])][int(row[1])][int(row[2])] = (
        #             [row[3], row[4]]
        #         )       

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
            self.logger.info("Case A")
            while nubar < self.particles[ID2].e_borders[0]:
                i += 1
                nubar = (self.particles[ID1].e_borders[i] *
                         self.geom.nuB)
            i = i - 1
            # Now  ph1.Eb[i] < ph2.Eb[0] <= ph1.Eb[i+1]
            for j in range(i, config['grid_' + ID1] - 1):
                self.particles[ID2].flux['0'][j - 1] = (
                    part1.flux['0'][j] *
                    (part1.e_borders[j + 1] * self.geom.nuB -
                     part2.e_borders[j - i]) /
                    (part1.e_diff[j] * self.geom.nuB) +
                    part1.flux['0'][j + 1] *
                    (part2.e_borders[j - i + 1] -
                     part1.e_borders[j + 1] * self.geom.nuB) /
                    (part1.e_diff[j + 1] * self.geom.nuB))
        # Case B : nubar > ph2.Eb[0]
        else:
            self.logger.info('Case B')
            i = np.where(nubar < part2.e_borders)[0][0] - 1
            # Now ph2.Eb[i] < ph1.Eb[0] * nuB < ph2.Eb[i+1]
            self.particles[ID2].flux['0'][i] = (
                part1.flux['0'][0] *
                (part2.e_borders[i + 1] -
                 part1.e_borders[0] * self.geom.nuB) /
                ((part1.e_diff[0]) * self.geom.nuB)
            )
            low = i + 1
            up = config['grid_' + ID1] - 1
            self.particles[ID2].flux['0'][low: up] = (
                part1.flux['0'][low - i: up - i] * (
                    part2.e_borders[low + 1: up + 1] -
                    part1.e_borders[low - i: up - i] * self.geom.nuB) /
                (part1.e_diff[low - i: up - i] * self.geom.nuB) +
                part1.flux['0'][low - i - 1: up - i - 1] *
                (-part2.e_borders[low: up] +
                 part1.e_borders[low - i: up - i] * self.geom.nuB) /
                ((part1.e_diff[low - i - 1: up - i - 1]) * self.geom.nuB)
            )
    
    def set_geometry(self, B, delta, R, d, z):
        """
        function: set_geometry
        Sets the geometry of the object
        Parameters:
            -float B:
                The B field
            -float delta:
                The delta parameter
            -float R:
                Radius of the object
            -float d:
                Distance to the object in MPC
            -float z:
                Redshift of the object
        Retruns:
            -None
        Note:
            d and z are both currently needed due to implementation.
            This should be changed to only one required parameter and the
            other should be calculated from it.
        """
        self.geom = geometry(B, delta, R, d, z)

    def solve_steady(self,
                     B=config['Bfield'], delta=config['delta'],
                     R=config['R'], d=config['d'],
                     z=config['z']):
        """
        function: solve_steady
        Solves for the steady state solution
        Parameters (optional):
            -float B:
                The B field
            -float delta:
                The delta parameter
            -float R:
                Radius of the object
            -float d:
                Distance to the object in MPC
            -float z:
                Redshift of the object
        Retruns:
            -None
        Note:
            d and z are both currently needed due to implementation.
            This should be changed to only one required parameter and the
            other should be calculated from it.
        """
        self.logger.info('Starting steady state solving...')
        self.logger.info('Setting the geometry...')
        self.set_geometry(B, delta, R, d, z)
        self.logger.info('Geometry set')
        self.logger.info('The initial electron spectrum...')
        self.__init_spectrum('11')
        self.logger.info('Set the initial spectrum')
        # Synchrotron
        self.logger.info('Synchrotron...')
        A = self.geom.A
        self.particles['22_local'].flux['0'] = np.zeros(config['grid_22_local'])
        # transition from electron_j to gamma_i
        self.particles['22_local'].flux['0'] = \
            A * np.dot(self.synch_rate,
                       self.particles['11'].flux['0'] /
                       self.particles['11'].e_diff) / self.particles['22'].e_diff
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
        self.logger.info('Inverse Compton...')
        self.particles['22_local'].flux['2'] = np.zeros(config['grid_22_local'])
        self.particles['22'].flux['0'] = np.zeros(config['grid_22'])
        self.particles['22'].flux['2'] = np.zeros(config['grid_22'])
        # Rescaling the photon grid according to the local version
        self.rescale('22_local', '22')
        summIC = np.zeros(config['grid_22'])
        for i in range(0, config['grid_22']):

            # This is now a matrix vector multiplication
            term = (np.dot(self.IC_rate[i, :, :, 0],  # * mask,
                    self.particles['11'].flux['0'] /
                    self.particles['11'].e_diff))
            summIC[i] = np.sum(self.particles['22'].flux['0'] * term)

        summICfinal = (
            0.75 * phys_const['c'] * phys_const['sigmaT'] * summIC /
            self.particles['22'].e_diff)

        # Multiply the emissivity by the volume.
        self.particles['22'].flux['2'] = (
            (4. * pi * self.geom.R**3 * self.geom.delta**3 /
             (self.geom.d**2 * unit_conv['Mpc_to_cm']**2 * 4. * pi * 3.)) *
            ((self.geom.R / phys_const['c']) * summICfinal) *
            phys_const['h'] * self.particles['22'].e_grid
        )
        self.particles['22'].flux['0'] = (
            (4. * pi * self.geom.R**3 * self.geom.delta**3 /
             (self.geom.d**2 * unit_conv['Mpc_to_cm']**2 * 4. * pi * 3.)) *
            self.particles['22'].flux['0'] * phys_const['h'] *
            self.particles['22'].e_grid)

        # Pair
        # Is not implemented
        # Normalizing pre-storage

        # Todo: Why? Think this is causes a problem when
        # changing parameters and refitting

        self.particles['22'].e_grid = (
            self.particles['22'].e_grid * self.geom.delta / (1. + self.geom.z)
        )
        # Storing results
        self.logger.info('Storing results...')
        saveString = '../data/results.pkl'
        with open(saveString, 'wb') as f:
                    cPickle.dump(self.particles,
                                 f, protocol=cPickle.HIGHEST_PROTOCOL)
        # Load with
        # with open(saveString, 'rb') as f:
        #             self.particles = cPickle.load(f)
        self.logger.info('Finished steady state solving...')
        self.logger.info('Results stored in particle fluxes...')
        handlers = self.logger.handlers[:]
        for handler in handlers:
            handler.close()
            self.logger.removeHandler(handler)
