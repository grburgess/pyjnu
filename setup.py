#!/usr/bin/env python

import os
import sys

import glob

from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install



class PostInstallCommand(install):
    """
    Post-installation for installation mode.
    """

    def run(self):
        from pyjnu.utils.create_rates import create_rates
        print('creating rate database')
        create_rates()
        install.run(self)

class PostDevelopCommand(develop):
    """
    Post-installation for installation mode.
    """

    def run(self):
        from pyjnu.utils.create_rates import create_rates
        print('creating rate database')
        create_rates()
        install.run(self)


# Get the version number
# with open("pyjnu/version.py") as f:
#     version_code = compile(f.read(), "pyjnu/version.py", 'exec')
#     exec(version_code)

# Now a global __version__ is available

# This dynamically loads a module and return it in a variable.
# Will use it for check optional dependencies

def is_module_available(module_name):

    # Fast path: see if the module has already been imported.

    try:

        exec('import %s' % module_name)

    except ImportError:

        return False

    else:

        return True


# Create list of data files
def find_data_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join('..', path, filename))

    return paths

extra_files = find_data_files('pyjnu/data')

# This list will contain the messages to print just before the end of the setup
# so that the user actually note them, instead of loosing them in the tons of
# messages of the build process

setup(

    name="pyjnu",

    packages=['pyjnu',
              'pyjnu/exceptions',
              'pyjnu/utils',
              'pyjnu/io',
#              'pyjnu/io/plotting',
#              'pyjnu/config',
#              'pyjnu/test'
              ],

    #version=__version__,

    description="pyjnu: Photon and Particle Emission from astrophysical sources",

    long_description="An emissivitiy solver blah balh",

    license='GPL',

    author='begue and burgess',

    author_email='dbeugue@mpe.mpg.de',

    url='https://github.com/grburgess/pyjnu',

#    download_url='https://github.com/grburgess/pyjnu/archive/%s' % __version__,

    keywords=['spectra', 'modeling', 'neutrino', 'photon', 'gamma-ray', 'x-ray', 'multi-wavelength'],

    classifiers=[],


    # NOTE: we use '' as package name because the extra_files already contain the full path from here

    package_data={'': extra_files, },
    include_package_data=True,

    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
    
    install_requires=[
        'numpy >= 1.6',
        'scipy >=0.18',
        'astropy>=1.3.3',
        'matplotlib',
        'pyyaml',
        'dill',
        'astromodels>=0.4.0',
        'pandas',
        'ipython<=5.9',
        'h5py'
    ],



    
    extras_require={
            'tests': [
                'pytest',],
            'docs': [
                'sphinx >= 1.4',
                'sphinx_rtd_theme',
                'nbsphinx']}

    ) # End of setup()



