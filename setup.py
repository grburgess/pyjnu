#!/usr/bin/env python

import os
import sys

import glob

from setuptools import setup

# Get the version number
with open("pyjnu/version.py") as f:
    version_code = compile(f.read(), "pyjnu/version.py", 'exec')
    exec(version_code)

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
              'pyjnu/io/plotting',
              'pyjnu/config',
              'pyjnu/test'
              ],

    version=__version__,

    description="pyjnu: Photon and Particle Emission from astrophysical sources",

    long_description="An emissivitiy solver blah balh",

    license='GPL',

    author='begue and burgess',

    author_email='dbeugue@mpe.mpg.de',

    url='https://github.com/grburgess/pyjnu',

    download_url='https://github.com/giacomov/3ML/archive/%s' % __version__,

    keywords=['spectra', 'modeling', 'neutrino', 'photon', 'gamma-ray', 'x-ray', 'multi-wavelength'],

    classifiers=[],

    # Install configuration file in user home and in the package repository

    # data_files=[(os.path.join(os.path.expanduser('~'), '.pyjnu'), ["pyjnu/config/pyjnu_config.yml"]),
    #             ('pyjnu/config', ["pyjnu/config/pyjnu_config.yml"])
    #             ],

    # NOTE: we use '' as package name because the extra_files already contain the full path from here

        package_data={'': extra_files, },
    include_package_data=True,

    install_requires=[
        'numpy >= 1.6',
        'scipy >=0.18',
        'astropy>=1.3.3',
        'matplotlib',
        'pyyaml',
        'dill',
        'astromodels>=0.4.0',
        'pandas',
        'ipython<=5.9'
    ],

    extras_require={
            'tests': [
                'pytest',],
            'docs': [
                'sphinx >= 1.4',
                'sphinx_rtd_theme',
                'nbsphinx']}

    ) # End of setup()

# Check for optional dependencies

# optional_dependencies = {'cpyjnu': [False,'needed by HAWC plugin'],
#                          'pymultinest': [False, 'provides the Multinest sampler for Bayesian analysis'],
#                          'pyOpt': [False, 'provides more optimizers'],
#                          'ROOT': [False, 'provides the ROOT optimizer'],
#                          'ipywidgets': [False, 'provides widget for jypyter (like the HTML progress bar)']}

# for dep_name in optional_dependencies:

#     optional_dependencies[dep_name][0] = is_module_available(dep_name)

# # Now print the final messages

# print("\n\n##################")
# print("OPTIONAL FEATURES:")
# print("##################\n\n")

# for dep_name in optional_dependencies:
    
#     if optional_dependencies[dep_name][0]:
        
#         status = 'available'
    
#     else:
        
#         status = '*NOT* available'
    
#     print(" * %s is %s (%s)\n" % (dep_name, status, optional_dependencies[dep_name][1]))
