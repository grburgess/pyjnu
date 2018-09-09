import h5py
import numpy as np

from pyjnu.io.package_data import get_path_of_data_file


def create_rates():

    # synchrotron

    # open the synchrotron txt file
    
    tmp_sync = np.genfromtxt(get_path_of_data_file('synchrotron_rate.txt'))

    # the first two columns are indices
    # we need to figure out how many elements there are
    nx = int(tmp_sync[:,0].max()) + 1
    ny = int(tmp_sync[:,1].max()) + 1

    # now create an appropriate size array
    synch_rate = np.zeros((nx, ny))

    # fill the matrix
    # TODO: there is a faster way to do this. 
    for i, (x, y) in enumerate(tmp_sync[:,[0,1]]):

        synch_rate[int(x), int(y)] = tmp_sync[i, 2]

    # tau

    # open the tau txt file
    
    tmp_tau = np.genfromtxt(get_path_of_data_file('tau.txt'))

    # the first two columns are indices
    # we need to figure out how many elements there are
    nx = int(tmp_tau[:,0].max()) + 1
    ny = int(tmp_tau[:,1].max()) + 1

    # now create an appropriate size array
    opactiy = np.zeros((nx, ny))

    # fill the matrix
    for i, (x, y) in enumerate(tmp_tau[:,[0,1]]):

        opactiy[int(x), int(y)] = tmp_tau[i, 2]
     
    # IC

    # open the IC txt file
    
    tmp_ic = np.genfromtxt(get_path_of_data_file('IC_rate.txt'))

    # the first two columns are indices
    # we need to figure out how many elements there are
    nx = int(tmp_ic[:,0].max()) + 1
    ny = int(tmp_ic[:,1].max()) + 1
    nz = int(tmp_ic[:,2].max()) + 1

    # now create an appropriate size array
    ic_rate = np.zeros((nx, ny, nz, 3))

    # fill the matrix
    for i, (x, y, z) in enumerate(tmp_ic[:,[0,1,2]]):

        ic_rate[int(x), int(y), int(z), : ] =  [ tmp_ic[i][3],tmp_ic[i][4] , 1 ]



    with h5py.File(get_path_of_data_file('rates.h5'), 'w') as f:

        f.create_dataset('synchrotron_rate', data=synch_rate, compression='lzf')

        f.create_dataset('opacity', data=opactiy, compression='lzf')

        f.create_dataset('ic_rate', data=ic_rate, compression='lzf')
