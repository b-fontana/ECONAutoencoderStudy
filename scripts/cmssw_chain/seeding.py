import os
import numpy as np
import h5py
import configuration as conf

# Event by event smoothing
storeIn  = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.SmoothingOut), mode='r')
storeOut = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.SeedingOut), mode='w')

for falgo in conf.FesAlgos:
    keys = [x for x in storeIn.keys() if falgo in x]

    for key in keys:
        energies, weighted_x, weighted_y = storeIn[key]

        # add top and bottom phi rows for boundary conditions
        # fill the rows with negative (unphysical) energy values
        phiPad = -1 * np.zeros((1,conf.NbinsPhi))
        energies = np.concatenate( (phiPad,energies,phiPad) )

        nrows = energies.shape[0]
        zeros = np.zeros((nrows,1))
        ones = np.ones((nrows,1))

        #remove padding
        s = slice(1,energies.shape[0]-1)

        south = np.roll(energies, shift=1,  axis=0)[s]
        north = np.roll(energies, shift=-1, axis=0)[s]
        east  = np.roll(energies, shift=-1, axis=1)[s]
        west  = np.roll(energies, shift=1,  axis=1)[s]
        northeast = np.roll(energies, shift=(-1,-1), axis=(0,1))[s]
        northwest = np.roll(energies, shift=(-1,1),  axis=(0,1))[s]
        southeast = np.roll(energies, shift=(1,-1),  axis=(0,1))[s]
        southwest = np.roll(energies, shift=(1,1),   axis=(0,1))[s]

        energies = energies[s]
        
        maxima = ( (energies > conf.histoThreshold ) &
                   (energies >= south) & (energies > north) & (energies >= east) & (energies > west) &
                   (energies >= northeast) & (energies > northwest) &
                   (energies > southeast) & (energies >= southwest) )

        seeds_idx = np.nonzero(maxima)

        res = (energies[seeds_idx], weighted_x[seeds_idx], weighted_y[seeds_idx])

        storeOut[key] = res
        storeOut[key].attrs['columns'] = ['energies', 'weighted_x', 'weighted_y']
        storeOut[key].attrs['doc'] = 'Smoothed energies and projected bin positions of seeds'

storeIn.close()
storeOut.close()    
