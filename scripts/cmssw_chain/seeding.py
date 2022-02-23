import os
import numpy as np
import h5py
import configuration as conf

# Event by event smoothing
storeIn  = h5py.File(conf.SmoothingOut, mode='r')
storeOut = h5py.File(conf.SeedingOut, mode='w')

def validation(outName, arr):
    with open(outName, 'w') as afile:
        afile.write('{}\t{}\t{}\n'.format(arr[0], arr[1], arr[2]))

for falgo in conf.FesAlgos:
    keys = [x for x in storeIn.keys() if falgo in x]

    for key in keys:
        energies, weighted_x, weighted_y = storeIn[key]

        # add unphysical top and bottom R/z rows for edge cases
        # fill the rows with negative (unphysical) energy values
        # boundary conditions on the phi axis are satisfied by 'np.roll'
        phiPad = -1 * np.ones((1,conf.NbinsPhi))
        energies = np.concatenate( (phiPad,energies,phiPad) )

        #remove padding
        slc = slice(1,energies.shape[0]-1)

        south = np.roll(energies, shift=1,  axis=0)[slc]
        north = np.roll(energies, shift=-1, axis=0)[slc]
        east  = np.roll(energies, shift=-1, axis=1)[slc]
        west  = np.roll(energies, shift=1,  axis=1)[slc]
        northeast = np.roll(energies, shift=(-1,-1), axis=(0,1))[slc]
        northwest = np.roll(energies, shift=(-1,1),  axis=(0,1))[slc]
        southeast = np.roll(energies, shift=(1,-1),  axis=(0,1))[slc]
        southwest = np.roll(energies, shift=(1,1),   axis=(0,1))[slc]

        energies = energies[slc]
        
        maxima = ( (energies > conf.histoThreshold ) &
                   (energies >= south) & (energies > north) & (energies >= east) & (energies > west) &
                   (energies >= northeast) & (energies > northwest) &
                   (energies >= southeast) & (energies > southwest) )

        seeds_idx = np.nonzero(maxima)

        res = (energies[seeds_idx], weighted_x[seeds_idx], weighted_y[seeds_idx])
        # if '4681' in key:
        #     validation('outLocalSeeding.txt', res)

        for r1,r2,r3 in zip(res[0],res[1],res[2]):
            print('E={}\tX={}\tY={}'.format(r1,r2,r3))

        storeOut[key] = res
        storeOut[key].attrs['columns'] = ['seedEn', 'seedX', 'seedY']
        storeOut[key].attrs['doc'] = 'Smoothed energies and projected bin positions of seeds'

storeIn.close()
storeOut.close()
