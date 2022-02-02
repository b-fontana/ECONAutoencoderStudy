import os
import numpy as np
import pandas as pd
import h5py
import configuration as conf
        
# Event by event smoothing
storeInSeeds  = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.SeedingOut), mode='r')
storeInTC  = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.FillingOut), mode='r')
storeOut = pd.HDFStore( os.path.join(os.environ['PWD'], conf.DataFolder, conf.ClusteringOut), mode='w')

for falgo in conf.FesAlgos:
    bin_keys = [x for x in storeInSeeds.keys() if falgo in x ]
    tc_keys  = [x for x in storeInTC.keys() if falgo in x and '_tc' in x]
    assert(len(bin_keys) == len(tc_keys))

    radiusCoeffB = conf.CoeffB

    for key1, key2 in zip(tc_keys, bin_keys):
        tc = storeInTC[key1]

        # check columns via `tc.attrs['columns']`
        radiusCoeffA = np.array( [conf.CoeffA[int(xi)-1] for xi in tc[:,5]] )
        minDist = radiusCoeffA + radiusCoeffB * (conf.MidRadius - np.abs(tc[:,4]))

        energies, weighted_x, weighted_y = storeInSeeds[key2]

        dRs = np.array([])
        for iseed, (en, wx, wy) in enumerate(zip(energies, weighted_x, weighted_y)):
            dR = np.sqrt( (tc[:,2]-wx)*(tc[:,2]-wx) + (tc[:,3]-wy)*(tc[:,3]-wy) )

            if dRs.shape == (0,):
                dRs = np.expand_dims(dR, axis=-1)
            else:
                dRs = np.concatenate((dRs, np.expand_dims(dR, axis=-1)), axis=1)

        #checks if each seeds has at least one seed which lies below the threshold
        pass_threshold = dRs < np.expand_dims(minDist, axis=-1)
        pass_threshold = np.logical_or.reduce(pass_threshold, axis=1)

        seeds_indexes = np.argmin(dRs, axis=1)
        seeds_energies = np.array( [energies[xi] for xi in seeds_indexes] )
        assert(tc[:].shape[0] == seeds_energies.shape[0]) # axis 0 stands for trigger cells

        seeds_indexes  = np.expand_dims( seeds_indexes[pass_threshold], axis=-1 )
        seeds_energies = np.expand_dims( seeds_energies[pass_threshold], axis=-1 )

        tc = tc[:][pass_threshold]
        res = np.concatenate((tc, seeds_indexes, seeds_energies), axis=1)

        key = key1.replace('_tc', '_clustered')
        cols = ['Rz_bin', 'tc_phi_bin', 'proj_x', 'proj_y',
                'tc_eta', 'tc_layer', 'tc_mipPt', 'seed_idx', 'seed_energy']
        assert(len(cols)==res.shape[1])
        storeOut[key] = pd.DataFrame(res, columns=cols)

storeInSeeds.close()
storeInTC.close()
storeOut.close()
