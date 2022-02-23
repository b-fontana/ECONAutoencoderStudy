import os
import re
import numpy as np
import pandas as pd
import h5py
import configuration as conf
        
# Event by event smoothing
def clustering():
    storeInSeeds  = h5py.File(conf.SeedingOut, mode='r')
    storeInTC  = h5py.File(conf.FillingOut, mode='r')
    storeOut = pd.HDFStore(conf.ClusteringOut, mode='w')

    for falgo in conf.FesAlgos:
        seed_keys = [x for x in storeInSeeds.keys() if falgo in x ]
        tc_keys  = [x for x in storeInTC.keys() if falgo in x and '_tc' in x]
        assert(len(seed_keys) == len(tc_keys))

        radiusCoeffB = conf.CoeffB

        print('CoeffA: ', conf.CoeffA)
        print('CoeffB: ', conf.CoeffB)
        print('kMidRadius: ', conf.MidRadius)

        for key1, key2 in zip(tc_keys, seed_keys):
            tc = storeInTC[key1]
            tc_cols = list(tc.attrs['columns'])

            projx = tc[:,2]/tc[:,4] #tc_x / tc_z
            projy = tc[:,3]/tc[:,4] #tc_y / tc_z

            # check columns via `tc.attrs['columns']`
            radiusCoeffA = np.array( [conf.CoeffA[int(xi)-1] for xi in tc[:,6]] )
            minDist = radiusCoeffA + radiusCoeffB * (conf.MidRadius - np.abs(tc[:,5]))
            print('minDist: ', minDist)
            
            seedEn, seedX, seedY = storeInSeeds[key2]

            dRs = np.array([])
            for iseed, (en, sx, sy) in enumerate(zip(seedEn, seedX, seedY)):
                dR = np.sqrt( (projx-sx)*(projx-sx) + (projy-sy)*(projy-sy) )
                print('d: ', dR)
                if dRs.shape == (0,):
                    dRs = np.expand_dims(dR, axis=-1)
                else:
                    dRs = np.concatenate((dRs, np.expand_dims(dR, axis=-1)), axis=1)

            #checks if each seeds has at least one seed which lies below the threshold
            pass_threshold = dRs < np.expand_dims(minDist, axis=-1)
            pass_threshold = np.logical_or.reduce(pass_threshold, axis=1)

            seeds_indexes = np.argmin(dRs, axis=1)
            seeds_energies = np.array( [seedEn[xi] for xi in seeds_indexes] )
            assert(tc[:].shape[0] == seeds_energies.shape[0]) # axis 0 stands for trigger cells

            seeds_indexes  = np.expand_dims( seeds_indexes[pass_threshold], axis=-1 )
            seeds_energies = np.expand_dims( seeds_energies[pass_threshold], axis=-1 )

            tc = tc[:][pass_threshold]

            res = np.concatenate((tc, seeds_indexes, seeds_energies), axis=1)

            key = key1.replace('_tc', '_cl')

            cols = tc_cols + [ 'seed_idx', 'seed_energy']
            assert(len(cols)==res.shape[1])
            df = pd.DataFrame(res, columns=cols)

            df['cl3d_pos_x'] = df['tc_x'] * df['tc_mipPt']
            df['cl3d_pos_y'] = df['tc_y'] * df['tc_mipPt']
            df['cl3d_pos_z'] = df['tc_z'] * df['tc_mipPt']

            cl3d = df.groupby(['seed_idx']).sum()[['cl3d_pos_x', 'cl3d_pos_y', 'cl3d_pos_z', 'tc_mipPt']]
            cl3d.rename(columns={'tc_mipPt': 'cl3d_en'}, inplace=True)
            cl3d['cl3d_pos_x'] /= cl3d['cl3d_en']
            cl3d['cl3d_pos_y'] /= cl3d['cl3d_en']
            cl3d['cl3d_pos_z'] /= cl3d['cl3d_en']

            cl3d['cl3d_Rz'] = np.sqrt(cl3d['cl3d_pos_x']*cl3d['cl3d_pos_x'] + cl3d['cl3d_pos_y']*cl3d['cl3d_pos_y']) / np.abs(cl3d['cl3d_pos_z'])
            cl3d['cl3d_phi'] = np.arctan2(cl3d['cl3d_pos_y'], cl3d['cl3d_pos_x'])

            event_number = re.search('Threshold_([0-9]{1,7})_tc', key1)
            if not event_number:
                raise ValueError('The event number was not extracted!')
            cl3d['event'] = event_number.group(1)
            storeOut[key] = cl3d

    storeInSeeds.close()
    storeInTC.close()
    storeOut.close()

def validation():
    storeInLocal = pd.HDFStore(conf.ClusteringOut, mode='r')
    storeInCMSSW = h5py.File(conf.FillingOut, mode='r')

    for falgo in conf.FesAlgos:
        local_keys = [x for x in storeInLocal.keys() if falgo in x]
        cmssw_keys = [x for x in storeInCMSSW.keys() if falgo in x and '_clpos' in x]
        assert(len(local_keys) == len(cmssw_keys))

        for key1, key2 in zip(local_keys, cmssw_keys):
            local = storeInLocal[key1]
            cmssw = storeInCMSSW[key2]

            event_number = re.search('Threshold_([0-9]{1,7})_cl', key1).group(1)
            print('Event: {}'.format(event_number))
            print('Custom: NClusters={}\tPhi={}\tRz={}\tEnergy={}'
                  .format(len(local['cl3d_phi'].to_numpy()),
                          local['cl3d_phi'].to_numpy(),
                          local['cl3d_Rz'].to_numpy(),
                          local['cl3d_en'].to_numpy()))
            print('CMSSW:  NClusters={}\tPhi={}\tRz={}\tEnergy={}'
                  .format(len(cmssw[:][0]), cmssw[:][0], cmssw[:][1], cmssw[:][2]))
            print()


clustering()
validation()
