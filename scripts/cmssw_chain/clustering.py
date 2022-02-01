import os
import numpy as np
import h5py
import configuration as conf
        
# Event by event smoothing
storeIn  = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.SeedingOut), mode='r')
storeOut = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.ClusteringOut), mode='w')

for falgo in conf.FesAlgos:
    keys = [x for x in storeIn.keys() if falgo in x]

    for key in keys:
        energies, weighted_x, weighted_y = storeIn[key]

        #storeOut[key] = res

storeIn.close()
storeOut.close()
