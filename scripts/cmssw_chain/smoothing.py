import os
import numpy as np
import h5py
import configuration as conf


def valid1(energies, infile, outfile):
    """
    compares all values of 2d histogram between local and CMSSW versions
    """
    flocal  = open('outLocalBeforeSmoothing.txt', 'w')
    fremote = open('outCMSSWBeforeSmoothing.txt', 'r')
    lines = fremote.readlines()

    for line in lines:
        l = line.split('\t')
        if l[0]=='\n' or '#' in l[0]:
            continue
        bin1 = int(l[0])
        bin2 = int(l[1])
        val_remote = float(l[2].replace('\n', ''))
        val_local = energies[bin1,bin2]
        if abs(val_remote-val_local)>0.001:
            print('Diff found! Bin1={}\t Bin2={}\tRemote={}\tLocal={}'.format(bin1, bin2, val_remote, val_local))

    for bin1 in range(conf.NbinsRz):
        for bin2 in range(conf.NbinsPhi):
            flocal.write('{}\t{}\t{}\n'.format(bin1, bin2, np.around(energies[bin1,bin2], 6)))
            
    flocal.close()
    fremote.close()


def smoothAlongRz(arr):
    """
    Smoothes the energy distribution of cluster energies deposited on trigger cells.
    Works along the Rz direction ("horizontal")
    """
    weights = 2 * np.ones_like(arr)
    weights[[0,conf.NbinsRz-1],:] = 1.5

    # add top and bottom phi rows for boundary conditions
    phiPad = np.zeros((1,conf.NbinsPhi))
    arr = np.concatenate( (phiPad,arr,phiPad) )

    arr_new = arr + ( np.roll(arr, shift=1, axis=0) + np.roll(arr, shift=-1, axis=0) ) * 0.5

    # remove top and bottom phi rows
    arr_new = arr_new[1:arr_new.shape[0]-1]
    assert(arr_new.shape[0] == conf.NbinsRz)
    return arr_new / weights

    # for bin1 in range(conf.NbinsRz):
    #     # Take into account edges with only one side up or down
    #     weight = 1.5 if (bin1 == 0 or bin1 == conf.NbinsRz - 1) else 2.

    #     for bin2 in range(conf.NbinsPhi):
    #         content = getBinContent(df, bin1, bin2, oldcol)
    #         contentDown = getBinContent(df, bin1 - 1, bin2, oldcol) if bin1 > 0 else 0.
    #         contentUp   = getBinContent(df, bin1 + 1, bin2, oldcol) if bin1 < (conf.NbinsRz - 1) else 0.
    #         content += (0.5 * contentDown + 0.5 * contentUp) / weight

    #         if content>0.:
    #             setBinContent(val=content,
    #                           df=df, bin1=bin1, bin2=bin2,
    #                           var=newcol, zerovar=oldcol)
    #         # bin.weighted_x = bin_orig.weighted_x;
    #         # bin.weighted_y = bin_orig.weighted_y;
    # return df.sort_index()

def smoothAlongPhi(arr):
    """
    Smoothes the energy distribution of cluster energies deposited on trigger cells.
    Works along the Phi direction ("horizontal")
    """
    arr_new = np.zeros_like(arr)

    nBinsSide = (np.array(conf.BinSums, dtype=np.int32) - 1) / 2; # one element per Rz bin
    assert(nBinsSide.shape[0] == conf.NbinsRz)
    area = (1 + 2.0 * (1 - 0.5**nBinsSide)) # one element per Rz bin

    if conf.SeedsNormByArea:
        R1 = conf.MinROverZ + bin1 * (conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz
        R2 = R1 + ((conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz)
        area = area * ((np.pi * (R2**2 - R1**2)) / conf.NbinsPhi);
    else:
        #compute quantities for non-normalised-by-area histoMax
        #The 0.1 factor in bin1_10pct is an attempt to keep the same rough scale for seeds.
        #The exact value is arbitrary.
        bin1_10pct = int(0.1) * conf.NbinsRz
        R1_10pct = conf.MinROverZ + bin1_10pct * (conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz
        R2_10pct = R1_10pct + ((conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz)
        area_10pct_ = ((np.pi * (R2_10pct**2 - R1_10pct**2)) / conf.NbinsPhi)
        area = area * area_10pct_;

    # print( "nBinsSide=", nBinsSide, "; ",
    #        "binSums[b1]=", conf.BinSums, "; ",
    #        "area=", area, "; ",
    #        "seeds_norm_by_area=", conf.SeedsNormByArea, "; ",
    #        "area_10pct_=", area_10pct_, "; ",
    #        "bin1_10pct=", bin1_10pct, "; ",
    #        "R1_10pct=", R1_10pct, "; ",
    #        "R2_10pct=", R2_10pct, "; ",
    #        "kROverZMin_=", conf.MinROverZ, "; ",
    #        "kROverZMax_=", conf.MaxROverZ, "; "
    #       )

    # loop per chunk of (equal Rz) rows with a common shift to speedup
    # unfortunately np.roll's 'shift' argument must be the same for different rows
    for idx in np.unique(nBinsSide):
        roll_indices = np.where(nBinsSide == idx)[0]
        arr_copy = arr[roll_indices,:]
        arr_smooth = arr[roll_indices,:]
        for nside in range(1, int(idx)+1):
            arr_smooth += ( (np.roll(arr_copy, shift=nside,  axis=1) +
                             np.roll(arr_copy, shift=-nside, axis=1))
                            / (2**nside) )
        arr_new[roll_indices,:] = arr_smooth / np.expand_dims(area[roll_indices], axis=-1)

    return arr_new * conf.areaPerTriggerCell

def printHistogram(arr):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if arr[i,j] == 0:
                print('-', end='|')
            else:
                print('X', end='|')
        print()

def validation(outName, arr):
    with open(outName, 'w') as afile:
        for bin1 in range(conf.NbinsRz):
            for bin2 in range(conf.NbinsPhi):
                afile.write('{}\t{}\t{}\n'.format(bin1, bin2, energies[bin1,bin2]))

def createHistogram(event):
    """
    Creates a 2D histogram with fixed (R/z vs Phi) size.
    The input event must be a 2D array where the inner axis encodes, in order:
    - 1: R/z bin index
    - 2: Phi bin index
    - 3: Value of interest ("counts" of the histogram, "z axis")
    """
    arr = np.zeros((conf.NbinsRz, conf.NbinsPhi))

    for ev in event[:]:
        assert(ev[0] >= 0)
        assert(ev[1] >= 0)
        rzbin = int(ev[0])
        phibin = int(ev[1])
        arr[rzbin,phibin] = ev[2]

    return arr

# Event by event smoothing
storeIn  = h5py.File(conf.FillingOut, mode='r')
storeOut = h5py.File(conf.SmoothingOut, mode='w')

for falgo in conf.FesAlgos:
    keys = [x for x in storeIn.keys() if falgo in x and '_group' in x]

    for key in keys:
        #print(key)
        energies   = createHistogram( storeIn[key][:,[0,1,2]] )            
        weighted_x = createHistogram( storeIn[key][:,[0,1,3]] )
        weighted_y = createHistogram( storeIn[key][:,[0,1,4]] )

        if '187544' in key:
            valid1(energies,
                   infile='outLocalBeforeSmoothing.txt',
                   outfile='outCMSSWBeforeSmoothing.txt')

            #validation('outLocalBeforeSmoothing.txt', energies)
        
        #printHistogram(ev)
        energies = smoothAlongPhi(energies)

        # if '187544' in key:
        #     valid1(energies, '187544',
        #            infile='outLocalHalfSmoothing.txt',
        #            outfile='outCMSSWHalfSmoothing.txt')

            #validation('outLocalHalfSmoothing.txt', energies)

        #printHistogram(ev)
        energies = smoothAlongRz(energies)
        #printHistogram(ev)
        storeOut[key] = (energies, weighted_x, weighted_y)
        storeOut[key].attrs['columns'] = ['energies', 'weighted_x', 'weighted_y']
        storeOut[key].attrs['doc'] = 'Smoothed energies and projected bin positions'
        
storeIn.close()
storeOut.close()
