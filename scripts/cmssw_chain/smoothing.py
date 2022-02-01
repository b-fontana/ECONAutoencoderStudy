import os
import numpy as np
import sys; np.set_printoptions(threshold=sys.maxsize)
import h5py
import configuration as conf

oldcol = 'sum_en'
newcol = 'sum_en_smooth'

def getBinContent(df, bin1, bin2, var):
    try:
        bin_orig = df.loc[(bin1,bin2)]
        content = bin_orig[var]
    except KeyError:
        content = 0.
    return content

def setBinContent(val, df, bin1, bin2, var, zerovar):
    """
    Sets a single cell (bin1,bin2) with a value for the smoothed energy.
    If the bin does not exist its "unsmoothed" energy should be first set to zero
    (otherwise it will be NaN, which will lead to issues during smoothing).
    """
    try:
        bin_orig = df.loc[(bin1,bin2)]
    except KeyError:
        df.at[ ( bin1,bin2 ), zerovar ] = 0.
    df.at[ ( bin1,bin2 ), var ] = val

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
        bin1_10pct = int(0.1 * conf.NbinsRz)
        R1_10pct = conf.MinROverZ + bin1_10pct * (conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz
        R2_10pct = R1_10pct + ((conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz)
        area_10pct_ = ((np.pi * (R2_10pct**2 - R1_10pct**2)) / conf.NbinsPhi)
        area = area * area_10pct_;

    # loop per chunk of (equal Rz) rows with a common shift to speedup
    # unforutnately np.roll's 'shift' argument must be the same for different rows
    for idx in np.unique(nBinsSide):
        roll_indices = np.where(nBinsSide == idx)[0]
        arr_copy = arr[roll_indices,:]
        arr_smooth = arr[roll_indices,:]
        for nside in range(1, int(idx)+1):
            arr_smooth += ( np.roll(arr_copy, shift=nside,  axis=1) / (2**nside) +
                            np.roll(arr_copy, shift=-nside, axis=1) / (2**nside) )
        arr_new[roll_indices,:] = arr_smooth / np.expand_dims(area[roll_indices], axis=-1)

    return arr_new * conf.areaPerTriggerCell

def removeOldColumn(df):
    #remove old energy sum and rename
    df[oldcol] = df[newcol]
    df = df.drop(columns=[newcol])
    return df

def printHistogram(arr):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if arr[i,j] == 0:
                print('-', end='|')
            else:
                print('X', end='|')
        print()

def createHistogram(event):
    arr = np.zeros((conf.NbinsRz, conf.NbinsPhi))

    for ev in event[:]:
        rzbin = int(ev[0]-1)
        phibin = int(ev[1]-1)
        arr[rzbin,phibin] = ev[2]

    return arr

# Event by event smoothing
storeIn  = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.FillingOut), mode='r')
storeOut = h5py.File( os.path.join(os.environ['PWD'], conf.DataFolder, conf.SmoothingOut), mode='w')

for falgo in conf.FesAlgos:
    keys = [x for x in storeIn.keys() if falgo in x]
    for key in keys:
        ev = createHistogram( storeIn[key] )
        #printHistogram(ev)
        ev = smoothAlongPhi(ev)
        #printHistogram(ev)
        ev = smoothAlongRz(ev)
        #printHistogram(ev)

        storeOut[key] = ev

storeIn.close()
storeOut.close()

# convert bin ids back to values (central values in the bin) 
# splittedClusters_tc['Rz'] = binConv(splittedClusters_tc['Rz_bin'], conf.BinDistRz, conf.MinROverZ)
# splittedClusters_tc['tc_phi'] = binConv(splittedClusters_tc['tc_phi_bin'], conf.BinDistPhi, -np.pi)
