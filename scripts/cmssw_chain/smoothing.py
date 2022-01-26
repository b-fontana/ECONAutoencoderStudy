import os
import pandas as pd
import configuration as conf

def smoothAlongPhi(df):
    # for bin1 in rannge(conf.NbinsRz):
    #     # Take into account edges with only one side up or down
    #     weight = 1.5 if (bin1 == 0 or bin1 == conf.NbinsRz - 1) else 2.

    #     for bin2 in range(conf.NbinsPhi):
    #         bin_orig = df[ (df.Rz_bin==bin1) & (df.tc_phi_bin==bin2) ]
    #         content = 0. if bin_orig.empty else bin_orig.sum_en
    #         contentDown = bin1 > 0 ? histoClusters.at(z_side, bin1 - 1, bin2).values[Bin::Content::Sum] : 0;
    #         contentUp = bin1 < (nBins1_ - 1) ? histoClusters.at(z_side, bin1 + 1, bin2).values[Bin::Content::Sum] : 0;

    #     auto& bin = histoSumRPhiClusters.at(z_side, bin1, bin2);
    #     bin.values[Bin::Content::Sum] = (content + 0.5 * contentDown + 0.5 * contentUp) / weight;
    #     bin.weighted_x = bin_orig.weighted_x;
    #     bin.weighted_y = bin_orig.weighted_y;
    pass

def smoothAlongRz(df):
    """
    Smoothes the energy distribution of cluster energies deposited on trigger cells.
    The input takes already into account the same binning as applied in CMSSW's chain.
    """
    for bin1 in rannge(conf.NbinsRz):
        nBinsSide = (conf.BinSums[bin1] - 1) / 2;
        # takes into account different area of bins in different R-rings + sum of quadratic weights used
        area = (1 + 2.0 * (1 - 0.5**nBinsSide))

        if conf.SeedsNormByArea:
            R1 = conf.MinROverZ + bin1 * (conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz
            R2 = R1 + ((conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz)
            area = area * ((math.pi * (R2**2 - R1**2)) / conf.NbinsPhi);
        else:
            #compute quantities for non-normalised-by-area histoMax
            #The 0.1 factor in bin1_10pct is an attempt to keep the same rough scale for seeds.
            #The exact value is arbitrary.
            bin1_10pct = int(0.1 * conf.NbinsRz)
            R1_10pct = conf.MinROverZ + bin1_10pct * (conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz
            R2_10pct = R1_10pct + ((conf.MaxROverZ - conf.MinROverZ) / conf.NbinsRz)
            area_10pct_ = ((math.pi * (R2_10pct**2 - R1_10pct**2)) / conf.NbinsPhi)
            area = area * area_10pct_;
          
        for bin2 in range(conf.NbinsPhi):
            bin_orig = df[ (df.Rz_bin==bin1) & (df.tc_phi_bin==bin2) ]
            content = 0. if bin_orig.empty else bin_orig.sum_en

            for (int bin22 = 1; bin22 <= nBinsSide; bin22++) {
          int binToSumLeft = bin2 - bin22;
          if (binToSumLeft < 0)
            binToSumLeft += nBins2_;
          unsigned binToSumRight = bin2 + bin22;
          if (binToSumRight >= nBins2_)
            binToSumRight -= nBins2_;

          content += histoClusters.at(z_side, bin1, binToSumLeft).values[Bin::Content::Sum] /
                     pow(2, bin22);  // quadratic kernel

          content += histoClusters.at(z_side, bin1, binToSumRight).values[Bin::Content::Sum] /
                     pow(2, bin22);  // quadratic kernel
        }



            
    
    
# Event by event smoothing
with pd.HDFStore( os.path.join(os.environ['PWD'], conf.DataFolder, conf.FillingOut), mode='r') as store:
    for falgo in conf.FesAlgos:
        print(store.keys())
        keys = [x for x in store.keys() if falgo in x]
        print(conf.FesAlgos, keys)
        for key in keys:
            ev = store[key]
            print(ev)
            ev = smoothAlongRz(ev)
            print(ev)
            ev = smoothAlongPhi(ev)
            print(ev)
            quit()


# convert bin ids back to values (central values in the bin) 
# splittedClusters_tc['Rz'] = binConv(splittedClusters_tc['Rz_bin'], conf.BinDistRz, conf.MinROverZ)
# splittedClusters_tc['tc_phi'] = binConv(splittedClusters_tc['tc_phi_bin'], conf.BinDistPhi, -np.pi)
    
