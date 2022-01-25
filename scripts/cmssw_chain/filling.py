import os
import random
random.seed(18) # fix seed for reproducibility
import argparse
import numpy as np
import pandas as pd
import uproot as up

from bokeh.io import export_png
from bokeh.plotting import figure

import configuration as conf
from utils import calculateRoverZfromEta

binConv = lambda vals,dist,amin : (vals*dist) + (dist/2) + amin

"""
Fills split clusters information according to the Stage2 FPGA fixed binning.
"""

### Data Extraction ####################################################
simDataPath = os.path.join(os.environ['PWD'], conf.DataFolder, conf.FillingIn)
simAlgoDFs, simAlgoFiles, simAlgoPlots = ({} for _ in range(3))
fes = ['Threshold']
for fe in fes:
    simAlgoFiles[fe] = [ os.path.join(simDataPath) ]

for fe,files in simAlgoFiles.items():
    name = fe
    dfs = []
    for file in files:
        with pd.HDFStore(file, mode='r') as store:
            dfs.append(store[name])
    simAlgoDFs[fe] = pd.concat(dfs)

simAlgoNames = sorted(simAlgoDFs.keys())
if conf.Debug:
    print('Input HDF5 keys:')
    print(simAlgoNames)

### Data Processing ######################################################
enrescuts = [-0.35]
assert(len(enrescuts)==len(fes))
for i,(fe,cut) in enumerate(zip(fes,enrescuts)):
    df = simAlgoDFs[fe]

    if conf.Debug:
        print('Cluster level information:')
        print(df.filter(regex='cl3d_*.'))

    df = df[ (df['genpart_exeta']>1.7) & (df['genpart_exeta']<2.8) ]
    df = df[ df['cl3d_eta']>0 ]
    df['enres'] = ( df['cl3d_energy']-df['genpart_energy'] ) / df['genpart_energy']

    #### Energy resolution histogram ####
    hist, edges = np.histogram(df['enres'], density=True, bins=150)

    p = figure( width=500, height=300, title='Energy Resolution: ' + fe,
                y_axis_type="log")
    virtualmin = 1e-4 #avoid log scale issues
    p.quad(top=hist, bottom=virtualmin, left=edges[:-1], right=edges[1:],
           fill_color="navy", line_color="white", alpha=0.7)
    p.line(x=[cut,cut], y=[virtualmin,max(hist)], line_color="#ff8888", line_width=4, alpha=0.9, legend_label="Cut")

    outname = os.path.join('out', 'filling_enrescut.png')
    export_png(p, filename=outname)
    ######################################

    nansel = pd.isna(df['enres']) 
    nandf = df[nansel]
    nandf['enres'] = 1.1
    df = df[~nansel]
    df = pd.concat([df,nandf], sort=False)

    # select events with splitted clusters
    splittedClusters = df[ df['enres'] < cut ]

    # random pick some events (fixing the seed for reproducibility)
    _events_remaining = list(splittedClusters.index.unique())
    _events_sample = random.sample(_events_remaining, conf.Nevents)
    splittedClusters = splittedClusters.loc[_events_sample]
    #splittedClusters.sample(n=NEVENTS, replace=False, random_state=8)

    if conf.Debug:
        print('SplitClusters Dataset: event random selection')
        print(splittedClusters)
        print(splittedClusters.columns)

    #splitting remaining data into cluster and tc to avoid tc data duplication
    _cl3d_vars = [x for x in splittedClusters.columns.to_list() if 'tc_' not in x]
    splittedClusters_3d = splittedClusters[_cl3d_vars]
    splittedClusters_3d = splittedClusters_3d.reset_index()
    _tc_vars = [x for x in splittedClusters.columns.to_list() if 'cl3d' not in x]

    #trigger cells info is repeated across clusters in the same event
    splittedClusters_tc = splittedClusters.groupby("event").head(1)[_tc_vars] #first() instead of head(1) also works

    _tc_vars = [x for x in _tc_vars if 'tc_' in x]
    splittedClusters_tc = splittedClusters_tc.explode( _tc_vars )

    for v in _tc_vars:
        splittedClusters_tc[v] = splittedClusters_tc[v].astype(np.float64)

    splittedClusters_tc['Rz'] = np.sqrt(splittedClusters_tc.tc_x*splittedClusters_tc.tc_x + splittedClusters_tc.tc_y*splittedClusters_tc.tc_y)  / abs(splittedClusters_tc.tc_z)
    splittedClusters_tc = splittedClusters_tc[ (splittedClusters_tc['Rz'] < conf.MaxROverZ) & (splittedClusters_tc.Rz > conf.MinROverZ) ]
    splittedClusters_tc = splittedClusters_tc.reset_index()
    splittedClusters_tc['Rz_bin'] = pd.cut( splittedClusters_tc.Rz, bins=conf.RzBinEdges, labels=False )
    splittedClusters_tc['tc_phi_bin'] = pd.cut( splittedClusters_tc['tc_phi'], bins=conf.PhiBinEdges, labels=False )

    #convert bin ids back to values (central values in the bin)
    splittedClusters_tc['Rz'] = binConv(splittedClusters_tc['Rz_bin'], conf.BinDistRz, conf.MinROverZ)
    splittedClusters_tc['tc_phi'] = binConv(splittedClusters_tc['tc_phi_bin'], conf.BinDistPhi, -np.pi)
    splittedClusters_tc.drop(['Rz_bin', 'tc_phi_bin'], axis=1)

    simAlgoPlots[fe] = (splittedClusters_3d, splittedClusters_tc)

### Event Processing ######################################################
with pd.HDFStore( os.path.join(os.environ['PWD'], conf.DataFolder, conf.FillingOut), mode='w') as store:

    for i,(_k,(df_3d,df_tc)) in enumerate(simAlgoPlots.items()):
        for ev in df_tc['event'].unique():
            ev_tc = df_tc[ df_tc.event == ev ]
            ev_3d = df_3d[ df_3d.event == ev ]
            _simCols_tc = ['tc_mipPt', 'tc_z', 'tc_phi', 'tc_eta',
                           'Rz', 'genpart_exeta', 'genpart_exphi']
            ev_tc = ev_tc.filter(items=_simCols_tc)
            ev_3d['cl3d_Roverz'] = calculateRoverZfromEta(ev_3d['cl3d_eta'])
            ev_3d['gen_Roverz']  = calculateRoverZfromEta(ev_3d['genpart_exeta'])

            cl3d_pos_rz, cl3d_pos_phi = ev_3d['cl3d_Roverz'].unique(), ev_3d['cl3d_phi'].unique()
            gen_pos_rz, gen_pos_phi = ev_3d['gen_Roverz'].unique(), ev_3d['genpart_exphi'].unique()
            ev_3d = ev_3d.drop(['cl3d_Roverz', 'cl3d_eta', 'cl3d_phi'], axis=1)
            assert( len(gen_pos_rz) == 1 and len(gen_pos_phi) == 1 )

            groupby = ev_tc.groupby(['Rz', 'tc_phi'], as_index=False)
            group = groupby.count()

            energy_sum = groupby.sum()['tc_mipPt']
            eta_mins = groupby.min()['tc_eta']
            eta_maxs = groupby.max()['tc_eta']

            group = group.rename(columns={'tc_z': 'nhits'}, errors='raise')
            group.insert(0, 'min_eta', eta_mins)
            group.insert(0, 'max_eta', eta_maxs)
            group.insert(0, 'sum_en', energy_sum)

            store[str(_k) + '_' + str(ev)] = group
