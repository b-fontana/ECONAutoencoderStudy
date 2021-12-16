import os
import argparse
import numpy as np
import pandas as pd
import uproot as up

from bokeh.io import output_file, show
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource,
                          LogColorMapper, LogTicker,
                          LinearColorMapper, BasicTicker,
                          PrintfTickFormatter,
                          Range1d,
                          Panel, Tabs)

from bokeh.plotting import figure
from bokeh.sampledata.unemployment1948 import data as testdata
from bokeh.transform import transform
from bokeh.palettes import viridis as _palette

import sys
sys.path.append( os.path.join(os.environ['HOME'], 'bokehplot') )
import bokehplot as bkp

import sys
sys.path.append( os.path.join(os.environ['HOME'], 'Documents/FPGAs') )
from utils.pandas import display_everything
#display_everything()

class dotDict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

#unicodes
UC = dotDict( dict(phi='\u03A6',) )

def set_figure_props(p, bincenters):
    p.axis.axis_line_color = 'black'
    p.axis.major_tick_line_color = 'black'
    p.axis.major_label_text_font_size = '10px'
    p.axis.major_label_standoff = 2
    #p.xaxis.major_label_orientation = 1.0
    p.xaxis.axis_label = r"$$\color{black} \phi$$"
    p.yaxis.axis_label = '$$R/z$$'
    p.yaxis.major_label_overrides = {dig: label for dig,label in enumerate(bincenters)}
    
    p.hover.tooltips = [
        ("#hits", "@{nhits}"),
    ]

#########################################################################
################### PARSER ############################
#########################################################################
parser = argparse.ArgumentParser(description='Plot trigger cells occupancy.')
parser.add_argument('-p', '--nphibins', help='number of uniform phi bins',
                    default=216, type=int)
parser.add_argument('-z', '--pos_endcap', help='Use only the positive endcap.',
                    default=True, type=bool)
parser.add_argument('--minROverZ', help='Minimum value of R/z, as defined in CMSSW.',
                    default=0.076, type=float)
parser.add_argument('--maxROverZ', help='Maximum value of R/z, as defined in CMSSW.',
                    default=0.58, type=float)
parser.add_argument('-d', '--debug', help='debug mode', action='store_true')
parser.add_argument('-l', '--log', help='use color log scale', action='store_true')

FLAGS = parser.parse_args()

#########################################################################
################### CONFIGURATION PARAMETERS ############################
#########################################################################
NEVENTS = 10
NBINSRZ = 42
NBINSPHI = FLAGS.nphibins

filename = 'gen_cl3d_tc.hdf5'
base_path = os.environ['PWD']
data_path = os.path.join(base_path, 'data')

vnames = dotDict( dict( RoverZ = 'Rz',
                        phi = 'tc_phi',
                        eta = 'tc_eta',
                        nhits = 'nhits',
                       )
                  )

algo_files = {}
fes = ['Threshold']
for fe in fes:
    algo_files[fe] = [ os.path.join(data_path, filename) ]

algos_dfs = {}
for fe,files in algo_files.items():
    name = fe
    dfs = []
    for file in files:
        with pd.HDFStore(file, mode='r') as store:
            dfs.append(store[name])
    algos_dfs[fe] = pd.concat(dfs)
algo_names = sorted(algos_dfs.keys())

b = bkp.BokehPlot( os.path.join(base_path, 'plots', 'clustersDistribution.html'),
                   nfigs=3, nframes=len(fes) )

#########################################################################
################### PLOT CONFIGURATION PARAMETERS #######################
#########################################################################
mypalette = _palette(40)    
shifth, shiftv = 1, 1
binwidth=1
title = ''
#title += '; Min(R/z)={} and Max(R/z)={})'.format(FLAGS.minROverZ, FLAGS.maxROverZ)

#########################################################################
################### DATA ANALYSIS #######################################
#########################################################################
enrescuts = [-0.3]
df_plots = {}
assert(len(enrescuts)==len(fes))
for i,(fe,cut) in enumerate(zip(fes,enrescuts)):
    df = algos_dfs[fe]
    #df = df[ (df['genpart_exeta']>1.7) & (df['genpart_exeta']<2.8) ]
    df = df[ df['cl3d_eta']>0 ]
    df['enres'] = ( df['cl3d_energy']-df['genpart_energy'] ) / df['genpart_energy']

    nansel = pd.isna(df['enres']) 
    nandf = df[nansel]
    nandf['enres'] = 1.1
    df = df[~nansel]
    df = pd.concat([df,nandf], sort=False)

    # select events with splitted clusters
    split = df[ df['enres'] < cut ]

    b.histogram(idx=0, iframe=i, data=np.histogram(df['enres'], bins=500), color='orange', style='v%0.8%red',
                fig_kwargs={'x.axis_label': r"\[ \frac{E_{Cl3d} - E_{Gen}}{E_{Gen}}\]",
                            'y.axis_label': 'Counts'})
    b.histogram(idx=1, iframe=i, data=np.histogram(split['enres'], bins=500), color='orange', style='v%0.8%red')
    b.histogram(idx=2, iframe=i, data=np.histogram(split['cl3d_eta'], bins=500), color='green', style='v%0.8%red')

    # random pick some events (fixing the seed for reproducibility)
    split = split.sample(n=NEVENTS, replace=False, random_state=9)

    #b.graph(idx=1, data=[np.arange(1,11),np.arange(1,20,2)], style='vbar%"%2%orange', line=True)
    #b.histogram(idx=2, data=np.histogram2d(arr[:,0],arr[:,1],bins=50), style='quad%Viridis')
    b.save_frame(show=False)

    if FLAGS.debug:
        print('Split Dataset: event random selection')
        print(split)
        print(split.columns)

    tc_vars = list(filter(lambda x: x.startswith('tc_'), split.columns.to_list()))
    split = split.explode( tc_vars )
    for v in tc_vars:
        split[v] = split[v].astype(np.float64)

    split[vnames.RoverZ] = np.sqrt(split.tc_x*split.tc_x + split.tc_y*split.tc_y) / abs(split.tc_z)

    split = split[ (split[vnames.RoverZ] < FLAGS.maxROverZ) & (split[vnames.RoverZ] > FLAGS.minROverZ) ]

    rzBinEdges = np.linspace( split[vnames.RoverZ].min(), split[vnames.RoverZ].max(), num=NBINSRZ )
    rzBinCenters = ['{:.2f}'.format(x) for x in ( rzBinEdges[1:] + rzBinEdges[:-1] ) / 2 ]

    split = split.reset_index()
    split[vnames.RoverZ] = pd.cut( split[vnames.RoverZ], bins=rzBinEdges, labels=False )
    split[vnames.phi] = pd.cut( split[vnames.phi], bins=NBINSPHI, labels=False )

    df_plots[fe] = split
    
#########################################################################
################### GROUPING AND PLOTTING ###############################
#########################################################################
tabs = []
for i,dfp in enumerate(df_plots.values()):
    for ev in dfp['event'].unique():
        df_tmp = dfp[ dfp.event == ev ]
        print(df_tmp.columns)
        quit()

        #CHANGE THE LINE BELOW TO GET A WEIGHTED HISTOGRAM
        group = df_tmp.groupby([vnames.RoverZ, vnames.phi], as_index=False).count()
        group = group.rename(columns={'tc_z': vnames.nhits})
        group = group[[vnames.phi, vnames.nhits, vnames.RoverZ]]

        source = ColumnDataSource(group)

        p = figure(width=1800, height=800, title=title,
                   x_range=Range1d(group[vnames.phi].min()-shifth,
                                   group[vnames.phi].max()+shifth),
                   y_range=Range1d(group[vnames.RoverZ].min()-shiftv, group[vnames.RoverZ].max()+shiftv),
                   tools="hover,box_select,box_zoom,reset,save", x_axis_location='below',
                   x_axis_type='linear', y_axis_type='linear',
                   )

        if FLAGS.log:
            mapper = LogColorMapper(palette=mypalette, low=group[vnames.nhits].min(), high=group[vnames.nhits].max())
        else:
            mapper = LinearColorMapper(palette=mypalette, low=group[vnames.nhits].min(), high=group[vnames.nhits].max())

        color_bar = ColorBar(color_mapper=mapper,
                             ticker= ( LogTicker(desired_num_ticks=len(mypalette))
                                       if FLAGS.log else BasicTicker(desired_num_ticks=int(len(mypalette)/4)) ),
                             formatter=PrintfTickFormatter(format="%d")
                             )
        p.add_layout(color_bar, 'right')

        p.rect( x=vnames.phi, y=vnames.RoverZ,
                source=source,
                width=binwidth, height=binwidth,
                width_units='data', height_units='data',
                line_color='black', fill_color=transform(vnames.nhits, mapper)
               )

        set_figure_props(p, rzBinCenters)
        
        tabs.append( Panel(child=p, title='Event {}'.format(ev)) )


if not FLAGS.debug:
    show( Tabs(tabs=tabs) )
