import os
import argparse
import numpy as np
import pandas as pd
import uproot as up

import bokehplot as bkp
from bokeh.io import output_file, show
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource,
                          LogColorMapper, PrintfTickFormatter,
                          LogTicker,
                          Range1d)
from bokeh.plotting import figure
from bokeh.sampledata.unemployment1948 import data as testdata
from bokeh.transform import transform
from bokeh.palettes import viridis as _palette

class dotDict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

parser = argparse.ArgumentParser(description='Plot trigger cells occupancy.')
parser.add_argument('-p', '--nphibins', help='number of uniform phi bins',
                    default=60, type=int)
parser.add_argument('-d', '--debug', help='debug mode', action='store_true')
FLAGS = parser.parse_args()
 
#########################################################################
################### CONFIGURATION PARAMETERS ############################
#########################################################################

NBINSRZ = 30
NBINSPHI = FLAGS.nphibins

data_path = '~/Downloads/hadd.root'
myfile = up.open(data_path)

folder = 'FloatingpointThresholdDummyHistomaxnoareath20Genclustersntuple'
tree_name = 'HGCalTriggerNtuple'
tree = myfile[ os.path.join(folder, tree_name) ]
if FLAGS.debug:
    for k in tree.keys():
        print(k)

#unicodes
UC = dotDict( dict(phi='\u03A6',
                   )
              )

variables = {'run', 'event', 'lumi',
             'gen_eta', 'gen_phi', 'gen_pt', 'gen_energy', 'gen_charge', 'gen_pdgid',
             'cl3d_pt', 'cl3d_energy', 'cl3d_eta', 'cl3d_phi',
             'vtx_z'
             }
assert(variables.issubset(tree.keys()))
variables = sorted(list(variables))

#########################################################################
################### DATA ANALYSIS #######################################
#########################################################################

df = tree.arrays(variables, library='pd')[0]

#df = df[ df.cl3d_eta > 0 ] #only look at positive endcap
# df = df.drop(['zside'], axis=1)
# variables.remove('zside')
print(df.head(20))
quit()
#########################################################################
################### PLOTTING ############################################
#########################################################################
b = bkp.BokehPlot('triggerCellsOccup.html',
                  nfigs=len(variables), nframes=1,
                  nwidgets=0,
                  fig_width=300, fig_height=200)

fkw = {'toolbar.logo': None,
       'y.axis_label': 'Counts',
       }
for i,var in enumerate(variables):
    print(data[var])
    fkw.update({'x.axis_label': var})
    b.histogram(idx=i, iframe=0,
                data=np.histogram(data[var].to_numpy(), bins=100),
                color='orange', style='v%0.8%red',
                fig_kwargs=fkw)

b.save_frame(nrows=4, ncols=5, show=True)
