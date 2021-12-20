import os
import argparse
import numpy as np
import pandas as pd
import uproot as up

import bokehplot as bkp
from bokeh.io import output_file, show
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource,
                          LogColorMapper, LogTicker,
                          LinearColorMapper, BasicTicker,
                          PrintfTickFormatter,
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
parser.add_argument('--nphibins', help='number of uniform phi bins',
                    default=216, type=int)
parser.add_argument('--nrzbins', help='number of uniform R/z bins',
                    default=42, type=int)
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
NBINSRZ = FLAGS.nrzbins
NBINSPHI = FLAGS.nphibins

data_path = '~/Downloads/test_triggergeom.root'
myfile = up.open(data_path)
folder = 'hgcaltriggergeomtester'
tree_name = 'TreeTriggerCells'
tree = myfile[ os.path.join(folder, tree_name) ]
if FLAGS.debug:
    print('Available keys in {}:'.format(tree_name) )
    for k in tree.keys():
        print('  ' + k)

vnames = dotDict( dict( RoverZ = 'Rz',
                        phi_calc = 'phi_calc',
                        phi = 'phi',
                        eta = 'eta',
                        nhits = 'nhits',
                        min_eta = 'min_eta',
                        max_eta = 'max_eta',
                       )
                  )
variables = {'zside', 'subdet', 'layer', vnames.phi, vnames.eta, 'x', 'y', 'z'}
assert(variables.issubset(tree.keys()))
variables = list(variables)

if FLAGS.debug:
    print('Input Tree:')
    print(tree.show())
    print('\n')

data = tree.arrays(variables, library='pd')
if FLAGS.debug:
    print( data.describe() )
    
selsections = { 'layer < 5': data.layer < 5,
                'layer >= 5 and layer < 10':  (data.layer >= 5)  & (data.layer < 10),
                'layer = 10': data.layer == 10,
                'layer > 10 and layer < 15': (data.layer > 10) & (data.layer < 15),
                'layer >= 15 and layer < 20': (data.layer >= 15) & (data.layer < 20),
                'layer >= 20 and layer < 25': (data.layer >= 20) & (data.layer < 25),
                'layer >= {}'.format(data.layer.max()): (data.layer < data.layer.max()) }

#########################################################################
################### DATA ANALYSIS #######################################
#########################################################################
if FLAGS.debug:
    subdetvals  = sorted(data.subdet.unique())
    layervals  = sorted(data.layer.unique())
    etavals  = sorted(data[vnames.eta].unique())
    phivals  = sorted(data[vnames.phi].unique())

    print('Values of Trigger Subdet ({}): {}'.format(len(subdetvals), subdetvals)) #ECAL (1), HCAL silicon (2) and HCAL scintillator (10)
    print('Values of Layers ({}): {}'.format(len(layervals), layervals))
    print('Values of Trigger Eta ({}): {}'.format(len(etavals), etavals))
    print('Values of Trigger Phi ({}): {}'.format(len(phivals), phivals))

if FLAGS.pos_endcap:
    data = data[ data.zside == 1 ] #only look at positive endcap
    data = data.drop(['zside'], axis=1)
    variables.remove('zside')

data[vnames.RoverZ] = np.sqrt(data.x*data.x + data.y*data.y) / abs(data.z)
data = data[ (data[vnames.RoverZ] < FLAGS.maxROverZ) & (data[vnames.RoverZ] > FLAGS.minROverZ) ]

#rzBinEdges = np.linspace( data[vnames.RoverZ].min(), data[vnames.RoverZ].max(), num=NBINSRZ )
rzBinEdges = np.linspace( FLAGS.minROverZ, FLAGS.maxROverZ, num=NBINSRZ )
rzBinCenters = ['{:.2f}'.format(x) for x in ( rzBinEdges[1:] + rzBinEdges[:-1] ) / 2 ]

phiBinEdges = np.linspace( -np.pi, np.pi, num=NBINSPHI )
phiBinCenters = ['{:.2f}'.format(x) for x in ( phiBinEdges[1:] + phiBinEdges[:-1] ) / 2 ]

data[vnames.RoverZ] = pd.cut( data[vnames.RoverZ], bins=rzBinEdges, labels=False )
data[vnames.phi] = pd.cut( data[vnames.phi], bins=phiBinEdges, labels=False )

groups = []
for _,v in selsections.items():
    groups.append( data[v] )
    groupby = groups[-1].groupby([vnames.RoverZ, vnames.phi], as_index=False)
    groups[-1] = groupby.count()
    eta_mins = groupby.min()[vnames.eta]
    eta_maxs = groupby.max()[vnames.eta]
    groups[-1].insert(0, vnames.min_eta, eta_mins)
    groups[-1].insert(0, vnames.max_eta, eta_maxs)
    groups[-1] = groups[-1].rename(columns={'z': vnames.nhits})
    groups[-1] = groups[-1][[vnames.phi, vnames.nhits, vnames.RoverZ, vnames.min_eta, vnames.max_eta]]

#########################################################################
################### PLOTTING ############################################
#########################################################################
for isel,(selk,_) in enumerate(selsections.items()):
    source = ColumnDataSource(groups[isel])

    mypalette = _palette(40)
    if FLAGS.log:
        mapper = LogColorMapper(palette=mypalette,
                                low=groups[isel][vnames.nhits].min(), high=groups[isel][vnames.nhits].max())
    else:
        mapper = LinearColorMapper(palette=mypalette,
                                   low=groups[isel][vnames.nhits].min(), high=groups[isel][vnames.nhits].max())

    shifth, shiftv = 1, 1
    title = r'Trigger Cells Occupancy ({} vs {} bins)'.format(NBINSPHI, NBINSRZ)
    if FLAGS.pos_endcap:
        title += '; Positive end-cap only'
    title += '; {}'.format(selk)
    title += '; Min(R/z)={} and Max(R/z)={}'.format(FLAGS.minROverZ, FLAGS.maxROverZ)
    p = figure(width=1800, height=800, title=title,
               x_range=Range1d(data.phi.min()-shifth, data.phi.max()+shifth),
               y_range=Range1d(data.Rz.min()-shiftv, data.Rz.max().max()+shiftv),
               tools="hover,box_select,box_zoom,undo,redo,reset,save", x_axis_location='below',
               x_axis_type='linear', y_axis_type='linear',
               )
    binwidth=1
    p.rect( x=vnames.phi, y=vnames.RoverZ,
            source=source,
            width=binwidth, height=binwidth,
            width_units='data', height_units='data',
            line_color='black', fill_color=transform(vnames.nhits, mapper)
           )

    color_bar = ColorBar(color_mapper=mapper,
                         ticker= ( LogTicker(desired_num_ticks=len(mypalette))
                                   if FLAGS.log else BasicTicker(desired_num_ticks=int(len(mypalette)/4)) ),
                         formatter=PrintfTickFormatter(format="%d")
                         )
    p.add_layout(color_bar, 'right')

    p.axis.axis_line_color = 'black'
    p.axis.major_tick_line_color = 'black'
    p.axis.major_label_text_font_size = '10px'
    p.axis.major_label_standoff = 2
    #p.xaxis.major_label_orientation = 1.0
    p.xaxis.axis_label = r"$$\color{black} \phi$$"
    p.xaxis.major_label_overrides = {dig: label for dig,label in enumerate(phiBinCenters)}
    p.yaxis.axis_label = '$$R/z$$'
    p.yaxis.major_label_overrides = {dig: label for dig,label in enumerate(rzBinCenters)}

    p.hover.tooltips = [
        ("#hits", "@{nhits}"),
        ("min(eta)", "@{min_eta}"),
        ("max(eta)", "@{max_eta}"),
    ]

    if not FLAGS.debug:
        output_file('triggerCellsOccup_sel{}.html'.format(isel))
        show(p)
