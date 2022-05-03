#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import tqdm
import sys
import uproot # uproot4
from datetime import date
import optparse
from itertools import chain

workdir=os.getcwd()

disconnectedTriggerLayers = [
    2,
    4,
    6,
    8,
    10,
    12,
    14,
    16,
    18,
    20,
    22,
    24,
    26,
    28
]

def deltar(df):
    df['deta']=df['cl3d_eta']-df['genpart_exeta']
    df['dphi']=np.abs(df['cl3d_phi']-df['genpart_exphi'])
    sel=df['dphi']>np.pi
    df['dphi']-=sel*(2*np.pi)
    return(np.sqrt(df['dphi']*df['dphi']+df['deta']*df['deta']))
    
def matching(event):
    if event.matches.sum()==0:
        return event.cl3d_pt==event.cl3d_pt.max()
    else:
        cond_a = event.matches==True
        cond_b = event.cl3d_pt==event[cond_a].cl3d_pt.max()
        return (cond_a&cond_b)

def create_dataframes(files, algo_trees, gen_tree, p):
    print('Input files: {}'.format(files))

    branches_gen = [ 'event', 'genpart_reachedEE', 'genpart_pid', 'genpart_gen',
                     'genpart_exphi', 'genpart_exeta', 'genpart_energy' ]
    branches_cl3d = [ 'event', 'cl3d_energy','cl3d_pt','cl3d_eta','cl3d_phi' ]
    branches_tc = [ 'event', 'tc_zside', 'tc_energy', 'tc_mipPt', 'tc_pt', 'tc_layer',
                    'tc_x', 'tc_y', 'tc_z', 'tc_phi', 'tc_eta', 'tc_id' ]

    batches_gen, batches_tc = ([] for _ in range(2))
    memsize_gen, memsize_tc = '128 MB', '64 MB'
    for filename in files:
        with uproot.open(filename + ':' + gen_tree) as data:
            #print( data.num_entries_for(memsize, expressions=branches_tc) )
            for ib,batch in enumerate(data.iterate(branches_gen, step_size=memsize_gen,
                                                   library='pd')):

                # reachedEE=2: photons that hit HGCAL
                batch = batch[ batch['genpart_reachedEE']==p.reachedEE ]
                batch = batch[ batch['genpart_gen']!=-1 ]
                batch = batch[ batch['genpart_pid']==22 ]
                #batch = batch.drop(columns=['genpart_reachedEE', 'genpart_gen', 'genpart_pid'])
                batch = batch.drop(columns=['genpart_gen', 'genpart_pid'])
                batch = batch[ batch['genpart_exeta']>0  ] #positive endcap only
                batch.set_index('event', inplace=True)

                batches_gen.append(batch)
                print('Step {}: +{} generated data processed.'.format(ib,memsize_gen))

            for ib,batch in enumerate(data.iterate(branches_tc, step_size=memsize_tc,
                                                   library='pd')):

                batch = batch[ batch['tc_zside']==1 ] #positive endcap
                batch = batch.drop(columns=['tc_zside'])
                #remove layers not read by trigger cells
                batch = batch[ ~batch['tc_layer'].isin(disconnectedTriggerLayers) ]
                #convert all the trigger cell hits in each event to a list
                batch = batch.groupby(by=['event']).aggregate(lambda x: list(x))
                batches_tc.append(batch)
                print('Step {}: +{} trigger cells data processed.'.format(ib,memsize_tc))

    df_gen = pd.concat(batches_gen)
    df_tc = pd.concat(batches_tc)
    
    """ This concatenation looks like:
       genpart_exphi  genpart_exeta  genpart_energy  tc_layer       tc_x        tc_y        tc_z
event                                                                                           
5           0.972961       -1.75371       47.971615      29.0  63.632179  112.073677 -367.699005
5           0.972961       -1.75371       47.971615      21.0 -56.870522 -106.310326  351.802765
5           0.972961       -1.75371       47.971615      19.0 -63.112049 -107.511497  348.832764
5           0.972961       -1.75371       47.971615       3.0 -67.323067  -95.413071  325.072754
(...)
6
(...)
    """

    df_algos = {}
    for algo_name, algo_tree in algo_trees.items():
        with uproot.open(filename)[algo_tree] as tree:
            df_algos[algo_name] = tree.arrays(branches_cl3d, library='pd')
            breakpoint()
            # Trick to read layers pTs, which is a vector of vector
            df_algos[algo_name]['cl3d_layer_pt'] = list(chain.from_iterable(tree.arrays(['cl3d_layer_pt'])[b'cl3d_layer_pt'].tolist()))

    return (df_gen, df_algos, df_tc)

def preprocessing(param):
    files            = param.files_photons
    threshold        = param.threshold
    gen_tree         = param.gen_tree
    algo_trees       = param.algo_trees
    output_file_name = param.output_file_name
    bestmatch_only   = param.bestmatch_only
    reachedEE        = param.reachedEE

    gen, algo, tc = create_dataframes(files, algo_trees, gen_tree, param)

    algo_clean={}

    # split df_gen_clean in two, one collection for each endcap
    #gen_neg = gen_clean[ gen_clean['genpart_exeta']<=0 ]
    #gen_pos = gen[ gen['genpart_exeta']>0  ]
    #df_gen = df_gen.join(pd.concat(batches_tc), how='left', rsuffix='_tc')
    
    for algo_name,df_algo in algo.items():
        # split clusters in two, one collection for each endcap
        algo_pos = df_algo[ df_algo['cl3d_eta']>0  ]
        
        #algo_neg = df_algo[ df_algo['cl3d_eta']<=0 ]

        #set the indices
        algo_pos.set_index('event', inplace=True)
        #algo_neg.set_index('event', inplace=True)

        #merging gen columns and cluster columns, keeping cluster duplicates (same event)
        algo_pos_merged=gen.join(algo_pos, how='right', rsuffix='_algo').dropna()

        #algo_neg_merged=gen_neg.join(algo_neg, how='left', rsuffix='_algo')

        # compute deltar
        algo_pos_merged['deltar']=deltar(algo_pos_merged)
        #algo_neg_merged['deltar']=deltar(algo_neg_merged)

        #could be better:
        algo_pos_merged['matches'] = algo_pos_merged.deltar<=threshold
        #algo_neg_merged['matches'] = algo_neg_merged.deltar<=threshold

        #matching
        # /!\ LP: but then, we want to remove only clusters that aren't "best match"
        #         best match could be:
        #              - Unmatched cluster with highest pT if no dr-matched cluster in evt
        #              - Matched cluster with highest pT *among dr-matched clusters*
        group=algo_pos_merged.groupby('event') # required when dealing with pile-up
        algo_pos_merged['best_match']=group.apply(matching).array
        #group=algo_neg_merged.groupby('event')
        #algo_neg_merged['best_match']=group.apply(matching).array

        #keep matched clusters only
        if bestmatch_only:
            sel=algo_pos_merged['best_match']==True
            algo_pos_merged=algo_pos_merged[sel]

            #sel=algo_neg_merged['best_match']==True
            #algo_neg_merged=algo_neg_merged[sel]

        #algo_clean[algo_name]=pd.concat([algo_neg_merged,algo_pos_merged], sort=False).sort_values('event')
        algo_clean[algo_name] = algo_pos_merged.sort_values('event')
        algo_clean[algo_name] = algo_clean[algo_name].join(tc, how='left', rsuffix='_tc')
        
        print(algo_name, algo_clean[algo_name].shape[0])


    #save files to savedir in HDF
    store = pd.HDFStore( os.path.join(param.output_dir, 'data', output_file_name), mode='w')
    for algo_name, df in algo_clean.items():
        store[algo_name] = df
    store.close()
        
#Run with: `python scripts/matching_v2.py --cfg scripts.custom_params`
if __name__=='__main__':
    parser = optparse.OptionParser()
    parser.add_option("--cfg",type="string", dest="params", help="select the path to the parameters file")
   
    (opt, args) = parser.parse_args()

    # Loading configuration parameters
    import importlib
    import sys
    current_dir = os.getcwd();
    sys.path.append(current_dir)
    param=importlib.import_module(opt.params)

    preprocessing(param)
