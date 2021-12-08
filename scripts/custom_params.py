from glob import glob
import itertools

local = True #local machine vs server machine
htcondor = False #whether to submit the script as multiple jobs to HTCondor
if htcondor and local:
    raise ValueError('No submission is possible from your local machine!')

# DeltaR matching threshold
threshold = 0.05

# Input files
if local:
    files_photons = ['/home/bruno/Downloads/hadd.root']
    files_electrons = []
    files_pions = []

else:
    files_photons = glob('/home/llr/cms/sauvan/DATA_UPG/HGCAL/Ntuples/study_autoencoder/3_22_1/SinglePhoton_PT2to200/GammaGun_Pt2_200_PU0_HLTWinter20_std_ae_xyseed/210430_091126/ntuple*.root')
    files_electrons = []
    files_pions = []
    if not htcondor:
        files_photons = files_photons[:1] if len(files_photons)>0 else []
        files_electrons = files_electrons[:1] if len(files_electrons)>0 else []
        files_pions = files_pions[:1] if len(files_pions)>0 else []
        
# Pick one of the different algos trees to retrieve the gen information
gen_tree = 'FloatingpointThresholdDummyHistomaxnoareath20Genclustersntuple'
# Store only information on the best match
bestmatch_only = True

if local:
    output_dir = '/home/bruno/Downloads/'
else:
    output_dir = '/home/llr/cms/sauvan/DATA_UPG/HGCAL/Dataframes/study_autoencoder/3_22_1/electron_photon_signaldriven/'
    file_per_batch_electrons = 5
    file_per_batch_pions = 2
    file_per_batch_photons = 2
    
algo_trees = {}
# List of ECON algorithms
fes = ['HGCalTriggerNtuple'
       , 'Threshold', 'Mixedbcstc',
       'AutoEncoderTelescopeMSE', 'AutoEncoderStride',
       'AutoEncoderQKerasTTbar', 'AutoEncoderQKerasEle',
]
ntuple_template = 'Floatingpoint{fe}Dummy{be}GenmatchGenclustersntuple/HGCalTriggerNtuple'
algo_trees = {}
for fe in fes:
    be = 'Histomaxxydr015'
    algo_trees[fe] = ntuple_template.format(fe=fe, be=be)
