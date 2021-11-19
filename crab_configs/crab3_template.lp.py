# TEMPLATE used for automatic script submission of multiple datasets

from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")

config.General.workArea = 'DefaultCrab3Area'

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True 
config.JobType.maxMemoryMB = 2500 


config.section_("Data")

config.Data.userInputFiles=open('phot0pu.txt').readlines()
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.totalUnits = -1 #number of event
config.Data.outLFNDirBase = '/store/user/lportale/HGCAL_L1'
config.Data.publication = False
config.Data.outputDatasetTag = 'HGCAL_L1'
config.Data.allowNonValidInputDataset = True

config.section_("Site")
config.Site.storageSite = 'T2_FR_GRIF_LLR'

