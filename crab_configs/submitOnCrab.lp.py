
import os
import sys
import re

# NB: you need to source crab environment before lauching this script:
# source /cvmfs/cms.cern.ch/crab3/crab.sh

###################################################################
#### Parameters to be changed for each production

datasetsFile = "samples.txt"

#PROCESS = ["ELECTRON_200PU"]
#tag = "ELE_200PU_TRUNC120_v4"

PROCESS = ["PHOTON_0PU"]
tag = "PHOT_0PU_TRUNC120_v4"

#PROCESS = ["MINBIAS_200PU"]
#tag = "MINB_200PU_TRUNC"

############################################################################################################

FastJobs = False # controls number of jobs 
VeryLong = False # controls time for each job
PublishDataset = False # publish dataset; set to false if producing ntuples

###################################################################
#### Automated script starting

# dataset block definition
comment = "#"
sectionBeginEnd = "==="

# check if file with dataset exist
if not os.path.isfile(datasetsFile):
    print "File %s not found!!!" % datasetsFile
    sys.exit()

#check if directory exists
crabJobsFolder = "crab3_" + tag
if os.path.isdir(crabJobsFolder):
    print "Folder %s already exists, please change tag name or delete it" % crabJobsFolder
    sys.exit()

# grep all datasets names, skip lines with # as a comment
# block between === * === are "sections" to be processed

currSection = ""
dtsetToLaunch = []

print " =========  Starting submission on CRAB ========"
print " Parameters: "
print " PROCESS: "
for pr in PROCESS: print "   * " , pr
print " tag: " , tag
print " Fast jobs?: " , FastJobs
print " Publish?: "   , PublishDataset

# READ INPUT FILE
with open(datasetsFile) as fIn:
    for line in fIn:
        line = line.strip() # remove newline at the end and leding/trailing whitespaces

        if not line: #skip empty lines
            continue

        if comment in line:
            continue

        #print line
        words = line.split()
        if len(words) >= 3:
            if words[0] == sectionBeginEnd and words[2] == sectionBeginEnd:
                currSection = words[1]
        else:
            if currSection in PROCESS:
                dtsetToLaunch.append(line)

# CREATE CRAB JOBS
os.system ("voms-proxy-init -voms cms")

for name in PROCESS: crabJobsFolder + "_" + name
print crabJobsFolder
os.system ("mkdir %s" % crabJobsFolder)

counter = 1 # appended to the request name to avoid overlaps between datasets with same name e.g. /DoubleEG/Run2015B-17Jul2015-v1/MINIAOD vs /DoubleEG/Run2015B-PromptReco-v1/MINIAOD
outlog = open ((crabJobsFolder + "/submissionLog.txt"), "w")
outlog.write (" =========  Starting submission on CRAB ========\n")
outlog.write (" Parameters: \n")
outlog.write (" PROCESS: \n")
for pr in PROCESS: outlog.write ("   * %s\n" % pr)
outlog.write (" tag: %s\n" % tag)
outlog.write (" Fast jobs?: %s\n" % str(FastJobs))
outlog.write (" Publish?: %s\n"   % str(PublishDataset))
outlog.write (" ===============================================\n\n\n")

print dtsetToLaunch
for dtset in dtsetToLaunch:
    dtsetNames = dtset
    if '/MINIAODSIM' in dtset:
        dtsetNames = dtset.replace('/MINIAODSIM', "")
    elif '/MINIAOD' in dtset:
        dtsetNames = dtset.replace('/MINIAOD', "")
    elif '/GEN-SIM-DIGI-RAW-MINIAOD' in dtset:
        dtsetNames = dtset.replace('/GEN-SIM-DIGI-RAW-MINIAOD', "")
    elif '/FEVT' in dtset:
        dtsetNames = dtset.replace('/FEVT', "")
        

    dtsetNames = dtsetNames.replace('/', "__")
    dtsetNames = dtsetNames.strip("__") # remove leading and trailing double __
    shortName = dtset.split('/')[1]

    if (len(shortName) > 95): # requestName not exceed 100 Characters!
        toRemove = len (shortName) - 95
        shortName = shortName[toRemove:]

    #dtSetName = dtsetNames[1]
    command = "crab submit -c crab3_template.lp.py"

    if 'MinBias' in dtset:
        command += " JobType.psetName=produce_ntuple_s1truncation_reduced_pt5_V11_cfg.py"
    else:
        command += " JobType.psetName=produce_ntuple_s1truncation_reduced_genmatch_V11_cfg.py"

    #if 'MinBias' in dtset:
    #    command += " JobType.psetName=produce_seedingnoarea_reduced_pt5_v11_cfg.py"
    #else:
    #    command += " JobType.psetName=produce_seedingnoarea_reduced_genmatch_V11_cfg.py"

    command += " General.requestName=%s" % (shortName + "_" + str(counter))
    command += " General.requestName=%s" % (shortName + "_" + str(counter))
    command += " General.workArea=%s" % crabJobsFolder
    #command += " Data.inputDataset=%s" % dtset
    
    #lst = []
    #for i in open('phot0pu.txt').readlines():
    #    lst.append(i.replace('\n',''))

    #command += " Data.userInputFiles=open('phot0pu.txt').readlines()"
    command += " Data.outLFNDirBase=/store/user/lportale/HGCAL_L1/%s" % (tag )
    command += " Data.outputDatasetTag=%s" % (shortName + "_" + tag + "_" + str(counter))
    if not PublishDataset : command += " Data.publication=False" # cannot publish flat root ntuples
    if VeryLong           : command += " JobType.maxJobRuntimeMin=2500" # 32 hours, default is 22 hours -- can do up to 2800 hrs
    print command ,  "\n"
    outlog.write(command + "\n\n")
    os.system (command)
    #outlog.write(command + "\n\n")
    counter = counter + 1
