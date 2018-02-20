from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsernameFromSiteDB

config = config()
config.General.requestName = 'Higgs_hzz_13TeV-madgraph-pythia8_dm_noMuCal_v4' 
config.General.workArea = 'crab_projects'
#config.General.transferOutputs = True
#config.General.transferLogs = True



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HiggsTozz_MiniAOD.py'
config.JobType.outputFiles = ['roottree_leptons.root']

config.Data.userInputFiles = open('hzz_mini.txt').readlines()

config.JobType.allowUndistributedCMSSW = True

config.Data.outputPrimaryDataset = 'MinBias'
config.Data.ignoreLocality = False
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits =49000
config.Data.outLFNDirBase = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'Higgs_hzz_13TeV-madgraph-pythia8_dm_noMuCal_v4' 
#config.Data.outputPrimaryDataset = 'CRAB_UserFiles'

config.Site.storageSite = 'T3_US_FNALLPC'

#config.Site.ignoreGlobalBlacklist = True
