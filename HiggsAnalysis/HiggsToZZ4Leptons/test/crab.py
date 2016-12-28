from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsernameFromSiteDB

config = config()
config.General.requestName = 'WW_TuneCUETP8M1_13TeV-pythia8_NtuplerRunII2015all_v3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HiggsTozz_MiniAOD.py'
config.JobType.outputFiles = ['roottree_leptons.root']

config.Data.inputDataset = '/WW_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'

config.JobType.allowUndistributedCMSSW = True


config.Data.inputDBS = 'global'
config.Data.ignoreLocality = False
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 2
config.Data.outLFNDirBase = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'WW_TuneCUETP8M1_13TeV-pythia8_NtuplerRunII2015all_v3'
#config.Data.outputPrimaryDataset = 'CRAB_UserFiles'

config.Site.storageSite = 'T3_US_FNALLPC'

#config.Site.ignoreGlobalBlacklist = True
