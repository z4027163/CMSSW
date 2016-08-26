from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsernameFromSiteDB

config = config()
config.General.requestName = 'DYToLLNtuplerRunIISpring2015_test'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HiggsTozz_MiniAOD.py'
config.JobType.outputFiles = ['roottree_leptons.root']

#config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/MINIAODSIM'

config.JobType.allowUndistributedCMSSW = True


config.Data.inputDBS = 'global'
config.Data.ignoreLocality = True
config.Site.whitelist = ["T2_US*"]
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 100
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'DYTOLL_13TeV_NTUPLER_v1_test'
config.Data.userInputFiles = list(open('DYLL.txt'))
#config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.section_('User')
config.section_('Site')

config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.blacklist = ['T3_US_Rutgers']
