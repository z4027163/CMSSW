from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsernameFromSiteDB

config = config()
config.General.requestName = 'DoubleMuon-Run2017D-17Nov2017-v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HiggsTozz_MiniAOD_data.py'
config.JobType.outputFiles = ['roottree_leptons.root']

config.Data.inputDataset = '/DoubleMuon/Run2017D-17Nov2017-v1/MINIAOD'
#'/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD'


config.JobType.allowUndistributedCMSSW = True

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.lumiMask = 'crab_projects/crab_DoubleEG-Run2016B-03Feb2017-noMuCal-part1/results/notFinishedLumis.json'
config.Data.inputDBS = 'global'
config.Data.ignoreLocality = False
#config.Site.whitelist = ["T2_US*"]
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 4900
config.Data.outLFNDirBase = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'DoubleMuon-Run2017D-17Nov2017-v1'
#config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.section_('User')
config.section_('Site')

config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.blacklist = ['T3_US_Rutgers']
