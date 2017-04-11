from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsernameFromSiteDB

config = config()
config.General.requestName = 'DoubleEG_Run2016B-23Sep2016-v3_part1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HiggsTozz_data_Mini.py'
config.JobType.outputFiles = ['roottree_leptons.root']

config.Data.inputDataset = '/DoubleEG/Run2016B-23Sep2016-v3/MINIAOD' 

config.JobType.allowUndistributedCMSSW = True

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
#config.Data.lumiMask = 'crab_projects/crab_DoubleMuon_Run2016D-23Sep2016-v1/results/notFinishedLumis.json'
config.Data.inputDBS = 'global'
config.Data.runRange = '273150-273600'# '273150-275376'
config.Data.ignoreLocality = False
config.Site.whitelist = ["T2_US*"]
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 10000
config.Data.outLFNDirBase = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'DoubleEG_Run2016B-23Sep2016-v3_part1'
#config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.section_('User')
config.section_('Site')

config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.blacklist = ['T3_US_Rutgers']
