from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsernameFromSiteDB

config = config()
config.General.requestName = 'DoubleMuon-Run2016H-03Feb2017-noMuCal-dm-part3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True



config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HiggsTozz_data_Mini.py'
config.JobType.outputFiles = ['roottree_leptons.root']

config.Data.inputDataset = '/DoubleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'
#'/DoubleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'
#'/DoubleEG/Run2016G-03Feb2017-v1/MINIAOD'
#'/DoubleMuon/Run2016B-03Feb2017_ver3-v2/MINIAOD'#'/DoubleEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'#'/DoubleMuon/Run2016D-03Feb2017-v1/MINIAOD'#'/DoubleMuon/Run2016C-03Feb2017-v1/MINIAOD'#'/DoubleEG/Run2016C-03Feb2017-v1/MINIAOD'

config.JobType.allowUndistributedCMSSW = True

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.lumiMask = 'crab_projects/crab_DoubleEG-Run2016B-03Feb2017-noMuCal-part1/results/notFinishedLumis.json'
config.Data.inputDBS = 'global'
#config.Data.runRange = '273150-274400' 
#B# '273150-275376' #'274401-274400
#H#'281613-284035''281613-283100' '283101-284035'
#G#'278820-279800''279801-280385'
#E# '276831-277420' #C#'275801-276283' #'275656-275800'# '275656-276283'
config.Data.ignoreLocality = False
config.Site.whitelist = ["T2_US*"]
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 3
#config.Data.totalUnits = 4900
config.Data.outLFNDirBase = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'DoubleMuon-Run2016H-03Feb2017-noMuCal-dm-part3'
#config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.section_('User')
config.section_('Site')

config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.blacklist = ['T3_US_Rutgers']
