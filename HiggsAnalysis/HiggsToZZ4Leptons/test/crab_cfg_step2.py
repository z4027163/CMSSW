from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'jj_13TeV_RECO_v11'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
#config.JobType.generator = 'lhe'
config.JobType.psetName = 'step2.py'

config.JobType.allowUndistributedCMSSW = True


config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Site.whitelist = ["T2_US*"]
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 1000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = True
config.Data.outputDatasetTag = 'jj_13TeV_RECO_v1'
config.Data.userInputFiles = list(open('jj1.txt'))
config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.section_('User')
config.section_('Site')

config.Site.storageSite = 'T3_US_FNALLPC'

