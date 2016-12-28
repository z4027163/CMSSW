from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'pp_jj_13TeV_vtest1'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'PrivateMC'
config.JobType.generator = 'lhe'
config.JobType.psetName = 'step1.py'
config.JobType.inputFiles = ['jj.lhe']

config.Data.primaryDataset = 'MinBias'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 50
NJOBS = 1
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.outLFN = '/store/user/wangz/' # or '/store/group/<subdir>'
config.Data.publication = True
config.Data.publishDataName = 'pp_jj_13TeV_vtest1'
config.JobType.allowUndistributedCMSSW = True

config.Site.storageSite = 'T3_US_FNALLPC'

