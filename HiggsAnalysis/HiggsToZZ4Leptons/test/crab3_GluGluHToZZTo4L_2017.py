from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = 'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8'
config.General.workArea = 'Fall17'
config.General.transferOutputs = True
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '/lustre/home/reham/TEST_ELE_cleaned/CMSSW_9_4_10/src/HiggsAnalysis/HiggsToZZ4Leptons/test/HiggsTozz_MiniAOD_mc.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['roottree_leptons.root']
config.JobType.maxMemoryMB = 8000
config.JobType.inputFiles = ['QGL_80X.db']
config.JobType.inputFiles = ['/lustre/home/reham/TEST_ELE_cleaned/CMSSW_9_4_10/src/roccor_Run2_v2/data/RoccoR2017.txt']
config.section_('Data')
config.Data.publication = False
config.Data.inputDataset = '/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM' 
# config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 1000
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'Automatic'
config.Data.outLFNDirBase = '/store/user/raly/MonoHiggs/MC2017processing'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
