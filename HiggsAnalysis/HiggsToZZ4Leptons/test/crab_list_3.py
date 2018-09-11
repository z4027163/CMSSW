from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

submitVersion = "QCD_inclusive"

mainOutputDir = '/store/user/wangz/%s' % submitVersion

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = 'HiggsTozz_MiniAOD.py'
config.Data.allowNonValidInputDataset = False
config.JobType.sendExternalFolder     = True
config.JobType.outputFiles = ['roottree_leptons.root']

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 

config.Site.storageSite = 'T3_US_FNALLPC'



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    ##### submit MC
    config.Data.outLFNDirBase = '%s' % mainOutputDir
    config.Data.splitting     = 'LumiBased'
    config.Data.unitsPerJob   = 20
#    config.Data.totalUnits =29000

    config.General.requestName  = 'WWTo2L2Nu_13TeV-powheg_v1' 
    config.Data.inputDataset    = '/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.Data.outputDatasetTag= 'WWTo2L2Nu_13TeV-powheg_v1' 
    submit(config)

 #   sys.exit(0)
    config.General.requestName  = 'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v1'
    config.Data.inputDataset    = '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.Data.outputDatasetTag = 'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v1'
    submit(config)

    config.General.requestName  = 'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v1'    
    config.Data.inputDataset    = '/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.Data.outputDatasetTag= 'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v1'
    submit(config)

 #   sys.exit(0)

