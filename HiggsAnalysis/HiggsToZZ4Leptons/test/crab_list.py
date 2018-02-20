from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

submitVersion = "GluGluTo4L"

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
    config.Data.unitsPerJob   = 10

    config.General.requestName  = 'GluGluToContinToZZTo4e_13TeV_DefaultShower_MCFM701_pythia8_v1'
    config.Data.inputDataset    = '/GluGluToContinToZZTo4e_13TeV_DefaultShower_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM'
    config.Data.outputDatasetTag= 'GluGluToContinToZZTo4e_13TeV_DefaultShower_MCFM701_pythia8_v1'
    submit(config)

 #   sys.exit(0)
    config.General.requestName  = 'GluGluToContinToZZTo4mu_13TeV_DefaultShower_MCFM701_pythia8_v1'
    config.Data.inputDataset    = '/GluGluToContinToZZTo4mu_13TeV_DefaultShower_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM'
    config.Data.outputDatasetTag = 'GluGluToContinToZZTo4mu_13TeV_DefaultShower_MCFM701_pythia8_v1'
    submit(config)

    config.General.requestName  = 'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_v1'
    config.Data.inputDataset    = '/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    config.Data.outputDatasetTag = 'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_v1' 
    submit(config)

    config.General.requestName  = 'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_v1'
    config.Data.inputDataset    = '/GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM' 
    config.Data.outputDatasetTag = 'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_v1'
    submit(config)

    config.General.requestName  = 'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_v1'
    config.Data.inputDataset    = '/GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM' 
    config.Data.outputDatasetTag = 'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_v1'
    submit(config)

    config.General.requestName  = 'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_v1'
    config.Data.inputDataset    = '/GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM' 
    config.Data.outputDatasetTag = 'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_v1'
    submit(config)

