from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

submitVersion = "QCD_inclusive_2017"

mainOutputDir = '/store/user/wangz/%s' % submitVersion

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = 'HiggsTozz_MiniAOD_mc.py'
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
    config.Data.totalUnits = 9990

    config.General.requestName  = 'QCD_Pt_30to50_TuneCP5_13TeV_pythia8_v1' 
    config.Data.inputDataset    = '/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM' 
    config.Data.outputDatasetTag= 'QCD_Pt_30to50_TuneCP5_13TeV_pythia8_v1'
    submit(config)

 #   sys.exit(0)
    config.General.requestName  = 'QCD_Pt_50to80_TuneCP5_13TeV_pythia8_v1'    
    config.Data.inputDataset    = '/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM' 
    config.Data.outputDatasetTag= 'QCD_Pt_50to80_TuneCP5_13TeV_pythia8_v1'
    submit(config)

 #   sys.exit(0)
    config.General.requestName  = 'QCD_Pt_80to120_TuneCP5_13TeV_pythia8_v1'
    config.Data.inputDataset    = '/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM' 
    config.Data.outputDatasetTag = 'QCD_Pt_80to120_TuneCP5_13TeV_pythia8_v1'
    submit(config)

    config.General.requestName  = 'QCD_Pt_120to170_TuneCP5_13TeV_pythia8_v1'    
    config.Data.inputDataset    = '/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' 
    config.Data.outputDatasetTag= 'QCD_Pt_120to170_TuneCP5_13TeV_pythia8_v1'
    submit(config)

    config.General.requestName  = 'QCD_Pt_170to300_TuneCP5_13TeV_pythia8_v1'                  
    config.Data.inputDataset    = '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM' 
    config.Data.outputDatasetTag= 'QCD_Pt_170to300_TuneCP5_13TeV_pythia8_v1'
    submit(config)

    config.General.requestName  = 'QCD_Pt_300to470_TuneCP5_13TeV_pythia8_v1'
    config.Data.inputDataset    = '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v3/MINIAODSIM' 
    config.Data.outputDatasetTag= 'QCD_Pt_300to470_TuneCP5_13TeV_pythia8_v1'
    submit(config)

