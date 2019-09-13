from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys
config = config()

submitVersion = "JetHT2017"

mainOutputDir = '/store/user/wangz/%s' % submitVersion

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = 'HiggsTozz_MiniAOD_data.py'
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

    config.General.requestName  = 'JetHT_Run2017B-31Mar2018-v1'
    config.Data.inputDataset    = '/JetHT/Run2017B-31Mar2018-v1/MINIAOD'
    config.Data.outputDatasetTag= 'JetHT_Run2017B-31Mar2018-v1'
    submit(config)

    config.General.requestName  = 'JetHT_Run2017C-31Mar2018-v1'
    config.Data.inputDataset    = '/JetHT/Run2017C-31Mar2018-v1/MINIAOD'
    config.Data.outputDatasetTag= 'JetHT_Run2017C-31Mar2018-v1'
    submit(config)

    config.General.requestName  = 'JetHT_Run2017D-31Mar2018-v1'
    config.Data.inputDataset    = '/JetHT/Run2017D-31Mar2018-v1/MINIAOD'
    config.Data.outputDatasetTag= 'JetHT_Run2017D-31Mar2018-v1'
    submit(config)

    config.General.requestName  = 'JetHT_Run2017E-31Mar2018-v1'
    config.Data.inputDataset    = '/JetHT/Run2017E-31Mar2018-v1/MINIAOD'
    config.Data.outputDatasetTag= 'JetHT_Run2017E-31Mar2018-v1'
    submit(config)

    config.General.requestName  = 'JetHT_Run2017F-31Mar2018-v1'
    config.Data.inputDataset    = '/JetHT/Run2017F-31Mar2018-v1/MINIAOD'
    config.Data.outputDatasetTag= 'JetHT_Run2017F-31Mar2018-v1'
    submit(config)

