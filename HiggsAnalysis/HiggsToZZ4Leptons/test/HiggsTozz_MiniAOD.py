
import FWCore.ParameterSet.Config as cms

process = cms.Process('MonoHiggs')

# Complete Preselection Sequence for 4l analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/Geometry/GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration/EventContent/EventContent_cff')


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')

# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.load('HiggsAnalysis.HiggsToZZ4Leptons.bunchSpacingProducer_cfi')
process.load('HiggsAnalysis.HiggsToZZ4Leptons.metFiltersMiniAOD_cff')
process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)
## process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)                                                                                                             
## process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)                                                                                       
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
## process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)   
## process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)                                                                                                         
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
## process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter) 
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
## process.Flag_trackingFailureFilter = cms.Path(process.goodVertices + process.trackingFailureFilter)                                                                                
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
## process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)                                                                                                           
## process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)                                                                                                                       
## process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)                                                                             
## proces..Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)                                                                                                             
## and the sub-filters                                                                                                                                                                
# process.Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)                                                                                                                  
# process.Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)                                                                                                            
# process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)                     

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMuonCalibrator_cfi')
process.hTozzTo4leptonsMuonCalibrator.isData = cms.bool(False) 

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.calibratedElectrons.isMC = cms.bool(True)

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 

process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
  #process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(False)    
process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)

process.genanalysis= cms.Sequence(
  process.hTozzTo4leptonsGenSequence                  *
  #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
  #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
  process.hTozzTo4leptonsMCDumper                     *                
  process.hTozzTo4leptonsMCCP                         )

process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.goodOfflinePrimaryVertices     *
  process.genanalysis *
  process.hTozzTo4leptonsSelectionSequenceData *
#  process.hTozzTo4leptonsMatchingSequence *
  process.hTozzTo4leptonsCommonRootTreePresel
  )


#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"

process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,
                                 process.Flag_HBHENoiseFilter,
                                 process.Flag_HBHENoiseIsoFilter,
          ###                       process.Flag_CSCTightHalo2015Filter,
                                 process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                 process.Flag_goodVertices,
                                 process.Flag_eeBadScFilter,
                                 process.hTozzTo4leptonsSelectionPath )


#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring(
'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/DY_mini.root'
#'root://cmsxrootd-site.fnal.gov//store/mc/RunIIFall15MiniAODv2/WW_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/0EE6E0B8-FFB7-E511-B0C6-0025904CDDFA.root'
#'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/MINIAOD_Spring15_test3.root'
#'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/MINIAOD_DY_2.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall15MiniAODv1/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/30000/00445FD9-BC9D-E511-A14D-003048D4DF6C.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2015C_25ns/DoubleMuon/MINIAOD/16Dec2015-v1/20000/081A3AE2-ABB5-E511-9A0D-7845C4FC368C.root'
#'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/MINIAODSIM/AVE30BX50_tsg_PHYS14_ST_V1-v1/30000/0080FDC4-5A8B-E411-AA4A-00259073E4F0.root',
  )

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

## # Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
