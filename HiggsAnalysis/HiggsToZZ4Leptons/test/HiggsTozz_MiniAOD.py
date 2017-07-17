
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
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                  calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(8675389),
                                                      engineName = cms.untracked.string('TRandom3'),
                                                      ),
                  calibratedPatPhotons    = cms.PSet( initialSeed = cms.untracked.uint32(8675389),
                                                      engineName = cms.untracked.string('TRandom3'),
                                                      ),   
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

process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
process.calibratedPatElectrons.isMC = cms.bool(True)
process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
process.calibratedPatPhotons.isMC = cms.bool(True)

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff')

process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
  #process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(True)    
process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)


process.genanalysis= cms.Sequence(
  process.hTozzTo4leptonsGenSequence                  *
  #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
  #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
  process.hTozzTo4leptonsMCDumper                     *                
  process.hTozzTo4leptonsMCCP                         )

#process.hTozzTo4leptonsSelectionSequenceData.remove(process.hTozzTo4leptonsHLTAnalysisFilter)

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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/00E54BE4-21E5-E611-BD4D-0025905A60B6.root'
#'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/DY_mini.root'
#'file:/eos/uscms/store/user/wangz/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DY_mini.root'
'file:/eos/uscms/store/user/wangz/mc/RunIISpring16MiniAODv1/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/ZZ_mini.root'
#'file:/eos/uscms/store/user/wangz/data/Spring16/TT_mini1.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/00E54BE4-21E5-E611-BD4D-0025905A60B6.root'
  )

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

## # Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
