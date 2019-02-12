import FWCore.ParameterSet.Config as cms

process = cms.Process('MonoHiggs')

# Complete Preselection Sequence for 4l analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/Geometry/GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')#reham
process.load('Configuration.StandardSequences.MagneticField_cff') #reham
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration/EventContent/EventContent_cff')


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '') # 

# Random generator 
#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#    calibratedPatElectrons = cms.PSet(
#        initialSeed = cms.untracked.uint32(1),
#        engineName = cms.untracked.string('TRandom3')
#    )
#)

process.load('HiggsAnalysis.HiggsToZZ4Leptons.bunchSpacingProducer_cfi')
#process.load('HiggsAnalysis.HiggsToZZ4Leptons.metFiltersMiniAOD_cff')

process.load('RecoMET.METFilters.metFilters_cff')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')

process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)

process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)                                                         
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.primaryVertexFilter.vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)                                    
#process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter) 
process.BadPFMuonFilter.muons  = cms.InputTag("slimmedMuons")
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.muons  = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates") #reham
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter) # Reham 




process.goodOfflinePrimaryVerticestwo = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMuonRochesterCalibrator_cfi')
process.hTozzTo4leptonsMuonRochesterCalibrator.isData = cms.bool(True)
process.hTozzTo4leptonsMuonRochesterCalibrator.MCTruth = cms.bool(False)



from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(

process,
jetSource = cms.InputTag('slimmedJets'),
labelName = 'UpdatedJEC',
jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')

)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)



from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
                           isData=True, #(or False),
                           postfix = "TEST"
                           )

#/////////////////////////////////////////////////////


from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       eleIDModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                                       'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                       runVID=True, #saves CPU time by not needlessly re-running VID
                      era='2017-Nov17ReReco')  


#/////////////////////////////////////////////////////


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 
#@#process.calibratedPatElectrons.isMC = cms.bool(False) 


process.hTozzTo4leptonsPFfsrPhoton.src = cms.InputTag("packedPFCandidates")
process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = False
process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
  #process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(False)    
process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)
#//@
#This variable isData to apply muon calibrator inside commonRooTree.h and get the error on muon pT
process.hTozzTo4leptonsCommonRootTreePresel.isData = cms.bool(True)

process.hTozzTo4leptonsCommonRootTreePresel.noiseFilterTag = cms.InputTag("TriggerResults","","RECO")

#for MC only but put need to run and not crash
process.hTozzTo4leptonsCommonRootTreePresel.LHEProduct = cms.InputTag("externalLHEProducer")# this inputTag depend on input mc sample 

process.genanalysis= cms.Sequence(
  # process.hTozzTo4leptonsGenSequence                  *
  #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
  #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
  process.hTozzTo4leptonsMCDumper                     
 # process.hTozzTo4leptonsMCCP                         
  )

process.hTozzTo4leptonsSelectionPath = cms.Path(
    process.goodOfflinePrimaryVerticestwo     *
    #  process.genanalysis * 
    process.jecSequence 
    process.fullPatMetSequenceTEST * 
    process.egammaPostRecoSeq * 
    process.hTozzTo4leptonsSelectionSequenceData *
    process.hTozzTo4leptonsCommonRootTreePresel 
    )

#///////////////////////////////////////////

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *   
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()  
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "Data_2017_DoubleMuon_RunB_hTozzToLeptons.root" 

process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew ) 
process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,                            
                                 #@#process.Flag_HBHENoiseFilter,
                                 #@#process.Flag_HBHENoiseIsoFilter,
                                 #@#process.Flag_globalSuperTightHalo2016Filter,
                                 #@#process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                 #@#process.Flag_goodVertices,
                                # process.Flag_eeBadScFilter,
                                 #@#process.Flag_BadPFMuonFilter,
                                 #@#process.Flag_BadChargedCandidateFilter,
                                 #process.Flag_ecalBadCalibFilter,  
                                 process.hTozzTo4leptonsSelectionPath,
                                 process.o)



process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source ("PoolSource",
                             
  fileNames = cms.untracked.vstring(
'file:data_DoubleMuon_2017_RunB_0852E0CB-E7D7-E711-B2DA-0025905C3DCE.root' 
#'file:Data_2017_DoubleMuon_RunB_hTozzToLeptons.root'
  )
)

## # Endpath
#process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew ) 
