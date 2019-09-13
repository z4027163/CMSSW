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
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v7', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17', '')

# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.load('HiggsAnalysis.HiggsToZZ4Leptons.bunchSpacingProducer_cfi')
#process.load('HiggsAnalysis.HiggsToZZ4Leptons.metFiltersMiniAOD_cff')

process.load('RecoMET.METFilters.metFilters_cff')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)

process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)                                                          
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.primaryVertexFilter.vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.BadPFMuonFilter.muons  = cms.InputTag("slimmedMuons")
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.muons  = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates") #
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter) # 

#ecalBadCalib
baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )


process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibReducedMINIAODFilter)


process.goodOfflinePrimaryVerticestwo = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        



from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(

process,
jetSource = cms.InputTag('slimmedJets'),
labelName = 'UpdatedJEC',
jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')

)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)


from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
                           isData=False, #(or False),
                           postfix = "TEST"
                           )



from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True,
                       eleIDModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                                       'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                      # runVID=True, #saves CPU time by not needlessly re-running VID
                      era='2017-Nov17ReReco')  



process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff')  

#@#process.calibratedPatElectrons.isMC = cms.bool(True)

process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(True)   
process.hTozzTo4leptonsCommonRootTreePresel.fillLHEinfo = cms.untracked.bool(True)
process.hTozzTo4leptonsCommonRootTreePresel.filljec = cms.untracked.bool(True)
process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)
#//@
#This variable isData to apply muon calibrator inside commonRooTree.h and get the error on muon pT
process.hTozzTo4leptonsCommonRootTreePresel.isData = cms.bool(False)

process.hTozzTo4leptonsCommonRootTreePresel.noiseFilterTag = cms.InputTag("TriggerResults","","PAT")

#for LHE informations for Jets
#process.hTozzTo4leptonsCommonRootTreePresel.LHEProduct = cms.InputTag("source") #when commented, default=externalLHEProducer
process.hTozzTo4leptonsCommonRootTreePresel.lheEventProduct = cms.InputTag("externalLHEProducer")# this inputTag depend on input mc sample 

process.genanalysis= cms.Sequence(
  process.hTozzTo4leptonsGenSequence                  *
  #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
  #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
  process.hTozzTo4leptonsMCDumper                     
 # process.hTozzTo4leptonsMCCP                         
  )


process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.goodOfflinePrimaryVerticestwo     *
  process.genanalysis *
  process.jecSequence* 
  process.fullPatMetSequenceTEST *
  process.egammaPostRecoSeq * 
  #process.hTozzTo4leptonsSelectionSequenceData * 
  process.hTozzTo4leptonsSelectionSequenceMC * 
  process.hTozzTo4leptonsMatchingSequence * # 
  process.hTozzTo4leptonsCommonRootTreePresel   
  )

#//////////////////////////////////////////////////////////////////

#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff') 
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *   
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()  
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "MC_Synch_2017hTozzToLeptons.root"  

#process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew )
process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,
                                 process.Flag_HBHENoiseFilter,
                                 process.Flag_HBHENoiseIsoFilter,
                                 process.Flag_globalSuperTightHalo2016Filter,
                                 process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                 process.Flag_goodVertices,
                                 #process.Flag_eeBadScFilter,
                                 process.Flag_BadPFMuonFilter, 
                                 #@#process.Flag_BadChargedCandidateFilter,
                                 process.Flag_ecalBadCalibFilter, 
                                 #process.ecalBadCalibReducedMINIAODFilter,
                                 process.hTozzTo4leptonsSelectionPath
                                 #process.o 
                                 )


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source ("PoolSource",
                             
  fileNames = cms.untracked.vstring(
'file:/eos/uscms/store/user/wangz/mc/RunIIFall17MiniAOD/ggZZ.root'
#'file:/eos/uscms/store/user/wangz/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DY_2017_mini.root'
#'file:/eos/uscms/store/user/wangz/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DY_mini.root'
 )#,
#skipEvents= cms.untracked.uint32(18000)
)

## # Endpath
#process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew ) 
