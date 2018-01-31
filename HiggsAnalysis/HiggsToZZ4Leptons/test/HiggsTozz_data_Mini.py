
import FWCore.ParameterSet.Config as cms

process = cms.Process('HZZdata')

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
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')

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
process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)

# filters
filters = []
# met filters
#if options.runMetFilter:
# run all the time and store result
print 'Preparing MET filters'
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
hltFilter = hltHighLevel.clone()
# PAT if miniaod by itself (MC) and RECO if at the same time as reco (data)
hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'PAT')
hltFilter.throw = cms.bool(True)
# ICHEP recommendation
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_8011_ICHEP_dataset
for flag in ['HBHENoiseFilter','HBHENoiseIsoFilter','EcalDeadCellTriggerPrimitiveFilter','goodVertices','eeBadScFilter','globalTightHalo2016Filter']:
    mod = hltFilter.clone(HLTPaths=cms.vstring('Flag_{0}'.format(flag)))
    modName = 'filter{0}'.format(flag)
    setattr(process,modName,mod)
    filters += [getattr(process,modName)]
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
filters += [process.BadChargedCandidateFilter]

process.metFilter = cms.Sequence()

for f in filters:
    process.metFilter += f

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 

process.calibratedPatElectrons.isMC = cms.bool(False)
process.calibratedPatPhotons.isMC = cms.bool(False)

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

process.hTozzTo4leptonsPFfsrPhoton.src = cms.InputTag("packedPFCandidates","","RECO")

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

jetCorr = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None')


updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrections = jetCorr,
)

process.hTozzTo4leptonsPFJetSelector.PFJetCollection = cms.InputTag("updatedPatJets")

process.load('PhysicsTools/PatAlgos/producersLayer1/jetUpdater_cff')
process.jecSequence = cms.Sequence( process.updatedPatJetCorrFactors * process.updatedPatJets)

## Following lines are for default MET for Type1 corrections.
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# If you only want to re-correct for JEC and get the proper uncertainties for the default MET
runMetCorAndUncFromMiniAOD(process,
                          isData=True,
                         )

# Now you are creating the bad muon corrected MET
process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)

from PhysicsTools.PatUtils.tools.muonRecoMitigation import muonRecoMitigation

muonRecoMitigation(
                       process = process,
                       pfCandCollection = "packedPFCandidates", #input PF Candidate Collection
                       runOnMiniAOD = True, #To determine if you are running on AOD or MiniAOD
                       selection="", #You can use a custom selection for your bad muons. Leave empty if you would like to use the bad muon recipe definition.
                       muonCollection="", #The muon collection name where your custom selection will be applied to. Leave empty if you would like to use the bad muon recipe definition.
                       cleanCollName="cleanMuonsPFCandidates", #output pf candidate collection ame
                       cleaningScheme="computeAllApplyClone", #Options are: "all", "computeAllApplyBad","computeAllApplyClone". Decides which (or both) bad muon collections to be used for MET cleaning coming from the bad muon recipe.
                       postfix="" #Use if you would like to add a post fix to your muon / pf collections
                       )

runMetCorAndUncFromMiniAOD(process,
                           isData=True,
                           pfCandColl="cleanMuonsPFCandidates",
                           recoMetFromPFCs=True,
                           postfix="MuClean"
                           )


process.mucorMET = cms.Sequence(                     
        process.badGlobalMuonTaggerMAOD *
        process.cloneGlobalMuonTaggerMAOD *
      #process.badMuons * # If you are using cleaning mode "all", uncomment this line
        process.cleanMuonsPFCandidates *
        process.fullPatMetSequenceMuClean
        )

 # Now you are creating the e/g corrected MET on top of the bad muon corrected MET (on re-miniaod)
from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
corMETFromMuonAndEG(process,
                    pfCandCollection="", #not needed         
                    electronCollection="slimmedElectronsBeforeGSFix",
                    photonCollection="slimmedPhotonsBeforeGSFix",
                    corElectronCollection="slimmedElectrons",
                    corPhotonCollection="slimmedPhotons",
                    allMETEGCorrected=True,
                    muCorrection=False,
                    eGCorrection=True,
                    runOnMiniAOD=True,
                    postfix="MuEGClean"
                    )
process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
del process.slimmedMETsMuEGClean.caloMET
 
     # If you are running in the scheduled mode:
process.egcorrMET = cms.Sequence(
        process.cleanedPhotonsMuEGClean+process.cleanedCorPhotonsMuEGClean+
        process.matchedPhotonsMuEGClean + process.matchedElectronsMuEGClean +
        process.corMETPhotonMuEGClean+process.corMETElectronMuEGClean+
        process.patPFMetT1MuEGClean+process.patPFMetRawMuEGClean+
        process.patPFMetT1SmearMuEGClean+process.patPFMetT1TxyMuEGClean+
        process.patPFMetTxyMuEGClean+process.patPFMetT1JetEnUpMuEGClean+
        process.patPFMetT1JetResUpMuEGClean+process.patPFMetT1SmearJetResUpMuEGClean+
        process.patPFMetT1ElectronEnUpMuEGClean+process.patPFMetT1PhotonEnUpMuEGClean+
        process.patPFMetT1MuonEnUpMuEGClean+process.patPFMetT1TauEnUpMuEGClean+
        process.patPFMetT1UnclusteredEnUpMuEGClean+process.patPFMetT1JetEnDownMuEGClean+
        process.patPFMetT1JetResDownMuEGClean+process.patPFMetT1SmearJetResDownMuEGClean+
        process.patPFMetT1ElectronEnDownMuEGClean+process.patPFMetT1PhotonEnDownMuEGClean+
        process.patPFMetT1MuonEnDownMuEGClean+process.patPFMetT1TauEnDownMuEGClean+
        process.patPFMetT1UnclusteredEnDownMuEGClean+process.slimmedMETsMuEGClean)

process.hTozzTo4leptonsCommonRootTreePresel.PfMETLabel = cms.InputTag("slimmedMETsMuEGClean","","HZZdata")

process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.goodOfflinePrimaryVertices     *
  process.metFilter *
  process.jecSequence *
  process.mucorMET  *
  process.fullPatMetSequence * # If you are re-correctign the default MET
  process.egcorrMET  *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsCommonRootTreePresel
  )


#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"

process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,
                                 process.hTozzTo4leptonsSelectionPath )


#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring(
#'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/DoubleMuon/RunII2016/data_H2.root'
#'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/DoubleEG/RunII2016/data_B.root'
'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/SingleMuon/SingleMuon_2016_data.root'
  )

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)


## # Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
