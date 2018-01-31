
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
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8', '')

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
process.calibratedPatElectrons.isMC = cms.bool(True)
process.calibratedPatPhotons.isMC = cms.bool(True)


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

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

jetCorr = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')


updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrections = jetCorr,
)

process.hTozzTo4leptonsPFJetSelector.PFJetCollection = cms.InputTag("updatedPatJets")

process.genanalysis= cms.Sequence(
  process.hTozzTo4leptonsGenSequence                  *
  #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
  #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
  process.hTozzTo4leptonsMCDumper                     *                
  process.hTozzTo4leptonsMCCP                         )

#process.hTozzTo4leptonsSelectionSequenceData.remove(process.hTozzTo4leptonsHLTAnalysisFilter)

process.load('PhysicsTools/PatAlgos/producersLayer1/jetUpdater_cff')
process.jecSequence = cms.Sequence( process.updatedPatJetCorrFactors * process.updatedPatJets)

## Following lines are for default MET for Type1 corrections.
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# If you only want to re-correct for JEC and get the proper uncertainties for the default MET
runMetCorAndUncFromMiniAOD(process,
                          isData=False,
                         )

process.hTozzTo4leptonsCommonRootTreePresel.PfMETLabel = cms.InputTag("slimmedMETs","","MonoHiggs")


process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.goodOfflinePrimaryVertices     *
  process.genanalysis *
  process.metFilter *
  process.jecSequence *
  process.fullPatMetSequence * # If you are re-correctign the default MET
  process.hTozzTo4leptonsSelectionSequenceData *
#  process.hTozzTo4leptonsMatchingSequence *
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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/00E54BE4-21E5-E611-BD4D-0025905A60B6.root'
#'file:/eos/uscms/store/user/wangz/MinBias/Higgs_hzz_13TeV_mini_Summer16_v1/171121_214250/0000/step3_999.root'
#'root://cmsxrootd.fnal.gov///store/user/wangz/MinBias/Higgs_hzz_13TeV_mini_Summer16_v1/171121_214250/0000/step3_1.root'
#'root://cmsxrootd-site.fnal.gov//store/user/wangz/data/DY_mini.root'
#'file:/eos/uscms/store/user/wangz/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ttz_mini.root'
#'file:/eos/uscms/store/user/wangz/mc/RunIISpring16MiniAODv1/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/TT.root'
'root://cmsxrootd.fnal.gov///store/user/wangz/MinBias/Higgs_hzz_13TeV_mini_Summer16_v1/171121_214250/0000/step3_1.root'
#'file:/eos/uscms/store/user/wangz/data/Spring16/TT_mini1.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/00E54BE4-21E5-E611-BD4D-0025905A60B6.root'
  )

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

## # Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
