import FWCore.ParameterSet.Config as cms

L1Validator = cms.EDAnalyzer('L1Validator',
  dirName=cms.string("L1T/L1TriggerVsGen"),
#  fileName=cms.string("L1Validation.root") #output file name
  GenSource=cms.InputTag("genParticles"),
  L1ExtraIsoEGSource=cms.InputTag("l1extraParticles", "Isolated"),
  L1ExtraNonIsoEGSource=cms.InputTag("l1extraParticles", "NonIsolated"),
  L1ExtraCenJetSource=cms.InputTag("l1extraParticles", "Central"),
  L1ExtraForJetSource=cms.InputTag("l1extraParticles", "Forward"),
  L1ExtraTauJetSource=cms.InputTag("l1extraParticles", "Tau"),
  srcToken = cms.InputTag("generator"),
  L1MuonBXSource=cms.InputTag("gmtStage2Digis", "Muon"),
  L1EGammaBXSource=cms.InputTag("caloStage2Digis", "EGamma"),
  L1TauBXSource=cms.InputTag("caloStage2Digis", "Tau"),
  L1JetBXSource=cms.InputTag("caloStage2Digis", "Jet"),
  L1ExtraMuonSource=cms.InputTag("l1extraParticles"),
  L1GenJetSource=cms.InputTag("ak4GenJets","")
)
