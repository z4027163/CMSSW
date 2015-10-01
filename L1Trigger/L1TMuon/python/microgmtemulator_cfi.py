import FWCore.ParameterSet.Config as cms

import os

l1tgmt_basedir = "L1Trigger/L1TMuon/"
lut_dir = os.path.join(l1tgmt_basedir, "data/microgmt_luts/")

microGMTEmulator = cms.EDProducer('l1t::MicroGMTEmulator',
    barrelTFInput = cms.InputTag("MicroGMTInputProducer", "BarrelTFMuons"),
    overlapTFInput = cms.InputTag("MicroGMTInputProducer", "OverlapTFMuons"),
    forwardTFInput = cms.InputTag("MicroGMTInputProducer", "ForwardTFMuons"),
    triggerTowerInput = cms.InputTag("MicroGMTInputProducer", "TriggerTowerSums"),
)

