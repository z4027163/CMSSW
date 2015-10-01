import FWCore.ParameterSet.Config as cms

import os

l1tgmt_basedir = "L1Trigger/L1TMuon/"
lut_dir = os.path.join(l1tgmt_basedir, "data/microgmt_luts/")

microGMTEmulator = cms.EDProducer('L1TMicroGMTProducer',
    barrelTFInput = cms.InputTag("L1TMicroGMTInputProducer", "BarrelTFMuons"),
    overlapTFInput = cms.InputTag("L1TMicroGMTInputProducer", "OverlapTFMuons"),
    forwardTFInput = cms.InputTag("L1TMicroGMTInputProducer", "ForwardTFMuons"),
    triggerTowerInput = cms.InputTag("L1TMicroGMTInputProducer", "TriggerTowerSums"),
)

