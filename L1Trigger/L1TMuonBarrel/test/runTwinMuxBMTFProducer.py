# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import commands

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(50)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.source = cms.Source(
    'PoolSource',
fileNames = cms.untracked.vstring('file:///afs/cern.ch/work/g/gflouris/public/SingleMuPt6180_noanti_50k_eta08.root')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

# PostLS1 geometry used
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
############################
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

####BMTF Emulator
process.load('L1Trigger.L1TMuonBarrel.bmtfDigis_cfi')
process.bmtfEmulator = cms.EDProducer("BMTrackFinder",
   DTDigi_Source = cms.InputTag("L1TTwinMuxProducer"),
   DTDigi_Theta_Source = cms.InputTag("simDtTriggerPrimitiveDigis"),
   Debug = cms.untracked.int32(0)
)
####TwinMux Emulator
process.load('L1Trigger.L1TMuonBarrel.L1TTwinMuxProducer_cfi')


process.L1TMuonSeq = cms.Sequence( process.L1TTwinMuxProducer +
                                   process.bmtfEmulator 
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

process.out = cms.OutputModule("PoolOutputModule", 
   fileName = cms.untracked.string("l1tbmtf_superprimitives.root"),
                               )

process.output_step = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.L1TMuonPath)
process.schedule.extend([process.output_step])
