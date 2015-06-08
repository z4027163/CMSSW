import FWCore.ParameterSet.Config as cms

process = cms.Process("L1MicroGMTEmulator")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))

process.dumper = cms.EDAnalyzer("l1t::uGMTLUTDumper",
    out_directory = cms.string("lut_dump"),
    SortRankLUTSettings = cms.PSet (
        pT_in_width = cms.int32(9), 
        qual_in_width = cms.int32(4), 
        out_width = cms.int32(10),
        filename = cms.string(""),
     )
)

process.dumpPath = cms.Path( process.dumper )
process.schedule = cms.Schedule(process.dumpPath)