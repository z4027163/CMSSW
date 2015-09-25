import FWCore.ParameterSet.Config as cms

MicroGMTInputProducer = cms.EDProducer('l1t::MicroGMTInputProducer',
    inputFileName = cms.string("../test/test.dat")
)
