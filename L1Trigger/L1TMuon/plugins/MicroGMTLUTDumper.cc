// -*- C++ -*-
//
// Package:    MicroGMTLUTDumper
// Class:      MicroGMTLUTDumper
//
/**\class MicroGMTLUTDumper MicroGMTLUTDumper.cc L1Trigger/L1TGlobalMuon/plugins/MicroGMTLUTDumper.cc

 Description: Takes txt-file input and produces barrel- / overlap- / forward TF muons

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joschka Philip Lingemann,40 3-B01,+41227671598,
//         Created:  Thu Oct  3 10:12:30 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTRankPtQualLUT.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTMatchQualLUT.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTLUTFactories.h"

#include "CondFormats/L1TObjects/interface/MicroGMTParams.h"
#include "CondFormats/DataRecord/interface/L1TMicroGMTParamsRcd.h"

#include <iostream>
//
// class declaration
//
namespace l1t {
class MicroGMTLUTDumper : public edm::EDAnalyzer {
   public:
      explicit MicroGMTLUTDumper(const edm::ParameterSet&);
      ~MicroGMTLUTDumper();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

   private:
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);

      void dumpLut(MicroGMTLUT*, const std::string&);

      // ----------member data ---------------------------
      MicroGMTParams* microGMTParams;
      std::string m_foldername;
      std::shared_ptr<MicroGMTRankPtQualLUT> m_rankLUT;

      std::shared_ptr<MicroGMTMatchQualLUT> m_boPosMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_boNegMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_foPosMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_foNegMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_brlSingleMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_ovlPosSingleMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_ovlNegSingleMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_fwdPosSingleMatchQualLUT;
      std::shared_ptr<MicroGMTMatchQualLUT> m_fwdNegSingleMatchQualLUT;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
MicroGMTLUTDumper::MicroGMTLUTDumper(const edm::ParameterSet& iConfig)
{
  //now do what ever other initialization is needed
  m_foldername = iConfig.getParameter<std::string> ("out_directory");
}


MicroGMTLUTDumper::~MicroGMTLUTDumper()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
void
MicroGMTLUTDumper::dumpLut(MicroGMTLUT* lut, const std::string& oName) {
  std::ofstream fStream(m_foldername+oName);
  lut->save(fStream);
  fStream.close();
}



// ------------ method called to produce the data  ------------
void
MicroGMTLUTDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  dumpLut(m_rankLUT.get(), std::string("/rank_lut.json"));
  dumpLut(m_boPosMatchQualLUT.get(), std::string("/boPosMatchQualLUT.json"));
  dumpLut(m_boNegMatchQualLUT.get(), std::string("/boNegMatchQualLUT.json"));
  dumpLut(m_foPosMatchQualLUT.get(), std::string("/foPosMatchQualLUT.json"));
  dumpLut(m_foNegMatchQualLUT.get(), std::string("/foNegMatchQualLUT.json"));
  dumpLut(m_brlSingleMatchQualLUT.get(), std::string("/brlSingleMatchQualLUT.json"));
  dumpLut(m_ovlPosSingleMatchQualLUT.get(), std::string("/ovlPosSingleMatchQualLUT.json"));
  dumpLut(m_ovlNegSingleMatchQualLUT.get(), std::string("/ovlNegSingleMatchQualLUT.json"));
  dumpLut(m_fwdPosSingleMatchQualLUT.get(), std::string("/fwdPosSingleMatchQualLUT.json"));
  dumpLut(m_fwdNegSingleMatchQualLUT.get(), std::string("/fwdNegSingleMatchQualLUT.json"));

}

// ------------ method called when starting to processes a run  ------------
void
MicroGMTLUTDumper::beginRun(edm::Run const& run, edm::EventSetup const& iSetup)
{
  const L1TMicroGMTParamsRcd& microGMTParamsRcd = iSetup.get<L1TMicroGMTParamsRcd>();
  edm::ESHandle<MicroGMTParams> microGMTParamsHandle;
  microGMTParamsRcd.get(microGMTParamsHandle);

  delete microGMTParams;
  microGMTParams = new (microGMTParams) MicroGMTParams(*microGMTParamsHandle.product());
  if (!microGMTParams) {
    edm::LogError("L1TMicroGMTLUTDumper") << "Could not retrieve parameters from Event Setup" << std::endl;
  }

  int fwVersion = microGMTParams->fwVersion();
  m_rankLUT = MicroGMTRankPtQualLUTFactory::create(microGMTParams->sortRankLUTParams()->filename(), fwVersion);
  m_boPosMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->bOPosMatchQualLUTParams()->filename(), cancel_t::omtf_bmtf_pos, fwVersion);
  m_boNegMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->bONegMatchQualLUTParams()->filename(), cancel_t::omtf_bmtf_neg, fwVersion);
  m_foPosMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->fOPosMatchQualLUTParams()->filename(), cancel_t::omtf_emtf_pos, fwVersion);
  m_foNegMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->fONegMatchQualLUTParams()->filename(), cancel_t::omtf_emtf_neg, fwVersion);
  m_brlSingleMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->brlSingleMatchQualLUTParams()->filename(), cancel_t::bmtf_bmtf, fwVersion);
  m_ovlPosSingleMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->ovlPosSingleMatchQualLUTParams()->filename(), cancel_t::omtf_omtf_pos, fwVersion);
  m_ovlNegSingleMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->ovlNegSingleMatchQualLUTParams()->filename(), cancel_t::omtf_omtf_neg, fwVersion);
  m_fwdPosSingleMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->fwdPosSingleMatchQualLUTParams()->filename(), cancel_t::emtf_emtf_pos, fwVersion);
  m_fwdNegSingleMatchQualLUT = MicroGMTMatchQualLUTFactory::create(microGMTParams->fwdNegSingleMatchQualLUTParams()->filename(), cancel_t::emtf_emtf_neg, fwVersion);
}

} // namespace l1t
//define this as a plug-in
DEFINE_FWK_MODULE(l1t::MicroGMTLUTDumper);
