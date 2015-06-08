// -*- C++ -*-
//
// Package:    uGMTLUTDumper
// Class:      uGMTLUTDumper
// 
/**\class uGMTLUTDumper uGMTLUTDumper.cc L1Trigger/L1TGlobalMuon/plugins/uGMTLUTDumper.cc

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
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTRankPtQualLUT.h"


#include <iostream>
//
// class declaration
//
namespace l1t {
class uGMTLUTDumper : public edm::EDAnalyzer {
   public:
      explicit uGMTLUTDumper(const edm::ParameterSet&);
      ~uGMTLUTDumper();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);  

   private:
      

      // ----------member data ---------------------------
      std::string m_foldername;
      MicroGMTRankPtQualLUT m_rankLUT;
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
uGMTLUTDumper::uGMTLUTDumper(const edm::ParameterSet& iConfig) : m_rankLUT(iConfig)
{
  //register your products

  //now do what ever other initialization is needed
  m_foldername = iConfig.getParameter<std::string> ("out_directory");
  

}


uGMTLUTDumper::~uGMTLUTDumper()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//



// ------------ method called to produce the data  ------------
void
uGMTLUTDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  std::ofstream rank_out(m_foldername+std::string("/rank_lut.json"));
  m_rankLUT.save(rank_out);

  // rank_out.write();
  rank_out.close();
}

} // namespace l1t
//define this as a plug-in
DEFINE_FWK_MODULE(l1t::uGMTLUTDumper);
