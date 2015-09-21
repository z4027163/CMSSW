// -*- C++ -*-
//
// Package:    MicroGMTInputProducer
// Class:      MicroGMTInputProducer
//
/**\class MicroGMTInputProducer MicroGMTInputProducer.cc L1Trigger/L1TGlobalMuon/plugins/MicroGMTInputProducer.cc

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
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/GMTInputCaloSumFwd.h"
#include "DataFormats/L1TMuon/interface/GMTInputCaloSum.h"

#include <iostream>
//
// class declaration
//
namespace l1t {
class MicroGMTInputProducer : public edm::EDProducer {
   public:
      explicit MicroGMTInputProducer(const edm::ParameterSet&);
      ~MicroGMTInputProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      void openFile();
      void skipHeader();
      int convertToInt(std::string &bitstr) const;
      static bool cmpProc(const l1t::RegionalMuonCand&, const l1t::RegionalMuonCand&);

      // ----------member data ---------------------------
      std::string m_fname;
      std::ifstream m_filestream;
      bool m_endOfBx;
      bool m_lastMuInBx;
      int m_currType;
      int m_currEvt;
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
MicroGMTInputProducer::MicroGMTInputProducer(const edm::ParameterSet& iConfig) :
  m_endOfBx(false), m_currType(0), m_currEvt(0)
{
  //register your products
  produces<RegionalMuonCandBxCollection>("BarrelTFMuons");
  produces<RegionalMuonCandBxCollection>("OverlapTFMuons");
  produces<RegionalMuonCandBxCollection>("ForwardTFMuons");
  produces<l1t::GMTInputCaloSumBxCollection>("TriggerTowerSums");

  //now do what ever other initialization is needed
  m_fname = iConfig.getParameter<std::string> ("inputFileName");

  openFile();
  skipHeader();
}


MicroGMTInputProducer::~MicroGMTInputProducer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  m_filestream.close();
}


//
// member functions
//
bool
MicroGMTInputProducer::cmpProc(const l1t::RegionalMuonCand& mu1, const l1t::RegionalMuonCand& mu2)
{
  return mu1.processor() < mu2.processor();
}

void
MicroGMTInputProducer::openFile()
{
  if (!m_filestream.is_open()) {
    m_filestream.open(m_fname.c_str());
    if (!m_filestream.good()) {
      cms::Exception("FileOpenError") << "Failed to open input file";
    }
  }
}

void
MicroGMTInputProducer::skipHeader()
{
  while (m_filestream.peek() == '#') {
    std::string tmp;
    getline(m_filestream, tmp);
  }
}

int
MicroGMTInputProducer::convertToInt(std::string &bitstr) const
{
  int num = 0;
  for (size_t cntr = 0; cntr < bitstr.size(); ++cntr) {
    char c = bitstr[cntr];
    num = (num << 1) |  // Shift the current set of bits to the left one bit
      (c - '0');        // Add in the current bit via a bitwise-or
  }
  return num;
}



// ------------ method called to produce the data  ------------
void
MicroGMTInputProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  std::auto_ptr<l1t::RegionalMuonCandBxCollection> barrelMuons (new l1t::RegionalMuonCandBxCollection());
  std::auto_ptr<l1t::RegionalMuonCandBxCollection> overlapMuons (new l1t::RegionalMuonCandBxCollection());
  std::auto_ptr<l1t::RegionalMuonCandBxCollection> endcapMuons (new l1t::RegionalMuonCandBxCollection());
  std::auto_ptr<l1t::GMTInputCaloSumBxCollection> towerSums (new l1t::GMTInputCaloSumBxCollection());

  l1t::RegionalMuonCand mu;
  l1t::GMTInputCaloSum tSum;
  m_endOfBx = false;
  int caloCounter = 0;
  std::vector<int> bar{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<int> ovl_neg{0, 0, 0, 0, 0, 0};
  std::vector<int> ovl_pos{0, 0, 0, 0, 0, 0};
  std::vector<int> fwd_neg{0, 0, 0, 0, 0, 0};
  std::vector<int> fwd_pos{0, 0, 0, 0, 0, 0};
  while(!m_endOfBx && !m_filestream.eof()) {
    std::string lineID;
    m_filestream >> lineID;
    std::string restOfLine;


    if (lineID == "BAR" || lineID == "OVL-" || lineID == "FWD-" || lineID == "OVL+" || lineID == "FWD+") {
      int tmp;
      m_filestream >> tmp; // cable no
      // if (lineID == "BAR") tmp += 12;
      // if (lineID == "OVL-") tmp = (tmp-6)+24;
      // if (lineID == "OVL+") tmp = tmp + 6;
      // if (lineID == "FWD-") tmp = (tmp-6)+30;

      // mu.setLink(tmp);
      m_filestream >> tmp;
      mu.setHwPt(tmp);

      m_filestream >> tmp;

      int globalPhi = int(tmp*0.560856864654333f); // correction from txt file producer!
      int globalWedgePhi = (globalPhi+24)%576; // this sets CMS phi = 0 to -15 deg
      int globalSectorPhi = (globalPhi-24); // this sets CMS phi = 0 to +15 deg
      if (globalSectorPhi < 0) {
        globalSectorPhi += 576;
      }


      // int globalMuonPhi = int(tmp*0.560856864654333f); // make sure scale is correct
      bool skip = false;
      if (lineID == "BAR") {
        int processor = globalWedgePhi / 48 + 1;
        int localPhi = globalWedgePhi%48;
        mu.setTFIdentifiers(processor, tftype::bmtf);
        mu.setHwPhi(localPhi);
        bar[processor-1]++;
        if (bar[processor-1] > 3) skip = true;
      }
      if (lineID == "OVL-") {
        int processor = globalSectorPhi / 96 + 1;
        int localPhi = globalSectorPhi%96;
        mu.setTFIdentifiers(processor, tftype::omtf_neg);
        mu.setHwPhi(localPhi);
        ovl_neg[processor-1]++;
        if (ovl_neg[processor-1] > 3) skip = true;
      }
      if (lineID == "OVL+") {
        int processor = globalSectorPhi / 96 + 1;
        int localPhi = globalSectorPhi%96;
        mu.setTFIdentifiers(processor, tftype::omtf_pos);
        mu.setHwPhi(localPhi);
        ovl_pos[processor-1]++;
        if (ovl_pos[processor-1] > 3) skip = true;
      }
      if (lineID == "FWD-") {
        int processor = globalSectorPhi / 96 + 1;
        int localPhi = globalSectorPhi%96;
        mu.setTFIdentifiers(processor, tftype::emtf_neg);
        mu.setHwPhi(localPhi);
        fwd_neg[processor-1]++;
        if (fwd_neg[processor-1] > 3) skip = true;
      }
      if (lineID == "FWD+") {
        int processor = globalSectorPhi / 96 + 1;
        int localPhi = globalSectorPhi%96;
        mu.setTFIdentifiers(processor, tftype::emtf_pos);
        mu.setHwPhi(localPhi);
        fwd_pos[processor-1]++;
        if (fwd_pos[processor-1] > 3) skip = true;
      }

      m_filestream >> tmp;
      tmp = int(tmp*0.9090909090f);
      mu.setHwEta(tmp);



      m_filestream >> tmp;
      mu.setHwSign(tmp);

      m_filestream >> tmp;
      mu.setHwSignValid(tmp);

      m_filestream >> tmp;
      mu.setHwQual(tmp);

      if (lineID == "BAR") m_currType = 0;
      if (lineID == "OVL-") m_currType = 1;
      if (lineID == "OVL+") m_currType = 2;
      if (lineID == "FWD-") m_currType = 3;
      if (lineID == "FWD+") m_currType = 4;

      if (m_currType == 0 && !skip)  barrelMuons->push_back(0, mu);
      if ((m_currType == 1 || m_currType == 2) && !skip) overlapMuons->push_back(0, mu);
      if ((m_currType == 3 || m_currType == 4) && !skip) endcapMuons->push_back(0, mu);
    }

    if (lineID == "EVT" && m_currEvt != 0) {
        m_endOfBx = true;
    } else if (lineID == "EVT") {
      m_currEvt++;
    }

    if (lineID == "CALO") {
      for (int i = 0; i < 28; ++i) {
        int ieta = i; //caloCounter%28;
        int iphi = caloCounter;
        int et;

        m_filestream >> et;
        tSum.setEtBits(et);
        tSum.setEtaBits(ieta);
        tSum.setPhiBits(iphi);
        tSum.setIndex(caloCounter*28+i);
        towerSums->push_back(0, tSum);
      }
      caloCounter++;
    }
    getline(m_filestream, restOfLine);
    //std::cout << restOfLine;
  }


  // std::sort(barrelMuons->begin(0), barrelMuons->end(0), MicroGMTInputProducer::cmpProc);
  // std::sort(overlapMuons->begin(0), overlapMuons->end(0), MicroGMTInputProducer::cmpProc);
  // std::sort(endcapMuons->begin(0), endcapMuons->end(0), MicroGMTInputProducer::cmpProc);

  iEvent.put(barrelMuons, "BarrelTFMuons");
  iEvent.put(overlapMuons, "OverlapTFMuons");
  iEvent.put(endcapMuons, "ForwardTFMuons");
  iEvent.put(towerSums, "TriggerTowerSums");
  m_currEvt++;

}

// ------------ method called once each job just before starting event loop  ------------
void
MicroGMTInputProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MicroGMTInputProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
MicroGMTInputProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
MicroGMTInputProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
MicroGMTInputProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MicroGMTInputProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MicroGMTInputProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
}
//define this as a plug-in
DEFINE_FWK_MODULE(l1t::MicroGMTInputProducer);
