// -*- C++ -*-
//
// Package:    MicroGMTEmulator
// Class:      MicroGMTEmulator
// 
/**\class MicroGMTEmulator MicroGMTEmulator.cc L1Trigger/L1TMuon/src/MicroGMTEmulator.cc

 Description: Takes txt-file input and produces barrel- / overlap- / forward TF muons 

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joschka Philip Lingemann,40 3-B01,+41227671598,
//         Created:  Thu Oct  3 16:31:34 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <fstream>
#include <sstream>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTRankPtQualLUT.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTIsolationUnit.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTCancelOutUnit.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidate.h"
#include "DataFormats/L1TMuon/interface/L1TGMTInternalMuon.h"

#include "TMath.h"
//
// class declaration
//
namespace l1t {
  class MicroGMTEmulator : public edm::EDProducer {
     public:
        explicit MicroGMTEmulator(const edm::ParameterSet&);
        ~MicroGMTEmulator();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

     private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        virtual void beginRun(edm::Run&, edm::EventSetup const&);
        virtual void endRun(edm::Run&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

        static bool compareMuons(const MicroGMTConfiguration::InterMuon&, const MicroGMTConfiguration::InterMuon&);

        void sortMuons(MicroGMTConfiguration::InterMuonList&, unsigned) const;
        void calculateRank(MicroGMTConfiguration::InterMuonList& muons) const;
        void splitAndConvertMuons(MicroGMTConfiguration::InputCollection const& in, MicroGMTConfiguration::InterMuonList& out_pos, MicroGMTConfiguration::InterMuonList& out_neg) const;
        void convertMuons(MicroGMTConfiguration::InputCollection const& in, MicroGMTConfiguration::InterMuonList& out) const;
        void addMuonsToCollections(MicroGMTConfiguration::InterMuonList& coll, MicroGMTConfiguration::InterMuonList& interout, std::auto_ptr<MuonBxCollection>& out) const;
        // ----------member data ---------------------------
        edm::InputTag m_barrelTfInputTag;
        edm::InputTag m_overlapTfInputTag;
        edm::InputTag m_forwardTfInputTag;
        edm::InputTag m_trigTowerTag;
        MicroGMTRankPtQualLUT m_rankPtQualityLUT;
        MicroGMTIsolationUnit m_isolationUnit;
        MicroGMTCancelOutUnit m_cancelOutUnit;
        std::ofstream m_debugOut;

  };
}
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
l1t::MicroGMTEmulator::MicroGMTEmulator(const edm::ParameterSet& iConfig) : m_rankPtQualityLUT(iConfig), m_isolationUnit(iConfig), m_cancelOutUnit(iConfig), m_debugOut("test/debug/iso_debug.dat")
{
  // edm::InputTag barrelTfInputTag = iConfig.getParameter<edm::InputTag>("barrelTFInput");
  // edm::InputTag overlapTfInputTag = iConfig.getParameter<edm::InputTag>("overlapTFInput");
  // edm::InputTag forwardTfInputTag = iConfig.getParameter<edm::InputTag>("forwardTFInput");

  // m_barrelTfInputToken = consumes<InputCollection>(barrelTfInputTag);
  // m_overlapTfInputToken = consumes<InputCollection>(overlapTfInputTag);
  // m_forwardTfInputToken = consumes<InputCollection>(forwardTfInputTag);
   //register your products
  produces<MuonBxCollection>();
  produces<MuonBxCollection>("imdMuonsBMTF");
  produces<MuonBxCollection>("imdMuonsEMTFPos");
  produces<MuonBxCollection>("imdMuonsEMTFNeg");
  produces<MuonBxCollection>("imdMuonsOMTFPos");
  produces<MuonBxCollection>("imdMuonsOMTFNeg");

  m_barrelTfInputTag = iConfig.getParameter<edm::InputTag>("barrelTFInput");
  m_overlapTfInputTag = iConfig.getParameter<edm::InputTag>("overlapTFInput");
  m_forwardTfInputTag = iConfig.getParameter<edm::InputTag>("forwardTFInput");
  m_trigTowerTag = iConfig.getParameter<edm::InputTag>("triggerTowerInput");

}

l1t::MicroGMTEmulator::~MicroGMTEmulator()
{
  m_debugOut.close();
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//



// ------------ method called to produce the data  ------------
void
l1t::MicroGMTEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  std::auto_ptr<MuonBxCollection> outMuons (new MuonBxCollection());
  std::auto_ptr<MuonBxCollection> imdMuonsBMTF (new MuonBxCollection());
  std::auto_ptr<MuonBxCollection> imdMuonsEMTFPos (new MuonBxCollection());
  std::auto_ptr<MuonBxCollection> imdMuonsEMTFNeg (new MuonBxCollection());
  std::auto_ptr<MuonBxCollection> imdMuonsOMTFPos (new MuonBxCollection());
  std::auto_ptr<MuonBxCollection> imdMuonsOMTFNeg (new MuonBxCollection());

  Handle<MicroGMTConfiguration::InputCollection> barrelMuons;
  Handle<MicroGMTConfiguration::InputCollection> forwardMuons;
  Handle<MicroGMTConfiguration::InputCollection> overlapMuons;
  Handle<MicroGMTConfiguration::CaloInputCollection> trigTowers;

  // iEvent.getByToken(m_barrelTfInputToken, barrelMuons);
  iEvent.getByLabel(m_barrelTfInputTag, barrelMuons);
  iEvent.getByLabel(m_forwardTfInputTag, forwardMuons);
  iEvent.getByLabel(m_overlapTfInputTag, overlapMuons);
  iEvent.getByLabel(m_trigTowerTag, trigTowers);
  
  m_isolationUnit.setTowerSums(*trigTowers);
  MicroGMTConfiguration::InterMuonList internalMuonsBarrel;
  MicroGMTConfiguration::InterMuonList internalMuonsEndcapPos;
  MicroGMTConfiguration::InterMuonList internalMuonsEndcapNeg;
  MicroGMTConfiguration::InterMuonList internalMuonsOverlapPos;
  MicroGMTConfiguration::InterMuonList internalMuonsOverlapNeg;

  // this converts the InputMuon type to the InternalMuon type
  // and splits them into positive / negative eta collections
  // necessary as LUTs may differ for pos / neg.
  convertMuons(*barrelMuons, internalMuonsBarrel);
  splitAndConvertMuons(*forwardMuons, internalMuonsEndcapPos, internalMuonsEndcapNeg);
  splitAndConvertMuons(*overlapMuons, internalMuonsOverlapPos, internalMuonsOverlapNeg);

  m_isolationUnit.extrapolateMuons(internalMuonsBarrel);
  m_isolationUnit.extrapolateMuons(internalMuonsEndcapNeg);
  m_isolationUnit.extrapolateMuons(internalMuonsEndcapPos);
  m_isolationUnit.extrapolateMuons(internalMuonsOverlapNeg);
  m_isolationUnit.extrapolateMuons(internalMuonsOverlapPos);

  calculateRank(internalMuonsBarrel);
  calculateRank(internalMuonsEndcapNeg);
  calculateRank(internalMuonsEndcapPos);
  calculateRank(internalMuonsOverlapNeg);
  calculateRank(internalMuonsOverlapPos);
  
  // The sort function both sorts and removes all but best "nSurvivors"
  sortMuons(internalMuonsBarrel, 8);
  sortMuons(internalMuonsOverlapPos, 4);
  sortMuons(internalMuonsOverlapNeg, 4);
  sortMuons(internalMuonsEndcapPos, 4);
  sortMuons(internalMuonsEndcapNeg, 4);


  // This combines the 5 streams into one InternalMuon collection for 
  // the final global sort.
  MicroGMTConfiguration::InterMuonList internalMuons;
  addMuonsToCollections(internalMuonsEndcapPos, internalMuons, imdMuonsEMTFPos);
  addMuonsToCollections(internalMuonsOverlapPos, internalMuons, imdMuonsOMTFPos);
  addMuonsToCollections(internalMuonsBarrel, internalMuons, imdMuonsBMTF);
  addMuonsToCollections(internalMuonsOverlapNeg, internalMuons, imdMuonsOMTFNeg);
  addMuonsToCollections(internalMuonsEndcapNeg, internalMuons, imdMuonsEMTFNeg);
  
  // sort internal muons and delete all but best 8
  sortMuons(internalMuons, 8);

  m_isolationUnit.isolatePreSummed(internalMuons);
  // copy muons to output collection...
  for (auto mu = internalMuons.begin(); mu != internalMuons.end(); ++mu) {
    if (mu->hwPt() > 0) {
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vec{};
      int iso = mu->hwAbsIso() + (mu->hwRelIso() << 1);
      // FIXME: once we debugged the change global -> local: Change hwLocalPhi -> hwGlobalPhi to test offsets
      Muon outMu{vec, mu->hwPt(), mu->hwEta(), mu->hwLocalPhi(), mu->hwQual(), mu->hwSign(), mu->hwSignValid(), iso, 0, true, mu->hwIsoSum(), mu->hwDPhi(), mu->hwDEta(), mu->hwRank()};
      m_debugOut << mu->hwCaloPhi() << " " << mu->hwCaloEta() << std::endl;
      outMuons->push_back(0, outMu);
    }
  }
  
  iEvent.put(outMuons);
  iEvent.put(imdMuonsBMTF, "imdMuonsBMTF");
  iEvent.put(imdMuonsEMTFPos, "imdMuonsEMTFPos");
  iEvent.put(imdMuonsEMTFNeg, "imdMuonsEMTFNeg");
  iEvent.put(imdMuonsOMTFPos, "imdMuonsOMTFPos");
  iEvent.put(imdMuonsOMTFNeg, "imdMuonsOMTFNeg");
}


bool 
l1t::MicroGMTEmulator::compareMuons(const MicroGMTConfiguration::InterMuon& mu1, const MicroGMTConfiguration::InterMuon& mu2) {
  return (mu1.hwWins() > mu2.hwWins());
}

void 
l1t::MicroGMTEmulator::sortMuons(MicroGMTConfiguration::InterMuonList& muons, unsigned nSurvivors) const {
  MicroGMTConfiguration::InterMuonList::iterator mu1;
  // reset from previous sort stage
  for (mu1 = muons.begin(); mu1 != muons.end(); ++mu1) {
    mu1->setHwWins(0);
  }
  
  for (mu1 = muons.begin(); mu1 != muons.end(); ++mu1) {
    auto mu2 = mu1;
    mu2++;
    for ( ; mu2 != muons.end(); ++mu2) {
      if (mu1->hwRank() >= mu2->hwRank() && mu1->hwCancelBit() != 1) {
        mu1->increaseWins();
      } else {
        mu2->increaseWins();
      }
    }
  }

  size_t nMuonsBefore = muons.size();
  mu1 = muons.begin();
  if (nMuonsBefore > nSurvivors) {
    while (mu1 != muons.end()) {
      if (mu1->hwWins() < (int)(nMuonsBefore-nSurvivors)) {
        muons.erase(mu1);
      }
      ++mu1;
    }
  }

  muons.sort(l1t::MicroGMTEmulator::compareMuons);
}



void 
l1t::MicroGMTEmulator::calculateRank(MicroGMTConfiguration::InterMuonList& muons) const 
{
  for (auto mu1 = muons.begin(); mu1 != muons.end(); ++mu1) {
    int rank = m_rankPtQualityLUT.lookup(mu1->hwPt(), mu1->hwQual());
    mu1->setHwRank(rank);
  }
}


void 
l1t::MicroGMTEmulator::addMuonsToCollections(MicroGMTConfiguration::InterMuonList& coll, MicroGMTConfiguration::InterMuonList& interout, std::auto_ptr<MuonBxCollection>& out) const {
  for (auto mu = coll.cbegin(); mu != coll.cend(); ++mu) { 
    interout.push_back(*mu);
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vec{};
    // FIXME: once we debugged the change global -> local: Change hwLocalPhi -> hwGlobalPhi to test offsets
    Muon outMu{vec, mu->hwPt(), mu->hwEta(), mu->hwLocalPhi(), mu->hwQual(), mu->hwSign(), mu->hwSignValid(), -1, 0, true, -1, mu->hwDPhi(), mu->hwDEta(), mu->hwRank()};
    out->push_back(0, outMu);
  }
}

void
l1t::MicroGMTEmulator::splitAndConvertMuons(const MicroGMTConfiguration::InputCollection& in, MicroGMTConfiguration::InterMuonList& out_pos, MicroGMTConfiguration::InterMuonList& out_neg) const
{
  for (auto mu = in.cbegin(); mu != in.cend(); ++mu) {
    if(mu->hwEta() > 0) {
      out_pos.emplace_back(*mu);
    } else {
      out_neg.emplace_back(*mu);
    }
  }
}
        
void 
l1t::MicroGMTEmulator::convertMuons(const MicroGMTConfiguration::InputCollection& in, MicroGMTConfiguration::InterMuonList& out) const
{
  for (auto mu = in.cbegin(); mu != in.cend(); ++mu) {
    out.emplace_back(*mu);
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
l1t::MicroGMTEmulator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
l1t::MicroGMTEmulator::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
l1t::MicroGMTEmulator::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
l1t::MicroGMTEmulator::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
l1t::MicroGMTEmulator::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
l1t::MicroGMTEmulator::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
l1t::MicroGMTEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(l1t::MicroGMTEmulator);
