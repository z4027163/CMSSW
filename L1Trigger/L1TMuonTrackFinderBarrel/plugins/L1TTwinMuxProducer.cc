#include "FWCore/Framework/interface/Event.h"
#include "L1Trigger/L1TMuonTrackFinderBarrel/src/Twinmux_v1/L1TTwinMuxAlgorithm.cc"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ConsumesCollector.h>
#include <FWCore/Framework/interface/one/EDProducer.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ConsumesCollector.h>


#include <iostream>
#include <iomanip>

using namespace std;

class L1TTwinMuxProducer: public edm::one::EDProducer<edm::one::SharedResources> {
public:
  L1TTwinMuxProducer(const edm::ParameterSet & pset);
  ~L1TTwinMuxProducer() {}
  void produce(edm::Event & e, const edm::EventSetup& c);
private:
  L1TTwinMuxAlgortithm * l1tma;

};




L1TTwinMuxProducer::L1TTwinMuxProducer(const edm::ParameterSet & pset) {
l1tma = new L1TTwinMuxAlgortithm();
consumes<L1MuDTChambPhContainer>(pset.getParameter<edm::InputTag>("DTDigi_Source"));
consumes<L1MuDTChambThContainer>(pset.getParameter<edm::InputTag>("DTThetaDigi_Source"));
consumes<RPCDigiCollection>(pset.getParameter<edm::InputTag>("RPC_Source"));
produces<L1MuDTChambPhContainer>();

}

void L1TTwinMuxProducer::produce(edm::Event& e, const edm::EventSetup& c) {

  edm::InputTag _src("simDtTriggerPrimitiveDigis");
  edm::Handle<L1MuDTChambPhContainer> phiDigis;
  edm::Handle<L1MuDTChambThContainer> thetaDigis;
  e.getByLabel(_src,phiDigis);
  e.getByLabel(_src,thetaDigis);
  edm::Handle<RPCDigiCollection> rpcDigis;
  edm::InputTag _srcrpc("simMuonRPCDigis");
  e.getByLabel(_srcrpc,rpcDigis);

  std::auto_ptr<L1MuDTChambPhContainer> l1ttmp = l1tma->produce(phiDigis, thetaDigis, rpcDigis,c);
  e.put(l1ttmp);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TTwinMuxProducer);
