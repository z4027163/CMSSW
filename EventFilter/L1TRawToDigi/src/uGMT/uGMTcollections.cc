#include "FWCore/Framework/interface/Event.h"

#include "uGMTcollections.h"

namespace l1t {
   UGMTcollections::~UGMTcollections()
   {
      event_.put(regionalMuonCandsBMTF_, "BMTF");
      event_.put(regionalMuonCandsOMTF_, "OMTF");
      event_.put(regionalMuonCandsEMTF_, "EMTF");
      event_.put(muons_);
   }
}
