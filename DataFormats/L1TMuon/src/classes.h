#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/L1TMuon/interface/GMTInputCaloSumFwd.h"
#include "DataFormats/L1TMuon/interface/GMTInputCaloSum.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

namespace {
  struct dictionary {
    l1t::GMTInputCaloSumBxCollection caloSum;
    edm::Wrapper<l1t::GMTInputCaloSumBxCollection> caloSumWrap;

    l1t::RegionalMuonCandBxCollection regCand;
    edm::Wrapper<l1t::RegionalMuonCandBxCollection> regCandWrap;
   
  };
}


