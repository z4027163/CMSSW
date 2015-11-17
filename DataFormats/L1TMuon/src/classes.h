#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/L1TMuon/interface/MuonRegionalTracksFwd.h"

#include "DataFormats/RPCDigi/interface/RPCDigiL1Link.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTTrackContainer.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"

#include "DataFormats/L1TMuon/interface/GMTInputCaloSumFwd.h"
#include "DataFormats/L1TMuon/interface/GMTInputCaloSum.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

namespace {
  using namespace L1TMuon;
  struct dictionary {
    l1t::GMTInputCaloSumBxCollection caloSum;
    edm::Wrapper<l1t::GMTInputCaloSumBxCollection> caloSumWrap;

    l1t::RegionalMuonCandBxCollection regCand;
    edm::Wrapper<l1t::RegionalMuonCandBxCollection> regCandWrap;

   
    RegionalCandBaseRef rcR2B;
    RegionalCandPtr     rcPtr;
    RegionalCandRef     rfRef;

    DTTrackCollection dtTrkColl;
    edm::Wrapper<DTTrackCollection> wdtTrkColl;
    DTTrackPtr dtTrkPtr;
    DTTrackRef dtTrackRef;

    CSCTrackCollection cscTrkColl;
    edm::Wrapper<CSCTrackCollection> wcscTrkColl;
    CSCTrackPtr cscTrkPtr;
    CSCTrackRef cscTrkRef;

    RPCL1LinkPtr prpcL1link;
    RPCL1LinkRef rrpcL1link;

    edm::reftobase::Holder<L1MuRegionalCand,RegionalCandRef>  r2rholder;
    edm::reftobase::Holder<L1MuRegionalCand,DTTrackRef>  r2dtholder;
    edm::reftobase::Holder<L1MuRegionalCand,CSCTrackRef>  r2cscholder;
    edm::reftobase::RefHolder<RegionalCandRef>  r2rrefholder;
    edm::reftobase::RefHolder<DTTrackRef>  r2dtrefholder;
    edm::reftobase::RefHolder<CSCTrackRef>  r2cscrefholder;
  };
}



#include <DataFormats/L1TMuon/interface/BMTrackCand.h>
#include <DataFormats/L1TMuon/interface/BMTrackContainer.h>
#include <DataFormats/Common/interface/Wrapper.h>

namespace DataFormats_L1BMTrackFinder {
  struct dictionary {
    BMTrackCand   tr_S;
    std::vector<BMTrackCand>   tr_V;

    BMTrackContainer   tr_K;

    edm::Wrapper<BMTrackContainer>   tr_W;
  };
}


