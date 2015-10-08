#ifndef RegionalMuonRawDigiTranslator_h
#define RegionalMuonRawDigiTranslator_h

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

namespace l1t {
  class RegionalMuonRawDigiTranslator {
    public:
      static void fillRegionalMuonCand(RegionalMuonCand&, uint32_t, uint32_t, int, tftype);
  };
}

#endif
