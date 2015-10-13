#ifndef RegionalMuonRawDigiTranslator_h
#define RegionalMuonRawDigiTranslator_h

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

namespace l1t {
  class RegionalMuonRawDigiTranslator {
    public:
      static void fillRegionalMuonCand(RegionalMuonCand&, uint32_t, uint32_t, int, tftype);
      static void fillRegionalMuonCand(RegionalMuonCand&, uint64_t, int, tftype);
      static void generatePackedDataWords(const RegionalMuonCand&, uint32_t&, uint32_t&);
      static uint64_t generate64bitDataWord(const RegionalMuonCand&);

    private:
      static const unsigned ptWidth_ = 0x1FF;
      static const unsigned ptShift_ = 0;
      static const unsigned qualWidth_ = 0xF;
      static const unsigned qualShift_ = 9;
      static const unsigned absEtaWidth_ = 0xFF;
      static const unsigned absEtaShift_ = 13;
      static const unsigned etaSignShift_ = 21;
      static const unsigned hfWidth_ = 0x1;
      static const unsigned hfShift_ = 22;
      static const unsigned phiWidth_ = 0xFF;
      static const unsigned phiShift_ = 23;
      static const unsigned signShift_ = 0;
      static const unsigned signValidShift_ = 1;
      static const unsigned trackAddressWidth_ = 0x1FFF;
      static const unsigned trackAddressShift_ = 4;
  };
}

#endif
