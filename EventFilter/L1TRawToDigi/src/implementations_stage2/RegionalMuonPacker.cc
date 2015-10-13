#include "FWCore/Framework/interface/Event.h"

#include "EventFilter/L1TRawToDigi/interface/Packer.h"

#include "L1Trigger/L1TMuon/interface/RegionalMuonRawDigiTranslator.h"
#include "GMTTokens.h"

namespace l1t {
   namespace stage2 {
      class RegionalMuonPacker : public Packer {
         public:
            virtual Blocks pack(const edm::Event&, const PackerTokens*) override;
         private:
            typedef std::map<unsigned int, std::vector<uint32_t>> LoadMap;
            void packTF(const edm::Event&, const edm::EDGetTokenT<RegionalMuonCandBxCollection>&, Blocks&);
      };
   }
}

// Implementation
namespace l1t {
   namespace stage2 {
      Blocks
      RegionalMuonPacker::pack(const edm::Event& event, const PackerTokens* toks)
      {
         auto bmtfToken = static_cast<const GMTTokens*>(toks)->getRegionalMuonCandTokenBMTF();
         auto omtfToken = static_cast<const GMTTokens*>(toks)->getRegionalMuonCandTokenOMTF();
         auto emtfToken = static_cast<const GMTTokens*>(toks)->getRegionalMuonCandTokenEMTF();

         Blocks blocks;

         // pack the muons for each TF in blocks
         packTF(event, bmtfToken, blocks);
         packTF(event, omtfToken, blocks);
         packTF(event, emtfToken, blocks);

         return blocks;
      }

      void
      RegionalMuonPacker::packTF(const edm::Event& event, const edm::EDGetTokenT<RegionalMuonCandBxCollection>& tfToken, Blocks &blocks)
      {
         edm::Handle<RegionalMuonCandBxCollection> muons;
         event.getByToken(tfToken, muons);
   
         LoadMap loadMap;
   
         for (int i = muons->getFirstBX(); i <= muons->getLastBX(); ++i) {
            if (muons->size(i) == 0)
               continue;
            for (auto mu = muons->begin(i); mu != muons->end(i); ++mu) {
               uint32_t msw = 0;
               uint32_t lsw = 0;

               RegionalMuonRawDigiTranslator::generatePackedDataWords(*mu, lsw, msw);

               loadMap[mu->link()*2].push_back(lsw);
               loadMap[mu->link()*2].push_back(msw);
            }

            // padding to 3 muons per block id (link) per BX
            for (auto &kv : loadMap) {
               while (kv.second.size()%6 != 0) {
                  kv.second.push_back(0);
               }
            }
         }

         // push everything in the blocks vector
         for (auto &kv : loadMap) {
            blocks.push_back(Block(kv.first, kv.second));
         }
      }
   }
}

DEFINE_L1T_PACKER(l1t::stage2::RegionalMuonPacker);
