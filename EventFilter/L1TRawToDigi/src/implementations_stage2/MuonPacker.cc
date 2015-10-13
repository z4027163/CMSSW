#include "FWCore/Framework/interface/Event.h"

#include "EventFilter/L1TRawToDigi/interface/Packer.h"

#include "L1Trigger/L1TMuon/interface/MuonRawDigiTranslator.h"
#include "GMTTokens.h"

namespace l1t {
   namespace stage2 {
      class MuonPacker : public Packer {
         public:
            virtual Blocks pack(const edm::Event&, const PackerTokens*) override;
         private:
            typedef std::map<unsigned int, std::vector<uint32_t>> LoadMap;
      };
   }
}

// Implementation
namespace l1t {
   namespace stage2 {
      Blocks
      MuonPacker::pack(const edm::Event& event, const PackerTokens* toks)
      {
         edm::Handle<MuonBxCollection> muons;
         event.getByToken(static_cast<const GMTTokens*>(toks)->getMuonToken(), muons);
   
         LoadMap loadMap;
   
         for (int i = muons->getFirstBX(); i <= muons->getLastBX(); ++i) {
            int muCtr = 0;
            int blkCtr = 1;
            if (muons->size(i) == 0)
               continue;
            for (auto mu = muons->begin(i); mu != muons->end(i); ++mu) {
               uint32_t msw = 0;
               uint32_t lsw = 0;

               MuonRawDigiTranslator::generatePackedDataWords(*mu, lsw, msw);

               // FIXME: need to define block id somehow. round robin for now
               loadMap[blkCtr].push_back(lsw);
               loadMap[blkCtr].push_back(msw);

               // FIXME: As long as the muons are not assigned to one link
               // skip every 3rd slot in block so that there are only 2 muons per block
               if ((muCtr+2)%3 == 0) {
                  ++muCtr;
                  blkCtr += 2;
               }
               ++muCtr;
            }
         }

         // padding empty muons to reach 3 muons per id (link)
         // and push everything in the blocks vector
         Blocks blocks;
         for (auto &kv : loadMap) {
            for (auto i = kv.second.size()-1; i < 5; ++i) {
               kv.second.push_back(0);
            }
            blocks.push_back(Block(kv.first, kv.second));
         }

         return blocks;
      }
   }
}

DEFINE_L1T_PACKER(l1t::stage2::MuonPacker);
