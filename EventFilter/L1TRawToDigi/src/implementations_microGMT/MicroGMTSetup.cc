#include "FWCore/Framework/interface/stream/EDProducerBase.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "EventFilter/L1TRawToDigi/interface/Packer.h"
#include "EventFilter/L1TRawToDigi/interface/Unpacker.h"

#include "EventFilter/L1TRawToDigi/interface/PackingSetup.h"

#include "MicroGMTCollections.h"
#include "MicroGMTTokens.h"

namespace l1t {
   class MicroGMTSetup : public PackingSetup {
      public:
         virtual std::unique_ptr<PackerTokens> registerConsumes(const edm::ParameterSet& cfg, edm::ConsumesCollector& cc) override {
            return std::unique_ptr<PackerTokens>(new MicroGMTTokens(cfg, cc));
         };

         virtual void fillDescription(edm::ParameterSetDescription& desc) override {};

         virtual PackerMap getPackers(int fed, unsigned int fw) override {
            PackerMap res;

            //res[{1, 1}] = {
            //   PackerFactory::get()->make("MicroGMTPacker"),
            //};

            return res;
         };

         virtual void registerProducts(edm::stream::EDProducerBase& prod) override {
            prod.produces<RegionalMuonCandBxCollection>("BMTF");
            prod.produces<RegionalMuonCandBxCollection>("OMTF");
            prod.produces<RegionalMuonCandBxCollection>("EMTF");
            prod.produces<MuonBxCollection>();
         };

         virtual std::unique_ptr<UnpackerCollections> getCollections(edm::Event& e) override {
            return std::unique_ptr<UnpackerCollections>(new MicroGMTCollections(e));
         };

         virtual UnpackerMap getUnpackers(int fed, int board, int amc, unsigned int fw) override {
            UnpackerMap res;

            auto microGMT_in_unp = UnpackerFactory::get()->make("MicroGMTInUnpacker");
            auto microGMT_out_unp = UnpackerFactory::get()->make("MicroGMTOutUnpacker");

            for (int iLink = 72; iLink < 144; iLink += 2)
                res[iLink] = microGMT_in_unp;
            for (int oLink = 1; oLink < 9; oLink += 2)
                res[oLink] = microGMT_out_unp;

            return res;
         };
   };
}

DEFINE_L1T_PACKING_SETUP(l1t::MicroGMTSetup);
