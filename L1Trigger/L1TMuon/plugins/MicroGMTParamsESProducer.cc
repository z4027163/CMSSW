// -*- C++ -*-
//
// Package:    L1Trigger/MicroGMTParamsESProducer
// Class:      MicroGMTParamsESProducer
// 
/**\class MicroGMTParamsESProducer MicroGMTParamsESProducer.h L1Trigger/MicroGMTParamsESProducer/plugins/MicroGMTParamsESProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Reis
//         Created:  Mon, 21 Sep 2015 13:28:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/ESProducts.h"

#include "CondFormats/L1TObjects/interface/MicroGMTParams.h"
#include "CondFormats/DataRecord/interface/L1TMicroGMTParamsRcd.h"

//
// class declaration
//

class MicroGMTParamsESProducer : public edm::ESProducer {
   public:
      MicroGMTParamsESProducer(const edm::ParameterSet&);
      ~MicroGMTParamsESProducer();

      typedef std::shared_ptr<l1t::MicroGMTParams> ReturnType;

      ReturnType produce(const MicroGMTParamsRcd&);
   private:
      l1t::MicroGMTParams m_params;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MicroGMTParamsESProducer::MicroGMTParamsESProducer(const edm::ParameterSet& iConfig)
{
   //the following line is needed to tell the framework what
   // data is being produced
   setWhatProduced(this);

   edm::ParameterSet AbsIsoCheckMemLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("AbsIsoCheckMemLUTSettings");
   m_params.absIsoCheckMemLUTParams()->setAreaSumInWidth(AbsIsoCheckMemLUTSettings_.getParameter<int>("areaSum_in_width"));
   m_params.absIsoCheckMemLUTParams()->setOutWidth(AbsIsoCheckMemLUTSettings_.getParameter<int>("out_width"));
   m_params.absIsoCheckMemLUTParams()->setFilename(AbsIsoCheckMemLUTSettings_.getParameter<std::string>("filename"));

   //edm::ParameterSet IdxSelMemPhiLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("IdxSelMemPhiLUTSettings");
   //edm::ParameterSet FwdPosSingleMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("FwdPosSingleMatchQualLUTSettings");
   //edm::ParameterSet BONegMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("BONegMatchQualLUTSettings");
   //edm::ParameterSet OvlNegSingleMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("OvlNegSingleMatchQualLUTSettings");
   //edm::ParameterSet IdxSelMemEtaLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("IdxSelMemEtaLUTSettings");
   //edm::ParameterSet FOPosMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("FOPosMatchQualLUTSettings");
   //edm::ParameterSet FwdNegSingleMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("FwdNegSingleMatchQualLUTSettings");
   //edm::ParameterSet BPhiExtrapolationLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("BPhiExtrapolationLUTSettings");
   //edm::ParameterSet BrlSingleMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("BrlSingleMatchQualLUTSettings");
   //edm::ParameterSet RelIsoCheckMemLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("RelIsoCheckMemLUTSettings");
   //edm::ParameterSet OPhiExtrapolationLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("OPhiExtrapolationLUTSettings");
   //edm::ParameterSet OvlPosSingleMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("OvlPosSingleMatchQualLUTSettings");
   //edm::ParameterSet FEtaExtrapolationLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("FEtaExtrapolationLUTSettings");
   //edm::ParameterSet BOPosMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("BOPosMatchQualLUTSettings");
   //edm::ParameterSet OEtaExtrapolationLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("OEtaExtrapolationLUTSettings");
   //edm::ParameterSet BEtaExtrapolationLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("BEtaExtrapolationLUTSettings");
   //edm::ParameterSet FPhiExtrapolationLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("FPhiExtrapolationLUTSettings");
   //edm::ParameterSet FONegMatchQualLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("FONegMatchQualLUTSettings");
   //edm::ParameterSet SortRankLUTSettings_ = iConfig.getParameter<edm::ParameterSet>("SortRankLUTSettings");
}


MicroGMTParamsESProducer::~MicroGMTParamsESProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
MicroGMTParamsESProducer::ReturnType
MicroGMTParamsESProducer::produce(const MicroGMTParamsRcd& iRecord)
{
   using namespace edm::es;
   std::shared_ptr<l1t::MicroGMTParams> pMicroGMTParams;

   pMicroGMTParams = std::shared_ptr<l1t::MicroGMTParams>(new l1t::MicroGMTParams(m_params));
   return products(pMicroGMTParams);
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(MicroGMTParamsESProducer);
