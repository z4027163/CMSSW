// -*- C++ -*-
//
// Package:    L1Trigger/L1TMicroGMTParamsESProducer
// Class:      L1TMicroGMTParamsESProducer
// 
/**\class L1TMicroGMTParamsESProducer L1TMicroGMTParamsESProducer.h L1Trigger/L1TMicroGMTParamsESProducer/plugins/L1TMicroGMTParamsESProducer.cc

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
#include "boost/shared_ptr.hpp"

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

using namespace l1t;

class L1TMicroGMTParamsESProducer : public edm::ESProducer {
   public:
      L1TMicroGMTParamsESProducer(const edm::ParameterSet&);
      ~L1TMicroGMTParamsESProducer();

      typedef boost::shared_ptr<MicroGMTParams> ReturnType;

      ReturnType produce(const L1TMicroGMTParamsRcd&);
   private:
      MicroGMTParams m_params;
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
L1TMicroGMTParamsESProducer::L1TMicroGMTParamsESProducer(const edm::ParameterSet& iConfig)
{
   //the following line is needed to tell the framework what
   // data is being produced
   setWhatProduced(this);

   // Firmware version
   m_params.setFwVersion(iConfig.getParameter<unsigned>("fwVersion"));

   // LUT paths
   m_params.setAbsIsoCheckMemLUTPath        (iConfig.getParameter<std::string>("AbsIsoCheckMemLUTPath"));
   m_params.setRelIsoCheckMemLUTPath        (iConfig.getParameter<std::string>("RelIsoCheckMemLUTPath"));
   m_params.setIdxSelMemPhiLUTPath          (iConfig.getParameter<std::string>("IdxSelMemPhiLUTPath"));
   m_params.setIdxSelMemEtaLUTPath          (iConfig.getParameter<std::string>("IdxSelMemEtaLUTPath"));
   m_params.setBrlSingleMatchQualLUTPath    (iConfig.getParameter<std::string>("BrlSingleMatchQualLUTPath"));
   m_params.setFwdPosSingleMatchQualLUTPath (iConfig.getParameter<std::string>("FwdPosSingleMatchQualLUTPath"));
   m_params.setFwdNegSingleMatchQualLUTPath (iConfig.getParameter<std::string>("FwdNegSingleMatchQualLUTPath"));
   m_params.setOvlPosSingleMatchQualLUTPath (iConfig.getParameter<std::string>("OvlPosSingleMatchQualLUTPath"));
   m_params.setOvlNegSingleMatchQualLUTPath (iConfig.getParameter<std::string>("OvlNegSingleMatchQualLUTPath"));
   m_params.setBOPosMatchQualLUTPath        (iConfig.getParameter<std::string>("BOPosMatchQualLUTPath"));
   m_params.setBONegMatchQualLUTPath        (iConfig.getParameter<std::string>("BONegMatchQualLUTPath"));
   m_params.setFOPosMatchQualLUTPath        (iConfig.getParameter<std::string>("FOPosMatchQualLUTPath"));
   m_params.setFONegMatchQualLUTPath        (iConfig.getParameter<std::string>("FONegMatchQualLUTPath"));
   m_params.setBPhiExtrapolationLUTPath     (iConfig.getParameter<std::string>("BPhiExtrapolationLUTPath"));
   m_params.setOPhiExtrapolationLUTPath     (iConfig.getParameter<std::string>("OPhiExtrapolationLUTPath"));
   m_params.setFPhiExtrapolationLUTPath     (iConfig.getParameter<std::string>("FPhiExtrapolationLUTPath"));
   m_params.setBEtaExtrapolationLUTPath     (iConfig.getParameter<std::string>("BEtaExtrapolationLUTPath"));
   m_params.setOEtaExtrapolationLUTPath     (iConfig.getParameter<std::string>("OEtaExtrapolationLUTPath"));
   m_params.setFEtaExtrapolationLUTPath     (iConfig.getParameter<std::string>("FEtaExtrapolationLUTPath"));
   m_params.setSortRankLUTPath              (iConfig.getParameter<std::string>("SortRankLUTPath"));
}


L1TMicroGMTParamsESProducer::~L1TMicroGMTParamsESProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
L1TMicroGMTParamsESProducer::ReturnType
L1TMicroGMTParamsESProducer::produce(const L1TMicroGMTParamsRcd& iRecord)
{
   using namespace edm::es;
   boost::shared_ptr<MicroGMTParams> pMicroGMTParams;

   pMicroGMTParams = boost::shared_ptr<MicroGMTParams>(new MicroGMTParams(m_params));
   return pMicroGMTParams;
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(L1TMicroGMTParamsESProducer);
