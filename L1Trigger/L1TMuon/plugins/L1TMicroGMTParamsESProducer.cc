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

   // LUT parameters
   edm::ParameterSet aicm_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("AbsIsoCheckMemLUTSettings");
   m_params.absIsoCheckMemLUTParams()->setAreaSumInWidth(aicm_lut_stgs_.getParameter<int>("areaSum_in_width"));
   m_params.absIsoCheckMemLUTParams()->setOutWidth(aicm_lut_stgs_.getParameter<int>("out_width"));
   m_params.absIsoCheckMemLUTParams()->setFilename(aicm_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet ismp_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("IdxSelMemPhiLUTSettings");
   m_params.idxSelMemPhiLUTParams()->setPhiInWidth(ismp_lut_stgs_.getParameter<int>("phi_in_width"));
   m_params.idxSelMemPhiLUTParams()->setOutWidth(ismp_lut_stgs_.getParameter<int>("out_width"));
   m_params.idxSelMemPhiLUTParams()->setFilename(ismp_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet fpsmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("FwdPosSingleMatchQualLUTSettings");
   m_params.fwdPosSingleMatchQualLUTParams()->setDeltaEtaRedInWidth(fpsmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.fwdPosSingleMatchQualLUTParams()->setDeltaPhiRedInWidth(fpsmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.fwdPosSingleMatchQualLUTParams()->setOutWidth(fpsmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.fwdPosSingleMatchQualLUTParams()->setFilename(fpsmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet bonmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("BONegMatchQualLUTSettings");
   m_params.bONegMatchQualLUTParams()->setDeltaEtaRedInWidth(bonmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.bONegMatchQualLUTParams()->setDeltaPhiRedInWidth(bonmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.bONegMatchQualLUTParams()->setOutWidth(bonmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.bONegMatchQualLUTParams()->setFilename(bonmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet onsmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("OvlNegSingleMatchQualLUTSettings");
   m_params.ovlNegSingleMatchQualLUTParams()->setDeltaEtaRedInWidth(onsmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.ovlNegSingleMatchQualLUTParams()->setDeltaPhiRedInWidth(onsmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.ovlNegSingleMatchQualLUTParams()->setOutWidth(onsmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.ovlNegSingleMatchQualLUTParams()->setFilename(onsmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet isme_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("IdxSelMemEtaLUTSettings");
   m_params.idxSelMemEtaLUTParams()->setEtaInWidth(isme_lut_stgs_.getParameter<int>("eta_in_width"));
   m_params.idxSelMemEtaLUTParams()->setOutWidth(isme_lut_stgs_.getParameter<int>("out_width"));
   m_params.idxSelMemEtaLUTParams()->setFilename(isme_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet fopmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("FOPosMatchQualLUTSettings");
   m_params.fOPosMatchQualLUTParams()->setDeltaEtaRedInWidth(fopmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.fOPosMatchQualLUTParams()->setDeltaPhiRedInWidth(fopmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.fOPosMatchQualLUTParams()->setOutWidth(fopmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.fOPosMatchQualLUTParams()->setFilename(fopmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet fnsmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("FwdNegSingleMatchQualLUTSettings");
   m_params.fwdNegSingleMatchQualLUTParams()->setDeltaEtaRedInWidth(fnsmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.fwdNegSingleMatchQualLUTParams()->setDeltaPhiRedInWidth(fnsmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.fwdNegSingleMatchQualLUTParams()->setOutWidth(fnsmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.fwdNegSingleMatchQualLUTParams()->setFilename(fnsmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet bpe_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("BPhiExtrapolationLUTSettings");
   m_params.bPhiExtrapolationLUTParams()->setEtaAbsRedInWidth(bpe_lut_stgs_.getParameter<int>("etaAbsRed_in_width"));
   m_params.bPhiExtrapolationLUTParams()->setPtRedInWidth(bpe_lut_stgs_.getParameter<int>("pTred_in_width"));
   m_params.bPhiExtrapolationLUTParams()->setOutWidth(bpe_lut_stgs_.getParameter<int>("out_width"));
   m_params.bPhiExtrapolationLUTParams()->setFilename(bpe_lut_stgs_.getParameter<std::string>("filename"));
   
   edm::ParameterSet bsmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("BrlSingleMatchQualLUTSettings");
   m_params.brlSingleMatchQualLUTParams()->setDeltaEtaRedInWidth(bsmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.brlSingleMatchQualLUTParams()->setDeltaPhiRedInWidth(bsmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.brlSingleMatchQualLUTParams()->setOutWidth(bsmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.brlSingleMatchQualLUTParams()->setFilename(bsmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet ricm_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("RelIsoCheckMemLUTSettings");
   m_params.relIsoCheckMemLUTParams()->setAreaSumInWidth(ricm_lut_stgs_.getParameter<int>("areaSum_in_width"));
   m_params.relIsoCheckMemLUTParams()->setPtInWidth(ricm_lut_stgs_.getParameter<int>("pT_in_width"));
   m_params.relIsoCheckMemLUTParams()->setOutWidth(ricm_lut_stgs_.getParameter<int>("out_width"));
   m_params.relIsoCheckMemLUTParams()->setFilename(ricm_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet ope_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("OPhiExtrapolationLUTSettings");
   m_params.oPhiExtrapolationLUTParams()->setEtaAbsRedInWidth(ope_lut_stgs_.getParameter<int>("etaAbsRed_in_width"));
   m_params.oPhiExtrapolationLUTParams()->setPtRedInWidth(ope_lut_stgs_.getParameter<int>("pTred_in_width"));
   m_params.oPhiExtrapolationLUTParams()->setOutWidth(ope_lut_stgs_.getParameter<int>("out_width"));
   m_params.oPhiExtrapolationLUTParams()->setFilename(ope_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet opsmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("OvlPosSingleMatchQualLUTSettings");
   m_params.ovlPosSingleMatchQualLUTParams()->setDeltaEtaRedInWidth(opsmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.ovlPosSingleMatchQualLUTParams()->setDeltaPhiRedInWidth(opsmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.ovlPosSingleMatchQualLUTParams()->setOutWidth(opsmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.ovlPosSingleMatchQualLUTParams()->setFilename(opsmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet fee_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("FEtaExtrapolationLUTSettings");
   m_params.fEtaExtrapolationLUTParams()->setEtaAbsRedInWidth(fee_lut_stgs_.getParameter<int>("etaAbsRed_in_width"));
   m_params.fEtaExtrapolationLUTParams()->setPtRedInWidth(fee_lut_stgs_.getParameter<int>("pTred_in_width"));
   m_params.fEtaExtrapolationLUTParams()->setOutWidth(fee_lut_stgs_.getParameter<int>("out_width"));
   m_params.fEtaExtrapolationLUTParams()->setFilename(fee_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet bopmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("BOPosMatchQualLUTSettings");
   m_params.bOPosMatchQualLUTParams()->setDeltaEtaRedInWidth(bopmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.bOPosMatchQualLUTParams()->setDeltaPhiRedInWidth(bopmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.bOPosMatchQualLUTParams()->setOutWidth(bopmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.bOPosMatchQualLUTParams()->setFilename(bopmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet oee_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("OEtaExtrapolationLUTSettings");
   m_params.oEtaExtrapolationLUTParams()->setEtaAbsRedInWidth(oee_lut_stgs_.getParameter<int>("etaAbsRed_in_width"));
   m_params.oEtaExtrapolationLUTParams()->setPtRedInWidth(oee_lut_stgs_.getParameter<int>("pTred_in_width"));
   m_params.oEtaExtrapolationLUTParams()->setOutWidth(oee_lut_stgs_.getParameter<int>("out_width"));
   m_params.oEtaExtrapolationLUTParams()->setFilename(oee_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet bee_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("BEtaExtrapolationLUTSettings");
   m_params.bEtaExtrapolationLUTParams()->setEtaAbsRedInWidth(bee_lut_stgs_.getParameter<int>("etaAbsRed_in_width"));
   m_params.bEtaExtrapolationLUTParams()->setPtRedInWidth(bee_lut_stgs_.getParameter<int>("pTred_in_width"));
   m_params.bEtaExtrapolationLUTParams()->setOutWidth(bee_lut_stgs_.getParameter<int>("out_width"));
   m_params.bEtaExtrapolationLUTParams()->setFilename(bee_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet fpe_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("FPhiExtrapolationLUTSettings");
   m_params.fPhiExtrapolationLUTParams()->setEtaAbsRedInWidth(fpe_lut_stgs_.getParameter<int>("etaAbsRed_in_width"));
   m_params.fPhiExtrapolationLUTParams()->setPtRedInWidth(fpe_lut_stgs_.getParameter<int>("pTred_in_width"));
   m_params.fPhiExtrapolationLUTParams()->setOutWidth(fpe_lut_stgs_.getParameter<int>("out_width"));
   m_params.fPhiExtrapolationLUTParams()->setFilename(fpe_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet fonmq_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("FONegMatchQualLUTSettings");
   m_params.fONegMatchQualLUTParams()->setDeltaEtaRedInWidth(fonmq_lut_stgs_.getParameter<int>("deltaEtaRed_in_width"));
   m_params.fONegMatchQualLUTParams()->setDeltaPhiRedInWidth(fonmq_lut_stgs_.getParameter<int>("deltaPhiRed_in_width"));
   m_params.fONegMatchQualLUTParams()->setOutWidth(fonmq_lut_stgs_.getParameter<int>("out_width"));
   m_params.fONegMatchQualLUTParams()->setFilename(fonmq_lut_stgs_.getParameter<std::string>("filename"));

   edm::ParameterSet sr_lut_stgs_ = iConfig.getParameter<edm::ParameterSet>("SortRankLUTSettings");
   m_params.sortRankLUTParams()->setPtInWidth(sr_lut_stgs_.getParameter<int>("pT_in_width"));
   m_params.sortRankLUTParams()->setQualInWidth(sr_lut_stgs_.getParameter<int>("qual_in_width"));
   m_params.sortRankLUTParams()->setOutWidth(sr_lut_stgs_.getParameter<int>("out_width"));
   m_params.sortRankLUTParams()->setFilename(sr_lut_stgs_.getParameter<std::string>("filename"));
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
