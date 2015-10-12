#include "sstream"

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"

#include "CondFormats/L1TObjects/interface/L1TMTFOverlapParams.h"
#include "CondFormats/DataRecord/interface/L1TMTFOverlapParamsRcd.h"

#include "L1Trigger/L1TMuonTrackFinderOverlap/interface/XMLConfigReader.h"
#include "L1Trigger/L1TMuonTrackFinderOverlap/plugins/L1TMTFOverlapParamsESProducer.h"

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
L1TMTFOverlapParamsESProducer::L1TMTFOverlapParamsESProducer(const edm::ParameterSet& theConfig){
   //the following line is needed to tell the framework what
   // data is being produced
   setWhatProduced(this);

   if ( !theConfig.exists("patternsXMLFiles") ) return;
   std::vector<std::string> fileNames;
   for(auto it: theConfig.getParameter<std::vector<edm::ParameterSet> >("patternsXMLFiles")){
     fileNames.push_back(it.getParameter<edm::FileInPath>("patternsXMLFile").fullPath());
   }  


   myOMTFConfig = new OMTFConfiguration(theConfig);
   
   XMLConfigReader myReader;
   for(auto it: fileNames){
     myReader.setPatternsFile(it);
     readXML(&myReader);
   }  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
L1TMTFOverlapParamsESProducer::~L1TMTFOverlapParamsESProducer() {}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
bool L1TMTFOverlapParamsESProducer::readXML(XMLConfigReader *aReader){

  

  std::cout<<"L1TMTFOverlapParamsESProducer::readXML BEGIN"<<std::endl;
  
  std::cout<<"L1TMTFOverlapParamsESProducer::readXML charge"<<std::endl;

  l1t::LUT chargeLUT;
  aReader->readLUT(&chargeLUT,"iCharge");
  m_params.setChargeLUT(chargeLUT);

  std::cout<<"L1TMTFOverlapParamsESProducer::readXML eta"<<std::endl;

  l1t::LUT etaLUT;
  aReader->readLUT(&etaLUT,"iEta");
  m_params.setEtaLUT(etaLUT);

  std::cout<<"L1TMTFOverlapParamsESProducer::readXML pt"<<std::endl;

  l1t::LUT ptLUT;
  aReader->readLUT(&ptLUT,"iPt");
  m_params.setPtLUT(ptLUT);

  std::cout<<"L1TMTFOverlapParamsESProducer::readXML meanDistPhi"<<std::endl;

  l1t::LUT meanDistPhiLUT;
  aReader->readLUT(&meanDistPhiLUT,"meanDistPhi");
  m_params.setMeanDistPhiLUT(meanDistPhiLUT);
  
  std::cout<<"L1TMTFOverlapParamsESProducer::readXML pdf"<<std::endl;

  l1t::LUT pdfLUT;
  aReader->readLUT(&pdfLUT,"pdf");
  m_params.setPdfLUT(pdfLUT);
  
  std::cout<<"L1TMTFOverlapParamsESProducer::readXML DONE"<<std::endl;
  
  return true;
  
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
L1TMTFOverlapParamsESProducer::ReturnType
L1TMTFOverlapParamsESProducer::produce(const L1TMTFOverlapParamsRcd& iRecord)
{
   using namespace edm::es;
   boost::shared_ptr<L1TMTFOverlapParams> aL1TMTFOverlapParams;

   aL1TMTFOverlapParams = boost::shared_ptr<L1TMTFOverlapParams>(new L1TMTFOverlapParams(m_params));
   return aL1TMTFOverlapParams;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

