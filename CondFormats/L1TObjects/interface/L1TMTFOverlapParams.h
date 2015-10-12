#ifndef L1TMTFOverlapParams_h
#define L1TMTFOverlapParams_h

#include <memory>
#include <iostream>
#include <vector>
#include <cmath>

#include "CondFormats/Serialization/interface/Serializable.h"
#include "CondFormats/L1TObjects/interface/LUT.h"

///////////////////////////////////////
///////////////////////////////////////
class L1TMTFOverlapParams {
  
 public:
  
  class Node {
  public:
    std::string type_;
    unsigned version_;
    l1t::LUT LUT_;
    std::vector<double> dparams_;
    std::vector<unsigned> uparams_;
    std::vector<int> iparams_;
    std::vector<std::string> sparams_;
    Node(){ type_="unspecified"; version_=0; }
    COND_SERIALIZABLE;
  };
  
  enum { Version = 1 };
  
  // DO NOT ADD ENTRIES ANYWHERE BUT DIRECTLY BEFORE "NUM_OMTFPARAMNODES"
  enum { CHARGE=0, ETA=1, PT=2, PDF=3, MEANDISTPHI=4, NUM_OMTFPARAMNODES=5};
  
  L1TMTFOverlapParams() { version_=Version; pnodes_.resize(NUM_OMTFPARAMNODES); }
  ~L1TMTFOverlapParams() {}

  l1t::LUT* chargeLUT()        { return &pnodes_[CHARGE].LUT_; }
  l1t::LUT* etaLUT()        { return &pnodes_[ETA].LUT_; }
  l1t::LUT* ptLUT()          { return &pnodes_[PT].LUT_; }
  l1t::LUT* pdfLUT()          { return &pnodes_[PDF].LUT_; }
  l1t::LUT* meanDistPhiLUT()          { return &pnodes_[MEANDISTPHI].LUT_; }

  void setChargeLUT (const l1t::LUT & lut) { pnodes_[CHARGE].type_ = "LUT"; pnodes_[CHARGE].LUT_ = lut; }
  void setEtaLUT (const l1t::LUT & lut) { pnodes_[ETA].type_ = "LUT"; pnodes_[ETA].LUT_ = lut; }
  void setPtLUT (const l1t::LUT & lut) { pnodes_[PT].type_ = "LUT"; pnodes_[PT].LUT_ = lut; }
  void setPdfLUT (const l1t::LUT & lut) { pnodes_[PDF].type_ = "LUT"; pnodes_[PDF].LUT_ = lut; }
  void setMeanDistPhiLUT (const l1t::LUT & lut) { pnodes_[MEANDISTPHI].type_ = "LUT"; pnodes_[MEANDISTPHI].LUT_ = lut; }
  
  
 private:
  unsigned version_;
    
  ///vector of LUT like parameters
  std::vector<Node> pnodes_;
  
  COND_SERIALIZABLE;
};    
#endif
