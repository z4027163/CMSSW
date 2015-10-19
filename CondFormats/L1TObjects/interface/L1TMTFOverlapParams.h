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
  enum { CHARGE=0, ETA=1, PT=2, PDF=3, MEANDISTPHI=4,
	 GENERAL = 5, HWLAYERS=6, BENDLAYERS = 7, CONNECTLAYERS=8,
	 NUM_OMTFPARAMNODES=9};

  // General configuration parameters indexes
  enum {GENERAL_ADDRBITS=0, GENERAL_VALBITS=1, GENERAL_HITSPERLAYER=2, GENERAL_PHIBITS=3, GENERAL_PHIBINS=4, GENERAL_NREFHITS=5, GENERAL_NTESTREFHITS=6};
	
  
  L1TMTFOverlapParams() { fwVersion_=Version; pnodes_.resize(NUM_OMTFPARAMNODES); }
  ~L1TMTFOverlapParams() {}

  //<GlobalData minPdfVal="0.001" nPdfAddrBits="7" nPdfValBits="6" nHitsPerLayer="6" nPhiBits="11" nPhiBins="5760" nRefHits="128" nTestRefHits="4"/>

  // Firmware version
  unsigned fwVersion() const { return fwVersion_; }
  void setFwVersion(unsigned fwVersion) { fwVersion_ = fwVersion; }

  ///General definitions
  l1t::LUT* generalLUT()        { return &pnodes_[GENERAL].LUT_; }
  void setGeneralLUT (const l1t::LUT & lut) { pnodes_[GENERAL].type_ = "INT"; pnodes_[GENERAL].LUT_ = lut; }

  ///Access to specific general settings.
  int nPdfAddrBits() { return pnodes_[GENERAL].iparams_[GENERAL_ADDRBITS];};

  int nPdfValBits() { return pnodes_[GENERAL].iparams_[GENERAL_VALBITS];};

  int nHitsPerLayer() { return pnodes_[GENERAL].iparams_[GENERAL_HITSPERLAYER];};

  int nPhiBits() { return pnodes_[GENERAL].iparams_[GENERAL_PHIBITS];};

  int nPhiBins() { return pnodes_[GENERAL].iparams_[GENERAL_PHIBINS];};

  int nRefHits() { return pnodes_[GENERAL].iparams_[GENERAL_NREFHITS];};
    
  int nTestRefHits() { return pnodes_[GENERAL].iparams_[GENERAL_NTESTREFHITS];};
  
  ///Connections definitions
  l1t::LUT* hwLayersLUT()        { return &pnodes_[HWLAYERS].LUT_; }
  void setHwLayersLUT (const l1t::LUT & lut) { pnodes_[HWLAYERS].type_ = "LUT"; pnodes_[HWLAYERS].LUT_ = lut; }

  l1t::LUT* bendLayersLUT()        { return &pnodes_[BENDLAYERS].LUT_; }
  void setBendLayersLUT (const l1t::LUT & lut) { pnodes_[BENDLAYERS].type_ = "LUT"; pnodes_[BENDLAYERS].LUT_ = lut; }

  l1t::LUT* connectLayersLUT()        { return &pnodes_[CONNECTLAYERS].LUT_; }
  void setConnectLayersLUT (const l1t::LUT & lut) { pnodes_[CONNECTLAYERS].type_ = "LUT"; pnodes_[CONNECTLAYERS].LUT_ = lut; }
  
  ///Golden Patterns definitions
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
  unsigned fwVersion_;
    
  ///vector of LUT like parameters
  std::vector<Node> pnodes_;
  
  COND_SERIALIZABLE;
};    
#endif
