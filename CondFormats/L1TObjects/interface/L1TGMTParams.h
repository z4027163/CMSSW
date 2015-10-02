///
/// \class L1TGMTParams
///
/// Description: Placeholder for MicroGMT parameters
///
/// Implementation:
///
/// \author: Thomas Reis
///

#ifndef L1TGMTParams_h
#define L1TGMTParams_h

#include <memory>
#include <iostream>
#include <vector>

#include "CondFormats/Serialization/interface/Serializable.h"

class L1TGMTParams {

public:
  enum { Version = 1 };

  L1TGMTParams() { version_=Version; }
  ~L1TGMTParams() {}

  // FW version
  unsigned fwVersion() const { return fwVersion_; }
  void setFwVersion(unsigned fwVersion) { fwVersion_ = fwVersion; }

  // LUT paths
  std::string absIsoCheckMemLUTPath()        { return aisocmlp_; }
  std::string relIsoCheckMemLUTPath()        { return risocmlp_; }
  std::string idxSelMemPhiLUTPath()          { return ismphilp_; }
  std::string idxSelMemEtaLUTPath()          { return ismetalp_; }
  std::string brlSingleMatchQualLUTPath()    { return bsinglemqlp_; }
  std::string fwdPosSingleMatchQualLUTPath() { return fposmqlp_; }
  std::string fwdNegSingleMatchQualLUTPath() { return fnegmqlp_; }
  std::string ovlPosSingleMatchQualLUTPath() { return oposmqlp_; }
  std::string ovlNegSingleMatchQualLUTPath() { return onegmqlp_; }
  std::string bOPosMatchQualLUTPath()        { return boposmqlp_; }
  std::string bONegMatchQualLUTPath()        { return bonegmqlp_; }
  std::string fOPosMatchQualLUTPath()        { return foposmqlp_; }
  std::string fONegMatchQualLUTPath()        { return fonegmqlp_; }
  std::string bPhiExtrapolationLUTPath()     { return bphieplp_; }
  std::string oPhiExtrapolationLUTPath()     { return ophieplp_; }
  std::string fPhiExtrapolationLUTPath()     { return fphieplp_; }
  std::string bEtaExtrapolationLUTPath()     { return betaeplp_; }
  std::string oEtaExtrapolationLUTPath()     { return oetaeplp_; }
  std::string fEtaExtrapolationLUTPath()     { return fetaeplp_; }
  std::string sortRankLUTPath()              { return srlp_; }
  void setAbsIsoCheckMemLUTPath        (std::string path) { aisocmlp_ = path; }
  void setRelIsoCheckMemLUTPath        (std::string path) { risocmlp_ = path; }
  void setIdxSelMemPhiLUTPath          (std::string path) { ismphilp_ = path; }
  void setIdxSelMemEtaLUTPath          (std::string path) { ismetalp_ = path; }
  void setBrlSingleMatchQualLUTPath    (std::string path) { bsinglemqlp_ = path; }
  void setFwdPosSingleMatchQualLUTPath (std::string path) { fposmqlp_ = path; }
  void setFwdNegSingleMatchQualLUTPath (std::string path) { fnegmqlp_ = path; }
  void setOvlPosSingleMatchQualLUTPath (std::string path) { oposmqlp_ = path; }
  void setOvlNegSingleMatchQualLUTPath (std::string path) { onegmqlp_ = path; }
  void setBOPosMatchQualLUTPath        (std::string path) { boposmqlp_ = path; }
  void setBONegMatchQualLUTPath        (std::string path) { bonegmqlp_ = path; }
  void setFOPosMatchQualLUTPath        (std::string path) { foposmqlp_ = path; }
  void setFONegMatchQualLUTPath        (std::string path) { fonegmqlp_ = path; }
  void setBPhiExtrapolationLUTPath     (std::string path) { bphieplp_ = path; }
  void setOPhiExtrapolationLUTPath     (std::string path) { ophieplp_ = path; }
  void setFPhiExtrapolationLUTPath     (std::string path) { fphieplp_ = path; }
  void setBEtaExtrapolationLUTPath     (std::string path) { betaeplp_ = path; }
  void setOEtaExtrapolationLUTPath     (std::string path) { oetaeplp_ = path; }
  void setFEtaExtrapolationLUTPath     (std::string path) { fetaeplp_ = path; }
  void setSortRankLUTPath              (std::string path) { srlp_ = path; }

  // print parameters to stream:
  void print(std::ostream&) const;
  friend std::ostream& operator<<(std::ostream& o, const L1TGMTParams & p) { p.print(o); return o; }

private:
  unsigned version_;
  unsigned fwVersion_;

  // LUT parameter objects
  std::string aisocmlp_;
  std::string risocmlp_;
  std::string ismetalp_;
  std::string ismphilp_;
  std::string bsinglemqlp_;
  std::string fposmqlp_;
  std::string fnegmqlp_;
  std::string oposmqlp_;
  std::string onegmqlp_;
  std::string boposmqlp_;
  std::string bonegmqlp_;
  std::string foposmqlp_;
  std::string fonegmqlp_;
  std::string bphieplp_;
  std::string ophieplp_;
  std::string fphieplp_;
  std::string betaeplp_;
  std::string oetaeplp_;
  std::string fetaeplp_;
  std::string srlp_;

  COND_SERIALIZABLE;
};
#endif
