///
/// \class l1t::MicroGMTParams
///
/// Description: Placeholder for MicroGMT parameters
///
/// Implementation:
///
/// \author: Thomas Reis
///

#ifndef MicroGMTParams_h
#define MicroGMTParams_h

#include <memory>
#include <iostream>
#include <vector>

#include "CondFormats/Serialization/interface/Serializable.h"

namespace l1t {
  class MicroGMTParams {

  public:
    enum { Version = 1 };

    // general LUT parameters
    class LUTParams {
    public:
      LUTParams() : out_width_(0), filename_("") {}
      
      int outWidth() const { return out_width_; }
      std::string filename() const { return filename_; }
      void setOutWidth(int width) { out_width_ = width; }
      void setFilename(std::string filename) { filename_ = filename; }

    private:
      int out_width_;
      std::string filename_;

      COND_SERIALIZABLE;
    };

    // MatchQual specific LUT parameters
    class MatchQualLUTParams : public LUTParams {
    public:
      MatchQualLUTParams() : deltaEtaRed_in_width_(0), deltaPhiRed_in_width_(0) {}

      int deltaEtaRedInWidth() const { return deltaEtaRed_in_width_; }
      int deltaPhiRedInWidth() const { return deltaPhiRed_in_width_; }
      void setDeltaEtaRedInWidth(int width) { deltaEtaRed_in_width_ = width; }
      void setDeltaPhiRedInWidth(int width) { deltaPhiRed_in_width_ = width; }

    private:
      int deltaEtaRed_in_width_;
      int deltaPhiRed_in_width_;

      COND_SERIALIZABLE;
    };

    // Extrapolation specific LUT parameters
    class ExtrapolationLUTParams : public LUTParams {
    public:
      ExtrapolationLUTParams() : etaAbsRed_in_width_(0), pTred_in_width_(0) {}

      int etaAbsRedInWidth() const { return etaAbsRed_in_width_; }
      int ptRedInWidth() const { return pTred_in_width_; }
      void setEtaAbsRedInWidth(int width) { etaAbsRed_in_width_ = width; }
      void setPtRedInWidth(int width) { pTred_in_width_ = width; }

    private:
      int etaAbsRed_in_width_;
      int pTred_in_width_;

      COND_SERIALIZABLE;
    };

    // IsoCheckMem specific LUT parameters
    class IsoCheckMemLUTParams : public LUTParams {
    public:
      IsoCheckMemLUTParams() : areaSum_in_width_(0) {}

      int areaSumInWidth() const { return areaSum_in_width_; }
      void setAreaSumInWidth(int width) { areaSum_in_width_ = width; }

    private:
      int areaSum_in_width_;

      COND_SERIALIZABLE;
    };

    // Rel. IsoCheckMem specific LUT parameters
    class RelIsoCheckMemLUTParams : public IsoCheckMemLUTParams {
    public:
      RelIsoCheckMemLUTParams() : pT_in_width_(0) {}

      int ptInWidth() const { return pT_in_width_; }
      void setPtInWidth(int width) { pT_in_width_ = width; }

    private:
      int pT_in_width_;

      COND_SERIALIZABLE;
    };

    // SelMemEta specific LUT parameters
    class SelMemEtaLUTParams : public LUTParams {
    public:
      SelMemEtaLUTParams() : eta_in_width_(0) {}

      int etaInWidth() const { return eta_in_width_; }
      void setEtaInWidth(int width) { eta_in_width_ = width; }

    private:
      int eta_in_width_;

      COND_SERIALIZABLE;
    };

    // SelMemPhi specific LUT parameters
    class SelMemPhiLUTParams : public LUTParams {
    public:
      SelMemPhiLUTParams() : phi_in_width_(0) {}

      int phiInWidth() const { return phi_in_width_; }
      void setPhiInWidth(int width) { phi_in_width_ = width; }

    private:
      int phi_in_width_;

      COND_SERIALIZABLE;
    };

    // SortRank specific LUT parameters
    class SortRankLUTParams : public LUTParams {
    public:
      SortRankLUTParams() : pT_in_width_(0), qual_in_width_(0) {}

      int ptInWidth() const { return pT_in_width_; }
      int qualInWidth() const { return qual_in_width_; }
      void setPtInWidth(int width) { pT_in_width_ = width; }
      void setQualInWidth(int width) { qual_in_width_ = width; }

    private:
      int pT_in_width_;
      int qual_in_width_;

      COND_SERIALIZABLE;
    };

    MicroGMTParams() { version_=Version; }
    ~MicroGMTParams() {}

    // FW version
    unsigned fwVersion() const { return fwVersion_; }
    void setFwVersion(unsigned fwVersion) { fwVersion_ = fwVersion; }

    // access to LUT parameter objects
    MatchQualLUTParams * brlSingleMatchQualLUTParams() { return &bsinglemqp_; }
    MatchQualLUTParams * fwdNegSingleMatchQualLUTParams() { return &fnegmqp_; }
    MatchQualLUTParams * fwdPosSingleMatchQualLUTParams() { return &fposmqp_; }
    MatchQualLUTParams * ovlNegSingleMatchQualLUTParams() { return &onegmqp_; }
    MatchQualLUTParams * ovlPosSingleMatchQualLUTParams() { return &oposmqp_; }
    MatchQualLUTParams * bONegMatchQualLUTParams() { return &bonegmqp_; }
    MatchQualLUTParams * bOPosMatchQualLUTParams() { return &boposmqp_; }
    MatchQualLUTParams * fONegMatchQualLUTParams() { return &fonegmqp_; }
    MatchQualLUTParams * fOPosMatchQualLUTParams() { return &foposmqp_; }
    ExtrapolationLUTParams * bPhiExtrapolationLUTParams() { return &bphiepp_; }
    ExtrapolationLUTParams * oPhiExtrapolationLUTParams() { return &ophiepp_; }
    ExtrapolationLUTParams * fPhiExtrapolationLUTParams() { return &fphiepp_; }
    ExtrapolationLUTParams * bEtaExtrapolationLUTParams() { return &betaepp_; }
    ExtrapolationLUTParams * oEtaExtrapolationLUTParams() { return &oetaepp_; }
    ExtrapolationLUTParams * fEtaExtrapolationLUTParams() { return &fetaepp_; }
    IsoCheckMemLUTParams * absIsoCheckMemLUTParams() { return &isocmp_; }
    RelIsoCheckMemLUTParams * relIsoCheckMemLUTParams() { return &risocmp_; }
    SelMemEtaLUTParams * idxSelMemEtaLUTParams() { return &smetap_; }
    SelMemPhiLUTParams * idxSelMemPhiLUTParams() { return &smphip_; }
    SortRankLUTParams * sortRankLUTParams() { return &srp_; }

    // print parameters to stream:
    void print(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream& o, const MicroGMTParams & p) { p.print(o); return o; }

  private:
    unsigned version_;
    unsigned fwVersion_;

    // LUT parameter objects
    MatchQualLUTParams bsinglemqp_;
    MatchQualLUTParams fnegmqp_;
    MatchQualLUTParams fposmqp_;
    MatchQualLUTParams onegmqp_;
    MatchQualLUTParams oposmqp_;
    MatchQualLUTParams bonegmqp_;
    MatchQualLUTParams boposmqp_;
    MatchQualLUTParams fonegmqp_;
    MatchQualLUTParams foposmqp_;
    ExtrapolationLUTParams bphiepp_;
    ExtrapolationLUTParams ophiepp_;
    ExtrapolationLUTParams fphiepp_;
    ExtrapolationLUTParams betaepp_;
    ExtrapolationLUTParams oetaepp_;
    ExtrapolationLUTParams fetaepp_;
    IsoCheckMemLUTParams isocmp_;
    RelIsoCheckMemLUTParams risocmp_;
    SelMemEtaLUTParams smetap_;
    SelMemPhiLUTParams smphip_;
    SortRankLUTParams srp_;

    COND_SERIALIZABLE;
  };

}// namespace
#endif
