#ifndef __l1t_regional_muon_candidate_h__
#define __l1t_regional_muon_candidate_h__

#include "L1TRegionalMuonCandidateFwd.h"
#include <iostream> 
namespace l1t {

class L1TRegionalMuonCandidate {
  public:
    L1TRegionalMuonCandidate() : 
      m_hwPt(0), m_hwPhi(0), m_hwEta(0), m_hwHF(false), m_hwSign(0), m_hwSignValid(0), m_hwQuality(0), 
      m_hwTrackAddress(0), m_link(0), m_processor(0), m_trackFinder(bmtf)
      {};

    L1TRegionalMuonCandidate(int pt, int phi, int eta, int sign, int signvalid, int quality, int processor, tftype trackFinder) : 
      m_hwPt(pt), m_hwPhi(phi), m_hwEta(eta), m_hwHF(false), m_hwSign(sign), m_hwSignValid(signvalid), m_hwQuality(quality), 
      m_hwTrackAddress(0), m_link(0), m_processor(processor), m_trackFinder(trackFinder)
      {};

    virtual ~L1TRegionalMuonCandidate() {};

    void setHwPt(int bits) { m_hwPt = bits; };
    void setHwPhi(int bits) { m_hwPhi = bits; };
    void setHwEta(int bits) { m_hwEta = bits; };
    void setHwSign(int bits) { m_hwSign = bits; };
    void setHwSignValid(int bits) { m_hwSignValid = bits; };
    void setHwQual(int bits) { m_hwQuality = bits; };
    void setHwHF(bool bit) { m_hwHF = bit; };
    void setHwTrackAddress(int bits) { m_hwTrackAddress = bits; };
    void setTFIdentifiers(int processor, tftype trackFinder);
    void setLink(int link);

    const int hwPt() const { return m_hwPt; };
    const int hwPhi() const { return m_hwPhi; };
    const int hwEta() const { return m_hwEta; };
    const int hwSign() const { return m_hwSign; };
    const int hwSignValid() const { return m_hwSignValid; };
    const int hwQual() const { return m_hwQuality; };
    const int hwTrackAddress() const { return m_hwTrackAddress; };
    const int link() const { return m_link; };
    const int processor() const { return m_processor; };
    const tftype trackFinderType() const { return m_trackFinder; };
  
  private:
    int m_hwPt;
    int m_hwPhi;
    int m_hwEta;
    bool m_hwHF;
    int m_hwSign;
    int m_hwSignValid;
    int m_hwQuality;
    int m_hwTrackAddress;
    int m_link;
    int m_processor;
    tftype m_trackFinder;
};

}

#endif /* define __l1t_regional_muon_candidate_h__ */