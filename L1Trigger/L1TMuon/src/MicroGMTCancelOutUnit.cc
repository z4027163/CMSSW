#include "../interface/MicroGMTCancelOutUnit.h"
#include "DataFormats/L1TMuon/interface/L1TGMTInternalMuon.h"

l1t::MicroGMTCancelOutUnit::MicroGMTCancelOutUnit (const edm::ParameterSet& iConfig) : 
    m_boPosMatchQualLUT(iConfig, "BOPos"),
    m_boNegMatchQualLUT(iConfig, "BONeg"),
    m_foPosMatchQualLUT(iConfig, "FOPos"),
    m_foNegMatchQualLUT(iConfig, "FONeg"),
    m_brlSingleMatchQualLUT(iConfig, "BrlSingle"),
    m_ovlPosSingleMatchQualLUT(iConfig, "OvlPosSingle"),
    m_ovlNegSingleMatchQualLUT(iConfig, "OvlNegSingle"),
    m_fwdPosSingleMatchQualLUT(iConfig, "FwdPosSingle"),
    m_fwdNegSingleMatchQualLUT(iConfig, "FwdNegSingle")
  {
    m_lutDict[tftype::bmtf+tftype::bmtf*5] = &m_brlSingleMatchQualLUT;
    m_lutDict[tftype::bmtf+tftype::omtf_neg*5] = &m_boNegMatchQualLUT;
    m_lutDict[tftype::bmtf+tftype::omtf_pos*5] = &m_boPosMatchQualLUT;
    m_lutDict[tftype::omtf_pos+tftype::omtf_pos*5] = &m_ovlPosSingleMatchQualLUT;
    m_lutDict[tftype::omtf_neg+tftype::omtf_neg*5] = &m_ovlNegSingleMatchQualLUT;
    m_lutDict[tftype::emtf_pos+tftype::emtf_pos*5] = &m_fwdPosSingleMatchQualLUT;
    m_lutDict[tftype::emtf_neg+tftype::emtf_neg*5] = &m_fwdNegSingleMatchQualLUT;
    m_lutDict[tftype::omtf_pos+tftype::emtf_pos*5] = &m_foPosMatchQualLUT;
    m_lutDict[tftype::omtf_neg+tftype::emtf_neg*5] = &m_foNegMatchQualLUT;


}

l1t::MicroGMTCancelOutUnit::~MicroGMTCancelOutUnit ()
{

}

void
l1t::MicroGMTCancelOutUnit::setCancelOutBits(MicroGMTConfiguration::InterMuonList& muons) 
{
  int cntr = 0;
  std::vector<MicroGMTConfiguration::InterMuonList::iterator> wedge1;
  wedge1.reserve(3);
  std::vector<MicroGMTConfiguration::InterMuonList::iterator> wedge2;
  MicroGMTConfiguration::InterMuonList::iterator mu;
  wedge2.reserve(3);
  for (mu = muons.begin(); mu != muons.end(); ++mu) {
    if ((cntr%6) < 3) {
      wedge1.push_back(mu);
    } else {
      wedge2.push_back(mu);
    }
    if (wedge1.size() == 3 && wedge2.size() == 3) {
      getCancelOutBits(wedge1, wedge2);
      wedge1.clear();
      wedge2.clear();
    }
    cntr++;
  }
}

void 
l1t::MicroGMTCancelOutUnit::getCancelOutBits( std::vector<MicroGMTConfiguration::InterMuonList::iterator> &wedge1, std::vector<MicroGMTConfiguration::InterMuonList::iterator> & wedge2)
{
  MicroGMTMatchQualLUT* matchLUT = m_lutDict[(*wedge1.begin())->trackFinderType()+(*wedge2.begin())->trackFinderType()*5];
  for (auto mu_w1 = wedge1.begin(); mu_w1 != wedge1.end(); ++mu_w1) {
    for (auto mu_w2 = wedge2.begin(); mu_w2 != wedge2.end(); ++mu_w2) {
      // phi coordinates shall be relative, do not have to worry about wrap around...
      int deltaPhi = std::abs((*mu_w1)->hwPhi() - (*mu_w2)->hwPhi()) >> (10 - matchLUT->getDeltaPhiWidth()); //diffbits = origwidth - widthweneed
      int deltaEta = std::abs((*mu_w1)->hwEta() - (*mu_w2)->hwEta()) >> (9 - matchLUT->getDeltaEtaWidth()); //diffbits = origwidth - widthweneed
      bool match = matchLUT->lookup(deltaEta, deltaPhi);
      if((*mu_w1)->hwQual() > (*mu_w2)->hwQual() && match) {
        (*mu_w2)->setHwCancelBit(1);
      } else {
        (*mu_w1)->setHwCancelBit(1);
      }
    }
  }
}