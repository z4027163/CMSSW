#include "../interface/MicroGMTCancelOutUnit.h"
#include "DataFormats/L1TMuon/interface/L1TGMTInternalMuon.h"

namespace l1t {
MicroGMTCancelOutUnit::MicroGMTCancelOutUnit (const edm::ParameterSet& iConfig) : 
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

MicroGMTCancelOutUnit::~MicroGMTCancelOutUnit ()
{

}

void
MicroGMTCancelOutUnit::setCancelOutBits(L1TGMTInternalWedges& wedges, tftype trackFinder, cancelmode mode) 
{ 
  std::vector<std::shared_ptr<L1TGMTInternalMuon>> coll1;
  coll1.reserve(3);
  std::vector<std::shared_ptr<L1TGMTInternalMuon>> coll2;
  coll2.reserve(3);
  int maxWedges = 12;
  if (trackFinder == bmtf) {
    maxWedges = 6;
  }

  for (int currentWedge = 1; currentWedge <= maxWedges; ++currentWedge) {
    for (auto mu : wedges[currentWedge]) {
      coll1.push_back(mu);
    }
    // handle wrap around: max "wedge" has to be compared to first "wedge"
    int neighbourWedge = ((currentWedge - 1) % maxWedges) + 1;
    for (auto mu : wedges[neighbourWedge]) {
      coll2.push_back(mu);
    }
    if (mode == cancelmode::coordinate) {
      getCoordinateCancelBits(coll1, coll2);
    } else {
      getTrackAddrCancelBits(coll1, coll2);
    }    coll1.clear();
    coll2.clear();
  }
}

void
MicroGMTCancelOutUnit::setCancelOutBitsOverlapBarrel(L1TGMTInternalWedges& omtfSectors, L1TGMTInternalWedges& bmtfWedges, cancelmode mode) 
{
  // overlap sector collection
  std::vector<std::shared_ptr<L1TGMTInternalMuon>> coll1;
  coll1.reserve(3);
  // barrel wedge collection with 4 wedges
  std::vector<std::shared_ptr<L1TGMTInternalMuon>> coll2;
  coll2.reserve(12);

  for (int currentSector = 1; currentSector <= 6; ++currentSector) {
    for (auto omtfMuon : omtfSectors[currentSector]) {
      coll1.push_back(omtfMuon);
    }
    // BMTF | 2  | 3  | 4  | 5  | 6  | 7  | 8  | 9  | 10 | 11 | 12 | 1  |
    // OMTF |    1    |    2    |    3    |    4    |    5    |    6    |
    // cancel OMTF sector x with corresponding BMTF wedge + the two on either side;
    // e.g. OMTF 1 with BMTF 1, 2, 3, 4, OMTF 2 with BMTF 3, 4, 5, 6 etc.
    for (int i = 0; i < 3; ++i) { 
      int currentWedge = currentSector * 2 - 1 + i;
      // handling the wrap-around: doing a shift by one for the modulo
      // as the wedge numbering starts at 1 instead of 0
      currentWedge = (currentWedge - 1) % 12 + 1;
      for (auto bmtfMuon : bmtfWedges[currentWedge]) {
        coll2.push_back(bmtfMuon);
      }
    }
    if (mode == cancelmode::coordinate) {
      getCoordinateCancelBits(coll1, coll2);
    } else {
      getTrackAddrCancelBits(coll1, coll2);
    }
    coll1.clear();
    coll2.clear();
  }
}

void
MicroGMTCancelOutUnit::setCancelOutBitsOverlapEndcap(L1TGMTInternalWedges& omtfSectors, L1TGMTInternalWedges& emtfSectors, cancelmode mode) 
{
  // overlap sector collection
  std::vector<std::shared_ptr<L1TGMTInternalMuon>> coll1;
  coll1.reserve(3);
  // endcap sector collection with 3 sectors
  std::vector<std::shared_ptr<L1TGMTInternalMuon>> coll2;
  coll2.reserve(9);

  for (int curOmtfSector = 1; curOmtfSector <= 6; ++curOmtfSector) {
    for (auto omtfMuon : omtfSectors[curOmtfSector]) {
      coll1.push_back(omtfMuon);
    }
    // OMTF |    1    |    2    |    3    |    4    |    5    |    6    |
    // EMTF |    1    |    2    |    3    |    4    |    5    |    6    |
    // cancel OMTF sector x with corresponding EMTF sector + the ones on either side;
    // e.g. OMTF 1 with EMTF 6, 1, 2; OMTF 2 with EMTF 1, 2, 3 etc.
    for (int i = 0; i < 2; ++i) {
      // handling the wrap around: doing shift by 6 (because 1 has to be compared to 6)
      // and the additional shift by one as above because of 1-indexed processor IDs
      int curEmtfSector = ((curOmtfSector + 6) - 1 + i) % 6 + 1;
      for (auto emtfMuon : emtfSectors[curEmtfSector]) {
        coll2.push_back(emtfMuon);
      }
    }
    if (mode == cancelmode::coordinate) {
      getCoordinateCancelBits(coll1, coll2);
    } else {
      getTrackAddrCancelBits(coll1, coll2);
    }    coll1.clear();
    coll2.clear();
  }

}

void 
MicroGMTCancelOutUnit::getCoordinateCancelBits(std::vector<std::shared_ptr<L1TGMTInternalMuon>>& coll1, std::vector<std::shared_ptr<L1TGMTInternalMuon>>& coll2)
{
  MicroGMTMatchQualLUT* matchLUT = m_lutDict[(*coll1.begin())->trackFinderType()+(*coll2.begin())->trackFinderType()*5];
  for (auto mu_w1 = coll1.begin(); mu_w1 != coll1.end(); ++mu_w1) {
    for (auto mu_w2 = coll2.begin(); mu_w2 != coll2.end(); ++mu_w2) {
      // phi coordinates shall be relative, do not have to worry about wrap around...
      int deltaPhi = std::abs((*mu_w1)->hwLocalPhi() - (*mu_w2)->hwLocalPhi()) >> (8 - matchLUT->getDeltaPhiWidth()); //diffbits = origwidth - widthweneed
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

void 
MicroGMTCancelOutUnit::getTrackAddrCancelBits(std::vector<std::shared_ptr<L1TGMTInternalMuon>>& coll1, std::vector<std::shared_ptr<L1TGMTInternalMuon>>& coll2)
{
  // not entirely clear how to do.. just a hook for now
}

} // namespace l1t