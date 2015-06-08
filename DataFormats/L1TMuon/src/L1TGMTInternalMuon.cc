#include "DataFormats/L1TMuon/interface/L1TGMTInternalMuon.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

namespace l1t {

L1TGMTInternalMuon::L1TGMTInternalMuon(const L1TRegionalMuonCandidate& regional) : 
	m_regional(regional), m_hwRank(0), m_hwCancelBit(0), m_hwIsoSum(0), m_hwDeltaEta(0), m_hwDeltaPhi(0), m_hwAbsIso(-1), m_hwRelIso(-1), m_link(-1)
{
  initLink();
}

void
L1TGMTInternalMuon::initLink() {
  switch (this->trackFinderType()) {
    case tftype::bmtf:
      m_link = this->processor() + 35;  // range 36...47
    case tftype::omtf_pos:
      m_link = this->processor() + 47;  // range 48...53
    case tftype::omtf_neg:
      m_link = this->processor() + 53;  // range 54...59
    case tftype::emtf_pos:
      m_link = this->processor() + 59;  // range 60...65
    case tftype::emtf_neg:
      m_link = this->processor() + 65;  // range 66...71
  }
}

} // namespace l1t