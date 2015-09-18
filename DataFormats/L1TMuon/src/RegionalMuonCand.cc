#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

namespace l1t {

void
RegionalMuonCand::setTFIdentifiers(int processor, tftype trackFinder) {
  int linkOffset = 36;
  m_trackFinder = trackFinder;
  m_processor = processor;

  switch (m_trackFinder) {
    case tftype::emtf_pos:
      m_link = m_processor + linkOffset;  // range 36...41
    case tftype::omtf_pos:
      m_link = m_processor + linkOffset;  // range 42...47
    case tftype::bmtf:
      m_link = m_processor + linkOffset;  // range 48...59
    case tftype::omtf_neg:
      m_link = m_processor + linkOffset;  // range 60...65
    case tftype::emtf_neg:
      m_link = m_processor + linkOffset;  // range 66...71
  }
}

} // namespace l1t
