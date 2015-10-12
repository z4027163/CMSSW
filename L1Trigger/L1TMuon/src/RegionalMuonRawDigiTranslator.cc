#include "L1Trigger/L1TMuon/interface/RegionalMuonRawDigiTranslator.h"

void
l1t::RegionalMuonRawDigiTranslator::fillRegionalMuonCand(RegionalMuonCand& mu, uint32_t raw_data_00_31, uint32_t raw_data_32_63, int proc, tftype tf)
{
  // translations as defined in DN-15-017
  mu.setHwPt((raw_data_00_31 >> ptShift_) & ptWidth_);
  mu.setHwQual((raw_data_00_31 >> qualShift_) & qualWidth_); 

  // eta is coded as two's complement
  int abs_eta = (raw_data_00_31 >> absEtaShift_) & absEtaWidth_;
  if ((raw_data_00_31 >> etaSignShift_) & etaSignWidth_) {
     mu.setHwEta(abs_eta - 256);
  } else {
     mu.setHwEta(abs_eta);
  }

  mu.setHwPhi((raw_data_00_31 >> phiShift_) & phiWidth_);
  // sign is coded as -1^signBit
  int signBit = (raw_data_32_63 >> signShift_) & signWidth_;
  mu.setHwSign(1 - 2*signBit);
  mu.setHwSignValid((raw_data_32_63 >> signValidShift_) & signValidWidth_);
  mu.setHwHF((raw_data_00_31 >> hfShift_) & hfWidth_);
  mu.setHwTrackAddress((raw_data_32_63 >> trackAddressShift_) & trackAddressWidth_);
  mu.setTFIdentifiers(proc, tf);
  mu.setDataword(raw_data_32_63, raw_data_00_31);
}

void
l1t::RegionalMuonRawDigiTranslator::fillRegionalMuonCand(RegionalMuonCand& mu, uint64_t dataword, int proc, tftype tf)
{
  fillRegionalMuonCand(mu, (uint32_t)(dataword & 0xFFFFFFFF), (uint32_t)((dataword >> 32) & 0xFFFFFFFF), proc, tf);
}

void
l1t::RegionalMuonRawDigiTranslator::generatePackedDataWords(const RegionalMuonCand& mu, uint32_t &raw_data_00_31, uint32_t &raw_data_32_63)
{
  raw_data_00_31 = (mu.hwPt() & ptWidth_) << ptShift_
                 | (mu.hwQual() & qualWidth_) << qualShift_
                 | (abs(mu.hwEta()) & absEtaWidth_) << absEtaShift_
                 | ((mu.hwEta() < 0) & etaSignWidth_) << etaSignShift_
                 | (mu.hwHF() & hfWidth_) << hfShift_
                 | (mu.hwPhi() & phiWidth_) << phiShift_;

  raw_data_32_63 = ((mu.hwSign() > 0) & signWidth_) << signShift_
                 | (mu.hwSignValid() & signValidWidth_) << signValidShift_
                 | (mu.hwTrackAddress() & trackAddressWidth_) << trackAddressShift_;
}

uint64_t 
l1t::RegionalMuonRawDigiTranslator::generate64bitDataWord(const RegionalMuonCand& mu)
{
  uint32_t lsw;
  uint32_t msw;

  generatePackedDataWords(mu, lsw, msw);
  return (((uint64_t)msw) << 32) + lsw;
}

