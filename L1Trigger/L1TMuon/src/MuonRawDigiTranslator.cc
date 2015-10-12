#include "L1Trigger/L1TMuon/interface/MuonRawDigiTranslator.h"

void
l1t::MuonRawDigiTranslator::fillMuon(Muon& mu, uint32_t raw_data_00_31, uint32_t raw_data_32_63)
{
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
  mu.setHwIso((raw_data_32_63 >> isoShift_) & isoWidth_); 
  // charge is coded as -1^chargeBit
  int chargeBit = (raw_data_32_63 >> chargeShift_) & chargeWidth_;
  mu.setHwCharge(1 - 2*chargeBit);
  mu.setHwChargeValid((raw_data_32_63 >> chargeValidShift_) & chargeValidWidth_);
}

void
l1t::MuonRawDigiTranslator::fillMuon(Muon& mu, uint64_t dataword)
{
  fillMuon(mu, (uint32_t)(dataword & 0xFFFFFFFF), (uint32_t)((dataword >> 32) & 0xFFFFFFFF));
}

void
l1t::MuonRawDigiTranslator::generatePackedDataWords(const Muon& mu, uint32_t &raw_data_00_31, uint32_t &raw_data_32_63)
{
  raw_data_00_31 = (mu.hwPt() & ptWidth_) << ptShift_
                 | (mu.hwQual() & qualWidth_) << qualShift_
                 | (abs(mu.hwEta()) & absEtaWidth_) << absEtaShift_
                 | ((mu.hwEta() < 0) & etaSignWidth_) << etaSignShift_
                 | (mu.hwPhi() & phiWidth_) << phiShift_;

  raw_data_32_63 = ((mu.hwCharge() > 0) & chargeWidth_) << chargeShift_
                 | (mu.hwChargeValid() & chargeValidWidth_) << chargeValidShift_
                 | (mu.hwIso() & isoWidth_) << isoShift_;
}

uint64_t 
l1t::MuonRawDigiTranslator::generate64bitDataWord(const Muon& mu)
{
  uint32_t lsw;
  uint32_t msw;

  generatePackedDataWords(mu, lsw, msw);
  return (((uint64_t)msw) << 32) + lsw;
}

