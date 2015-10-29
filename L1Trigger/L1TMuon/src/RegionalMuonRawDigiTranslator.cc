#include "L1Trigger/L1TMuon/interface/RegionalMuonRawDigiTranslator.h"

void
l1t::RegionalMuonRawDigiTranslator::fillRegionalMuonCand(RegionalMuonCand& mu, uint32_t raw_data_00_31, uint32_t raw_data_32_63, int proc, tftype tf)
{
  // translations as defined in DN-15-017
  mu.setHwPt((raw_data_00_31 >> ptShift_) & ptWidth_);
  mu.setHwQual((raw_data_00_31 >> qualShift_) & qualWidth_); 

  // eta is coded as two's complement
  int abs_eta = (raw_data_00_31 >> absEtaShift_) & absEtaWidth_;
  if ((raw_data_00_31 >> etaSignShift_) & 0x1) {
     mu.setHwEta(abs_eta - (1 << (etaSignShift_ - absEtaShift_)));
  } else {
     mu.setHwEta(abs_eta);
  }

  mu.setHwPhi((raw_data_00_31 >> phiShift_) & phiWidth_);
  // sign is coded as -1^signBit
  int signBit = (raw_data_32_63 >> signShift_) & 0x1;
  mu.setHwSign(1 - 2*signBit);
  mu.setHwSignValid((raw_data_32_63 >> signValidShift_) & 0x1);
  mu.setHwHF((raw_data_00_31 >> hfShift_) & hfWidth_);

  // set track address with subaddresses
  int rawTrackAddress = (raw_data_32_63 >> trackAddressShift_) & trackAddressWidth_;
  if (tf == bmtf) {
    int wheel = (rawTrackAddress >> bmtfTrAddrWheelShift_) & bmtfTrAddrWheelWidth_;
    int statAddr1 = (rawTrackAddress >> bmtfTrAddrStat1Shift_) & bmtfTrAddrStat1Width_;
    int statAddr2 = (rawTrackAddress >> bmtfTrAddrStat2Shift_) & bmtfTrAddrStat2Width_;
    int statAddr3 = (rawTrackAddress >> bmtfTrAddrStat3Shift_) & bmtfTrAddrStat3Width_;
    int statAddr4 = (rawTrackAddress >> bmtfTrAddrStat4Shift_) & bmtfTrAddrStat4Width_;
    mu.setTrackSubAddress(RegionalMuonCand::kWheel, wheel);
    mu.setTrackSubAddress(RegionalMuonCand::kStat1, statAddr1);
    mu.setTrackSubAddress(RegionalMuonCand::kStat2, statAddr2);
    mu.setTrackSubAddress(RegionalMuonCand::kStat3, statAddr3);
    mu.setTrackSubAddress(RegionalMuonCand::kStat4, statAddr4);
  } else if (tf == emtf_neg || tf == emtf_pos) {
    int me12 = (rawTrackAddress >> emtfTrAddrMe12Shift_) & emtfTrAddrMe12Width_;
    int me22 = (rawTrackAddress >> emtfTrAddrMe22Shift_) & emtfTrAddrMe22Width_;
    mu.setTrackSubAddress(RegionalMuonCand::kME12, me12);
    mu.setTrackSubAddress(RegionalMuonCand::kME22, me22);
  } else {
    std::map<int, int> trackAddr;
    trackAddr[0] = rawTrackAddress;
    mu.setTrackAddress(trackAddr);
  }

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
  int abs_eta = mu.hwEta();
  if (abs_eta < 0) {
    abs_eta += (1 << (etaSignShift_ - absEtaShift_));
  }
  raw_data_00_31 = (mu.hwPt() & ptWidth_) << ptShift_
                 | (mu.hwQual() & qualWidth_) << qualShift_
                 | (abs_eta & absEtaWidth_) << absEtaShift_
                 | (mu.hwEta() < 0) << etaSignShift_
                 | (mu.hwHF() & hfWidth_) << hfShift_
                 | (mu.hwPhi() & phiWidth_) << phiShift_;

  int tf = mu.trackFinderType();
  int rawTrkAddr = 0;
  if (tf == bmtf) {
    rawTrkAddr = (mu.trackSubAddress(RegionalMuonCand::kWheel) & bmtfTrAddrWheelWidth_) << bmtfTrAddrWheelShift_
               | (mu.trackSubAddress(RegionalMuonCand::kStat1) & bmtfTrAddrStat1Width_) << bmtfTrAddrStat1Shift_
               | (mu.trackSubAddress(RegionalMuonCand::kStat2) & bmtfTrAddrStat2Width_) << bmtfTrAddrStat2Shift_
               | (mu.trackSubAddress(RegionalMuonCand::kStat3) & bmtfTrAddrStat3Width_) << bmtfTrAddrStat3Shift_
               | (mu.trackSubAddress(RegionalMuonCand::kStat4) & bmtfTrAddrStat4Width_) << bmtfTrAddrStat4Shift_;
  } else if (tf == emtf_neg || tf == emtf_pos) {
    rawTrkAddr = (mu.trackSubAddress(RegionalMuonCand::kME12) & emtfTrAddrMe12Width_) << emtfTrAddrMe12Shift_
               | (mu.trackSubAddress(RegionalMuonCand::kME22) & emtfTrAddrMe22Width_) << emtfTrAddrMe22Shift_;
  } else {
    rawTrkAddr = mu.trackAddress().at(0);
  }

  raw_data_32_63 = (mu.hwSign() < 0) << signShift_
                 | mu.hwSignValid() << signValidShift_
                 | (rawTrkAddr & trackAddressWidth_) << trackAddressShift_;
}

uint64_t 
l1t::RegionalMuonRawDigiTranslator::generate64bitDataWord(const RegionalMuonCand& mu)
{
  uint32_t lsw;
  uint32_t msw;

  generatePackedDataWords(mu, lsw, msw);
  return (((uint64_t)msw) << 32) + lsw;
}

