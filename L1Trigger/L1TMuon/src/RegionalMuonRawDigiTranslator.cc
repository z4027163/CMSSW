#include "L1Trigger/L1TMuon/interface/RegionalMuonRawDigiTranslator.h"

void
l1t::RegionalMuonRawDigiTranslator::fillRegionalMuonCand(RegionalMuonCand& mu, uint32_t raw_data_00_31, uint32_t raw_data_32_63, int proc, tftype tf)
{
  // translations as defined in DN-15-017
  mu.setHwPt((raw_data_00_31 >> 0) & 0x1FF);
  mu.setHwQual((raw_data_00_31 >> 9) & 0xF); 

  // eta is coded as two's complement
  int abs_eta = (raw_data_00_31 >> 13) & 0xFF;
  if ((raw_data_00_31 >> 21) & 0x1) {
     mu.setHwEta(abs_eta - 256);
  } else {
     mu.setHwEta(abs_eta);
  }

  mu.setHwPhi((raw_data_00_31 >> 23) & 0xFF);
  // sign is coded as -1^signBit
  int signBit = (raw_data_32_63 >> 0) & 0x1;
  mu.setHwSign(1 - 2*signBit);
  mu.setHwSignValid((raw_data_32_63 >> 1) & 0x1);
  mu.setHwHF((raw_data_00_31 >> 22) & 0x1);
  mu.setHwTrackAddress((raw_data_32_63 >> 4) & 0x1FFF);
  mu.setTFIdentifiers(proc, tf);
  mu.setDataword(raw_data_32_63, raw_data_00_31);
}

