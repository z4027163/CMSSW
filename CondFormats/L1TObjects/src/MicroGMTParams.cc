#include "CondFormats/L1TObjects/interface/MicroGMTParams.h"

using namespace l1t;

void MicroGMTParams::print(std::ostream& out) const {

  out << "L1 MicroGMT Parameters" << std::endl;

  out << "Firmware version: " << fwVersion_ << std::endl;

  out << "LUT paths" << std::endl;
  out << " Abs isolation checkMem LUT path: "        << aisocmlp_ << std::endl;
  out << " Rel isolation checkMem LUT path: "        << risocmlp_ << std::endl;
  out << " Index selMem phi LUT path: "              << ismphilp_ << std::endl;
  out << " Index selMem eta LUT path: "              << ismetalp_ << std::endl;
  out << " Barrel Single MatchQual LUT path: "       << bsinglemqlp_ << std::endl;
  out << " Forward pos MatchQual LUT path: "         << fposmqlp_ << std::endl;
  out << " Forward neg MatchQual LUT path: "         << fnegmqlp_ << std::endl;
  out << " Overlap pos MatchQual LUT path: "         << oposmqlp_ << std::endl;
  out << " Overlap neg MatchQual LUT path: "         << onegmqlp_ << std::endl;
  out << " Barrel-Overlap pos MatchQual LUT path: "  << boposmqlp_ << std::endl;
  out << " Barrel-Overlap neg MatchQual LUT path: "  << bonegmqlp_ << std::endl;
  out << " Forward-Overlap pos MatchQual LUT path: " << foposmqlp_ << std::endl;
  out << " Forward-Overlap neg MatchQual LUT path: " << fonegmqlp_ << std::endl;
  out << " Barrel phi extrapolation LUT path: "      << bphieplp_ << std::endl;
  out << " Overlap phi extrapolation LUT path: "     << ophieplp_ << std::endl;
  out << " Forward phi extrapolation LUT path: "     << fphieplp_ << std::endl;
  out << " Barrel eta extrapolation LUT path: "      << betaeplp_ << std::endl;
  out << " Overlap eta extrapolation LUT path: "     << oetaeplp_ << std::endl;
  out << " Forward eta extrapolation LUT path: "     << fetaeplp_ << std::endl;
  out << " Sort rank LUT path: "                     << srlp_ << std::endl;
}
