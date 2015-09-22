#include "CondFormats/L1TObjects/interface/MicroGMTParams.h"

using namespace l1t;

void MicroGMTParams::print(std::ostream& out) const {

  out << "L1 MicroGMT Parameters" << std::endl;

  out << "Firmware version: " << fwVersion_ << std::endl;

  out << "LUT configurations" << std::endl;
  out << " Barrel Single MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << bsinglemqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << bsinglemqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << bsinglemqp_.outWidth() << std::endl;
  out << "  Filename                     : " << bsinglemqp_.filename() << std::endl;

  out << " Forward neg MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << fnegmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << fnegmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << fnegmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << fnegmqp_.filename() << std::endl;

  out << " Forward pos MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << fposmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << fposmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << fposmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << fposmqp_.filename() << std::endl;

  out << " Overlap neg MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << onegmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << onegmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << onegmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << onegmqp_.filename() << std::endl;

  out << " Overlap pos MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << oposmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << oposmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << oposmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << oposmqp_.filename() << std::endl;

  out << " Barrel-Overlap neg MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << bonegmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << bonegmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << bonegmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << bonegmqp_.filename() << std::endl;

  out << " Barrel-Overlap pos MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << boposmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << boposmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << boposmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << boposmqp_.filename() << std::endl;

  out << " Forward-Overlap neg MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << fonegmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << fonegmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << fonegmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << fonegmqp_.filename() << std::endl;

  out << " Forward-Overlap pos MatchQual" << std::endl;
  out << "  Delta eta reduced input width: " << foposmqp_.deltaEtaRedInWidth() << std::endl;
  out << "  Delta phi reduced input width: " << foposmqp_.deltaPhiRedInWidth() << std::endl;
  out << "  Output width                 : " << foposmqp_.outWidth() << std::endl;
  out << "  Filename                     : " << foposmqp_.filename() << std::endl;

  out << " Barrel phi extrapolation" << std::endl;
  out << "  Abs eta reduced input width  : " << betaepp_.etaAbsRedInWidth() << std::endl;
  out << "  pT reduced input width       : " << betaepp_.ptRedInWidth() << std::endl;
  out << "  Output width                 : " << betaepp_.outWidth() << std::endl;
  out << "  Filename                     : " << betaepp_.filename() << std::endl;

  out << " Overlap eta extrapolation" << std::endl;
  out << "  Abs eta reduced input width  : " << oetaepp_.etaAbsRedInWidth() << std::endl;
  out << "  pT reduced input width       : " << oetaepp_.ptRedInWidth() << std::endl;
  out << "  Output width                 : " << oetaepp_.outWidth() << std::endl;
  out << "  Filename                     : " << oetaepp_.filename() << std::endl;

  out << " Forward eta extrapolation" << std::endl;
  out << "  Abs eta reduced input width  : " << fetaepp_.etaAbsRedInWidth() << std::endl;
  out << "  pT reduced input width       : " << fetaepp_.ptRedInWidth() << std::endl;
  out << "  Output width                 : " << fetaepp_.outWidth() << std::endl;
  out << "  Filename                     : " << fetaepp_.filename() << std::endl;

  out << " Barrel eta extrapolation" << std::endl;
  out << "  Abs eta reduced input width  : " << bphiepp_.etaAbsRedInWidth() << std::endl;
  out << "  pT reduced input width       : " << bphiepp_.ptRedInWidth() << std::endl;
  out << "  Output width                 : " << bphiepp_.outWidth() << std::endl;
  out << "  Filename                     : " << bphiepp_.filename() << std::endl;

  out << " Overlap phi extrapolation" << std::endl;
  out << "  Abs eta reduced input width  : " << ophiepp_.etaAbsRedInWidth() << std::endl;
  out << "  pT reduced input width       : " << ophiepp_.ptRedInWidth() << std::endl;
  out << "  Output width                 : " << ophiepp_.outWidth() << std::endl;
  out << "  Filename                     : " << ophiepp_.filename() << std::endl;

  out << " Forward phi extrapolation" << std::endl;
  out << "  Abs eta reduced input width  : " << fphiepp_.etaAbsRedInWidth() << std::endl;
  out << "  pT reduced input width       : " << fphiepp_.ptRedInWidth() << std::endl;
  out << "  Output width                 : " << fphiepp_.outWidth() << std::endl;
  out << "  Filename                     : " << fphiepp_.filename() << std::endl;

  out << " Abs isolation checkMem" << std::endl;
  out << "  Area sum input width         : " << isocmp_.areaSumInWidth() << std::endl;
  out << "  Output width                 : " << isocmp_.outWidth() << std::endl;
  out << "  Filename                     : " << isocmp_.filename() << std::endl;

  out << " Rel isolation checkMem" << std::endl;
  out << "  Area sum input width         : " << risocmp_.areaSumInWidth() << std::endl;
  out << "  pT input width               : " << risocmp_.ptInWidth() << std::endl;
  out << "  Output width                 : " << risocmp_.outWidth() << std::endl;
  out << "  Filename                     : " << risocmp_.filename() << std::endl;

  out << " Index selMem eta" << std::endl;
  out << "  Eta input width              : " << smetap_.etaInWidth() << std::endl;
  out << "  Output width                 : " << smetap_.outWidth() << std::endl;
  out << "  Filename                     : " << smetap_.filename() << std::endl;

  out << " Index selMem phi" << std::endl;
  out << "  Phi input width              : " << smphip_.phiInWidth() << std::endl;
  out << "  Output width                 : " << smphip_.outWidth() << std::endl;
  out << "  Filename                     : " << smphip_.filename() << std::endl;

  out << " Sort rank" << std::endl;
  out << "  pT input width               : " << srp_.ptInWidth() << std::endl;
  out << "  Quality input width          : " << srp_.qualInWidth() << std::endl;
  out << "  Output width                 : " << srp_.outWidth() << std::endl;
  out << "  Filename                     : " << srp_.filename() << std::endl;
}
