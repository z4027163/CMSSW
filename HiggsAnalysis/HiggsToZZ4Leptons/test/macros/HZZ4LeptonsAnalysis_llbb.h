//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  7 12:14:04 2012 by ROOT version 5.32/00
// from TTree HZZ4LeptonsAnalysis/HZZ4Leptons Analysis Tree
// found on file: roottree_leptons_Fall11_0706.root
//////////////////////////////////////////////////////////

#ifndef HZZ4LeptonsAnalysis_h
#define HZZ4LeptonsAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>

#define nn 50
#define nnn 200
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class HZZ4LeptonsAnalysis {
public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          Run;
   UInt_t          Event;
   UInt_t          LumiSection;
   Float_t         Avginstlumi;
   Int_t           num_PU_vertices;
   Int_t           PU_BunchCrossing;
   Float_t         MC_weighting;
   Float_t         MC_weighting_un[9];
   Float_t         PDF_weighting_un;
   Int_t           RECO_nMuHLTMatch;
   Float_t         RECOMU_PT_MuHLTMatch[nn];
   Float_t         RECOMU_ETA_MuHLTMatch[nn];
   Bool_t          RECOMU_sm_MuHLTMatch[nn];
   Int_t           RECOMU_dm_MuHLTMatch[nn]; 
   Float_t         RECOELE_PT_EleHLTMatch[nn];
   Float_t         RECOELE_ETA_EleHLTMatch[nn];
   Bool_t          RECOELE_se_EleHLTMatch[nn];
   Int_t           RECOELE_de_EleHLTMatch[nn];
   Int_t           RECOBOT_MatchingMCTruth[nn];
   Bool_t          dm_trig;
   Bool_t          sm_trig;
   Bool_t          de_trig;
   Bool_t          se_trig;
   Bool_t          tri_trig;
/*
   Float_t         MC_E[7];
   Float_t         MC_PT[7];
   Float_t         MC_ETA[7];
   Float_t         MC_THETA[7];
   Float_t         MC_PHI[7];
   Float_t         MC_MASS[7];
   Float_t         MC_PDGID[7];
   Float_t         MC_LEPT_PT[4];
   Float_t         MC_LEPT_ETA[4];
   Float_t         MC_LEPT_PHI[4];
   Float_t         MC_LEPT_THETA[4];
   Float_t         MC_LEPT_PDGID[4];
*/
   Float_t         MC_Z_PT[2][5];
   Float_t         MC_Z_ETA[2][5];
   Float_t         MC_Z_PHI[2][5];
   Float_t         MC_Z_THETA[2][5];
   Float_t         MC_Z_MASS[2][5];
   Float_t         MC_Z_PDGID[2][5];

   Float_t         MC_GENJET_PT[nn];
   Float_t         MC_GENJET_ETA[nn];
   Float_t         MC_GENJET_PHI[nn];
   Float_t         MC_GENMET;
   Float_t         RECOELE_E[nn];
   Float_t         RECOELE_PT[nn];
   Float_t         RECOELE_PTError[nn];
   Float_t         RECOELE_PT_uncorr[nn];
   Float_t         RECOELE_P[nn];
   Float_t         RECOELE_ETA[nn];
   Float_t         RECOELE_THETA[nn];
   Float_t         RECOELE_PHI[nn];
   Float_t         RECOELE_MASS[nn];
   Float_t         RECOELE_CHARGE[nn];
   UChar_t         RECOELE_isEcalDriven[nn];
   UChar_t         RECOELE_isTrackerDriven[nn];
/*
   Float_t         RECOELE_gsftrack_NPixHits[nn];
   Float_t         RECOELE_gsftrack_NStripHits[nn];
   Float_t         RECOELE_gsftrack_chi2[nn];
   Float_t         RECOELE_gsftrack_dxyB[nn];
*/
   Float_t         RECOELE_gsftrack_dxy[nn];
   Float_t         RECOELE_gsftrack_dxyError[nn];
   Float_t         RECOELE_gsftrack_dzB[nn];
   Float_t         RECOELE_gsftrack_dz[nn];
   Float_t         RECOELE_gsftrack_dzError[nn];
   Float_t         RECOELE_scl_E[nn];
   Float_t         RECOELE_scl_Et[nn];
   Float_t         RECOELE_scl_Eta[nn];
   Float_t         RECOELE_scl_Phi[nn];
/*
   Float_t         RECOELE_sigmaIetaIeta[nn];
   Float_t         RECOELE_sigmaEtaEta[nn];
   Float_t         RECOELE_he[nn];
   Float_t         RECOELE_r9[nn];
   Float_t         RECOELE_mva[nn];
   Float_t         RECOELE_fbrem[nn];
   Int_t           RECOELE_nbrems[nn];
   Int_t           RECOELE_Class[nn];
   Double_t        RECOELE_fbrem_mode[nn];
   Double_t        RECOELE_fbrem_mean[nn];
   Float_t         RECOELE_EGMTRACKISO[nn];
   Float_t         RECOELE_EGMHCALISO[nn];
   Float_t         RECOELE_EGMECALISO[nn];
   Float_t         RECOELE_EGMX[nn];
*/
   Double_t        RECOELE_PFchAllPart[nn];
   Double_t        RECOELE_PFchHad[nn];
   Double_t        RECOELE_PFneuHad[nn];
   Double_t        RECOELE_PFphoton[nn];
   Double_t        RECOELE_PFPUchAllPart[nn];
   Double_t        RECOELE_PFX_dB[nn];
   Double_t        RECOELE_PFX_rho[nn];
   Double_t        RECOELE_PFX_rho_new[nn];
   Double_t        RECOELE_regEnergy[nn];
   Double_t        RECOELE_regEnergyError[nn];
   Float_t         RECOELE_SIP[nn];
   Double_t        RECOELE_sclRawE[nn];
   Double_t        RECOELE_sclX[nn];
   Double_t        RECOELE_sclY[nn];
   Double_t        RECOELE_sclZ[nn];
   Double_t        RECOELE_mvaTrigV0[nn];
   Double_t        RECOELE_mvaNonTrigV0[nn];
   Double_t        RECOELE_COV[nn][3][3];
   Float_t         RECOELE_energyScaleUp[nn];
   Float_t         RECOELE_energyScaleDown[nn];
   Float_t         RECOELE_energySigmaUp[nn];
   Float_t         RECOELE_energySigmaDown[nn];
   UChar_t         RECOMU_isPFMu[nn];
   UChar_t         RECOMU_isMedium[nn];
   Float_t         RECOMU_muInnertrkvalidFraction[nn];
   UChar_t         RECOMU_isGlobalMu[nn];
   UChar_t         RECOMU_isStandAloneMu[nn];
   UChar_t         RECOMU_isTrackerMu[nn];
   UChar_t         RECOMU_isCaloMu[nn];
   UChar_t         RECOMU_isTrackerHighPtMu[nn];
   Float_t         RECOMU_E[nn];
   Float_t         RECOMU_PT[nn];
   Float_t         RECOMU_P[nn];
   Float_t         RECOMU_ETA[nn];
   Float_t         RECOMU_THETA[nn];
   Float_t         RECOMU_PHI[nn];
   Float_t         RECOMU_MASS[nn];
   Float_t         RECOMU_CHARGE[nn];
   Double_t        RECOMU_COV[nn][3][3];
   Float_t         RECOMU_TRACKISO[nn];
   Float_t         RECOMU_TRACKISO_SUMPT[nn];
   Float_t         RECOMU_HCALISO[nn];
   Float_t         RECOMU_ECALISO[nn];
   Float_t         RECOMU_X[nn];
   Double_t        RECOMU_PFchHad[nn];
   Double_t        RECOMU_PFneuHad[nn];
   Double_t        RECOMU_PFphoton[nn];
   Double_t        RECOMU_PFPUchAllPart[nn];
   Double_t        RECOMU_PFX_dB[nn];
   Double_t        RECOMU_PFX_dB_new[nn];
   Double_t        RECOMU_PFX_rho[nn];
   Double_t        RECOPFPHOT_PFchHad[20];
   Double_t        RECOPFPHOT_PFneuHad[20];
   Double_t        RECOPFPHOT_PFphoton[20];
   Double_t        RECOPFPHOT_PFPUchAllPart[20];
   Double_t        RECOPFPHOT_PFX_rho[20];
   Float_t         RECOMU_mutrkPT[nn];
   Float_t         RECOMU_mutrkPTError[nn];
   Float_t         RECOMU_mutrkDxy[nn];
   Float_t         RECOMU_mutrkDxyError[nn];
   Float_t         RECOMU_mutrkDxyB[nn];
   Float_t         RECOMU_mutrkDz[nn];
   Float_t         RECOMU_mutrkDzError[nn];
   Float_t         RECOMU_mutrkDzB[nn];
   Float_t         RECOMU_mutrkChi2PerNdof[nn];
   Float_t         RECOMU_mutrkCharge[nn];
   Float_t         RECOMU_mutrkNHits[nn];
   Float_t         RECOMU_mutrkNStripHits[nn];
   Float_t         RECOMU_mutrkNPixHits[nn];
   Float_t         RECOMU_mutrkNMuonHits[nn];
   Float_t         RECOMU_mutrktrackerLayersWithMeasurement[nn];
   Int_t           RECOMU_mubesttrkType[nn];
   Float_t         RECOMU_mubesttrkDxy[nn];
   Float_t         RECOMU_mubesttrkDxyError[nn];
   Float_t         RECOMU_mubesttrkDz[nn];
   Float_t         RECOMU_mubesttrkDzError[nn];
   UChar_t         RECOMU_MatchingMCTruth[nn];
   Float_t         RECOMU_MatchingMCpT[nn];
   Float_t         RECOMU_MatchingMCEta[nn];
   Float_t         RECOMU_MatchingMCPhi[nn];
   UChar_t         RECOELE_MatchingMCTruth[nn];
   Float_t         RECOELE_MatchingMCpT[nn];
   Float_t         RECOELE_MatchingMCEta[nn];
   Float_t         RECOELE_MatchingMCPhi[nn];
   UChar_t         RECOPHOT_MatchingMCTruth[nn];
   Float_t         RECOPHOT_MatchingMCpT[nn];
   Float_t         RECOPHOT_MatchingMCEta[nn];
   Float_t         RECOPHOT_MatchingMCPhi[nn];
   Int_t           RECO_NMU;
   Int_t           RECO_NELE;
   Int_t           RECO_NTRACK;
   Int_t           RECO_NPHOT;
   Float_t         RECOPHOT_PT[20];
   Float_t         RECOPHOT_ETA[20];
   Float_t         RECOPHOT_PHI[20];
   Float_t         RECOPHOT_THETA[20];
   Int_t           RECO_NPFPHOT;
   Float_t         RECOPFPHOT_PT[20];
   Float_t         RECOPFPHOT_PTError[20];
   Float_t         RECOPFPHOT_ETA[20];
   Float_t         RECOPFPHOT_PHI[20];
   Float_t         RECOPFPHOT_THETA[20];
   Double_t        BeamSpot_X;
   Double_t        BeamSpot_Y;
   Double_t        BeamSpot_Z;
   Int_t           RECO_NVTX;
   Int_t           RECO_PFJET_N;
   Int_t           RECO_PFJET_CHARGE[nn];
   Float_t         RECO_PFJET_ET[nn];
   Float_t         RECO_PFJET_PT[nn];
   Float_t         RECO_PFJET_PT_UP[nn];
   Float_t         RECO_PFJET_PT_DOW[nn];
   Float_t         RECO_PFJET_ETA[nn];
   Float_t         RECO_PFJET_PHI[nn];
   Int_t           RECO_PFJET_PUID[nn];
   Float_t         RECO_PFJET_PUID_MVA[nn];
   Int_t           RECO_PFJET_nconstituents[nn];
   Int_t           RECO_PFJET_NCH[nn];
   Float_t         RECO_PFJET_NHF[nn];
   Float_t         RECO_PFJET_NEF[nn];
   Float_t         RECO_PFJET_CHF[nn];
   Float_t         RECO_PFJET_CEF[nn];
   Float_t         RECO_PFJET_MUF[nn];
   Double_t        RHO_ele;
   Double_t        RHO_mu;
   Float_t         RECO_CALOMET;
   Float_t         RECO_PFMET;
   Float_t         RECO_PFMET_X;
   Float_t         RECO_PFMET_Y;
   Float_t         RECO_PFMET_PHI;
   Float_t         RECO_PFMET_THETA;
   Float_t         RECO_TCMET;
   Float_t         RECO_CORMETMUONS;
   Float_t         tCHighEff_BTagJet_PT[nn];
   Float_t         tCHighPur_BTagJet_PT[nn];
   Float_t         cSV_BTagJet_PT[nn];
   Float_t         tCHighEff_BTagJet_ETA[nn];
   Float_t         tCHighPur_BTagJet_ETA[nn];
   Float_t         cSV_BTagJet_ETA[nn];
   Float_t         tCHighEff_BTagJet_PHI[nn];
   Float_t         tCHighPur_BTagJet_PHI[nn];
   Float_t         cSV_BTagJet_PHI[nn];
   Float_t         tCHighEff_BTagJet_DISCR[nn];
   Float_t         tCHighPur_BTagJet_DISCR[nn];
   Float_t         cSV_BTagJet_DISCR[nn];
   Float_t         cSV_BTagJet_ET[nn];

   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_ils;   //!
   TBranch        *b_Avginstlumi;   //!
   TBranch        *b_num_PU_vertices;   //!
   TBranch        *b_PU_BunchCrossing;   //!
   TBranch        *b_MC_weighting;   //!
   TBranch        *b_MC_weighting_un;   //!
   TBranch        *b_PDF_weighting_un;  //!
   TBranch        *b_RECO_nMuHLTMatch;   //!
   TBranch        *b_RECOMU_PT_MuHLTMatch;   //!
   TBranch        *b_RECOMU_ETA_MuHLTMatch;  //!
   TBranch        *b_RECOMU_sm_MuHLTMatch;   //!
   TBranch        *b_RECOMU_dm_MuHLTMatch;   //!
   TBranch        *b_RECOELE_PT_EleHLTMatch;   //!
   TBranch        *b_RECOELE_ETA_EleHLTMatch;  //!
   TBranch        *b_RECOELE_se_EleHLTMatch;   //!
   TBranch        *b_RECOELE_de_EleHLTMatch;   //!
   TBranch        *b_RECOBOT_MatchingMCTruth;   //!
   TBranch        *b_dm_trig;  //!
   TBranch        *b_sm_trig;  //!
   TBranch        *b_de_trig;  //!
   TBranch        *b_se_trig;   //!
   TBranch        *b_tri_trig; //!
/*   TBranch        *b_MC_E;   //!
   TBranch        *b_MC_PT;   //!
   TBranch        *b_MC_ETA;   //!
   TBranch        *b_MC_THETA;   //!
   TBranch        *b_MC_PHI;   //!
   TBranch        *b_MC_MASS;   //!
   TBranch        *b_MC_PDGID;   //!
   TBranch        *b_MC_LEPT_PT;   //!
   TBranch        *b_MC_LEPT_ETA;   //!
   TBranch        *b_MC_LEPT_PHI;   //!
   TBranch        *b_MC_LEPT_THETA;   //!
   TBranch        *b_MC_LEPT_PDGID;   //!
*/
   TBranch        *b_MC_Z_PT;   //!
   TBranch        *b_MC_Z_ETA;   //!
   TBranch        *b_MC_Z_PHI;   //!
   TBranch        *b_MC_Z_THETA;   //!
   TBranch        *b_MC_Z_MASS;   //!
   TBranch        *b_MC_Z_PDGID;   //!

   TBranch        *b_MC_GENJET_PT;  //!
   TBranch        *b_MC_GENJET_ETA; //!
   TBranch        *b_MC_GENJET_PHI;  //!
   TBranch        *b_MC_GENMET;   //!
   TBranch        *b_RECOELE_E;   //!
   TBranch        *b_RECOELE_energyScaleUp;
   TBranch        *b_RECOELE_energyScaleDown;
   TBranch        *b_RECOELE_energySigmaUp; 
   TBranch        *b_RECOELE_energySigmaDown;
   TBranch        *b_RECOELE_PT;   //!
   TBranch        *b_RECOELE_PTError;   //!
   TBranch        *b_RECOELE_PT_uncorr; //!
   TBranch        *b_RECOELE_P;   //!
   TBranch        *b_RECOELE_ETA;   //!
   TBranch        *b_RECOELE_THETA;   //!
   TBranch        *b_RECOELE_PHI;   //!
   TBranch        *b_RECOELE_MASS;   //!
   TBranch        *b_RECOELE_CHARGE;   //!
   TBranch        *b_RECOELE_isEcalDriven;   //!
   TBranch        *b_RECOELE_isTrackerDriven;   //!
/*
   TBranch        *b_RECOELE_gsftrack_NPixHits;   //!
   TBranch        *b_RECOELE_gsftrack_NStripHits;   //!
   TBranch        *b_RECOELE_gsftrack_chi2;   //!
   TBranch        *b_RECOELE_gsftrack_dxyB;   //!
*/
   TBranch        *b_RECOELE_gsftrack_dxy;   //!
   TBranch        *b_RECOELE_gsftrack_dxyError;   //!
   TBranch        *b_RECOELE_gsftrack_dzB;   //!
   TBranch        *b_RECOELE_gsftrack_dz;   //!
   TBranch        *b_RECOELE_gsftrack_dzError;   //!
   TBranch        *b_RECOELE_scl_E;   //!
   TBranch        *b_RECOELE_scl_Et;   //!
   TBranch        *b_RECOELE_scl_Eta;   //!
   TBranch        *b_RECOELE_scl_Phi;   //!
/*
   TBranch        *b_RECOELE_sigmaIetaIeta;   //!
   TBranch        *b_RECOELE_sigmaEtaEta;   //!
   TBranch        *b_RECOELE_he;   //!
   TBranch        *b_RECOELE_r9;   //!
   TBranch        *b_RECOELE_mva;   //!
   TBranch        *b_RECOELE_fbrem;   //!
   TBranch        *b_RECOELE_nbrems;   //!
   TBranch        *b_RECOELE_Class;   //!
   TBranch        *b_RECOELE_fbrem_mode;   //!
   TBranch        *b_RECOELE_fbrem_mean;   //!
   TBranch        *b_RECOELE_EGMTRACKISO;   //!
   TBranch        *b_RECOELE_EGMHCALISO;   //!
   TBranch        *b_RECOELE_EGMECALISO;   //!
   TBranch        *b_RECOELE_EGMX;   //!
*/
   TBranch        *b_RECOELE_PFchAllPart;   //!
   TBranch        *b_RECOELE_PFchHad;   //!
   TBranch        *b_RECOELE_PFneuHad;   //!
   TBranch        *b_RECOELE_PFphoton;   //!
   TBranch        *b_RECOELE_PFPUchAllPart;   //!
   TBranch        *b_RECOELE_PFX_dB;   //!
   TBranch        *b_RECOELE_PFX_rho;   //!
   TBranch        *b_RECOELE_regEnergy;   //!
   TBranch        *b_RECOELE_regEnergyError;   //!
   TBranch        *b_RECOELE_SIP;   //!
   TBranch        *b_RECOELE_sclRawE;   //!
   TBranch        *b_RECOELE_sclX;   //!
   TBranch        *b_RECOELE_sclY;   //!
   TBranch        *b_RECOELE_sclZ;   //!
   TBranch        *b_RECOELE_mvaTrigV0;   //!
   TBranch        *b_RECOELE_mvaNonTrigV0;   //!
   TBranch        *b_RECOELE_COV;   //!
   TBranch        *b_RECOMU_isPFMu;   //!
   TBranch        *b_RECOMU_isMedium; //!
   TBranch        *b_RECOMU_muInnertrkvalidFraction; //!
   TBranch        *b_RECOMU_isGlobalMu;   //!
   TBranch        *b_RECOMU_isStandAloneMu;   //!
   TBranch        *b_RECOMU_isTrackerMu;   //!
   TBranch        *b_RECOMU_isCaloMu;   //!
   TBranch        *b_RECOMU_isTrackerHighPtMu;   //!
   TBranch        *b_RECOMU_E;   //!
   TBranch        *b_RECOMU_PT;   //!
   TBranch        *b_RECOMU_P;   //!
   TBranch        *b_RECOMU_ETA;   //!
   TBranch        *b_RECOMU_THETA;   //!
   TBranch        *b_RECOMU_PHI;   //!
   TBranch        *b_RECOMU_MASS;   //!
   TBranch        *b_RECOMU_CHARGE;   //!
   TBranch        *b_RECOMU_COV;   //!
   TBranch        *b_RECOMU_TRACKISO;   //!
   TBranch        *b_RECOMU_TRACKISO_SUMPT;   //!
   TBranch        *b_RECOMU_HCALISO;   //!
   TBranch        *b_RECOMU_ECALISO;   //!
   TBranch        *b_RECOMU_X;   //!
   TBranch        *b_RECOMU_PFchHad;   //!
   TBranch        *b_RECOMU_PFneuHad;   //!
   TBranch        *b_RECOMU_PFphoton;   //!
   TBranch        *b_RECOMU_PFPUchAllPart;   //!
   TBranch        *b_RECOMU_PFX_dB;   //!
   TBranch        *b_RECOMU_PFX_rho;   //!
   TBranch        *b_RECOPFPHOT_PFchHad;   //!
   TBranch        *b_RECOPFPHOT_PFneuHad;   //!
   TBranch        *b_RECOPFPHOT_PFphoton;   //!
   TBranch        *b_RECOPFPHOT_PFPUchAllPart;   //!
   TBranch        *b_RECOPFPHOT_PFX_rho;   //!
   TBranch        *b_RECOMU_mutrkPT;   //!
   TBranch        *b_RECOMU_mutrkPTError;   //!
   TBranch        *b_RECOMU_mutrkDxy;   //!
   TBranch        *b_RECOMU_mutrkDxyError;   //!
   TBranch        *b_RECOMU_mutrkDxyB;   //!
   TBranch        *b_RECOMU_mutrkDz;   //!
   TBranch        *b_RECOMU_mutrkDzError;   //!
   TBranch        *b_RECOMU_mutrkDzB;   //!
   TBranch        *b_RECOMU_mutrkChi2PerNdof;   //!
   TBranch        *b_RECOMU_mutrkCharge;   //!
   TBranch        *b_RECOMU_mutrkNHits;   //!
   TBranch        *b_RECOMU_mutrkNStripHits;   //!
   TBranch        *b_RECOMU_mutrkNPixHits;   //!
   TBranch        *b_RECOMU_mutrkNMuonHits;   //!
   TBranch        *b_RECOMU_mutrktrackerLayersWithMeasurement;   //!
   TBranch        *b_RECOMU_mubesttrkType;   //!
   TBranch        *b_RECOMU_mubesttrkDxy;   //!
   TBranch        *b_RECOMU_mubesttrkDxyError;   //!
   TBranch        *b_RECOMU_mubesttrkDz;   //!
   TBranch        *b_RECOMU_mubesttrkDzError;   //!
   TBranch        *b_RECOMU_MatchingMCTruth;   //!
   TBranch        *b_RECOMU_MatchingMCpT;   //!
   TBranch        *b_RECOMU_MatchingMCEta;   //!
   TBranch        *b_RECOMU_MatchingMCPhi;   //!
   TBranch        *b_RECOELE_MatchingMCTruth;   //!
   TBranch        *b_RECOELE_MatchingMCpT;   //!
   TBranch        *b_RECOELE_MatchingMCEta;   //!
   TBranch        *b_RECOELE_MatchingMCPhi;   //!
   TBranch        *b_RECOPHOT_MatchingMCTruth;   //!
   TBranch        *b_RECOPHOT_MatchingMCpT;   //!
   TBranch        *b_RECOPHOT_MatchingMCEta;   //!
   TBranch        *b_RECOPHOT_MatchingMCPhi;   //!
   TBranch        *b_RECO_NMU;   //!
   TBranch        *b_RECO_NELE;   //!
   TBranch        *b_RECO_NTRACK;   //!
   TBranch        *b_RECO_NPHOT;   //!
   TBranch        *b_RECOPHOT_PT;   //!
   TBranch        *b_RECOPHOT_ETA;   //!
   TBranch        *b_RECOPHOT_PHI;   //!
   TBranch        *b_RECOPHOT_THETA;   //!
   TBranch        *b_RECO_NPFPHOT;   //!
   TBranch        *b_RECOPFPHOT_PT;   //!
   TBranch        *b_RECOPFPHOT_PTError;   //!
   TBranch        *b_RECOPFPHOT_ETA;   //!
   TBranch        *b_RECOPFPHOT_PHI;   //!
   TBranch        *b_RECOPFPHOT_THETA;   //!
   TBranch        *b_BeamSpot_X;   //!
   TBranch        *b_BeamSpot_Y;   //!
   TBranch        *b_BeamSpot_Z;   //!
   TBranch        *b_RECO_NVTX;   //!
   TBranch        *b_RECO_PFJET_N;   //!
   TBranch        *b_RECO_PFJET_CHARGE;   //!
   TBranch        *b_RECO_PFJET_ET;   //!
   TBranch        *b_RECO_PFJET_PT;   //!
   TBranch        *b_RECO_PFJET_PT_UP;   //!
   TBranch        *b_RECO_PFJET_PT_DOW;   //!
   TBranch        *b_RECO_PFJET_ETA;   //!
   TBranch        *b_RECO_PFJET_PHI;   //!
   TBranch        *b_RECO_PFJET_PUID;   //!
   TBranch        *b_RECO_PFJET_PUID_MVA; //!
   TBranch        *b_RECO_PFJET_NCH;   //!
   TBranch        *b_RECO_PFJET_nconstituents;   //!
   TBranch        *b_RECO_PFJET_NHF;   //!
   TBranch        *b_RECO_PFJET_NEF;   //!
   TBranch        *b_RECO_PFJET_CHF;   //!
   TBranch        *b_RECO_PFJET_CEF;   //!
   TBranch        *b_RECO_PFJET_MUF;
   TBranch        *b_RHO_ele;   //!
   TBranch        *b_RHO_mu;   //!
   TBranch        *b_RECO_CALOMET;   //!
   TBranch        *b_RECO_PFMET;   //!
   TBranch        *b_RECO_PFMET_X;   //!
   TBranch        *b_RECO_PFMET_Y;   //!
   TBranch        *b_RECO_PFMET_PHI;   //!
   TBranch        *b_RECO_PFMET_THETA;   //!
   TBranch        *b_RECO_TCMET;   //!
   TBranch        *b_RECO_CORMETMUONS;   //!
   TBranch        *b_tCHighEff_BTagJet_PT;   //!
   TBranch        *b_tCHighPur_BTagJet_PT;   //!
   TBranch        *b_cSV_BTagJet_PT;   //!
   TBranch        *b_tCHighEff_BTagJet_ETA;   //!
   TBranch        *b_tCHighPur_BTagJet_ETA;   //!
   TBranch        *b_cSV_BTagJet_ETA;   //!
   TBranch        *b_tCHighEff_BTagJet_PHI;   //!
   TBranch        *b_tCHighPur_BTagJet_PHI;   //!
   TBranch        *b_cSV_BTagJet_PHI;   //!
   TBranch        *b_tCHighEff_BTagJet_DISCR;   //!
   TBranch        *b_tCHighPur_BTagJet_DISCR;   //!
   TBranch        *b_cSV_BTagJet_DISCR;   //!
   TBranch        *b_cSV_BTagJet_ET;   //!


   HZZ4LeptonsAnalysis(TTree *tree=0,Double_t weight_=1.,string DATA_type_="DATA",string MC_type_="MC");
   virtual ~HZZ4LeptonsAnalysis();
   Double_t weight;
   string DATA_type,MC_type;
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Char_t *name);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   ofstream bnn_file;
   double EAele(int ,bool );
   double masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr );
   void printelebnn(int i);
   void printmubnn(int i);
   float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState);
   float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState);
   float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState);
   float kfactor_ggZZ(float GENmassZZ, int finalState);
};

#endif

#ifdef HZZ4LeptonsAnalysis_cxx
HZZ4LeptonsAnalysis::HZZ4LeptonsAnalysis(TTree *tree,Double_t weight_, string DATA_type_, string MC_type_) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   weight = weight_;
   DATA_type = DATA_type_;
   MC_type = MC_type_;

   if (tree == 0) {
      TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
      chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root_a");
      chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_1.root_b");

      tree = chain;

   }
   Init(tree);
}

HZZ4LeptonsAnalysis::~HZZ4LeptonsAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HZZ4LeptonsAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HZZ4LeptonsAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HZZ4LeptonsAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);


   fChain->SetBranchAddress("Run", &Run, &b_irun);
   fChain->SetBranchAddress("Event", &Event, &b_ievt);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_ils);
   fChain->SetBranchAddress("Avginstlumi", &Avginstlumi, &b_Avginstlumi);
   fChain->SetBranchAddress("num_PU_vertices", &num_PU_vertices, &b_num_PU_vertices);
   fChain->SetBranchAddress("PU_BunchCrossing", &PU_BunchCrossing, &b_PU_BunchCrossing);
   fChain->SetBranchAddress("MC_weighting", &MC_weighting, &b_MC_weighting);
   fChain->SetBranchAddress("MC_weighting_un", &MC_weighting_un, &b_MC_weighting_un);
   fChain->SetBranchAddress("PDF_weighting_un", &PDF_weighting_un, &b_PDF_weighting_un);
   fChain->SetBranchAddress("RECO_nMuHLTMatch", &RECO_nMuHLTMatch, &b_RECO_nMuHLTMatch);
   fChain->SetBranchAddress("RECOMU_PT_MuHLTMatch", RECOMU_PT_MuHLTMatch, &b_RECOMU_PT_MuHLTMatch);
   fChain->SetBranchAddress("RECOMU_ETA_MuHLTMatch", RECOMU_ETA_MuHLTMatch, &b_RECOMU_ETA_MuHLTMatch);
   fChain->SetBranchAddress("RECOMU_dm_MuHLTMatch", RECOMU_dm_MuHLTMatch, &b_RECOMU_dm_MuHLTMatch);
   fChain->SetBranchAddress("RECOMU_sm_MuHLTMatch", RECOMU_sm_MuHLTMatch, &b_RECOMU_sm_MuHLTMatch);
   fChain->SetBranchAddress("RECOELE_PT_EleHLTMatch", RECOELE_PT_EleHLTMatch, &b_RECOELE_PT_EleHLTMatch);
   fChain->SetBranchAddress("RECOELE_ETA_EleHLTMatch", RECOELE_ETA_EleHLTMatch, &b_RECOELE_ETA_EleHLTMatch);
   fChain->SetBranchAddress("RECOELE_de_MuHLTMatch", RECOELE_de_EleHLTMatch, &b_RECOELE_de_EleHLTMatch);
   fChain->SetBranchAddress("RECOELE_se_MuHLTMatch", RECOELE_se_EleHLTMatch, &b_RECOELE_se_EleHLTMatch);
   fChain->SetBranchAddress("RECOBOT_MatchingMCTruth",RECOBOT_MatchingMCTruth,&b_RECOBOT_MatchingMCTruth);
   fChain->SetBranchAddress("dm_trig",&dm_trig,&b_dm_trig);
   fChain->SetBranchAddress("sm_trig",&sm_trig,&b_sm_trig);
   fChain->SetBranchAddress("de_trig",&de_trig,&b_de_trig);
   fChain->SetBranchAddress("se_trig",&se_trig,&b_se_trig);
   fChain->SetBranchAddress("tri_trig",&tri_trig,&b_tri_trig);
/*
   fChain->SetBranchAddress("MC_E", MC_E, &b_MC_E);
   fChain->SetBranchAddress("MC_PT", MC_PT, &b_MC_PT);
   fChain->SetBranchAddress("MC_ETA", MC_ETA, &b_MC_ETA);
   fChain->SetBranchAddress("MC_THETA", MC_THETA, &b_MC_THETA);
   fChain->SetBranchAddress("MC_PHI", MC_PHI, &b_MC_PHI);
   fChain->SetBranchAddress("MC_MASS", MC_MASS, &b_MC_MASS);
   fChain->SetBranchAddress("MC_PDGID", MC_PDGID, &b_MC_PDGID);
   fChain->SetBranchAddress("MC_LEPT_PT", MC_LEPT_PT, &b_MC_LEPT_PT);
   fChain->SetBranchAddress("MC_LEPT_ETA", MC_LEPT_ETA, &b_MC_LEPT_ETA);
   fChain->SetBranchAddress("MC_LEPT_PHI", MC_LEPT_PHI, &b_MC_LEPT_PHI);
   fChain->SetBranchAddress("MC_LEPT_THETA", MC_LEPT_THETA, &b_MC_LEPT_THETA);
   fChain->SetBranchAddress("MC_LEPT_PDGID", MC_LEPT_PDGID, &b_MC_LEPT_PDGID);
*/
   fChain->SetBranchAddress("MC_Z_PT", MC_Z_PT, &b_MC_Z_PT);
   fChain->SetBranchAddress("MC_Z_ETA", MC_Z_ETA, &b_MC_Z_ETA);
   fChain->SetBranchAddress("MC_Z_PHI", MC_Z_PHI, &b_MC_Z_PHI);
   fChain->SetBranchAddress("MC_Z_THETA", MC_Z_THETA, &b_MC_Z_THETA);
   fChain->SetBranchAddress("MC_Z_MASS", MC_Z_MASS, &b_MC_Z_MASS);
   fChain->SetBranchAddress("MC_Z_PDGID", MC_Z_PDGID, &b_MC_Z_PDGID);

   fChain->SetBranchAddress("MC_GENJET_PT", MC_GENJET_PT, &b_MC_GENJET_PT);
   fChain->SetBranchAddress("MC_GENJET_ETA", MC_GENJET_ETA, &b_MC_GENJET_ETA);
   fChain->SetBranchAddress("MC_GENJET_PHI", MC_GENJET_PHI, &b_MC_GENJET_PHI);
   fChain->SetBranchAddress("MC_GENMET", &MC_GENMET, &b_MC_GENMET);
   fChain->SetBranchAddress("RECOELE_E", RECOELE_E, &b_RECOELE_E);
   fChain->SetBranchAddress("RECOELE_energyScaleUp",RECOELE_energyScaleUp,&b_RECOELE_energyScaleUp);
   fChain->SetBranchAddress("RECOELE_energyScaleDown",RECOELE_energyScaleDown,&b_RECOELE_energyScaleDown);
   fChain->SetBranchAddress("RECOELE_energySigmaUp",RECOELE_energySigmaUp,&b_RECOELE_energySigmaUp);
   fChain->SetBranchAddress("RECOELE_energySigmaDown",RECOELE_energySigmaDown,&b_RECOELE_energySigmaDown);

   fChain->SetBranchAddress("RECOELE_PT", RECOELE_PT, &b_RECOELE_PT);
   fChain->SetBranchAddress("RECOELE_PTError", RECOELE_PTError, &b_RECOELE_PTError);
   fChain->SetBranchAddress("RECOELE_PT_uncorr",RECOELE_PT_uncorr,&b_RECOELE_PT_uncorr);
   fChain->SetBranchAddress("RECOELE_P", RECOELE_P, &b_RECOELE_P);
   fChain->SetBranchAddress("RECOELE_ETA", RECOELE_ETA, &b_RECOELE_ETA);
   fChain->SetBranchAddress("RECOELE_THETA", RECOELE_THETA, &b_RECOELE_THETA);
   fChain->SetBranchAddress("RECOELE_PHI", RECOELE_PHI, &b_RECOELE_PHI);
   fChain->SetBranchAddress("RECOELE_MASS", RECOELE_MASS, &b_RECOELE_MASS);
   fChain->SetBranchAddress("RECOELE_CHARGE", RECOELE_CHARGE, &b_RECOELE_CHARGE);
   fChain->SetBranchAddress("RECOELE_isEcalDriven", RECOELE_isEcalDriven, &b_RECOELE_isEcalDriven);
   fChain->SetBranchAddress("RECOELE_isTrackerDriven", RECOELE_isTrackerDriven, &b_RECOELE_isTrackerDriven);
/*
   fChain->SetBranchAddress("RECOELE_gsftrack_NPixHits", RECOELE_gsftrack_NPixHits, &b_RECOELE_gsftrack_NPixHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_NStripHits", RECOELE_gsftrack_NStripHits, &b_RECOELE_gsftrack_NStripHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_chi2", RECOELE_gsftrack_chi2, &b_RECOELE_gsftrack_chi2);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyB", RECOELE_gsftrack_dxyB, &b_RECOELE_gsftrack_dxyB);
*/
   fChain->SetBranchAddress("RECOELE_gsftrack_dxy", RECOELE_gsftrack_dxy, &b_RECOELE_gsftrack_dxy);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyError", RECOELE_gsftrack_dxyError, &b_RECOELE_gsftrack_dxyError);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzB", RECOELE_gsftrack_dzB, &b_RECOELE_gsftrack_dzB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dz", RECOELE_gsftrack_dz, &b_RECOELE_gsftrack_dz);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzError", RECOELE_gsftrack_dzError, &b_RECOELE_gsftrack_dzError);
   fChain->SetBranchAddress("RECOELE_scl_E", RECOELE_scl_E, &b_RECOELE_scl_E);
   fChain->SetBranchAddress("RECOELE_scl_Et", RECOELE_scl_Et, &b_RECOELE_scl_Et);
   fChain->SetBranchAddress("RECOELE_scl_Eta", RECOELE_scl_Eta, &b_RECOELE_scl_Eta);
   fChain->SetBranchAddress("RECOELE_scl_Phi", RECOELE_scl_Phi, &b_RECOELE_scl_Phi);
/*
   fChain->SetBranchAddress("RECOELE_sigmaIetaIeta", RECOELE_sigmaIetaIeta, &b_RECOELE_sigmaIetaIeta);
   fChain->SetBranchAddress("RECOELE_sigmaEtaEta", RECOELE_sigmaEtaEta, &b_RECOELE_sigmaEtaEta);
   fChain->SetBranchAddress("RECOELE_he", RECOELE_he, &b_RECOELE_he);
   fChain->SetBranchAddress("RECOELE_r9", RECOELE_r9, &b_RECOELE_r9);
   fChain->SetBranchAddress("RECOELE_mva", RECOELE_mva, &b_RECOELE_mva);
   fChain->SetBranchAddress("RECOELE_fbrem", RECOELE_fbrem, &b_RECOELE_fbrem);
   fChain->SetBranchAddress("RECOELE_nbrems", RECOELE_nbrems, &b_RECOELE_nbrems);
   fChain->SetBranchAddress("RECOELE_Class", RECOELE_Class, &b_RECOELE_Class);
   fChain->SetBranchAddress("RECOELE_fbrem_mode", RECOELE_fbrem_mode, &b_RECOELE_fbrem_mode);
   fChain->SetBranchAddress("RECOELE_fbrem_mean", RECOELE_fbrem_mean, &b_RECOELE_fbrem_mean);
   fChain->SetBranchAddress("RECOELE_EGMTRACKISO", RECOELE_EGMTRACKISO, &b_RECOELE_EGMTRACKISO);
   fChain->SetBranchAddress("RECOELE_EGMHCALISO", RECOELE_EGMHCALISO, &b_RECOELE_EGMHCALISO);
   fChain->SetBranchAddress("RECOELE_EGMECALISO", RECOELE_EGMECALISO, &b_RECOELE_EGMECALISO);
   fChain->SetBranchAddress("RECOELE_EGMX", RECOELE_EGMX, &b_RECOELE_EGMX);
*/
   fChain->SetBranchAddress("RECOELE_PFchAllPart", RECOELE_PFchAllPart, &b_RECOELE_PFchAllPart);
   fChain->SetBranchAddress("RECOELE_PFchHad", RECOELE_PFchHad, &b_RECOELE_PFchHad);
   fChain->SetBranchAddress("RECOELE_PFneuHad", RECOELE_PFneuHad, &b_RECOELE_PFneuHad);
   fChain->SetBranchAddress("RECOELE_PFphoton", RECOELE_PFphoton, &b_RECOELE_PFphoton);
   fChain->SetBranchAddress("RECOELE_PFPUchAllPart", RECOELE_PFPUchAllPart, &b_RECOELE_PFPUchAllPart);
   fChain->SetBranchAddress("RECOELE_PFX_dB", RECOELE_PFX_dB, &b_RECOELE_PFX_dB);
   fChain->SetBranchAddress("RECOELE_PFX_rho", RECOELE_PFX_rho, &b_RECOELE_PFX_rho);
   fChain->SetBranchAddress("RECOELE_regEnergy", RECOELE_regEnergy, &b_RECOELE_regEnergy);
   fChain->SetBranchAddress("RECOELE_regEnergyError", RECOELE_regEnergyError, &b_RECOELE_regEnergyError);
   fChain->SetBranchAddress("RECOELE_SIP", RECOELE_SIP, &b_RECOELE_SIP);
   fChain->SetBranchAddress("RECOELE_sclRawE", RECOELE_sclRawE, &b_RECOELE_sclRawE);
   fChain->SetBranchAddress("RECOELE_sclX", RECOELE_sclX, &b_RECOELE_sclX);
   fChain->SetBranchAddress("RECOELE_sclY", RECOELE_sclY, &b_RECOELE_sclY);
   fChain->SetBranchAddress("RECOELE_sclZ", RECOELE_sclZ, &b_RECOELE_sclZ);
   fChain->SetBranchAddress("RECOELE_mvaTrigV0", RECOELE_mvaTrigV0, &b_RECOELE_mvaTrigV0);
   fChain->SetBranchAddress("RECOELE_mvaNonTrigV0", RECOELE_mvaNonTrigV0, &b_RECOELE_mvaNonTrigV0);
   fChain->SetBranchAddress("RECOELE_COV", RECOELE_COV, &b_RECOELE_COV);
   fChain->SetBranchAddress("RECOMU_isPFMu", RECOMU_isPFMu, &b_RECOMU_isPFMu);
   fChain->SetBranchAddress("RECOMU_isMedium", RECOMU_isMedium, &b_RECOMU_isMedium);
   fChain->SetBranchAddress("RECOMU_muInnertrkvalidFraction", RECOMU_muInnertrkvalidFraction, &b_RECOMU_muInnertrkvalidFraction);
   fChain->SetBranchAddress("RECOMU_isGlobalMu", RECOMU_isGlobalMu, &b_RECOMU_isGlobalMu);
   fChain->SetBranchAddress("RECOMU_isStandAloneMu", RECOMU_isStandAloneMu, &b_RECOMU_isStandAloneMu);
   fChain->SetBranchAddress("RECOMU_isTrackerMu", RECOMU_isTrackerMu, &b_RECOMU_isTrackerMu);
   fChain->SetBranchAddress("RECOMU_isCaloMu", RECOMU_isCaloMu, &b_RECOMU_isCaloMu);
   fChain->SetBranchAddress("RECOMU_isTrackerHighPtMu", RECOMU_isTrackerHighPtMu, &b_RECOMU_isTrackerHighPtMu);
   fChain->SetBranchAddress("RECOMU_E", RECOMU_E, &b_RECOMU_E);
   fChain->SetBranchAddress("RECOMU_PT", RECOMU_PT, &b_RECOMU_PT);
   fChain->SetBranchAddress("RECOMU_P", RECOMU_P, &b_RECOMU_P);
   fChain->SetBranchAddress("RECOMU_ETA", RECOMU_ETA, &b_RECOMU_ETA);
   fChain->SetBranchAddress("RECOMU_THETA", RECOMU_THETA, &b_RECOMU_THETA);
   fChain->SetBranchAddress("RECOMU_PHI", RECOMU_PHI, &b_RECOMU_PHI);
   fChain->SetBranchAddress("RECOMU_MASS", RECOMU_MASS, &b_RECOMU_MASS);
   fChain->SetBranchAddress("RECOMU_CHARGE", RECOMU_CHARGE, &b_RECOMU_CHARGE);
   fChain->SetBranchAddress("RECOMU_COV", RECOMU_COV, &b_RECOMU_COV);
   fChain->SetBranchAddress("RECOMU_TRACKISO", RECOMU_TRACKISO, &b_RECOMU_TRACKISO);
   fChain->SetBranchAddress("RECOMU_TRACKISO_SUMPT", RECOMU_TRACKISO_SUMPT, &b_RECOMU_TRACKISO_SUMPT);
   fChain->SetBranchAddress("RECOMU_HCALISO", RECOMU_HCALISO, &b_RECOMU_HCALISO);
   fChain->SetBranchAddress("RECOMU_ECALISO", RECOMU_ECALISO, &b_RECOMU_ECALISO);
   fChain->SetBranchAddress("RECOMU_X", RECOMU_X, &b_RECOMU_X);
   fChain->SetBranchAddress("RECOMU_PFchHad", RECOMU_PFchHad, &b_RECOMU_PFchHad);
   fChain->SetBranchAddress("RECOMU_PFneuHad", RECOMU_PFneuHad, &b_RECOMU_PFneuHad);
   fChain->SetBranchAddress("RECOMU_PFphoton", RECOMU_PFphoton, &b_RECOMU_PFphoton);
   fChain->SetBranchAddress("RECOMU_PFPUchAllPart", RECOMU_PFPUchAllPart, &b_RECOMU_PFPUchAllPart);
   fChain->SetBranchAddress("RECOMU_PFX_dB", RECOMU_PFX_dB, &b_RECOMU_PFX_dB);
   fChain->SetBranchAddress("RECOMU_PFX_rho", RECOMU_PFX_rho, &b_RECOMU_PFX_rho);
   fChain->SetBranchAddress("RECOPFPHOT_PFchHad", RECOPFPHOT_PFchHad, &b_RECOPFPHOT_PFchHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFneuHad", RECOPFPHOT_PFneuHad, &b_RECOPFPHOT_PFneuHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFphoton", RECOPFPHOT_PFphoton, &b_RECOPFPHOT_PFphoton);
   fChain->SetBranchAddress("RECOPFPHOT_PFPUchAllPart", RECOPFPHOT_PFPUchAllPart, &b_RECOPFPHOT_PFPUchAllPart);
   fChain->SetBranchAddress("RECOPFPHOT_PFX_rho", RECOPFPHOT_PFX_rho, &b_RECOPFPHOT_PFX_rho);
   fChain->SetBranchAddress("RECOMU_mutrkPT", RECOMU_mutrkPT, &b_RECOMU_mutrkPT);
   fChain->SetBranchAddress("RECOMU_mutrkPTError", RECOMU_mutrkPTError, &b_RECOMU_mutrkPTError);
   fChain->SetBranchAddress("RECOMU_mutrkDxy", RECOMU_mutrkDxy, &b_RECOMU_mutrkDxy);
   fChain->SetBranchAddress("RECOMU_mutrkDxyError", RECOMU_mutrkDxyError, &b_RECOMU_mutrkDxyError);
   fChain->SetBranchAddress("RECOMU_mutrkDxyB", RECOMU_mutrkDxyB, &b_RECOMU_mutrkDxyB);
   fChain->SetBranchAddress("RECOMU_mutrkDz", RECOMU_mutrkDz, &b_RECOMU_mutrkDz);
   fChain->SetBranchAddress("RECOMU_mutrkDzError", RECOMU_mutrkDzError, &b_RECOMU_mutrkDzError);
   fChain->SetBranchAddress("RECOMU_mutrkDzB", RECOMU_mutrkDzB, &b_RECOMU_mutrkDzB);
   fChain->SetBranchAddress("RECOMU_mutrkChi2PerNdof", RECOMU_mutrkChi2PerNdof, &b_RECOMU_mutrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_mutrkCharge", RECOMU_mutrkCharge, &b_RECOMU_mutrkCharge);
   fChain->SetBranchAddress("RECOMU_mutrkNHits", RECOMU_mutrkNHits, &b_RECOMU_mutrkNHits);
   fChain->SetBranchAddress("RECOMU_mutrkNStripHits", RECOMU_mutrkNStripHits, &b_RECOMU_mutrkNStripHits);
   fChain->SetBranchAddress("RECOMU_mutrkNPixHits", RECOMU_mutrkNPixHits, &b_RECOMU_mutrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mutrkNMuonHits", RECOMU_mutrkNMuonHits, &b_RECOMU_mutrkNMuonHits);
   fChain->SetBranchAddress("RECOMU_mutrktrackerLayersWithMeasurement", RECOMU_mutrktrackerLayersWithMeasurement, &b_RECOMU_mutrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_mubesttrkType", RECOMU_mubesttrkType, &b_RECOMU_mubesttrkType);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxy", RECOMU_mubesttrkDxy, &b_RECOMU_mubesttrkDxy);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxyError", RECOMU_mubesttrkDxyError, &b_RECOMU_mubesttrkDxyError);
   fChain->SetBranchAddress("RECOMU_mubesttrkDz", RECOMU_mubesttrkDz, &b_RECOMU_mubesttrkDz);
   fChain->SetBranchAddress("RECOMU_mubesttrkDzError", RECOMU_mubesttrkDzError, &b_RECOMU_mubesttrkDzError);
   fChain->SetBranchAddress("RECOMU_MatchingMCTruth", RECOMU_MatchingMCTruth, &b_RECOMU_MatchingMCTruth);
   fChain->SetBranchAddress("RECOMU_MatchingMCpT", RECOMU_MatchingMCpT, &b_RECOMU_MatchingMCpT);
   fChain->SetBranchAddress("RECOMU_MatchingMCEta", RECOMU_MatchingMCEta, &b_RECOMU_MatchingMCEta);
   fChain->SetBranchAddress("RECOMU_MatchingMCPhi", RECOMU_MatchingMCPhi, &b_RECOMU_MatchingMCPhi);
   fChain->SetBranchAddress("RECOELE_MatchingMCTruth", RECOELE_MatchingMCTruth, &b_RECOELE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOELE_MatchingMCpT", RECOELE_MatchingMCpT, &b_RECOELE_MatchingMCpT);
   fChain->SetBranchAddress("RECOELE_MatchingMCEta", RECOELE_MatchingMCEta, &b_RECOELE_MatchingMCEta);
   fChain->SetBranchAddress("RECOELE_MatchingMCPhi", RECOELE_MatchingMCPhi, &b_RECOELE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCTruth", RECOPHOT_MatchingMCTruth, &b_RECOPHOT_MatchingMCTruth);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCpT", RECOPHOT_MatchingMCpT, &b_RECOPHOT_MatchingMCpT);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCEta", RECOPHOT_MatchingMCEta, &b_RECOPHOT_MatchingMCEta);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCPhi", RECOPHOT_MatchingMCPhi, &b_RECOPHOT_MatchingMCPhi);
   fChain->SetBranchAddress("RECO_NMU", &RECO_NMU, &b_RECO_NMU);
   fChain->SetBranchAddress("RECO_NELE", &RECO_NELE, &b_RECO_NELE);
   fChain->SetBranchAddress("RECO_NTRACK", &RECO_NTRACK, &b_RECO_NTRACK);
   fChain->SetBranchAddress("RECO_NPHOT", &RECO_NPHOT, &b_RECO_NPHOT);
   fChain->SetBranchAddress("RECOPHOT_PT", RECOPHOT_PT, &b_RECOPHOT_PT);
   fChain->SetBranchAddress("RECOPHOT_ETA", RECOPHOT_ETA, &b_RECOPHOT_ETA);
   fChain->SetBranchAddress("RECOPHOT_PHI", RECOPHOT_PHI, &b_RECOPHOT_PHI);
   fChain->SetBranchAddress("RECOPHOT_THETA", RECOPHOT_THETA, &b_RECOPHOT_THETA);
   fChain->SetBranchAddress("RECO_NPFPHOT", &RECO_NPFPHOT, &b_RECO_NPFPHOT);
   fChain->SetBranchAddress("RECOPFPHOT_PT", RECOPFPHOT_PT, &b_RECOPFPHOT_PT);
   fChain->SetBranchAddress("RECOPFPHOT_PTError", RECOPFPHOT_PTError, &b_RECOPFPHOT_PTError);
   fChain->SetBranchAddress("RECOPFPHOT_ETA", RECOPFPHOT_ETA, &b_RECOPFPHOT_ETA);
   fChain->SetBranchAddress("RECOPFPHOT_PHI", RECOPFPHOT_PHI, &b_RECOPFPHOT_PHI);
   fChain->SetBranchAddress("RECOPFPHOT_THETA", RECOPFPHOT_THETA, &b_RECOPFPHOT_THETA);
   fChain->SetBranchAddress("BeamSpot_X", &BeamSpot_X, &b_BeamSpot_X);
   fChain->SetBranchAddress("BeamSpot_Y", &BeamSpot_Y, &b_BeamSpot_Y);
   fChain->SetBranchAddress("BeamSpot_Z", &BeamSpot_Z, &b_BeamSpot_Z);
   fChain->SetBranchAddress("RECO_NVTX", &RECO_NVTX, &b_RECO_NVTX);
   fChain->SetBranchAddress("RECO_PFJET_N", &RECO_PFJET_N, &b_RECO_PFJET_N);
   fChain->SetBranchAddress("RECO_PFJET_CHARGE", RECO_PFJET_CHARGE, &b_RECO_PFJET_CHARGE);
   fChain->SetBranchAddress("RECO_PFJET_ET", RECO_PFJET_ET, &b_RECO_PFJET_ET);
   fChain->SetBranchAddress("RECO_PFJET_PT", RECO_PFJET_PT, &b_RECO_PFJET_PT);
   fChain->SetBranchAddress("RECO_PFJET_PT_UP", RECO_PFJET_PT_UP, &b_RECO_PFJET_PT_UP);
   fChain->SetBranchAddress("RECO_PFJET_PT_DOW", RECO_PFJET_PT_DOW, &b_RECO_PFJET_PT_DOW);
   fChain->SetBranchAddress("RECO_PFJET_ETA", RECO_PFJET_ETA, &b_RECO_PFJET_ETA);
   fChain->SetBranchAddress("RECO_PFJET_PHI", RECO_PFJET_PHI, &b_RECO_PFJET_PHI);
   fChain->SetBranchAddress("RECO_PFJET_PUID", RECO_PFJET_PUID, &b_RECO_PFJET_PUID);
   fChain->SetBranchAddress("RECO_PFJET_PUID_MVA", RECO_PFJET_PUID_MVA, &b_RECO_PFJET_PUID_MVA);

   fChain->SetBranchAddress( "RECO_PFJET_nconstituents",  RECO_PFJET_nconstituents,  &b_RECO_PFJET_nconstituents);
   fChain->SetBranchAddress( "RECO_PFJET_NCH",  RECO_PFJET_NCH,  &b_RECO_PFJET_NCH);
   fChain->SetBranchAddress( "RECO_PFJET_NHF",  RECO_PFJET_NHF,  &b_RECO_PFJET_NHF);
   fChain->SetBranchAddress( "RECO_PFJET_NEF",  RECO_PFJET_NEF, &b_RECO_PFJET_NEF);
   fChain->SetBranchAddress( "RECO_PFJET_CHF", RECO_PFJET_CHF, &b_RECO_PFJET_CHF);
   fChain->SetBranchAddress( "RECO_PFJET_CEF", RECO_PFJET_CEF, &b_RECO_PFJET_CEF);
   fChain->SetBranchAddress( "RECO_PFJET_MUF", RECO_PFJET_MUF, &b_RECO_PFJET_MUF);

   fChain->SetBranchAddress("RHO_ele", &RHO_ele, &b_RHO_ele);
   fChain->SetBranchAddress("RHO_mu", &RHO_mu, &b_RHO_mu);
   fChain->SetBranchAddress("RECO_CALOMET", &RECO_CALOMET, &b_RECO_CALOMET);
   fChain->SetBranchAddress("RECO_PFMET", &RECO_PFMET, &b_RECO_PFMET);
   fChain->SetBranchAddress("RECO_PFMET_X", &RECO_PFMET_X, &b_RECO_PFMET_X);
   fChain->SetBranchAddress("RECO_PFMET_Y", &RECO_PFMET_Y, &b_RECO_PFMET_Y);
   fChain->SetBranchAddress("RECO_PFMET_PHI", &RECO_PFMET_PHI, &b_RECO_PFMET_PHI);
   fChain->SetBranchAddress("RECO_PFMET_THETA", &RECO_PFMET_THETA, &b_RECO_PFMET_THETA);
   fChain->SetBranchAddress("RECO_TCMET", &RECO_TCMET, &b_RECO_TCMET);
   fChain->SetBranchAddress("RECO_CORMETMUONS", &RECO_CORMETMUONS, &b_RECO_CORMETMUONS);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PT", tCHighEff_BTagJet_PT, &b_tCHighEff_BTagJet_PT);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PT", tCHighPur_BTagJet_PT, &b_tCHighPur_BTagJet_PT);
   fChain->SetBranchAddress("cSV_BTagJet_PT", cSV_BTagJet_PT, &b_cSV_BTagJet_PT);
   fChain->SetBranchAddress("tCHighEff_BTagJet_ETA", tCHighEff_BTagJet_ETA, &b_tCHighEff_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighPur_BTagJet_ETA", tCHighPur_BTagJet_ETA, &b_tCHighPur_BTagJet_ETA);
   fChain->SetBranchAddress("cSV_BTagJet_ETA", cSV_BTagJet_ETA, &b_cSV_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PHI", tCHighEff_BTagJet_PHI, &b_tCHighEff_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PHI", tCHighPur_BTagJet_PHI, &b_tCHighPur_BTagJet_PHI);
   fChain->SetBranchAddress("cSV_BTagJet_PHI", cSV_BTagJet_PHI, &b_cSV_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighEff_BTagJet_DISCR", tCHighEff_BTagJet_DISCR, &b_tCHighEff_BTagJet_DISCR);
   fChain->SetBranchAddress("tCHighPur_BTagJet_DISCR", tCHighPur_BTagJet_DISCR, &b_tCHighPur_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_DISCR", cSV_BTagJet_DISCR, &b_cSV_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_ET", cSV_BTagJet_ET, &b_cSV_BTagJet_ET);
 
   Notify();
}

Bool_t HZZ4LeptonsAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HZZ4LeptonsAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HZZ4LeptonsAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void HZZ4LeptonsAnalysis::printelebnn(int i){    
  bnn_file 
                << RECOELE_PT[i]  << " "  
                << RECOELE_ETA[i] << " "  
                << RECOELE_PHI[i] << " "  
                << RECOELE_CHARGE[i] << " "
                << RECOELE_PFX_rho[i] << " "
                << RECOELE_SIP[i] << " ";
}


void HZZ4LeptonsAnalysis::printmubnn(int i){    
  bnn_file 
                << RECOMU_PT[i]  << " "  
                << RECOMU_ETA[i] << " "  
                << RECOMU_PHI[i] << " "  
                << RECOMU_CHARGE[i] << " "
                << RECOMU_PFX_dB[i] << " ";
}

#endif // #ifdef HZZ4LeptonsAnalysis_cxx
