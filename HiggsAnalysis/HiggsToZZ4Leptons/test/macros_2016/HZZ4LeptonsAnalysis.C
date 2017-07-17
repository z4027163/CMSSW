#define HZZ4LeptonsAnalysis_cxx
#include "HZZ4LeptonsAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
//#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TNtuple.h>

#include <TMath.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <libgen.h>


using namespace std;
// using namespace RooFit;
// using namespace meMCFM;

const double Zmass = 91.188; // nominal Z boson mass
const double mPI = 3.141592654; 

void HZZ4LeptonsAnalysis::Loop(Char_t *output)
{

   if (fChain == 0) return;

   // Declare MEM class
   
   // MuonCalibrator
      
   // BNN
   Char_t datasetChar[500],bnnOUT[500],eventsOUT[500];
  
   cout << "The output file is " << output << endl;
   TString out = output;
   TString datasetName=out.ReplaceAll(".root","");
   sprintf(datasetChar,"%s",datasetName.Data());
   sprintf(bnnOUT,"%s_bnn.txt",datasetName.Data());
   sprintf(eventsOUT,"%s_bnn.root",datasetName.Data());
   cout << "bnnOUT= " << bnnOUT << endl;
   bnn_file.open(bnnOUT);

     
   // Book txt file for candidate events
   Char_t txtOUT[500];
   sprintf(txtOUT,"%s_txt.txt",datasetName.Data());
   cout << "Opening a txt file with candidate events " << txtOUT << endl;
   ofstream output_txt;
   output_txt.open(txtOUT);
   // isSignal DM
   
   // Pileup reweighting 2016 data vs Spring16 MC in 80x

   TFile *_filePU;
   _filePU= TFile::Open("pileup_MC_80x_271036-276811_69200.root");
   TH1D *puweight = (TH1D*)_filePU->Get("puweight");

   /////////////Lepton Efficiency Scale Factrons/////////////
   // Load histograms
   //

   // Book root file (for output):
   TFile * theFile = new TFile(output,"RECREATE");


   double DELTAPHI( double , double ) ; //call the function  
   double invmass (float M1, float PT1, float ETA1, float PHI1, float M2, float PT2, float ETA2, float PHI2 );
   
   // Book relevant variables -- counters:

   int N_0 = 0;	  // MC truth & acceptance
   int N_01 = 0;
   int N_02 = 0;

   int N_1 = 0;	  // Skim
   int N_2 = 0;

   int N_3a = 0;
   int N_3_FSR = 0;
   int N_3b = 0;

   int N_4a = 0;
   int N_4b = 0;
   int N_4c = 0;
   int N_4d = 0;

   int N_5 = 0;
   int N_6 = 0;
   int N_7 = 0;
   int N_8 = 0;
   int N_8_a = 0;
   int N_9 = 0;

   // counter weighted
   double N_0_w = 0;	  // MC truth & acceptance
   double N_01_w = 0;
   double N_02_w = 0;

   double N_1_w = 0;	  // Skim
   double N_2_w = 0;

   double N_3a_w = 0;
   double N_3_FSR_w = 0;
   double N_3b_w = 0;

   double N_4a_w = 0;
   double N_4b_w = 0;
   double N_4c_w = 0;
   double N_4d_w = 0;

   double N_5_w = 0;
   double N_6_w = 0;


   // Book Histos ***
   TH1D *nEvent_4l_w = new TH1D("nEvent_4l_w", "nEventComplete Weightd", 20, 0., 20.);
   TH1D *nEvent_4l = new TH1D("nEvent_4l", "nEventComplete", 20, 0., 20.);

   TH1F *Gen_H_MASS              = new TH1F("Gen_H_MASS", "Gen_H_MASS",8000,0.,2000.);  
   TH1F *Gen_H_MASS_ReWeighted   = new TH1F("Gen_H_MASS_ReWeighted", "Gen_H_MASS_ReWeighted",8000,0.,2000.);  
   TH1F *Gen_H_MASS_ReWeightedP  = new TH1F("Gen_H_MASS_ReWeightedP", "Gen_H_MASS_ReWeightedP",8000,0.,2000.);  
   TH1F *Gen_H_MASS_ReWeightedM  = new TH1F("Gen_H_MASS_ReWeightedM", "Gen_H_MASS_ReWeightedM",8000,0.,2000.);  

   // Pileup reweighting
   TH1F *hPUvertices             = new TH1F("hPUvertices", "hPUvertices",70,0.,70.);  
   TH1F *hPUvertices_ReWeighted  = new TH1F("hPUvertices_ReWeighted", "hPUvertices_ReWeighted",70,0.,70.);  

   //step 3
   TH1F * hMZ1_3 = new TH1F("hMZ1_3", "Mass of Z1 after selection step 3", 200 , -0.5 , 199.5 );
   hMZ1_3->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_3 = new TH1F("hPtZ1_3", "Pt of Z1 after selection step 3", 200 , -0.5 , 199.5 );
   hPtZ1_3->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_3 = new TH1F("hYZ1_3", "Y of Z1 after selection step 3", 500 , -5. , 5.);
   hYZ1_3->SetXTitle("Y_Z1");
   
   TH1F * hPtLep1_3 = new TH1F("hPtLep1_3", "Pt of Lep1 after selection step 3", 200 , -0.5 , 199.5 );
   hPtLep1_3->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_3 = new TH1F("hEtaLep1_3", "Y of Lep1 after selection step 3", 500 , -5. , 5. );
   hYLep1_3->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_3 = new TH1F("hIsoLep1_3", "Isolation of Lep1 after selection step 3", 2000 , -10., 10.);
   hIsoLep1_3->SetXTitle("Iso");
   TH1F * hSipLep1_3 = new TH1F("hSipLep1_3", "Sip of Lep1 after selection step 3",  1000 , -20. , 40. );
   hSipLep1_3->SetXTitle("Sip");
   TH1F * hIpLep1_3 = new TH1F("hIpLep1_3", "Ip of Lep1 after selection step 3",  1000 , -20. , 40. );
   hIpLep1_3->SetXTitle("Ip");
   TH1F * hIpErLep1_3 = new TH1F("hIpErLep1_3", "IpEr of Lep1 after selection step 3",  1000 , 0. , 10. );
   hIpErLep1_3->SetXTitle("IpEr");

   TH1F * hPtLep2_3 = new TH1F("hPtLep2_3", "Pt of Lep2 after selection step 3", 200 , -0.5 , 199.5 );
   hPtLep2_3->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_3 = new TH1F("hEtaLep2_3", "Y of Lep2 after selection step 3", 500 , -5. , 5. );
   hYLep2_3->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_3 = new TH1F("hIsoLep2_3", "Isolation of Lep2 after selection step 3", 2000 , -10., 10. );
   hIsoLep2_3->SetXTitle("Iso");
   TH1F * hSipLep2_3 = new TH1F("hSipLep2_3", "Sip of Lep2 after selection step 3",  1000 , -20. , 40. );
   hSipLep2_3->SetXTitle("Sip");
   TH1F * hIpLep2_3 = new TH1F("hIpLep2_3", "Ip of Lep2 after selection step 3",  1000 , -20. , 40. );
   hIpLep2_3->SetXTitle("Ip");
   TH1F * hIpErLep2_3 = new TH1F("hIpErLep2_3", "IpEr of Lep2 after selection step 3",  1000 , 0. , 10. );
   hIpErLep2_3->SetXTitle("IpEr");

   TH1F * hIso_3 = new TH1F("hIso_3", "Isolation maxima after selection step 3", 2000 , -10. , 10. );
   hIso_3->SetXTitle("Iso");
   TH1F * hSip_3 = new TH1F("hSip_3", "Sip maxima after selection step 3",  1000 , -20. , 40. );
   hSip_3->SetXTitle("Sip");
   TH1F * hIp_3 = new TH1F("hIp_3", "Ip maxima after selection step 3",  1000 , -20. , 40. );
   hIp_3->SetXTitle("Ip");

   TH1F * hDjj_3 = new TH1F("hDjj_3", "Delta jets vbf selection step 1", 200, -19.5, 19.5 );
   hDjj_3->SetXTitle("Delta jets");
   TH1F * hMjj_3 = new TH1F("hMjj_3", "Mass jets vbf selection step 1", 200, -0.5, 499.5 );
   hMjj_3->SetXTitle("Mass jets");
   TH1F * hVD_3 = new TH1F("hVD_3", "Discriminant vbf selection step 1", 200, -0.5, 9.5 );
   hMjj_3->SetXTitle("Discriminant");
   
   TH1F * hPFMET_3 = new TH1F("hPFMET_3", "PF MET after selection step 3", 1000 , 0., 1000.);

   //step 5
   TH1F * hM4l_5 = new TH1F("hM4l_5", "Mass of four leptons after selection step 5", 1200, 4.5,1204.5 );
   hM4l_5->SetXTitle("4 lepton mass  (GeV)");


   TH1F * hMZ1_5 = new TH1F("hMZ1_5", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5 = new TH1F("hPtZ1_5", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5 = new TH1F("hYZ1_5", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5->SetXTitle("Y_Z1");

   TH1F * hMZ2_5 = new TH1F("hMZ2_5", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5 = new TH1F("hPtZ2_5", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5 = new TH1F("hYZ2_5", "Y of Z2 after selection step 5", 500 , -5. , 5. );
   hYZ2_5->SetXTitle("Y_Z2");

   TH1F * hPtLep1_5 = new TH1F("hPtLep1_5", "Pt of Lep1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep1_5->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_5 = new TH1F("hEtaLep1_5", "Y of Lep1 after selection step 5", 500 , -5. , 5. );
   hYLep1_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_5 = new TH1F("hIsoLep1_5", "Isolation of Lep1 after selection step 5", 2000 , -10. , 10. );
   hIsoLep1_5->SetXTitle("Iso");
   TH1F * hSipLep1_5 = new TH1F("hSipLep1_5", "Sip of Lep1 after selection step 5",  1000 , -20. , 40. );
   hSipLep1_5->SetXTitle("Sip");
   TH1F * hIpLep1_5 = new TH1F("hIpLep1_5", "Ip of Lep1 after selection step 5",  1000 , -20. , 40. );
   hIpLep1_5->SetXTitle("Ip");
   TH1F * hIpErLep1_5 = new TH1F("hIpErLep1_5", "IpEr of Lep1 after selection step 5",  1000 , 0. , 10. );
   hIpErLep1_5->SetXTitle("IpEr");

   TH1F * hPtLep2_5 = new TH1F("hPtLep2_5", "Pt of Lep2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep2_5->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_5 = new TH1F("hEtaLep2_5", "Y of Lep2 after selection step 5", 500 , -5. , 5. );
   hYLep2_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_5 = new TH1F("hIsoLep2_5", "Isolation of Lep2 after selection step 5", 2000 , -10. , 10. );
   hIsoLep2_5->SetXTitle("Iso");
   TH1F * hSipLep2_5 = new TH1F("hSipLep2_5", "Sip of Lep2 after selection step 5",  1000 , -20. , 40. );
   hSipLep2_5->SetXTitle("Sip");
   TH1F * hIpLep2_5 = new TH1F("hIpLep2_5", "Ip of Lep2 after selection step 5",  1000 , -20. , 40. );
   hIpLep2_5->SetXTitle("Ip");
   TH1F * hIpErLep2_5 = new TH1F("hIpErLep2_5", "IpEr of Lep2 after selection step 5",  1000 , 0. , 10. );
   hIpErLep2_5->SetXTitle("IpEr");
   
   TH1F * hPtLep3_5 = new TH1F("hPtLep3_5", "Pt of Lep3 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep3_5->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_5 = new TH1F("hEtaLep3_5", "Y of Lep3 after selection step 5", 500 , -5. , 5. );
   hYLep3_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_5 = new TH1F("hIsoLep3_5", "Isolation of Lep3 after selection step 5", 2000 , -10. , 10. );
   hIsoLep3_5->SetXTitle("Iso");
   TH1F * hSipLep3_5 = new TH1F("hSipLep3_5", "Sip of Lep3 after selection step 5",  1000 , -20. , 40. );
   hSipLep3_5->SetXTitle("Sip");
   TH1F * hIpLep3_5 = new TH1F("hIpLep3_5", "Ip of Lep3 after selection step 5",  1000 , -20. , 40. );
   hIpLep3_5->SetXTitle("Ip");
   TH1F * hIpErLep3_5 = new TH1F("hIpErLep3_5", "IpEr of Lep3 after selection step 5",  1000 , 0. , 10. );
   hIpErLep3_5->SetXTitle("IpEr");

   TH1F * hPtLep4_5 = new TH1F("hPtLep4_5", "Pt of Lep4 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep4_5->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_5 = new TH1F("hEtaLep4_5", "Y of Lep4 after selection step 5", 50 , -5. , 5. );
   hYLep4_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_5 = new TH1F("hIsoLep4_5", "Isolation of Lep4 after selection step 5", 2000 , -10. , 10. );
   hIsoLep4_5->SetXTitle("Iso");
   TH1F * hSipLep4_5 = new TH1F("hSipLep4_5", "Sip of Lep4 after selection step 5",  1000 , -20. , 40. );
   hSipLep4_5->SetXTitle("Sip");
   TH1F * hIpLep4_5 = new TH1F("hIpLep4_5", "Ip of Lep4 after selection step 5",  1000 , -20. , 40. );
   hIpLep4_5->SetXTitle("Ip");
   TH1F * hIpErLep4_5 = new TH1F("hIpErLep4_5", "IpEr of Lep4 after selection step 5",  1000 , 0. , 10. );
   hIpErLep4_5->SetXTitle("IpEr");


   TH1F * hIso_5 = new TH1F("hIso_5", "Isolation maxima after selection step 5", 2000 , -10. , 10. );
   hIso_5->SetXTitle("Iso");
   TH1F * hSip_5 = new TH1F("hSip_5", "Sip maxima after selection step 5",  1000 , -20. , 40. );
   hSip_5->SetXTitle("Sip");
   TH1F * hIp_5 = new TH1F("hIp_5", "Ip maxima after selection step 5",  1000 , -20. , 40. );
   hIp_5->SetXTitle("Ip");
   

   //step 6
   
   TH1F * hminMll_6 = new TH1F("hminMll_6", "minMll at selection step 6", 400 , 0. , 200.);
   hminMll_6->SetXTitle("minMll  (GeV)");
   

   //step 7
   TH1F * hM4l_7 = new TH1F("hM4l_7", "Mass of four leptons after selection step 7", 1200, 4.5,1204.5 );
   hM4l_7->SetXTitle("4 lepton mass  (GeV)");

   TH1F * hMZ1_7 = new TH1F("hMZ1_7", "Mass of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hMZ1_7->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_7 = new TH1F("hPtZ1_7", "Pt of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ1_7->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_7 = new TH1F("hYZ1_7", "Y of Z1 after selection step 7", 500 , -5. , 5.);
   hYZ1_7->SetXTitle("Y_Z1");

   TH1F * hMZ2_7 = new TH1F("hMZ2_7", "Mass of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hMZ2_7->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7 = new TH1F("hPtZ2_7", "Pt of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ2_7->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7 = new TH1F("hYZ2_7", "Y of Z2 after selection step 7", 500 , -5. , 5.);
   hYZ2_7->SetXTitle("Y_Z2");

   TH1F * hPtLep1_7 = new TH1F("hPtLep1_7", "Pt of Lep1 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep1_7->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_7 = new TH1F("hEtaLep1_7", "Y of Lep1 after selection step 7", 500 , -5. , 5. );
   hYLep1_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_7 = new TH1F("hIsoLep1_7", "Isolation of Lep1 after selection step 7", 2000 , -10. , 10. );
   hIsoLep1_7->SetXTitle("Iso");
   TH1F * hSipLep1_7 = new TH1F("hSipLep1_7", "Sip of Lep1 after selection step 7",  1000 , -20. , 40. );
   hSipLep1_7->SetXTitle("Sip");
   TH1F * hIpLep1_7 = new TH1F("hIpLep1_7", "Ip of Lep1 after selection step 7",  1000 , -20. , 40. );
   hIpLep1_7->SetXTitle("Ip");
   TH1F * hIpErLep1_7 = new TH1F("hIpErLep1_7", "IpEr of Lep1 after selection step 7",  1000 , 0. , 10. );
   hIpErLep1_7->SetXTitle("IpEr");


   TH1F * hPtLep2_7 = new TH1F("hPtLep2_7", "Pt of Lep2 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep2_7->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_7 = new TH1F("hEtaLep2_7", "Y of Lep2 after selection step 7", 500 , -5. , 5. );
   hYLep2_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_7 = new TH1F("hIsoLep2_7", "Isolation of Lep2 after selection step 7", 2000 , -10. , 10. );
   hIsoLep2_7->SetXTitle("Iso");
   TH1F * hSipLep2_7 = new TH1F("hSipLep2_7", "Sip of Lep2 after selection step 7",  1000 , -20. , 40. );
   hSipLep2_7->SetXTitle("Sip");
   TH1F * hIpLep2_7 = new TH1F("hIpLep2_7", "Ip of Lep2 after selection step 7",  1000 , -20. , 40. );
   hIpLep2_7->SetXTitle("Ip");
   TH1F * hIpErLep2_7 = new TH1F("hIpErLep2_7", "IpEr of Lep2 after selection step 7",  1000 , 0. , 10. );
   hIpErLep2_7->SetXTitle("IpEr");

   TH1F * hPtLep3_7 = new TH1F("hPtLep3_7", "Pt of Lep3 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep3_7->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_7 = new TH1F("hEtaLep3_7", "Y of Lep3 after selection step 7", 500 , -5. , 5. );
   hYLep3_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_7 = new TH1F("hIsoLep3_7", "Isolation of Lep3 after selection step 7", 2000 , -10. , 10. );
   hIsoLep3_7->SetXTitle("Iso");
   TH1F * hSipLep3_7 = new TH1F("hSipLep3_7", "Sip of Lep3 after selection step 7",  1000 , -20. , 40. );
   hSipLep3_7->SetXTitle("Sip");
   TH1F * hIpLep3_7 = new TH1F("hIpLep3_7", "Ip of Lep3 after selection step 7",  1000 , -20. , 40. );
   hIpLep3_7->SetXTitle("Ip");
   TH1F * hIpErLep3_7 = new TH1F("hIpErLep3_7", "IpEr of Lep3 after selection step 7",  1000 , 0. , 10. );
   hIpErLep3_7->SetXTitle("IpEr");

   TH1F * hPtLep4_7 = new TH1F("hPtLep4_7", "Pt of Lep4 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep4_7->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_7 = new TH1F("hEtaLep4_7", "Y of Lep4 after selection step 7", 500 , -5. , 5. );
   hYLep4_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_7 = new TH1F("hIsoLep4_7", "Isolation of Lep4 after selection step 7", 2000 , -10. , 10. );
   hIsoLep4_7->SetXTitle("Iso");
   TH1F * hSipLep4_7 = new TH1F("hSipLep4_7", "Sip of Lep4 after selection step 7",  1000 , -20. , 40. );
   hSipLep4_7->SetXTitle("Sip");
   TH1F * hIpLep4_7 = new TH1F("hIpLep4_7", "Ip of Lep4 after selection step 7",  1000 , -20. , 40. );
   hIpLep4_7->SetXTitle("Ip");
   TH1F * hIpErLep4_7 = new TH1F("hIpErLep4_7", "IpEr of Lep4 after selection step 7",  1000 , 0. , 10. );
   hIpErLep4_7->SetXTitle("IpEr");

   TH1F * hIso_7 = new TH1F("hIso_7", "Isolation maxima after selection step 7", 2000 , -10. , 10. );
   hIso_7->SetXTitle("Iso");
   TH1F * hSip_7 = new TH1F("hSip_7", "Sip maxima after selection step 7",  1000 , -20. , 40. );
   hSip_7->SetXTitle("Sip");
   TH1F * hIp_7 = new TH1F("hIp_7", "Ip maxima after selection step 7",  1000 , -20. , 40. );
   hIp_7->SetXTitle("Ip");
   
   //step 8

   // Step 9 with PFMET cut
   

   //global histos (during step 2..)
   
   TH1F * hN_loose_mu = new TH1F("hN_loose_mu", "N_loose_mu", 30 , 0. , 30. );
   hN_loose_mu->SetXTitle("N_loose_mu");
   TH1F * hN_loose_e = new TH1F("hN_loose_e", "N_loose_e", 30 , 0. , 30. );
   hN_loose_e->SetXTitle("N_loose_e");

   TH1F * hIso_loose_mu = new TH1F("hIso_loose_mu", "Isolation maxima after loose selection ", 2000 , -10. , 10. );
   hIso_loose_mu->SetXTitle("Iso");
   TH1F * hSip_loose_mu = new TH1F("hSip_loose_mu", "Sip maxima after loose selection ",  1000 , -20. , 40. );
   hSip_loose_mu->SetXTitle("Sip");
   TH1F * hIp_loose_mu = new TH1F("hIp_loose_mu", "Ip maxima after loose selection ",  1000 , -20. , 40. );
   hIp_loose_mu->SetXTitle("Ip");

   TH1F * hIso_loose_e = new TH1F("hIso_loose_e", "Isolation maxima after loose selection ", 2000 , -10. , 10. );
   hIso_loose_e->SetXTitle("Iso");
   TH1F * hSip_loose_e = new TH1F("hSip_loose_e", "Sip maxima after loose selection ",  1000 , -20. , 40. );
   hSip_loose_e->SetXTitle("Sip");
   TH1F * hIp_loose_e = new TH1F("hIp_loose_e", "Ip maxima after loose selection ",  1000 , -20. , 40. );
   hIp_loose_e->SetXTitle("Ip");


   TH1F * hN_good_lep = new TH1F("hN_good_lep", "N_good_lep", 30 , 0. , 30. );
   hN_good_lep->SetXTitle("N_good_lep");

   TH1F * hN_good_mu = new TH1F("hN_good_mu", "N_good_mu", 30 , 0. , 30. );

   hN_good_mu->SetXTitle("N_good_mu");
   TH1F * hN_good_ele = new TH1F("hN_good_ele", "N_good_ele", 30 , 0. , 30. );
   hN_good_ele->SetXTitle("N_good_ele");
   TH1F * hN_good_phot = new TH1F("hN_good_phot", "N_good_phot", 30 , 0. , 30. );
   hN_good_phot->SetXTitle("N_good_phot");

   TH1F * hMELA_8 = new TH1F("hMELA_8", "MELA after selection step 8", 300,-0.00166,1.00166 );
   hMELA_8->SetXTitle("MELA discriminant (4mu)");
   TH2F * hMELA_vs_M4l_8 = new TH2F("hMELA_vs_M4l_8", "MELA after selection step 8",1200, 4.5,1204.5,300,-0.00166,1.00166 );
   hMELA_vs_M4l_8->SetXTitle("MELA discriminant (4mu)");

   TH1F * hMELA_9 = new TH1F("hMELA_9", "MELA after selection step 9", 300,-0.00166,1.00166 );
   hMELA_9->SetXTitle("MELA discriminant (4mu)");
   TH2F * hMELA_vs_M4l_9 = new TH2F("hMELA_vs_M4l_9", "MELA after selection step 9",1200, 4.5,1204.5,300,-0.00166,1.00166 );
   hMELA_vs_M4l_9->SetXTitle("MELA discriminant (4mu)");
   
   TH1I * hVBF_PUID = new TH1I("hVBF_PUID", "PUID after step 8", 6,-1,4);
   hVBF_PUID->SetXTitle("PUID VBF");



   //PFJET Plots
   /* TH1I * hN_PFJET_6 = new TH1I("hN_PFJET_6", "Number of PFJets after step 6", 102, 0, 102);
   hN_PFJET_6->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_6 = new TH1F("hChg_PFJET_6", "Charge of PFJets after step 6", 50, -10., 10.);
   hChg_PFJET_6->SetXTitle("Charge");
   TH1F * hEt_PFJET_6 = new TH1F("hEt_PFJET_6", "Et of PFJets after step 6", 201, -0.5, 400.5 );
   hEt_PFJET_6->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_6 = new TH1F("hPt_PFJET_6", "Pt of PFJets after step 6", 201, -0.5, 400.5 );
   hPt_PFJET_6->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_6 = new TH1F("hEta_PFJET_6", "Eta of PFJets after step 6", 50, -10., 10.);
   hEta_PFJET_6->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_6 = new TH1F("hPhi_PFJET_6", "Phi of PFJets after step 6", 50, -10., 10.);
   hPhi_PFJET_6->SetXTitle("phi_PFJET");

   TH1I * hN_PFJET_8 = new TH1I("hN_PFJET_8", "Number of PFJets after step 8", 102, 0, 102);
   hN_PFJET_8->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_8 = new TH1F("hChg_PFJET_8", "Charge of PFJets after step 8", 50, -10., 10.);
   hChg_PFJET_8->SetXTitle("Charge");
   TH1F * hEt_PFJET_8 = new TH1F("hEt_PFJET_8", "Et of PFJets after step 8", 201, -0.5, 400.5 );
   hEt_PFJET_8->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_8 = new TH1F("hPt_PFJET_8", "Pt of PFJets after step 8", 201, -0.5, 400.5 );
   hPt_PFJET_8->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_8 = new TH1F("hEta_PFJET_8", "Eta of PFJets after step 8", 50, -10., 10.);
   hEta_PFJET_8->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_8 = new TH1F("hPhi_PFJET_8", "Phi of PFJets after step 8", 50, -10., 10.);
   hPhi_PFJET_8->SetXTitle("phi_PFJET");

   TH1I * hN_PFJET_VBF = new TH1I("hN_PFJET_VBF", "Number of PFJets after VBF selection", 102, 0, 102);
   hN_PFJET_VBF->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_VBF = new TH1F("hChg_PFJET_VBF", "Charge of PFJets after VBF selection", 50, -10., 10.);
   hChg_PFJET_VBF->SetXTitle("Charge");
   TH1F * hEt_PFJET_VBF = new TH1F("hEt_PFJET_VBF", "Et of PFJets after step VBF selection", 201, -0.5, 400.5 );
   hEt_PFJET_VBF->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_VBF = new TH1F("hPt_PFJET_VBF", "Pt of PFJets after step VBF selection", 201, -0.5, 400.5 );
   hPt_PFJET_VBF->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_VBF = new TH1F("hEta_PFJET_VBF", "Eta of PFJets after step VBF selection", 50, -10., 10.);
   hEta_PFJET_VBF->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_VBF = new TH1F("hPhi_PFJET_VBF", "Phi of PFJets after step VBF selection", 50, -10., 10.);
   hPhi_PFJET_VBF->SetXTitle("phi_PFJET");
   */


   // end book histo ***
      
   TTree *newtree = new TTree("HZZ4LeptonsAnalysisReduced", "reduced ttree");
  
   // Add branches to output rootuple 
   Float_t f_weight, f_int_weight, f_pu_weight, f_eff_weight, f_lept1_pt, f_lept1_eta, f_lept1_phi, f_lept1_charge, f_lept1_pfx, f_lept1_sip, f_lept1_mvaid, f_lept2_pt, f_lept2_eta, f_lept2_phi, f_lept2_charge, f_lept2_pfx, f_lept2_sip, f_lept2_mvaid, f_lept3_pt, f_lept3_eta, f_lept3_phi, f_lept3_charge, f_lept3_pfx, f_lept3_sip, f_lept3_mvaid, f_lept4_pt, f_lept4_eta, f_lept4_phi, f_lept4_charge, f_lept4_pfx, f_lept4_sip, f_lept4_mvaid, f_iso_max, f_sip_max, f_Z1mass, f_Z2mass, f_angle_costhetastar, f_angle_costheta1, f_angle_costheta2, f_angle_phi, f_angle_phistar1, f_eta4l, f_pt4l, f_mass4l, f_mass4lErr, f_njets_pass, f_deltajj, f_massjj, f_D_jet, f_jet1_pt, f_jet1_eta, f_jet1_phi, f_jet1_e, f_jet2_pt, f_jet2_eta, f_jet2_phi, f_jet2_e;
   Float_t f_D_bkg_kin,f_D_bkg,f_D_gg,f_D_g4,f_Djet_VAJHU; 
   Float_t f_pfmet;
 
   Int_t f_run, f_lumi, f_event;
   
   TBranch *b_run= newtree->Branch("f_run", &f_run,"f_run/I");
   TBranch *b_lumi= newtree->Branch("f_lumi", &f_lumi,"f_lumi/I");    
   TBranch *b_event= newtree->Branch("f_event", &f_event,"f_event/I");    
   
   TBranch *b_weight= newtree->Branch("f_weight", &f_weight,"f_weight/F");
   TBranch *b_int_weight= newtree->Branch("f_int_weight", &f_int_weight,"f_int_weight/F");
   TBranch *b_pu_weight= newtree->Branch("f_pu_weight", &f_pu_weight,"f_pu_weight/F");
   TBranch *b_eff_weight= newtree->Branch("f_eff_weight", &f_eff_weight,"f_eff_weight/F");
   TBranch *b_lept1_pt= newtree->Branch("f_lept1_pt", &f_lept1_pt,"f_lept1_pt/F");
   TBranch *b_lept1_eta= newtree->Branch("f_lept1_eta", &f_lept1_eta,"f_lept1_eta/F");
   TBranch *b_lept1_phi= newtree->Branch("f_lept1_phi", &f_lept1_phi,"f_lept1_phi/F");
   TBranch *b_lept1_charge= newtree->Branch("f_lept1_charge", &f_lept1_charge,"f_lept1_charge/F");
   TBranch *b_lept1_pfx= newtree->Branch("f_lept1_pfx", &f_lept1_pfx,"f_lept1_pfx/F");
   TBranch *b_lept1_sip= newtree->Branch("f_lept1_sip", &f_lept1_sip,"f_lept1_sip/F");
   TBranch *b_lept2_pt= newtree->Branch("f_lept2_pt", &f_lept2_pt,"f_lept2_pt/F");
   TBranch *b_lept2_eta= newtree->Branch("f_lept2_eta", &f_lept2_eta,"f_lept2_eta/F");
   TBranch *b_lept2_phi= newtree->Branch("f_lept2_phi", &f_lept2_phi,"f_lept2_phi/F");
   TBranch *b_lept2_charge= newtree->Branch("f_lept2_charge", &f_lept2_charge,"f_lept2_charge/F");
   TBranch *b_lept2_pfx= newtree->Branch("f_lept2_pfx", &f_lept2_pfx,"f_lept2_pfx/F");
   TBranch *b_lept2_sip= newtree->Branch("f_lept2_sip", &f_lept2_sip,"f_lept2_sip/F");
   TBranch *b_lept3_pt= newtree->Branch("f_lept3_pt", &f_lept3_pt,"f_lept3_pt/F");
   TBranch *b_lept3_eta= newtree->Branch("f_lept3_eta", &f_lept3_eta,"f_lept3_eta/F");
   TBranch *b_lept3_phi= newtree->Branch("f_lept3_phi", &f_lept3_phi,"f_lept3_phi/F");
   TBranch *b_lept3_charge= newtree->Branch("f_lept3_charge", &f_lept3_charge,"f_lept3_charge/F");
   TBranch *b_lept3_pfx= newtree->Branch("f_lept3_pfx", &f_lept3_pfx,"f_lept3_pfx/F");
   TBranch *b_lept3_sip= newtree->Branch("f_lept3_sip", &f_lept3_sip,"f_lept3_sip/F");
   TBranch *b_lept4_pt= newtree->Branch("f_lept4_pt", &f_lept4_pt,"f_lept4_pt/F");
   TBranch *b_lept4_eta= newtree->Branch("f_lept4_eta", &f_lept4_eta,"f_lept4_eta/F");
   TBranch *b_lept4_phi= newtree->Branch("f_lept4_phi", &f_lept4_phi,"f_lept4_phi/F");
   TBranch *b_lept4_charge= newtree->Branch("f_lept4_charge", &f_lept4_charge,"f_lept4_charge/F");
   TBranch *b_lept4_pfx= newtree->Branch("f_lept4_pfx", &f_lept4_pfx,"f_lept4_pfx/F");
   TBranch *b_lept4_sip= newtree->Branch("f_lept4_sip", &f_lept4_sip,"f_lept4_sip/F");
   TBranch *b_iso_max= newtree->Branch("f_iso_max", &f_iso_max,"f_iso_max/F");
   TBranch *b_sip_max= newtree->Branch("f_sip_max", &f_sip_max,"f_sip_max/F");
   TBranch *b_Z1mass= newtree->Branch("f_Z1mass", &f_Z1mass,"f_Z1mass/F");
   TBranch *b_Z2mass= newtree->Branch("f_Z2mass", &f_Z2mass,"f_Z2mass/F");
   TBranch *b_angle_costhetastar= newtree->Branch("f_angle_costhetastar", &f_angle_costhetastar,"f_angle_costhetastar/F");
   TBranch *b_angle_costheta1= newtree->Branch("f_angle_costheta1", &f_angle_costheta1,"f_angle_costheta1/F");
   TBranch *b_angle_costheta2= newtree->Branch("f_angle_costheta2", &f_angle_costheta2,"f_angle_costheta2/F");
   TBranch *b_angle_phi= newtree->Branch("f_angle_phi", &f_angle_phi,"f_angle_phi/F");
   TBranch *b_angle_phistar1= newtree->Branch("f_angle_phistar1", &f_angle_phistar1,"f_angle_phistar1/F");
   TBranch *b_pt4l= newtree->Branch("f_pt4l", &f_pt4l,"f_pt4l/F");
   TBranch *b_eta4l= newtree->Branch("f_eta4l", &f_eta4l,"f_eta4l/F");
   TBranch *b_mass4l= newtree->Branch("f_mass4l", &f_mass4l,"f_mass4l/F");
   TBranch *b_mass4lErr= newtree->Branch("f_mass4lErr", &f_mass4lErr,"f_mass4lErr/F");
   TBranch *b_njets_pass= newtree->Branch("f_njets_pass", &f_njets_pass,"f_njets_pass/F");
   TBranch *b_deltajj= newtree->Branch("f_deltajj", &f_deltajj,"f_deltajj/F");
   TBranch *b_massjj= newtree->Branch("f_massjj", &f_massjj,"f_massjj/F");
   TBranch *b_D_jet= newtree->Branch("f_D_jet", &f_D_jet,"f_D_jet/F");
   TBranch *b_jet1_pt= newtree->Branch("f_jet1_pt", &f_jet1_pt,"f_jet1_pt/F");
   TBranch *b_jet1_eta= newtree->Branch("f_jet1_eta", &f_jet1_eta,"f_jet1_eta/F");
   TBranch *b_jet1_phi= newtree->Branch("f_jet1_phi", &f_jet1_phi,"f_jet1_phi/F");
   TBranch *b_jet1_e= newtree->Branch("f_jet1_e", &f_jet1_e,"f_jet1_e/F");
   TBranch *b_jet2_pt= newtree->Branch("f_jet2_pt", &f_jet2_pt,"f_jet2_pt/F");
   TBranch *b_jet2_eta= newtree->Branch("f_jet2_eta", &f_jet2_eta,"f_jet2_eta/F");
   TBranch *b_jet2_phi= newtree->Branch("f_jet2_phi", &f_jet2_phi,"f_jet2_phi/F");
   TBranch *b_jet2_e= newtree->Branch("f_jet2_e", &f_jet2_e,"f_jet2_e/F");
   TBranch *b_D_bkg_kin= newtree->Branch("f_D_bkg_kin", &f_D_bkg_kin,"f_D_bkg_kin/F");
   TBranch *b_D_bkg= newtree->Branch("f_D_bkg", &f_D_bkg,"f_D_bkg/F");
   TBranch *b_D_gg= newtree->Branch("f_D_gg", &f_D_gg,"f_D_gg/F");
   TBranch *b_D_g4= newtree->Branch("f_D_g4", &f_D_g4,"f_D_g4/F");
   TBranch *b_Djet_VAJHU= newtree->Branch("f_Djet_VAJHU", &f_Djet_VAJHU,"f_Djet_VAJHU/F");
   TBranch *b_pfmet= newtree->Branch("f_pfmet", &f_pfmet,"f_pfmet/F");
    
   float newweight=1.;
   
   // New tree with clone of events passing the final selection
   TFile *skimfile = new TFile(eventsOUT,"recreate");
   // Clone tree for final events
   TTree *finaltree = fChain->CloneTree(0);

   // loop on entries
   
   Long64_t nentries = fChain->GetEntries();

   cout << "\n****************************"  <<endl;
   cout << "Analyzing " << nentries << " entries"  <<endl;     

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      if(jentry%1 == 5000) cout << "Analyzing entry: " << jentry << endl;
      

      if( RECO_NMU > 100 ) RECO_NMU = 100;
      if( RECO_NELE > 100 ) RECO_NELE = 100;
      if( RECO_NPFPHOT > 20 ) RECO_NPFPHOT = 20;
      
      bool debug=true;  //debug flag  -- default false

      newweight=weight;
      cout << "Starting weight= " << newweight << endl;

      // pileup reweighting 2015
      hPUvertices->Fill(num_PU_vertices,weight);

      double pu_weight=1.;
      if (MC_type == "Fall15"){
	Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
	cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;	
	pu_weight=double(puweight->GetBinContent(binx));
	
      }      
       
      hPUvertices_ReWeighted->Fill(num_PU_vertices,weight*pu_weight);
      cout << "Pileup interations and weight is= " << num_PU_vertices << " " << " and weight= " << pu_weight << endl;  
      
      // Changing the weight for pileup
      newweight=weight*pu_weight;
      cout << "Starting weight + pileup = " << newweight << endl;
           
      
      float pFill[11];for(int pf=0;pf<11;pf++)pFill[11]=-999.;

      // ** Step 0:
      // simply number of entries...
      if( debug ) cout << "\n** Step 0: \nAnalyzing entry: " << jentry << " Run: " << Run << " Event: " << Event << " LumiSection: " << LumiSection << endl ;
      ++N_0 ;  // fill counter
      N_0_w=N_0_w+newweight;
      
     
      ++N_1 ;  // fill counter
      N_1_w=N_1_w+newweight;
      
      
      
      // Loose lepton identification
      
      int N_loose_mu = 0;
      int N_loose_e = 0;
      double max_Iso_loose_mu = -1 ;
      double max_Sip_loose_mu = -1 ;
      double max_Ip_loose_mu = -1 ;
      double max_Iso_loose_e = -1 ;
      double max_Sip_loose_e = -1 ;
      double max_Ip_loose_e = -1 ;
      
      int* arraysize_mu = new int[1];
      arraysize_mu[0] = RECO_NMU;
      int iL_loose_mu[arraysize_mu[0]];
      delete [] arraysize_mu;

      for( int i = 0; i < RECO_NMU; ++i ){
	iL_loose_mu[i]=-999.;
      }
 
      for( int i = 0; i < RECO_NMU; ++i ){

        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
			 << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
			 << "\nRECOMU_PT[i] " << RECOMU_PT[i]
			 << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
			 << "\nfabs( RECOMU_mubestrkDxy[i] ) " << fabs( RECOMU_mubesttrkDxy[i] )
			 << "\nfabs( RECOMU_mubesttrkDz[i] ) " << fabs( RECOMU_mubesttrkDz[i] )
			 << endl ;
       	
 	if( ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0) )
	    && RECOMU_mubesttrkType[i]!=2
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.
	    ){ 
	  iL_loose_mu[N_loose_mu]=i;
	  ++N_loose_mu ;
	  if( RECOMU_PFX_dB[i] > max_Iso_loose_mu ) max_Iso_loose_mu = RECOMU_PFX_dB[i] ;
	  if( fabs( RECOMU_SIP[i] ) > max_Sip_loose_mu ) max_Sip_loose_mu = fabs( RECOMU_SIP[i] ) ;
	  if( fabs( RECOMU_IP[i] ) > max_Ip_loose_mu ) max_Ip_loose_mu = fabs( RECOMU_IP[i] ) ;
	}
	
      } // end loop on muons

      
      int* arraysize_e = new int[1];
      arraysize_e[0] = RECO_NELE;
      int iL_loose_e[arraysize_e[0]];
      delete [] arraysize_e;

      for( int i = 0; i < RECO_NELE; ++i ){
	iL_loose_e[i]=-999.;
      }

      for( int i = 0; i < RECO_NELE; ++i ){
	
        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOELE_PT[i] " << RECOELE_PT[i]
			 << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
			 << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
			 << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
			 << endl ;
       	
 	if( RECOELE_PT[i] > 7. 
	    && fabs(RECOELE_ETA[i]) < 2.5 
	    // && RECOELE_gsftrack_expected_inner_hits[i]<=1  not used anymore
	    && fabs(RECOELE_gsftrack_dxy[i]) < .5 
	    && fabs(RECOELE_gsftrack_dz[i]) < 1. 
	    ) {	  
	  iL_loose_e[N_loose_e]=i;
	  ++N_loose_e ;
	  if( RECOELE_PFX_rho[i] > max_Iso_loose_e ) max_Iso_loose_e = RECOELE_PFX_rho[i] ;
	  if( fabs( RECOELE_SIP[i] ) > max_Sip_loose_e ) max_Sip_loose_e = fabs( RECOELE_SIP[i] ) ;
	  if( fabs( RECOELE_IP[i] ) > max_Ip_loose_e ) max_Ip_loose_e = fabs( RECOELE_IP[i] ) ;
	}
	
      }// end loop on electrons
      
      hN_loose_mu->Fill( N_loose_mu,newweight );
      hN_loose_e->Fill( N_loose_e,newweight );
      hIso_loose_mu->Fill( max_Iso_loose_mu,newweight);
      hSip_loose_mu->Fill( max_Sip_loose_mu,newweight );
      hIp_loose_mu->Fill( max_Ip_loose_mu,newweight );
      hIso_loose_e->Fill( max_Iso_loose_e,newweight );
      hSip_loose_e->Fill( max_Sip_loose_e,newweight );
      hIp_loose_e->Fill( max_Ip_loose_e,newweight );
      
      // Electron Cross Cleaning  -- eles separated from muons (deltaR > 0.05)
      
      for(int e = 0; e < RECO_NELE; ++e)
      	for(int mu = 0; mu < RECO_NMU; ++mu){
	  
	  if( RECOMU_isPFMu[mu] 
	      && (RECOMU_isGlobalMu[mu] || (RECOMU_isTrackerMu[mu] && RECOMU_numberOfMatches[mu]>0))
	      && RECOMU_mubesttrkType[mu]!=2
	      && RECOMU_PT[mu] > 5. 
	      && fabs(RECOMU_ETA[mu]) < 2.4 
	      && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. 
	      && fabs(RECOMU_SIP[mu])<4. // TightID + SIP cut
	      );
	  else continue;
	  
	  double deltaR = sqrt( pow( DELTAPHI( RECOMU_PHI[mu] , RECOELE_PHI[e] ),2) + pow(RECOMU_ETA[mu] - RECOELE_ETA[e],2) );
	  
	  if( deltaR <= 0.05 ){
	    
	    if( debug )cout << "Electrom not passing the cross cleaning" << endl;
	    
	    RECOELE_PT[e]  = -0.01;
	    RECOELE_ETA[e] = -99.;
	    RECOELE_PHI[e] = -99.;
	    RECOELE_SIP[e] = -99.;
	  }
	}
      
                  
      // Lepton identification -- no iso
      
      int iL[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1};
      
      int N_good = 0 ;
      
      for( int i = 0; i < RECO_NMU; ++i ){
	
        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOMU_isPFMu[i] " << int(RECOMU_isPFMu[i])
			 << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
			 << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
			 << "\nRECOMU_PT[i] " << RECOMU_PT[i]
			 << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
			 << "\nRECOMU_PFX_dB[i] " << RECOMU_PFX_dB[i]
			 << "\nfabs( RECOMU_SIP[i] ) " << fabs( RECOMU_SIP[i] )
			 << "\nfabs( RECOMU_mubesttrkDxy[i] ) " << fabs( RECOMU_mubesttrkDxy[i] )
			 << "\nfabs( RECOMU_mubesttrkDz[i] ) " << fabs( RECOMU_mubesttrkDz[i] )
			 << endl ;
	
       	// Tight muons
 	if( RECOMU_isPFMu[i] 
	    && ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0))
	    && RECOMU_mubesttrkType[i]!=2	 
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.
	    ){
	  
	  iL[ N_good ] = i ;
	  ++N_good ;	  
	}
      } // end loop on muons
      
      hN_good_mu->Fill( N_good,newweight );
      
      if( debug ) cout << "\nLeptons' indeces: "
		       << "\niL[0]: " << iL[0]
		       << "\niL[1]: " << iL[1]
		       << "\niL[2]: " << iL[2]
		       << "\niL[3]: " << iL[3]
		       << "\niL[4]: " << iL[4]
		       << "\niL[5]: " << iL[5]
		       << "\niL[6]: " << iL[6]
		       << "\niL[7]: " << iL[7]
		       << "\nNumber of good muons: " << N_good
		       << endl ;
      
      ++N_2 ;  // fill counter
      N_2_w=N_2_w+newweight;

    
      
      /// *** FSR
      // Photon identification & cleaning
      // ele identification is also needed
      // Effective AREA
      bool tag_2011=false;
      if (DATA_type=="2010" || DATA_type=="2011" || MC_type=="Fall11"){
        tag_2011=true;
      } 
      //FSR photon identifications, will be used with MELA later
      int FSR_Z1_photid=-1;
      int FSR_Z2_photid=-1;
      int FSR_Z1_lepid=-1;
      int FSR_Z2_lepid=-1;

      //electrons:
      
      int Ne_good = 0 ;
      int iLe[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1}; //electrons

      for( int i = 0; i < RECO_NELE; ++i ){

        if( debug ) cout << "\n Electron i="<< i <<" properties: "
      		  << "\nRECOELE_PT[i] " << RECOELE_PT[i]
      		  << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
		  << "\nfabs(RECOELE_scl_Eta[i]) " << fabs(RECOELE_scl_Eta[i])
      		  << "\nRECOELE_PFX_rho[i] " << RECOELE_PFX_rho[i]
      		  << "\nfabs( RECOELE_SIP[i] ) " << fabs( RECOELE_SIP[i] )
      		  << "\nRECOELE_mvaNonTrigV0[i] " << RECOELE_mvaNonTrigV0[i]
		  << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
      		  << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
		  << endl ;
       	
 	if( RECOELE_PT[i] > 7. && fabs(RECOELE_ETA[i]) < 2.5 );
	  // && RECOELE_gsftrack_expected_inner_hits[i]<=1 ) /* ok */ ;
	else continue ;
	
	bool BDT_ok = 0; // Fall15 with CMSSW_7_6_x
	if( RECOELE_PT[i] > 7. &&  RECOELE_PT[i] <= 10. ){
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.265 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) < 1.479 )
						 && RECOELE_mvaNonTrigV0[i] > -0.556 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) >= 1.479 && RECOELE_mvaNonTrigV0[i] > -0.551 ) BDT_ok = 1 ;
	}
	else { 
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.072 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
						 && RECOELE_mvaNonTrigV0[i] > -0.286 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > -0.267 ) BDT_ok = 1 ;
	}	
	if( !BDT_ok ) continue ;
	
	if( fabs(RECOELE_gsftrack_dxy[i]) < .5 
	 && fabs(RECOELE_gsftrack_dz[i])  < 1. ) /* ok */ ;
	else continue ; 
		
	iLe[ Ne_good ] = i ;
	++Ne_good ;

      }// end loop on electrons

      hN_good_ele->Fill( Ne_good,newweight );
      hN_good_lep->Fill( N_good + Ne_good,newweight );


      if( debug ) cout << "\n Electrons' indeces: "
 		  << "\niLe[0]: " << iLe[0]
  		  << "\niLe[1]: " << iLe[1]
 		  << "\niLe[2]: " << iLe[2]
		  << "\niLe[3]: " << iLe[3]
 		  << "\niLe[4]: " << iLe[4]
  		  << "\niLe[5]: " << iLe[5]
 		  << "\niLe[6]: " << iLe[6]
		  << "\niLe[7]: " << iLe[7]
		  << "\nNumber of good eles: " << Ne_good
		  << endl ;
            
      // Define a new isolation array to allocate the contribution of photons
      // float RECOMU_PFX_dB_new[100],RECOELE_PFX_rho_new[100];
      for (int i=0;i<100;i++){
	RECOMU_PFX_dB_new[i]=RECOMU_PFX_dB[i];
	RECOELE_PFX_rho_new[i]=RECOELE_PFX_rho[i];	
      }
      //


      // photon definition & cleaning:
      int iLp[30];
      	for( int i = 0 ; i < 30 ; ++i )iLp[i] = -1;
      
      int Nphotons = 0;

      for( int i = 0; i < RECO_NPFPHOT; ++i ){
	
        if( debug ) cout << "\n Photon i="<< i <<" properties: "
			 << "\n RECOPFPHOT_PT[i] " << RECOPFPHOT_PT[i]
			 << "\n fabs(RECOPFPHOT_ETA[i]) " << fabs(RECOPFPHOT_ETA[i])
			 << "\n RECOPFPHOT_PHI[i] " << RECOPFPHOT_PHI[i]
			 << "\n RECOPFPHOT_PFX_rho[i] " << RECOPFPHOT_PFX_rho[i]
			 << endl ;
	
	if ( RECOPFPHOT_PT[i] > 2. && fabs(RECOPFPHOT_ETA[i]) < 2.4 && RECOPFPHOT_PFX_rho[i]<1.8) {
	  
	  bool is_clean = 1;
	  
	  // cleaning
	  for(int e = 0; e < N_loose_e; ++e){
  	    if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;  // loose ID + SIP cut	    
	    double deltaPhi = DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ) ;
	    double deltaEta = fabs( RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]] );
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	    cout << "debug: " << RECOELE_PT[iL_loose_e[e]] << " " << deltaPhi << " " << deltaEta << " " << deltaR << endl;
	    if( ( fabs(deltaPhi) < 2 && fabs(deltaEta) < 0.05 ) || deltaR <= 0.15 ){		  
	      if( debug )cout << "Photon not passing the electron cleaning" << endl;	
	      is_clean = 0;	  
	      
	    }
	  } // end loop on eles		             	
	  
	  if( !is_clean ) continue ;
	  
	  
	  iLp[ Nphotons ] = i ;
	  ++Nphotons ;
	  
	}
      }// end loop on photons
      
      hN_good_phot->Fill( Nphotons,newweight );
      
      if( debug ) cout << "Photons' indeces: "
		       << "\niLp[0]: " << iLp[0]
		       << "\niLp[1]: " << iLp[1]
		       << "\niLp[2]: " << iLp[2]
		       << "\niLp[3]: " << iLp[3]
		       << "\niLp[4]: " << iLp[4]
		       << "\niLp[5]: " << iLp[5]
		       << "\niLp[6]: " << iLp[6]
		       << "\niLp[7]: " << iLp[7]
		       << "\nNumber of good photons: " << Nphotons
		       << endl ;
      
      
      // assign to each photon the closest lepton
      int iLp_l[30];
      for( int i = 0 ; i < 30 ; ++i )iLp_l[i] = -1;
      int iLp_tagEM[30];
      for( int i = 0 ; i < 30 ; ++i )iLp_tagEM[i] = -1;  // tag  0: mu  1: ele
      
      float RECOPFPHOT_DR[30];
      for( int i = 0 ; i < 30 ; ++i ) RECOPFPHOT_DR[i] = -999; 

      for( int i = 0; i < Nphotons; ++i ){
	
	double min_deltaR = 1000;
	int  l_min_deltaR = -1;
	int  tag_min_deltaR = -1;   // 0: mu  1: ele
	
	for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	  if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	  cout << "DeltaR= " << deltaR << " " <<  deltaR/pow(RECOPFPHOT_PT[iLp[i]],2) << endl;
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to a muon with pT= " << RECOMU_PT[iL_loose_mu[l]]<< endl;
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 0;
	  }
	  
	}//end loop on muons  
	
	for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
	  if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[l]],2) );
	  cout << "DeltaR= " << deltaR << " " <<  deltaR/pow(RECOPFPHOT_PT[iLp[i]],2) << endl;
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to an electron with pT= " << RECOELE_PT[iL_loose_e[l]]<< endl;
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 1;
	  }
	  
	}//end loop on electrons  

	
	if( min_deltaR < 0.5 ){
	  if (tag_min_deltaR==0) iLp_l[ i ] = iL_loose_mu[l_min_deltaR];
	  if (tag_min_deltaR==1) iLp_l[ i ] = iL_loose_e[l_min_deltaR];
	  iLp_tagEM[ i ] = tag_min_deltaR;
	  RECOPFPHOT_DR[iLp[i]]=min_deltaR; 	
	}

      }
      
      if( debug ) cout << "Indeces of loose leptons associated to photons: "
		       << "\niLp_l[0]: " << iLp_l[0]
		       << "\niLp_l[1]: " << iLp_l[1]
		       << "\niLp_l[2]: " << iLp_l[2]
		       << "\niLp_l[3]: " << iLp_l[3]
		       << "\niLp_l[4]: " << iLp_l[4]
		       << "\niLp_l[5]: " << iLp_l[5]
		       << "\niLp_l[6]: " << iLp_l[6]
		       << "\niLp_l[7]: " << iLp_l[7]
		       << endl ;
      
      if( debug ) cout << "Tag of leptons associated to photons: (0: mu , 1:ele)"
 		  << "\niLp_tagEM[0]: " << iLp_tagEM[0]
  		  << "\niLp_tagEM[1]: " << iLp_tagEM[1]
 		  << "\niLp_tagEM[2]: " << iLp_tagEM[2]
		  << "\niLp_tagEM[3]: " << iLp_tagEM[3]
 		  << "\niLp_tagEM[4]: " << iLp_tagEM[4]
  		  << "\niLp_tagEM[5]: " << iLp_tagEM[5]
 		  << "\niLp_tagEM[6]: " << iLp_tagEM[6]
		  << "\niLp_tagEM[7]: " << iLp_tagEM[7]
		  << endl ;


      // Multiple photons associated to the same lepton: the lowest-ΔR(γ,l)/ETγ2 has to be selected.
      double min_deltaR_ET2=1000;
      int p_min_deltaR_ET2=-1;

      for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	    cout <<  "index muon" << iL_loose_mu[l] << endl;
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    if (deltaR_ET2<min_deltaR_ET2) {
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
	    }
	  }
	}
	
	if (p_min_deltaR_ET2!=-1){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	      if (p!=p_min_deltaR_ET2){
		iLp_l[ p ] = -1;
		iLp_tagEM[ p ] = -1;
	      }
	    }
	  }
	}
	
      }

   
      
      //
      min_deltaR_ET2=1000;
      p_min_deltaR_ET2=-1;
      
      for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
	if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	    cout <<  "index electron" << iL_loose_e[l] << endl;
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOELE_ETA[iL_loose_e[l]],2));
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    cout << " deltaR_ET2= " << deltaR_ET2 <<endl;
	    if (deltaR_ET2<min_deltaR_ET2){
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
	      cout << " p_min_deltaR_ET2= " << p_min_deltaR_ET2 <<endl;
	    }
	  }	  
	}
	
	if (p_min_deltaR_ET2!=-1){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	      if (p!=p_min_deltaR_ET2){
		iLp_l[ p ] = -1;
		iLp_tagEM[ p ] = -1;
	      }
	    }
	  }	  
	}
	
      }	
     
      
      if( debug ) cout << "Indeces of loose leptons associated to the photon with lowest DeltaR/ET2: "
		       << "\niLp_l[0]: " << iLp_l[0]
		       << "\niLp_l[1]: " << iLp_l[1]
		       << "\niLp_l[2]: " << iLp_l[2]
		       << "\niLp_l[3]: " << iLp_l[3]
		       << "\niLp_l[4]: " << iLp_l[4]
		       << "\niLp_l[5]: " << iLp_l[5]
		       << "\niLp_l[6]: " << iLp_l[6]
		       << "\niLp_l[7]: " << iLp_l[7]
		       << endl ;
      
      if( debug ) cout << "Tag of leptons associated to the photon with lowest DetaR/ET2: (0: mu , 1:ele)"
		       << "\niLp_tagEM[0]: " << iLp_tagEM[0]
		       << "\niLp_tagEM[1]: " << iLp_tagEM[1]
		       << "\niLp_tagEM[2]: " << iLp_tagEM[2]
		       << "\niLp_tagEM[3]: " << iLp_tagEM[3]
		       << "\niLp_tagEM[4]: " << iLp_tagEM[4]
		       << "\niLp_tagEM[5]: " << iLp_tagEM[5]
		       << "\niLp_tagEM[6]: " << iLp_tagEM[6]
		       << "\niLp_tagEM[7]: " << iLp_tagEM[7]
		       << endl ;

      
      for(int i=0.;i<Nphotons;i++) {
	if (iLp_l[i]!=-1 && iLp_tagEM[i]==0) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[i]] << endl;
	if (iLp_l[i]!=-1 && iLp_tagEM[i]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to an electron with pT= " << RECOELE_PT[iLp_l[i]] << endl;
      };

       // Exclude that photon from the isolation cone all leptons in the event passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons
      
      cout << "Rho for electron pileup isolation correction is= " << RHO_ele << endl;
      double EffectiveArea=-9999.;

	    
      for(int i=0.;i<Nphotons;i++) {
	if (iLp_l[i]==-1) continue;
	
	for(int e = 0; e < N_loose_e; ++e){
	  //if(!( iLp_l[i] == iL_loose_e[e] && iLp_tagEM[i] == 1 ) ) continue;
	  if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;
	  //double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[e]],2) );
	  cout << "deltaR for photon subtraction= " << deltaR << endl;
	  if( deltaR<=0.3 && (RECOELE_scl_Eta[iL_loose_e[e]]< 1.479 || deltaR>0.08) ){ // 0.3 in 76x              
	    if( debug )cout << "Subtracting the photon isolation from the electron isolation value " << endl;
	    
	    EffectiveArea=EAele(iL_loose_e[e],tag_2011);
	    RECOELE_PFX_rho_new[iL_loose_e[e]]=
              (RECOELE_PFchHad[iL_loose_e[e]]+
               max(0.,RECOELE_PFneuHad[iL_loose_e[e]]+
                   (RECOELE_PFphoton[iL_loose_e[e]]-RECOPFPHOT_PT[iLp[i]] )-
                   max(RHO_ele,0.0)*(EffectiveArea)))/RECOELE_PT[iL_loose_e[e]];	    
	  }
	} // end loop on ele
	
	for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	  //if(!( iLp_l[i] == iL_loose_mu[l] && iLp_tagEM[i] == 0 ) ) continue;
          if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;
          double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	  cout << "deltaR for photon subtraction= " << deltaR << endl;
	  if( deltaR<=0.3 && deltaR>0.01){ // 0.3 is the isolation cone for muons in 76x
	    if( debug )cout << "Subtracting the photon isolation from the muon isolation value " << endl;
	    
	    RECOMU_PFX_dB_new[iL_loose_mu[l]]=
              (RECOMU_PFchHad[iL_loose_mu[l]]+
               max(0.,RECOMU_PFneuHad[iL_loose_mu[l]]+
                   (RECOMU_PFphoton[iL_loose_mu[l]]-RECOPFPHOT_PT[iLp[i]] )-
                   0.5*RECOMU_PFPUchAllPart[iL_loose_mu[l]]))/RECOMU_PT[iL_loose_mu[l]];
	    
	  }
	} // end loop on mu
	
	
	
      }	
      
      
      // *** end FSR
      
      
      // **** Step 3:
      struct candidateZ {
	float massvalue;
	int ilept1;
	float pt1;
	float isol1;
	bool ilept1_FSR;
	float eta1;
	float phi1;
	int charge1;
	int charge2;
	int ilept2;
	float pt2;
	float isol2;
	bool ilept2_FSR;
	float eta2;
	float phi2;
	float pxZ;
	float pyZ;
	float pzZ;
	float EZ;
	bool withFSR;
	float ptFSR;
	int tag;
      };

      vector<candidateZ> Zcandvector;
      Zcandvector.clear();
      vector<candidateZ> Zcandisolvector;
      Zcandisolvector.clear();

      // a) pair #1: mass closest to Z1
      // b) mLL in ] 40,120 [
      if( debug ) cout  << "\nStep 3: Number of good leptons: " << N_good+Ne_good << endl;
 
      if( N_good + Ne_good < 2 ) continue ; 	

      int Zxx_tag = 0;    // 1: Zmumu  ,  2: Zee

      int i1 = -1; //index of the first lepton (from Z1)
      int j1 = -1; //index of the second lepton (from Z1)
      int pi1 = -1; 
      int pj1 = -1;
      
      bool has_FSR_Z1 = 0;
      TLorentzVector Lepton1,Lepton2,DiLepton,LeptonCorrection;

      for(int i = 0; i < N_good; ++i){
        for(int j = i + 1; j < N_good; ++j){
	  if (fabs(RECOMU_SIP[iL[i]])>=4.) continue; // SIP cut
	  if (fabs(RECOMU_SIP[iL[j]])>=4.) continue;
	  if (fabs(RECOMU_PFX_dB_new[iL[i]])>=0.35) continue; // Isolation
	  if (fabs(RECOMU_PFX_dB_new[iL[j]])>=0.35) continue;
	  
	  if(RECOMU_CHARGE[ iL[j] ] == RECOMU_CHARGE[ iL[i] ]) continue; // opposite charge

	  cout << "\n Pairing muons with pT= " << RECOMU_PT[ iL[i] ] << " and " <<  RECOMU_PT[ iL[j] ] << endl;
		  
	  // evaluate the mass &
	  double pxZ, pyZ, pzZ;
	  double EZ;
	  double massZ;
	  double massZ_noFSR = 0;
	  
	  int tempphotid=-1;
	  int templepid=-1;

	  float pTphot=-999.;
	  Lepton1.SetPtEtaPhiM(RECOMU_PT[iL[i]], RECOMU_ETA[iL[i]], RECOMU_PHI[iL[i]], 0.105);
	  Lepton2.SetPtEtaPhiM(RECOMU_PT[iL[j]], RECOMU_ETA[iL[j]], RECOMU_PHI[iL[j]], 0.105);
	  DiLepton=Lepton1+Lepton2;	  
	  massZ = DiLepton.M();	  
	  massZ_noFSR = massZ;
	  if (debug) cout << "Mass Z= " << massZ << endl;
	  pxZ=DiLepton.Px();
	  pyZ=DiLepton.Py();
	  pzZ=DiLepton.Pz();
	  EZ=DiLepton.E();	  
	 
	  Zxx_tag=1;	 

	  // ** Association of FSR to Z
	  if( debug ) cout  << "Step Z+FSR  " << endl;
	  
	  bool has_FSR_Z = 0;
	  int N_FSR_Z = 0;
	  double max_pt_FSR_Z = -1.;
	  int pi = -1; 
	  int pj = -1;
	  
	  
	  for( int p = 0; p < Nphotons; ++p ){
	    
	    if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )  {  // exists a photon associated to a lepton mu
	      cout << "Attaching a photon to muon and then to the Z" << endl;
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton1=Lepton1+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();	      	    

	      //cout << mllp << " " << Zmass << " " << massZ << endl;
	      pi = p; 
	      has_FSR_Z = 1;
	      ++N_FSR_Z;
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      cout << "Mass Z with FSR= "<< massZ << endl;

	    }
	    if( iLp_l[ p ] == iL[j] && iLp_tagEM[ p ] == 0 )  { 
	      cout << "Attaching a photon to muon and then to the Z" << endl;
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton2=Lepton2+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();

	      //cout << mllp << " " << Zmass << " " << massZ << endl;
	      pj = p;
	      has_FSR_Z = 1;
	      ++N_FSR_Z; 
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      cout << "Mass Z with FSR= "<< massZ << endl;

	    }
	  } // end loop on FSR photons

	 
	  
	  
	  if( debug && has_FSR_Z) {
	    cout  << " Z has FSR! " << endl;
	    cout  << "  N_FSR_Z " << N_FSR_Z << endl;
	    cout  << "  max_pt of photon FSR_Z " << max_pt_FSR_Z << endl;
	    if( pi > -1 ) cout  << "  pi " << pi << " --> index photon: " << iLp[pi] << " associated lepton: " << iLp_l[pi] << " (= "<< iL[i]<<" ? )  tag: " << iLp_tagEM[pi] << endl;
	    if( pj > -1 ) cout  << "  pj " << pj << " --> index photon: " << iLp[pj] << " associated lepton: " << iLp_l[pj] << " (= "<< iL[j]<<" ? )  tag: " << iLp_tagEM[pj] << endl;
	  }
	  else {
	    cout << "No FSR photon attached" << endl;
	  }
	  
	  
	  if( has_FSR_Z ){ // if Z has FSR
	    
	    ++N_3_FSR; // fill the counter
	    N_3_FSR_w=N_3_FSR_w+newweight;
	    	    
	    // do not recompute isolation here
	    if( debug ) cout  << "Z Isolation (not corrected for photon): "
                              << "\n RECOMU_PFX_dB[ iL[i] ] " << RECOMU_PFX_dB[ iL[i] ]
                              << "\n RECOMU_PFX_dB[ iL[j] ] " << RECOMU_PFX_dB[ iL[j] ]
                              << endl;
	    	   
	    if( pi != -1 ){
	      //double deltaR_i = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pi]] , RECOMU_PHI[iL[i]] ),2) + pow(RECOPFPHOT_ETA[iLp[pi]] - RECOMU_ETA[iL[i]],2) );
              //if( deltaR_i < 0.4 && deltaR_i > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pi]];	      
	    }
	    else if( pj != -1 ){	      
              //double deltaR_j = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pj]] , RECOMU_PHI[iL[j]] ),2) + pow(RECOPFPHOT_ETA[iLp[pj]] - RECOMU_ETA[iL[j]],2) );
              //if( deltaR_j < 0.4 && deltaR_j > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pj]];	      
	    }

 
	  } // end if has FSR
	  else{
	    if( debug ) cout  << "Z Isolation: "  
			      << "\n RECOMU_PFX_dB_new[ iL[i] ] " << RECOMU_PFX_dB_new[ iL[i] ]
			      << "\n RECOMU_PFX_dB_new[ iL[j] ] " << RECOMU_PFX_dB_new[ iL[j] ]
			      << endl;	    
	  }
	  // ** end association of FSR to Z
	  
	  //if( massZ == 0 || i1 == -1 || j1 == -1) continue;
	  
	  cout << "2mu2e: " << Zxx_tag << endl; 
	  
	  cout << "Filling a struct for Z" << endl; 
	  candidateZ *Z = new candidateZ;
	  Z->massvalue=massZ;
	  Z->ilept1=iL[i];
	  Z->ilept2=iL[j];
	  Z->pt1=RECOMU_PT[iL[i]];
	  Z->pt2=RECOMU_PT[iL[j]];
	  Z->eta1=RECOMU_ETA[iL[i]];
	  Z->eta2=RECOMU_ETA[iL[j]];
	  Z->phi1=RECOMU_PHI[iL[i]];
	  Z->phi2=RECOMU_PHI[iL[j]];
	  Z->charge1=RECOMU_CHARGE[iL[i]];
	  Z->charge2=RECOMU_CHARGE[iL[j]];
	  Z->isol1=RECOMU_PFX_dB_new[ iL[i] ];
	  Z->isol2=RECOMU_PFX_dB_new[ iL[j] ];
	  if( pi != -1 ) Z->ilept1_FSR=true;
	  if( pj != -1 ) Z->ilept2_FSR=true;
	  Z->pxZ=pxZ;
	  Z->pyZ=pyZ;
	  Z->pzZ=pzZ;
	  Z->EZ=EZ;
	  if( has_FSR_Z ) {
	    Z->withFSR=1;
	    Z->ptFSR=pTphot;	    
	  }	      
	  else {
	    Z->withFSR=0;
	    Z->ptFSR=0.;
	  }
	  Z->tag=Zxx_tag;
	  	 	  
	  Zcandvector.push_back(*Z);	  
	  
	}
      } // end loop on pairs

      // 2mu2e
      for(int i = 0; i < Ne_good; ++i){
        for(int j = i + 1; j < Ne_good; ++j){
	  if (fabs(RECOELE_SIP[iLe[i]])>=4.) continue; // SIP cut
	  if (fabs(RECOELE_SIP[iLe[j]])>=4.) continue;
	  if (fabs(RECOELE_PFX_rho_new[iLe[i]])>=0.35) continue; // Isolation cut
	  if (fabs(RECOELE_PFX_rho_new[iLe[j]])>=0.35) continue;
	  
	  if(RECOELE_CHARGE[ iLe[j] ] == RECOELE_CHARGE[ iLe[i] ]) continue; // opposite charge

	  cout << "\n Pairing electrons with pT= " << RECOELE_PT[ iLe[i] ] << " and " <<  RECOELE_PT[ iLe[j] ] << endl;
	  
	  // evaluate the mass &
	  double pxZ, pyZ, pzZ;
	  double EZ;
	  double massZ;
	  double massZ_noFSR = 0;
	  
	 
	  int tempphotid=-1;
	  int templepid=-1;
	  
	  float pTphot=-999.;
	  Lepton1.SetPtEtaPhiM(RECOELE_PT[iLe[i]], RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]], 0.000511);
	  Lepton2.SetPtEtaPhiM(RECOELE_PT[iLe[j]], RECOELE_ETA[iLe[j]], RECOELE_PHI[iLe[j]], 0.000511);
	  DiLepton=Lepton1+Lepton2;	  
	  massZ = DiLepton.M();	  
	  massZ_noFSR = massZ;
	  if (debug) cout << "Mass Z= " << massZ << endl;
	  pxZ=DiLepton.Px();
	  pyZ=DiLepton.Py();
	  pzZ=DiLepton.Pz();
	  EZ=DiLepton.E();
	  
	  Zxx_tag=2;
	
	  // ** Association of FSR to Z
	  if( debug ) cout  << "Step Z+FSR  " << endl;
	  
	  bool has_FSR_Z = 0;
	  int N_FSR_Z = 0;
	  double max_pt_FSR_Z = -1.;
	  int pi = -1; 
	  int pj = -1;
	  
	  
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {  // exit a photon associated to a lepton electron
	      
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton1=Lepton1+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();	      	      	   
	      
	      has_FSR_Z = 1; 
	      pi = p; 
	      ++N_FSR_Z;
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      cout << "Mass Z with FSR= "<< massZ << endl;
		      
	    }
	    
	    if( iLp_l[ p ] == iLe[j] && iLp_tagEM[ p ] == 1 )  { 
	      
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton2=Lepton2+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();
	      	    
	      pj = p;
	      has_FSR_Z = 1;
	      ++N_FSR_Z; 
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      cout << "Mass Z with FSR= "<< massZ << endl;

	    }
	  } // end loop on FSR photons
	  
	  
	  
	  //if( has_FSR_Z ) debug = 1;
	  
	  if( debug && has_FSR_Z) {
	    cout  << " Z has FSR! " << endl;
	    cout  << "  N_FSR_Z " << N_FSR_Z << endl;
	    cout  << "  max_pt of photon FSR_Z " << max_pt_FSR_Z << endl;
	    if( pi > -1 ) cout  << "  pi " << pi << " --> index photon: " << iLp[pi] << " associated lepton: " << iLp_l[pi] << " (= "<< iLe[i]<<" ? )  tag: " << iLp_tagEM[pi] << endl;
	    if( pj > -1 ) cout  << "  pj " << pj << " --> index photon: " << iLp[pj] << " associated lepton: " << iLp_l[pj] << " (= "<< iLe[j]<<" ? )  tag: " << iLp_tagEM[pj] << endl;
	  }
	  else {
	    cout << "No FSR photon attached" << endl;
	  }
	  
	  
	  if( has_FSR_Z ){ // if Z has FSR
	    
	    ++N_3_FSR; // fill the counter
	    N_3_FSR_w=N_3_FSR_w+newweight;	    
	    
	    // do not recompute isolation here
	    if( debug ) cout  << "Z Isolation (not corrected for photon): "
                              << "\n RECOELE_PFX_rho[ iLe[i] ] " << RECOELE_PFX_rho[ iLe[i] ]
                              << "\n RECOELE_PFX_rho[ iLe[j] ] " << RECOELE_PFX_rho[ iLe[j] ]
                              << endl;
	    
	    if( pi != -1 ){
	      //double deltaR_i = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pi]] , RECOELE_PHI[iLe[i]] ),2) + pow(RECOPFPHOT_ETA[iLp[pi]] - RECOELE_ETA[iLe[i]],2) );
              //if( deltaR_i < 0.4 && deltaR_i > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pi]];	      
	    }
	    else if( pj != -1 ){	      
              //double deltaR_j = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pj]] , RECOELE_PHI[iLe[j]] ),2) + pow(RECOPFPHOT_ETA[iLp[pj]] - RECOELE_ETA[iLe[j]],2) );
              //if( deltaR_j < 0.4 && deltaR_j > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pj]];	      
	    }
	    
	    
	  } // end if has FSR
	  else{
	    
	    if( debug ) cout  << "Z Isolation: "  
			      << "\n RECOELE_PFX_rho_new[ iLe[i] ] " << RECOELE_PFX_rho_new[ iLe[i] ]
			      << "\n RECOELE_PFX_rho_new[ iLe[j] ] " << RECOELE_PFX_rho_new[ iLe[j] ]
			      << endl;	    
	  }
	  // ** end association of FSR to Z
	  
	  //if( massZ == 0 || i1 == -1 || j1 == -1) continue;
	  
	  cout << "2e2mu: " << Zxx_tag << endl;
	   
	  cout << "Filling a struct for Z" << endl; 
	  candidateZ *Z = new candidateZ;
	  Z->massvalue=massZ;
	  Z->ilept1=iLe[i];
	  Z->ilept2=iLe[j];
	  Z->pt1=RECOELE_PT[iLe[i]];
	  Z->pt2=RECOELE_PT[iLe[j]];
	  Z->eta1=RECOELE_ETA[iLe[i]];
	  Z->eta2=RECOELE_ETA[iLe[j]];
	  Z->phi1=RECOELE_PHI[iLe[i]];
	  Z->phi2=RECOELE_PHI[iLe[j]];
	  Z->charge1=RECOELE_CHARGE[iLe[i]];
	  Z->charge2=RECOELE_CHARGE[iLe[j]];
	  Z->isol1=RECOELE_PFX_rho_new[ iLe[i] ];
	  Z->isol2=RECOELE_PFX_rho_new[ iLe[j] ];
	  if( pi != -1 ) Z->ilept1_FSR=true;
	  if( pj != -1 ) Z->ilept2_FSR=true;
	  Z->pxZ=pxZ;
	  Z->pyZ=pyZ;
	  Z->pzZ=pzZ;
	  Z->EZ=EZ;
	  if( has_FSR_Z ) {
	    Z->withFSR=1;
	    Z->ptFSR=pTphot;
	  }	      
	  else {
	    Z->withFSR=0;
	    Z->ptFSR=0.;
	  }
	  	 
	  Z->tag=Zxx_tag;
	  Zcandvector.push_back(*Z);	  
	  
	}
      } // end loop on couples
      
      
      if (Zcandvector.size()<1) {
	cout << "Less than one Z pairs with isolated leptons...exiting" << endl;
	continue; 
      }
      

      ++N_3a ;  // fill counter
      N_3a_w=N_3a_w+newweight;
      

      // Mass cut on Z
      vector<candidateZ> Zcandisolmassvector;
      Zcandisolmassvector.clear();

      for (int index=0; index<Zcandvector.size();index++){
	if (!(Zcandvector.at(index).massvalue > 12 && Zcandvector.at(index).massvalue < 240)) continue;
	cout << "Z passing the 12 < mll < 240 cut with mass= " << Zcandvector.at(index).massvalue<< endl;
	Zcandisolmassvector.push_back(Zcandvector.at(index));
      };
      
      if (Zcandisolmassvector.size()<1) {
	cout << "No Z passing the mass cut"<< endl;
	continue;
      }

      cout << "Number of Z passing the isolation and the 12 << mll < 240 cut is= " << Zcandisolmassvector.size() << endl;

      ++N_3b ;  // fill counter
      N_3b_w=N_3b_w+newweight;
      

      cout << "Starting weight + pileup + LineShape + efficiency= " << newweight << endl;

      hPFMET_3->Fill(RECO_PFMET,newweight);
      
      // **** Step 4:
       // a) 4 leptons
      // b) pair #2
      // c) highest pt
      // d) mZ2 in ] 4,120 [

      int issamesign = 0;

      int N_Z2_pairs = 0;

      int i2 = -1; //index of the first lepton (from Z1)
      int j2 = -1; //index of the second lepton (from Z1)
      int pi2 = -1; 
      int pj2 = -1; 
      
      /*
      double pxZ2 = 0;  //Z2 kinematics
      double pyZ2 = 0;
      double pzZ2 = 0;
      double ptZ2 = 0;
      double EZ2 = 0;
      double Y_Z2 = -9;
      double massZ2 = 0;
      double massZ2_noFSR = 0;
      double sum_ptZ2 = 0.;
      */
      bool has_FSR_Z2 = 0;
      

      // if( N_good < 4 ) continue ; 

      // ++N_4a ;  // fill counter
      // N_4a_w=N_4a_w+newweight;

//mark1
/*
      // ZZ object pairs
      array<int,2> isolgoodZs;
      vector<std::array<int, 2> > isolgoodZsv;
        
      
      // Ghost removal delta R > 0.02 
      vector<candidateZ> goodZ;

      for (int i=0;i<Zcandisolmassvector.size();i++){
	cout << "Checking mass i= " << Zcandisolmassvector.at(i).massvalue << endl;
	cout << "Ghost removal check 0: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(i).phi2 ),2) 
		  + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(i).eta2,2) ) << endl;
	if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(i).phi2 ),2) 
		  + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(i).eta2,2) ) <= 0.02) continue;
	
	
	for (int j=i+1;j<Zcandisolmassvector.size();j++){
	  if (Zcandisolmassvector.at(i).tag==Zcandisolmassvector.at(j).tag) continue;	// select only 2e2mu pairs 
	  cout << "Checking mass j= " << Zcandisolmassvector.at(j).massvalue << endl;
	  cout <<Zcandisolmassvector.at(i).pt1 << " " << Zcandisolmassvector.at(j).pt1 << " " << Zcandisolmassvector.at(j).pt2 << endl;
	  cout <<Zcandisolmassvector.at(i).pt2 << " " << Zcandisolmassvector.at(j).pt1 << " " << Zcandisolmassvector.at(j).pt2 << endl;
	  if (Zcandisolmassvector.at(i).pt1==Zcandisolmassvector.at(j).pt1 || Zcandisolmassvector.at(i).pt1==Zcandisolmassvector.at(j).pt2) continue;
	  if (Zcandisolmassvector.at(i).pt2==Zcandisolmassvector.at(j).pt1 || Zcandisolmassvector.at(i).pt2==Zcandisolmassvector.at(j).pt2) continue;
	 	  
	  cout << "Ghost removal check 1: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(j).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(j).eta1-Zcandisolmassvector.at(j).eta2,2) )<< endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(j).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(j).eta1-Zcandisolmassvector.at(j).eta2,2) ) <= 0.02) continue;
	  
	  cout << "Ghost removal check 2: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta1,2) )<< endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta1,2) ) <= 0.02) continue;

	  cout << "Ghost removal check 3: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta2,2) ) << endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta2,2) ) <= 0.02) continue;

	  cout << "Ghost removal check 4: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta1,2) ) << endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta1,2) ) <= 0.02) continue;

	  cout << "Ghost removal check 5: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta2,2) ) << endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta2,2) ) <= 0.02) continue;

	  cout << "There is a set of 4 leptons passing the ghost removal (deltaR > 0.02)" << endl;
	  goodZ.push_back(Zcandisolmassvector.at(i));
	  goodZ.push_back(Zcandisolmassvector.at(j));
	  cout << "Filling isolgoodZsv vector " << endl;
	  isolgoodZs={i,j};
	  isolgoodZsv.push_back(isolgoodZs);
	}
      }

      //cout << "Debug ZZ= " << Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ).massvalue << " " << Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1)).massvalue << endl;
      if (isolgoodZsv.size()==0) {
	cout << "No ZZ (2e2mu) combination passing the cuts  ...exiting " << endl;
	continue;
      }
 */ //mark2     
      
      // PT,20/10 for any di-lepton
      vector<candidateZ> firstpTcleanedgoodZ;    
      vector<float> leptonspTcleaned;
      
//      array<int,2> ileptonspTcleanedgoodZs;
//      vector<std::array<int, 2> > ileptonspTcleanedgoodZsv;
/*      
      // PT,20/10 for any di-lepton      
      for (int l=0;l<isolgoodZsv.size();l++){
	leptonspTcleaned.clear();	
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ).pt1);
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ).pt2);
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1) ).pt1);
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1) ).pt2);
	std::sort(leptonspTcleaned.rbegin(),leptonspTcleaned.rend());
	
	if (leptonspTcleaned.at(0)>20. && leptonspTcleaned.at(1)>10.) {
	  firstpTcleanedgoodZ.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ));
	  firstpTcleanedgoodZ.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1) ));
	  ileptonspTcleanedgoodZs={ (isolgoodZsv.at(l)).at(0), (isolgoodZsv.at(l)).at(1)};
	  ileptonspTcleanedgoodZsv.push_back(ileptonspTcleanedgoodZs);
	}	 	
      }
      
      cout << "Cleaned Good ZZ passing pT cuts are " << ileptonspTcleanedgoodZsv.size() << endl;            
      if (ileptonspTcleanedgoodZsv.size()==0) continue;
*/

      for (int l=0;l<Zcandisolmassvector.size();l++){
        leptonspTcleaned.clear();       
        leptonspTcleaned.push_back((Zcandisolmassvector.at(l)).pt1);
        leptonspTcleaned.push_back((Zcandisolmassvector.at(l)).pt1);
        std::sort(leptonspTcleaned.rbegin(),leptonspTcleaned.rend());

        if (leptonspTcleaned.at(0)>20. && leptonspTcleaned.at(1)>10.) {
          firstpTcleanedgoodZ.push_back(Zcandisolmassvector.at(l));
        }               
      }

      // **** Step 6:
      // QCD suppression: mll>4 GeV cut on all OS- any flavour pairs (4/4)
//mark1
      vector<candidateZ> pTcleanedgoodZ;
      pTcleanedgoodZ=firstpTcleanedgoodZ; 
/*     
      double min_mass_2L = 10000.;
      TLorentzVector Lepton1qcd,Lepton2qcd,DiLeptonQCD;
      vector<int> ileptonsmumu,ileptonsee,ileptonsemu;

      array<int,2> iqcdcleanedgoodZs;
      vector<std::array<int, 2> > iqcdcleanedgoodZsv;
      vector<candidateZ> pTcleanedgoodZ; 
      
      for (int l=0;l<ileptonspTcleanedgoodZsv.size();l++){
	cout << "checking masses for QCD suppression= " 
	     << Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).massvalue << " and " 
	     << Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).massvalue
	     << endl;
	min_mass_2L = 10000.;
	ileptonsmumu.clear();
	ileptonsee.clear();
	ileptonsemu.clear();
	
	if (Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).tag==1){	  
	  ileptonsmumu.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).ilept1);
	  ileptonsmumu.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).ilept2);
	  ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).ilept1);
	  ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).ilept2);
	}
	else if (goodZ.at(ileptonspTcleanedgoodZsv.at(l).at(0)).tag==2){
	  cout << "Checking ee pairs for QCD rejection" << endl;
	  ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).ilept1);
	  ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).ilept2);
	  ileptonsmumu.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).ilept1);
	  ileptonsmumu.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).ilept2);
	}

	double mass;
	
	// mumu
	cout << "Checking mumu pairs for QCD rejection" << endl;
	Lepton1qcd.SetPtEtaPhiM(RECOMU_PT[ileptonsmumu.at(0)], RECOMU_ETA[ileptonsmumu.at(0)], RECOMU_PHI[ileptonsmumu.at(0)], 0.105);
	Lepton2qcd.SetPtEtaPhiM(RECOMU_PT[ileptonsmumu.at(1)], RECOMU_ETA[ileptonsmumu.at(1)], RECOMU_PHI[ileptonsmumu.at(1)], 0.105);
	DiLeptonQCD=Lepton1qcd+Lepton2qcd;       
	mass = DiLeptonQCD.M();
	cout << "Mass of mumu for QCD rejection= " << mass << endl;
	if( mass < min_mass_2L ) min_mass_2L = mass ;	

	// ee	
	cout << "Checking ee pairs for QCD rejection" << endl;       
	Lepton1qcd.SetPtEtaPhiM(RECOELE_PT[ileptonsee.at(0)], RECOELE_ETA[ileptonsee.at(0)], RECOELE_PHI[ileptonsee.at(0)], 0.000511);
	Lepton2qcd.SetPtEtaPhiM(RECOELE_PT[ileptonsee.at(1)], RECOELE_ETA[ileptonsee.at(1)], RECOELE_PHI[ileptonsee.at(1)], 0.000511);
	DiLeptonQCD=Lepton1qcd+Lepton2qcd;       
	mass = DiLeptonQCD.M();
	cout << "Mass of ee for QCD rejection= " << mass << endl;
	if( mass < min_mass_2L ) min_mass_2L = mass ;	  
	

	//emu
	cout << "Checking emu pairs for QCD rejection" << endl;

	for (int i=0;i<2;i++){ 
	  for (int j=0;j<2;j++){
	    
	    if ( RECOMU_CHARGE[ileptonsmumu.at(i)] == RECOELE_CHARGE[ileptonsee.at(j)] ) continue; // shoud be OS	  
	    // evaluate the mass
	    double mass;
	    
	    Lepton1qcd.SetPtEtaPhiM(RECOMU_PT[ileptonsmumu.at(i)], RECOMU_ETA[ileptonsmumu.at(i)], RECOMU_PHI[ileptonsmumu.at(i)], 0.105);
	    Lepton2qcd.SetPtEtaPhiM(RECOELE_PT[ileptonsee.at(j)], RECOELE_ETA[ileptonsee.at(j)], RECOELE_PHI[ileptonsee.at(j)], 0.000511);
	    DiLeptonQCD=Lepton1qcd+Lepton2qcd;       
	    mass = DiLeptonQCD.M();
	    cout << "Mass of emu for QCD rejection= " << mass << endl;
	    if( mass < min_mass_2L ) min_mass_2L = mass ;
	  }
	}

	hminMll_6->Fill( min_mass_2L,newweight );
	if( min_mass_2L <= 4 ) { 
	  cout << "Not passing the mll>4 cut" << endl;
	  continue ;
	}
	else {
	  pTcleanedgoodZ.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)));
	  pTcleanedgoodZ.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)));
	  
	  iqcdcleanedgoodZs={ileptonspTcleanedgoodZsv.at(l).at(0),ileptonspTcleanedgoodZsv.at(l).at(1)};
	  iqcdcleanedgoodZsv.push_back(iqcdcleanedgoodZs);
	}
	
      }
//mark2
*/            
      
      // Z1 selection
      double pxZ1 = 0;  //Z1 kinematics
      double pyZ1 = 0;
      double pzZ1 = 0;
      double ptZ1 = 0;
      double EZ1 = 0;
      double Y_Z1 = -9;
      double massZ1 = 0;
      double massZ1_noFSR = 0;
      double sum_ptZ1 = 0.;
      int indexlep1Z1 = -1;
      int indexlep2Z1 = -1;
      int indexZ1= -1;
      int Z1tag=-999;
      
      // Choice of Z1 as the closest to the Z mass
      for (int i=0;i<pTcleanedgoodZ.size();++i){
	
	if( fabs(pTcleanedgoodZ.at(i).massvalue - Zmass) < fabs(massZ1 - Zmass) ){
	  
	  massZ1 = pTcleanedgoodZ.at(i).massvalue;
	  indexZ1=i;
	  
	  pxZ1 = pTcleanedgoodZ.at(i).pxZ;
	  pyZ1 = pTcleanedgoodZ.at(i).pyZ;
	  pzZ1 = pTcleanedgoodZ.at(i).pzZ;
	  EZ1  = pTcleanedgoodZ.at(i).EZ;
	  
	  ptZ1 = sqrt( pxZ1*pxZ1 + pyZ1*pyZ1 );
	  sum_ptZ1 = pTcleanedgoodZ.at(i).pt1+pTcleanedgoodZ.at(i).pt2;
	  
	  Y_Z1 = 0.5 * log ( (EZ1 + pzZ1)/(EZ1 - pzZ1) );
	  indexlep1Z1=pTcleanedgoodZ.at(i).ilept1;
	  indexlep2Z1=pTcleanedgoodZ.at(i).ilept2;
	  Z1tag=pTcleanedgoodZ.at(i).tag;
	}
      }
      
      if (massZ1 < 40.) {
	cout << "The mass of Z1 is < 40 GeV...exiting" << endl;
	continue;
      } 
      
      if( debug ) cout  << "\n Final Z1 properties: "
			<< "\n pxZ1 " << pxZ1
			<< "\n pyZ1 " << pyZ1
			<< "\n pzZ1 " << pzZ1
			<< "\n ptZ1 " << ptZ1
			<< "\n EZ1 "  << EZ1
			<< "\n Y_Z1 " << Y_Z1
			<< "\n massZ1 " << massZ1
			<< "\n indexlep1 " << indexlep1Z1
			<< "\n indexlep2 " << indexlep2Z1
			<< "\n indexZ1 " << indexZ1 
			<< "\n Z1 tag (1 for 2mu and 2 for 2e) " << Z1tag
		    
			<< endl;
      
  /*    
      // ZZ objects with alternative pairing
      array<int,2> icleanedgoodZs;
      vector<std::array<int, 2> > icleanedgoodZsv;
      
      for (int i=0;i<pTcleanedgoodZ.size();i++){
	for (int j=i+1;j<pTcleanedgoodZ.size();j++){
	  if (pTcleanedgoodZ.at(j).tag==pTcleanedgoodZ.at(i).tag) continue; // we want 2e2mu or 2mu2e
	  float massZa=-999.,massZb=-999.;
	  cout << "Masses= " << pTcleanedgoodZ.at(i).massvalue << " " << pTcleanedgoodZ.at(j).massvalue << endl;
	  if (fabs(pTcleanedgoodZ.at(i).massvalue-Zmass) < fabs(pTcleanedgoodZ.at(j).massvalue-Zmass)) {
	    massZa=pTcleanedgoodZ.at(i).massvalue;
	    massZb=pTcleanedgoodZ.at(j).massvalue;
	  }
	  else {
	    massZa=pTcleanedgoodZ.at(j).massvalue;
	    massZb=pTcleanedgoodZ.at(i).massvalue;
	  }	  
	  
	  //  if ( fabs(massZa-Zmass) < fabs(massZ1-Zmass) && massZb<12) continue; // exclude those pairs   smart cut non needed for 2e2mu
	  cout << "mass Za and b= " << massZa << " " << massZb << endl;
	  //	 cout << "i j " << i << " " << j << endl;
	  icleanedgoodZs={i,j};
	  icleanedgoodZsv.push_back(icleanedgoodZs);
	}
      }
      
      cout << "How many ZZ objects are present? " << icleanedgoodZsv.size() << endl;
      if (icleanedgoodZsv.size()==0) continue;
      
      // // sort Z by mass value                                                                                                                
      // struct SortCandByClosestToZ {                                                                                                          
      //   bool operator()( candidateZ c1, candidateZ c2) {                                                                                     
      //          return (fabs(c1.massvalue - Zmass) < fabs(c2.massvalue - Zmass));                                                             
      //   }                                                                                                                                    
      // };                                                                                                                                     
      // std::sort(vcandZ2.begin(),vcandZ2.end(),SortCandByClosestToZ());                                                                       
      
      double pxZ2 = 0;  //Z2 kinematics                                                                                                         
      double pyZ2 = 0;
      double pzZ2 = 0;
      double ptZ2 = 0;
      double EZ2 = 0;
      double Y_Z2 = -9;
      double massZ2 = 0;
      double massZ2_noFSR = 0;
      double sum_ptZ2 = 0.;
      int indexlep1Z2=-1;
      int indexlep2Z2=-1;
      int indexZ2=-1;
      int Z2tag=-999;
      
      // identify pairs  Z1+Z2s
      vector<int> vindexZZ;
      for (int l=0;l<icleanedgoodZsv.size();l++){
	cout << "Z+Z pairs with masses/indices: " 
	     << pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).massvalue << " and index " << icleanedgoodZsv.at(l).at(0) << " " 
	     << pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).massvalue << " and index " << icleanedgoodZsv.at(l).at(1) << endl;
	
	if ( (icleanedgoodZsv.at(l).at(0)==indexZ1 
	      //&&  
	      //!(pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept1==indexlep1Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept2==indexlep1Z1 || 
	      //  pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept1==indexlep2Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept2==indexlep2Z1)
	      ) 
	     || 
	     (icleanedgoodZsv.at(l).at(1)==indexZ1 
	      //&& 
	      //!(pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept1==indexlep1Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept2==indexlep1Z1 ||
	      //  pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept1==indexlep2Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept2==indexlep2Z1)
	      )
	     ) 
	  {
	    cout << "Found a Z1 + X  pair with no lepton shared with masses= " << pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).massvalue << " " <<  pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).massvalue << endl;
	    vindexZZ.push_back(l);
	  }
      }
      
      cout << "How many Z1+Z2s objects are present? " << vindexZZ.size() << endl;
      
      // If more than one Z2 (couple to the same Z1), choose the one with the highest-pT Z2 leptons
      if (vindexZZ.size()==1){
	cout << "Just one Z1+Z2 pair" << endl;
	if (icleanedgoodZsv.at(vindexZZ.at(0)).at(0)==indexZ1) indexZ2=icleanedgoodZsv.at(vindexZZ.at(0)).at(1); 
	if (icleanedgoodZsv.at(vindexZZ.at(0)).at(1)==indexZ1) indexZ2=icleanedgoodZsv.at(vindexZZ.at(0)).at(0);
      }
      else {
	cout << "more than one Z2 (couple to the same Z1), choose the one with the highest-pT Z2 leptons" << endl;
	int indexZ2tmp=-999;
	float sumpT=-999.,tmpsumpT=-999.;
	for (int ll=0;ll<vindexZZ.size();ll++){
	  if (icleanedgoodZsv.at(vindexZZ.at(ll)).at(0)==indexZ1) indexZ2tmp=icleanedgoodZsv.at(vindexZZ.at(ll)).at(1);
	  if (icleanedgoodZsv.at(vindexZZ.at(ll)).at(1)==indexZ1) indexZ2tmp=icleanedgoodZsv.at(vindexZZ.at(ll)).at(0);
	  if (indexZ2tmp>=0){
	    
	    TLorentzVector LorentzZ1,LorentzZ2,LorentzZZ;
	    LorentzZ1.SetPxPyPzE(pTcleanedgoodZ.at(indexZ1).pxZ, pTcleanedgoodZ.at(indexZ1).pyZ, pTcleanedgoodZ.at(indexZ1).pzZ, pTcleanedgoodZ.at(indexZ1).EZ);
	    LorentzZ2.SetPxPyPzE(pTcleanedgoodZ.at(indexZ2tmp).pxZ, pTcleanedgoodZ.at(indexZ2tmp).pyZ, pTcleanedgoodZ.at(indexZ2tmp).pzZ, pTcleanedgoodZ.at(indexZ2tmp).EZ);
	    LorentzZZ=LorentzZ1+LorentzZ2;
	    if (LorentzZZ.M()<70.) continue; // cut m4l>70
	    cout << "Passed m4l>70. cut"<< endl;
	    
	    tmpsumpT=
	      pTcleanedgoodZ.at(indexZ2tmp).pt1+
	      pTcleanedgoodZ.at(indexZ2tmp).pt2;
	    if (tmpsumpT>sumpT) {
	      sumpT=tmpsumpT;
	      indexZ2=indexZ2tmp;
	    }
	  }       
	}
      }
      
      // Z2
      if (indexZ2<0) continue;
      cout << "The highest pT leptons Z2 has mass= " <<  pTcleanedgoodZ.at(indexZ2).massvalue << endl;
      
      massZ2 = pTcleanedgoodZ.at(indexZ2).massvalue;	 
      pxZ2 = pTcleanedgoodZ.at(indexZ2).pxZ;
      pyZ2 = pTcleanedgoodZ.at(indexZ2).pyZ;
      pzZ2 = pTcleanedgoodZ.at(indexZ2).pzZ;
      EZ2  = pTcleanedgoodZ.at(indexZ2).EZ;	 
      ptZ2 = sqrt( pxZ2*pxZ2 + pyZ2*pyZ2 );
      sum_ptZ2 = pTcleanedgoodZ.at(indexZ2).pt1+pTcleanedgoodZ.at(indexZ2).pt2;	 
      Y_Z2 = 0.5 * log ( (EZ2 + pzZ2)/(EZ2 - pzZ2) );
      indexlep1Z2=pTcleanedgoodZ.at(indexZ2).ilept1;
      indexlep2Z2=pTcleanedgoodZ.at(indexZ2).ilept2;
      Z2tag=pTcleanedgoodZ.at(indexZ2).tag;
      
      // Z1 and Z2 final 
      cout << "Z1 has index= " << indexZ1 << "  Z2 has index= " << indexZ2 << endl;
      
      if (Z1tag==1) cout << "PTs= " << RECOMU_PT[indexlep1Z1] << " " << RECOMU_PT[indexlep2Z1] << " " <<  RECOELE_PT[indexlep1Z2] << " " << RECOELE_PT[indexlep2Z2]<< endl; 
      if (Z1tag==2) cout << "PTs= " << RECOELE_PT[indexlep1Z1] << " " << RECOELE_PT[indexlep2Z1] << " " <<  RECOMU_PT[indexlep1Z2] << " " << RECOMU_PT[indexlep2Z2]<< endl; 
      
      if (std::isnan(massZ2)) {
	cout << "No Z2 found" << endl;
	continue; 
      }
      else {
	if( debug ) cout  << "\n Final Z2 properties: "  
			  << "\n pxZ2 " << pxZ2
			  << "\n pyZ2 " << pyZ2
			  << "\n pzZ2 " << pzZ2
			  << "\n ptZ2 " << ptZ2
			  << "\n EZ2 "  << EZ2
			  << "\n Y_Z2 " << Y_Z2
			  << "\n massZ2 " << massZ2
			  << "\n indexlep1_Z2 " << indexlep1Z2
			  << "\n indexlep2_Z2 " << indexlep2Z2
			  << endl;
      }
      
      if( debug && has_FSR_Z2) {
      	cout  << " Z2 has FSR! " << endl;
     	if (Z1tag==1) {
	  cout  << "  pi " << pi2 << " --> index: " << iLp[pi2] << " associated lepton: " << iLp_l[pi2] << " (= "<< iL[i2]<<" ? )  tag: " << iLp_tagEM[pi2] << endl;
	  cout  << "  pj " << pj2 << " --> index: " << iLp[pj2] << " associated lepton: " << iLp_l[pj2] << " (= "<< iL[j2]<<" ? )  tag: " << iLp_tagEM[pj2] << endl;
	}
	else if (Z1tag==0){
	  cout  << "  pi " << pi2 << " --> index: " << iLp[pi2] << " associated lepton: " << iLp_l[pi2] << " (= "<< iLe[i2]<<" ? )  tag: " << iLp_tagEM[pi2] << endl;
	  cout  << "  pj " << pj2 << " --> index: " << iLp[pj2] << " associated lepton: " << iLp_l[pj2] << " (= "<< iLe[j2]<<" ? )  tag: " << iLp_tagEM[pj2] << endl;
	}
      }
      
      ++N_4b ;  // fill counter
      N_4b_w=N_4b_w+newweight;

      if( N_Z2_pairs > 1 ) {
	++N_4c ;  // fill counter
	N_4c_w=N_4c_w+newweight;
      }
      
      
      ++N_4d ;  // fill counter
      N_4d_w=N_4d_w+newweight;
    */  
      // **** Step 5:

      
      // int iL_4L[4] = {iL[i1],iL[j1],iL[i2],iL[j2]};
      
      // if( debug ) cout << "\n 4L Leptons' indeces: "
      // 		       << "\niL_4L[0]: " << iL_4L[0]
      // 		       << "\niL_4L[1]: " << iL_4L[1]
      // 		       << "\niL_4L[2]: " << iL_4L[2]
      // 		       << "\niL_4L[3]: " << iL_4L[3]
      // 		       << endl ;
      

      // // Execute Efficiency Reweighting
      Double_t eff_weight = 1.;
      
 
      TLorentzVector hP4,Z1P4,Z2P4;
      Z1P4.SetPxPyPzE(pxZ1,pyZ1,pzZ1,EZ1);
//      Z2P4.SetPxPyPzE(pxZ2,pyZ2,pzZ2,EZ2);
//      hP4 = Z1P4 + Z2P4;
//      cout << "Mass of 4l best candidate= " << hP4.M() << endl;
//      double mass4l=hP4.M();      

//      hM4l_5->Fill( mass4l,newweight );
  
//only fill double electron,remember to change it back----- qier
      if(Z1tag==2){ 
      hMZ1_5->Fill( massZ1,newweight );
      hPtZ1_5->Fill( ptZ1,newweight );
      hYZ1_5->Fill( Y_Z1,newweight );
      }
//      hMZ2_5->Fill( massZ2,newweight );
//      hPtZ2_5->Fill( ptZ2,newweight );
//      hYZ2_5->Fill( Y_Z2,newweight );

//mark1   
/*   
      // sort index by pt (kinematics not corrected for FSR)
      int ipt[4] ;
      double tmp_pt[4];
      if (Z1tag==1) cout << "PTs= " << RECOMU_PT[indexlep1Z1] << " " << RECOMU_PT[indexlep2Z1] << " " <<  RECOELE_PT[indexlep1Z2] << " " << RECOELE_PT[indexlep2Z2]<< endl;
      if (Z1tag==2) cout << "PTs= " << RECOELE_PT[indexlep1Z1] << " " << RECOELE_PT[indexlep2Z1] << " " <<  RECOMU_PT[indexlep1Z2] << " " << RECOMU_PT[indexlep2Z2]<< endl;

      int indexleptonfinal[4]={indexlep1Z1,indexlep2Z1,indexlep1Z2,indexlep2Z2};
      //cout << "PTs= " << RECOMU_PT[indexleptonfinal[0]] << " " << RECOMU_PT[indexleptonfinal[1]] << " " <<  RECOMU_PT[indexleptonfinal[2]] << " " << RECOMU_PT[indexleptonfinal[3]]<< endl;

      for(int i = 0; i < 2; ++i){ 
	if (Z1tag==1) tmp_pt[i] =  RECOMU_PT[indexleptonfinal[i]];
	if (Z1tag==2) tmp_pt[i] =  RECOELE_PT[indexleptonfinal[i]];
	cout << tmp_pt[i] << endl;
      }
      for(int i = 2; i < 4; ++i){ 
	if (Z2tag==1) tmp_pt[i] =  RECOMU_PT[indexleptonfinal[i]];
	if (Z2tag==2) tmp_pt[i] =  RECOELE_PT[indexleptonfinal[i]];
	cout << tmp_pt[i] << endl;
      }
      
      float sortedpT[4];

      for(int i = 0; i < 4; ++i){		
        double tmp_max_pt = 0;
      	int jj = i;
        for(int j = 0; j < 4; ++j){
      		if( tmp_pt[j] > tmp_max_pt ){
      			tmp_max_pt = tmp_pt[j];
      			jj  = j;
      		}
      	}
	sortedpT[i]=tmp_max_pt;
      	ipt[i] = indexleptonfinal[jj];
      	tmp_pt[jj] = 0;	
      }
      if(debug)
	for(int i = 0; i < 4; ++i)
	  cout << "\n ipt[" << i << "] = " << ipt[i] << "\t pt = " << sortedpT[i] << endl; 
      //end sorting
     
*/
      // Format lepton syncronization                                                                                                                                      
      // {run}:{lumi}:{event}:{pdgId}:{pT:.2f}:{eta:.2f}:{phi:.2f}{SIP:.2f}:{PFChargedHadIso:.2f}:{PFNeutralHadIso:.2f}:{PFPhotonIso:.2f}:{PUCorr:.2f}:{combRelIsoPF:.3f}:{eleBDT:.3f}:{photpT:.2f}:{photDR:.2f}:{photRelIso:.2f}          
  //mark1
/*    
      Char_t leptformat[20000];

      for(int i = 0; i < 4; ++i){      
	bool ismuon=false,iselectron=false;
	float dummy=0.;
	int flagFSR_tag=-999;
	int pfsr=-999;
	
	for (int j=0;j<RECO_NMU;j++) {	    
	  if (RECOMU_PT[j]==sortedpT[i]) {
	    ismuon=true;
	    break;
	  }
	  else {
	    for (int j=0;j<RECO_NELE;j++) {	    
	      if (RECOELE_PT[j]==sortedpT[i]) {
		iselectron=true;
		break;
	      }
	    }
	  }
	}
	
	
	cout << "isMuon= " << ismuon << " and isElectron= " << iselectron << endl;

	if (ismuon){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == ipt[i] && iLp_tagEM[ p ] == 0 )  {
	      cout << "Muon with pT= " << RECOMU_PT[ipt[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	      flagFSR_tag=0;
	      pfsr=p;
	      break;
	    }
	  }
	}
	    
	else if (iselectron){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == ipt[i] && iLp_tagEM[ p ] == 1 )  {
	      cout << "Electron with pT= " << RECOELE_PT[ipt[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	      flagFSR_tag=1;
	      pfsr=p;
	      break;
	    }
	  }
	}
	

	if (ismuon && flagFSR_tag==0){
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f:%.2f:%.2f:%.2f",
		   Run,LumiSection,Event,
		   int(-13*RECOMU_CHARGE[ipt[i]]),
		   RECOMU_PT[ipt[i]],RECOMU_ETA[ipt[i]],RECOMU_PHI[ipt[i]],RECOMU_SIP[ipt[i]],
		   RECOMU_PFchHad[ipt[i]],RECOMU_PFneuHad[ipt[i]],RECOMU_PFphoton[ipt[i]],RECOMU_PFPUchAllPart[ipt[i]],RECOMU_PFX_dB[ipt[i]],dummy,
		   RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_DR[iLp[pfsr]],RECOPFPHOT_PFX_rho[iLp[pfsr]]
		   );
	}
	else if (iselectron && flagFSR_tag==1){
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f:%.2f:%.2f:%.2f",
		   Run,LumiSection,Event,
		   int(-11*RECOELE_CHARGE[ipt[i]]),
		   RECOELE_PT[ipt[i]],RECOELE_ETA[ipt[i]],RECOELE_PHI[ipt[i]],RECOELE_SIP[ipt[i]],
		   RECOELE_PFchHad[ipt[i]],RECOELE_PFneuHad[ipt[i]],RECOELE_PFphoton[ipt[i]],RHO_ele,RECOELE_PFX_rho[ipt[i]],RECOELE_mvaNonTrigV0[ipt[i]],
		   RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_DR[iLp[pfsr]],RECOPFPHOT_PFX_rho[iLp[pfsr]]
		   );
	}	
	else if (ismuon){	  	  	      
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f",
		   Run,LumiSection,Event,
		   int(-13*RECOMU_CHARGE[ipt[i]]),
		   RECOMU_PT[ipt[i]],RECOMU_ETA[ipt[i]],RECOMU_PHI[ipt[i]],RECOMU_SIP[ipt[i]],
		   RECOMU_PFchHad[ipt[i]],RECOMU_PFneuHad[ipt[i]],RECOMU_PFphoton[ipt[i]],RECOMU_PFPUchAllPart[ipt[i]],RECOMU_PFX_dB[ipt[i]],dummy
		   );
	}
	else if (iselectron){
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f",
		   Run,LumiSection,Event,
		   int(-11*RECOELE_CHARGE[ipt[i]]),
		   RECOELE_PT[ipt[i]],RECOELE_ETA[ipt[i]],RECOELE_PHI[ipt[i]],RECOELE_SIP[ipt[i]],
		   RECOELE_PFchHad[ipt[i]],RECOELE_PFneuHad[ipt[i]],RECOELE_PFphoton[ipt[i]],RHO_ele,RECOELE_PFX_rho[ipt[i]],RECOELE_mvaNonTrigV0[ipt[i]]
		   );
	}	  
	  
	output_txt  << leptformat << endl;


      // N.B. Do NOT Update the Isolation values and correct the 4 momenta of leptons for FSR
      for(int i = 0; i < N_good; ++i){
	int flagFSR=0;
	int pfsr=-999;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if (iLp[p]==-1) continue;
	  if (iLp_l[p]==-1) continue;
	  
	  cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << iL[i] << endl;
	  if( iLp_l[ p ] == iL[i] && iLp_tagEM[ p ] == 0 )  {
	    cout << "Muon with pT= " << RECOMU_PT[iL[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    
	    flagFSR=1;
	    pfsr=p;
	  }
	}
	
	
	if (flagFSR==1){
	  cout << "Before correcting for FSR; muon pT= " << RECOMU_PT[iL[i]] << " Eta= " << RECOMU_ETA[iL[i]] << " Phi= " << RECOMU_PHI[iL[i]] << endl;
	  TLorentzVector Lept,LeptCorrection;
	  Lept.SetPtEtaPhiM(RECOMU_PT[iL[i]], RECOMU_ETA[iL[i]], RECOMU_PHI[iL[i]], 0.105);
	  LeptCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_ETA[iLp[pfsr]],RECOPFPHOT_PHI[iLp[pfsr]],0);
	  Lept+=LeptCorrection;
	  RECOMU_PT[iL[i]]=Lept.Pt();
	  RECOMU_ETA[iL[i]]=Lept.Eta();
	  RECOMU_PHI[iL[i]]=Lept.Phi();
	  cout << "After correcting for FSR; muon pT= " << RECOMU_PT[iL[i]] << " Eta= " << RECOMU_ETA[iL[i]] << " Phi= " << RECOMU_PHI[iL[i]] << endl;
	}
      }
      
      // N.B. DO NOT Update the Isolation values and correct the 4 momenta of leptons for FSR
      for(int i = 0; i < Ne_good; ++i){
	int flagFSR=0;
	int pfsr=-999;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if (iLp[p]==-1) continue;
	  if (iLp_l[p]==-1) continue;
	  
	  cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << iLe[i] << endl;
	  if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {
	    cout << "Electron with pT= " << RECOELE_PT[iLe[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    flagFSR=1;
	    pfsr=p;
	  }
	}
	
	if (flagFSR==1){
	  cout << "Before correcting for FSR; electron pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	  TLorentzVector Lept,LeptCorrection;
	  Lept.SetPtEtaPhiM(RECOELE_PT[iLe[i]], RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]], 0.105);
	  LeptCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_ETA[iLp[pfsr]],RECOPFPHOT_PHI[iLp[pfsr]],0);
	  Lept+=LeptCorrection;
	  RECOELE_PT[iLe[i]]=Lept.Pt();
	  RECOELE_ETA[iLe[i]]=Lept.Eta();
	  RECOELE_PHI[iLe[i]]=Lept.Phi();
	  cout << "After correcting for FSR; muon pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	}
      }
      
      cout << "Kinematics of leptons corrected for FSR photons (if existing)" << endl;
      

     //Basic cuts to jets AND delta R section
     int njets_pass=0;
     TLorentzVector JET1,JET2;
     int jet1=-999,jet2=-999;      
     int jetfail[100];

     for(int i=0;i<100;i++) jetfail[i]=0;
     
     for(int i=0;i<RECO_PFJET_N;i++){
       cout<<i<<" Jet with pt= "<<RECO_PFJET_PT[i]<<" ETA "<<RECO_PFJET_ETA[i]<<" PUID "<<RECO_PFJET_PUID[i] << " PUID_MVA "<< RECO_PFJET_PUID_MVA[i]<<endl;
       
       if(RECO_PFJET_PT[i]>30. && fabs(RECO_PFJET_ETA[i])<4.7 ){
       
      	 //Check that jet has deltaR>0.4 away from any tight lepton corrected for FSR
	 for(int mu = 0; mu < N_good; ++mu){
	   if (fabs(RECOMU_SIP[iL[mu]])>=4.) continue;
      	   if (RECOMU_PFX_dB_new[iL[mu]]>=0.35) continue;
	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOMU_PHI[iL[mu]]),2) + pow(RECO_PFJET_ETA[i] - RECOMU_ETA[iL[mu]],2));
	   cout << "1st lepton muon: " << " pT=" << RECOMU_PT[iL[mu]] <<" deltaR "<< deltaR <<endl;	   
	   if (deltaR<0.4){
	     jetfail[i]=1;
     	     cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
      	 for(int ele = 0; ele < Ne_good; ++ele){
      	   if (fabs(RECOELE_SIP[iLe[ele]])>=4.) continue;
	   if (RECOELE_PFX_rho_new[iLe[ele]]>=0.35) continue;
      	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOELE_PHI[iLe[ele]]),2) + pow(RECO_PFJET_ETA[i] - RECOELE_ETA[iLe[ele]],2));
     	   cout << "1st lepton electron: " << " pT=" << RECOELE_PT[iLe[ele]] <<" deltaR "<< deltaR <<endl;
	   if (deltaR<0.4){
     	     jetfail[i]=1;
     	     cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
	 // cleaning w.r.t FSR photons attached to leptons
	 for(int j=0.;j<Nphotons;j++) {
           if (iLp_l[j]!=-1 && (iLp_tagEM[j]==0 || iLp_tagEM[j]==1) ) {
	     if (iLp_tagEM[j]==0) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[j]] << endl;
	     if (iLp_tagEM[j]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a electron with pT= " << RECOELE_PT[iLp_l[j]] << endl;
	     double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOPFPHOT_PHI[iLp[j]]),2) + pow(RECO_PFJET_ETA[i] - RECOPFPHOT_ETA[iLp[j]],2));
	     if (deltaR<0.4){
	       jetfail[i]=1;
	       cout << " jetfail " << jetfail[i] <<endl;
	       break;
	     }
	   }
         }
	 // 


	 if (jetfail[i]==0){
	   cout<< " PASS jet " <<i<<" PT= "<<RECO_PFJET_PT[i]<<" ETA= "<<RECO_PFJET_ETA[i]<<" PUID= "<<RECO_PFJET_PUID[i]<<endl;
	   njets_pass++;
	   if (njets_pass==1){
	     jet1=i;
	     JET1.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	   }
	   if (njets_pass==2){
	     jet2=i;
	     JET2.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	   }
	 }

       }
       else{
      	 jetfail[i]=1;
       }
       //cout<<" JETFAIL "<<jetfail[i]<<endl;
     }
     
     // b-tagged jets - ordered in pT
     int n_bjets=0;
     int index_bjets[2]={-999,-999};
     
     for (int i=0;i<50;i++){
	 if (cSV_BTagJet_DISCR[i] > 0.8){ // 76x
	 if(cSV_BTagJet_PT[i]>30. && fabs(cSV_BTagJet_ETA[i])<4.7 ) cout << "Found a bjet (pT>30 and |eta|<2.4) with pT= " << cSV_BTagJet_PT[i] << endl;	 
	 n_bjets++;
	 if (n_bjets==1) index_bjets[0]=i; 
	 if (n_bjets==2) index_bjets[1]=i;	   
       }
     }
     cout << "Number of b-jets passing the selection= " << n_bjets << endl;
*/     
     // VBF-tagged category - category 2
     // filling branches in the reduced tree
     f_weight = newweight;
     
     f_int_weight = -1;
     
     f_pu_weight = pu_weight;
     f_eff_weight = eff_weight;

     f_run = Run;
     f_event = Event;
     f_lumi = LumiSection;
//mark1
     if (Z1tag==1){
       f_lept1_pt  = RECOMU_PT[indexlep1Z1] ;
       f_lept1_eta = RECOMU_ETA[indexlep1Z1] ;
       f_lept1_phi = RECOMU_PHI[indexlep1Z1];
       f_lept1_charge = RECOMU_CHARGE[indexlep1Z1];
       f_lept1_pfx = RECOMU_PFX_dB_new[indexlep1Z1];
       f_lept1_sip = RECOMU_SIP[indexlep1Z1];

       f_lept2_pt  = RECOMU_PT[indexlep2Z1] ;
       f_lept2_eta = RECOMU_ETA[indexlep2Z1] ;
       f_lept2_phi = RECOMU_PHI[indexlep2Z1];
       f_lept2_charge = RECOMU_CHARGE[indexlep2Z1];
       f_lept2_pfx = RECOMU_PFX_dB_new[indexlep2Z1];
       f_lept2_sip = RECOMU_SIP[indexlep2Z1];
     }
     else if (Z1tag==2) {
       f_lept1_pt = RECOELE_PT[indexlep1Z1] ;
       f_lept1_eta = RECOELE_ETA[indexlep1Z1] ;
       f_lept1_phi = RECOELE_PHI[indexlep1Z1];
       f_lept1_charge = RECOELE_CHARGE[indexlep1Z1];
       f_lept1_pfx = RECOELE_PFX_rho_new[indexlep1Z1];
       f_lept1_sip = RECOELE_SIP[indexlep1Z1];
       f_lept1_mvaid = RECOELE_mvaNonTrigV0[indexlep1Z1];
       f_lept2_pt = RECOELE_PT[indexlep2Z1] ;
       f_lept2_eta = RECOELE_ETA[indexlep2Z1] ;
       f_lept2_phi = RECOELE_PHI[indexlep2Z1];
       f_lept2_charge = RECOELE_CHARGE[indexlep2Z1];
       f_lept2_pfx = RECOELE_PFX_rho_new[indexlep2Z1];
       f_lept2_sip = RECOELE_SIP[indexlep2Z1];
       f_lept2_mvaid = RECOELE_mvaNonTrigV0[indexlep2Z1];
}
/*
     if (Z1tag==1){
       f_lept1_pt  = RECOMU_PT[indexleptonfinal[0]] ;
       f_lept1_eta = RECOMU_ETA[indexleptonfinal[0]] ;
       f_lept1_phi = RECOMU_PHI[indexleptonfinal[0]];
       f_lept1_charge = RECOMU_CHARGE[indexleptonfinal[0]];
       f_lept1_pfx = RECOMU_PFX_dB_new[indexleptonfinal[0]];
       f_lept1_sip = RECOMU_SIP[indexleptonfinal[0]];
       //    f_lept1_mvaid = RECOMU_mvaNonTrigV0[indexleptonfinal[0]];
       f_lept2_pt  = RECOMU_PT[indexleptonfinal[1]] ;
       f_lept2_eta = RECOMU_ETA[indexleptonfinal[1]] ;
       f_lept2_phi = RECOMU_PHI[indexleptonfinal[1]];
       f_lept2_charge = RECOMU_CHARGE[indexleptonfinal[1]];
       f_lept2_pfx = RECOMU_PFX_dB_new[indexleptonfinal[1]];
       f_lept2_sip = RECOMU_SIP[indexleptonfinal[1]];
       //    f_lept2_mvaid = RECOMU_mvaNonTrigV0[indexleptonfinal[1]];
       f_lept3_pt  = RECOELE_PT[indexleptonfinal[2]] ;
       f_lept3_eta = RECOELE_ETA[indexleptonfinal[2]] ;
       f_lept3_phi = RECOELE_PHI[indexleptonfinal[2]];
       f_lept3_charge = RECOELE_CHARGE[indexleptonfinal[2]];
       f_lept3_pfx = RECOELE_PFX_rho_new[indexleptonfinal[2]];
       f_lept3_sip = RECOELE_SIP[indexleptonfinal[2]];
       f_lept3_mvaid = RECOELE_mvaNonTrigV0[indexleptonfinal[2]];

       f_lept4_pt  = RECOELE_PT[indexleptonfinal[3]] ;
       f_lept4_eta = RECOELE_ETA[indexleptonfinal[3]] ;
       f_lept4_phi = RECOELE_PHI[indexleptonfinal[3]];
       f_lept4_charge = RECOELE_CHARGE[indexleptonfinal[3]];
       f_lept4_pfx = RECOELE_PFX_rho_new[indexleptonfinal[3]];
       f_lept4_sip = RECOELE_SIP[indexleptonfinal[3]];
       f_lept4_mvaid = RECOELE_mvaNonTrigV0[indexleptonfinal[3]];
     }
     else if (Z1tag==2) {
       f_lept1_pt = RECOELE_PT[indexleptonfinal[0]] ;
       f_lept1_eta = RECOELE_ETA[indexleptonfinal[0]] ;
       f_lept1_phi = RECOELE_PHI[indexleptonfinal[0]];
       f_lept1_charge = RECOELE_CHARGE[indexleptonfinal[0]];
       f_lept1_pfx = RECOELE_PFX_rho_new[indexleptonfinal[0]];
       f_lept1_sip = RECOELE_SIP[indexleptonfinal[0]];
       f_lept1_mvaid = RECOELE_mvaNonTrigV0[indexleptonfinal[0]];
       f_lept2_pt = RECOELE_PT[indexleptonfinal[1]] ;
       f_lept2_eta = RECOELE_ETA[indexleptonfinal[1]] ;
       f_lept2_phi = RECOELE_PHI[indexleptonfinal[1]];
       f_lept2_charge = RECOELE_CHARGE[indexleptonfinal[1]];
       f_lept2_pfx = RECOELE_PFX_rho_new[indexleptonfinal[1]];
       f_lept2_sip = RECOELE_SIP[indexleptonfinal[1]];
       f_lept2_mvaid = RECOELE_mvaNonTrigV0[indexleptonfinal[1]];

       f_lept3_pt = RECOMU_PT[indexleptonfinal[2]] ;
       f_lept3_eta = RECOMU_ETA[indexleptonfinal[2]] ;
       f_lept3_phi = RECOMU_PHI[indexleptonfinal[2]];
       f_lept3_charge = RECOMU_CHARGE[indexleptonfinal[2]];
       f_lept3_pfx = RECOMU_PFX_dB_new[indexleptonfinal[2]];
       f_lept3_sip = RECOMU_SIP[indexleptonfinal[2]];
       //    f_lept3_mvaid = RECOMU_mvaNonTrigV0[indexleptonfinal[2]];
       f_lept4_pt = RECOMU_PT[indexleptonfinal[3]] ;
       f_lept4_eta = RECOMU_ETA[indexleptonfinal[3]] ;
       f_lept4_phi = RECOMU_PHI[indexleptonfinal[3]];
       f_lept4_charge = RECOMU_CHARGE[indexleptonfinal[3]];
       f_lept4_pfx = RECOMU_PFX_dB_new[indexleptonfinal[3]];
       f_lept4_sip = RECOMU_SIP[indexleptonfinal[3]];
       //    f_lept4_mvaid = RECOMU_mvaNonTrigV0[indexleptonfinal[3]];
     }
*/     
     //f_iso_max = Iso_max;
     //f_sip_max = Sip_max;
     f_Z1mass = massZ1;
//     f_Z2mass = massZ2;
     f_angle_costhetastar = -99.;
     f_angle_costheta1 = -99.;
     f_angle_costheta2 = -99.;
     f_angle_phi = -99.;
     f_angle_phistar1 = -99.;
//   f_pt4l = -99.;
//   f_eta4l = -99;
//     f_mass4l = mass4l;
//     f_mass4lErr = massErr;
  //   f_njets_pass = njets_pass;
     f_deltajj = -999.;
     f_massjj = -999.;
     f_D_jet = -999.;
//mark
/*     
     if (njets_pass<1) {
       f_jet1_pt = -999;
       f_jet1_eta = -999;
       f_jet1_phi = -999;
       f_jet1_e = -999;
       f_jet2_pt = -999;
       f_jet2_eta = -999;
       f_jet2_phi = -999;
       f_jet2_e = -999;
     }
     else if (njets_pass==1) { // store the highest pT jets passing the criteria
       f_jet1_pt = RECO_PFJET_PT[jet1];
       f_jet1_eta = RECO_PFJET_ETA[jet1];
       f_jet1_phi = RECO_PFJET_PHI[jet1];
       f_jet1_e = RECO_PFJET_ET[jet1];
     }
     else if (njets_pass>=2) { // store the second highest pT jets passing the criteria
       f_jet2_pt = RECO_PFJET_PT[jet2];
       f_jet2_eta = RECO_PFJET_ETA[jet2];
       f_jet2_phi = RECO_PFJET_PHI[jet2];
       f_jet2_e = RECO_PFJET_ET[jet2];

       TLorentzVector JET1new,JET2new,mJJnew;
       JET1new.SetPtEtaPhiE(RECO_PFJET_PT[jet1],RECO_PFJET_ETA[jet1],RECO_PFJET_PHI[jet1],RECO_PFJET_ET[jet1]*TMath::CosH(RECO_PFJET_ETA[jet1]));
       JET2new.SetPtEtaPhiE(RECO_PFJET_PT[jet2],RECO_PFJET_ETA[jet2],RECO_PFJET_PHI[jet2],RECO_PFJET_ET[jet2]*TMath::CosH(RECO_PFJET_ETA[jet2]));
       mJJnew=JET1new+JET2new;
       f_deltajj = fabs(JET1new.Eta()-JET2new.Eta());
       f_massjj = mJJnew.M();
       f_D_jet =0.18*fabs(JET1new.Eta()-JET2new.Eta())+1.92E-4*mJJnew.M();
     }
 */      
     f_pfmet=RECO_PFMET;

     
     //MELA discriminant - Lepton kinematics already corrected for FSR
     TLorentzVector L11P4,L12P4,L21P4,L22P4;
     float angle_costheta1,angle_costheta2,angle_Phi,angle_costhetastar,angle_Phi1;

 //    cout << "Final: pT of final leptons from Z1 and Z2 =";
//     if (Z1tag==1) cout << RECOMU_PT[indexleptonfinal[0]] << " " << RECOMU_PT[indexleptonfinal[1]] << " " << RECOELE_PT[indexleptonfinal[2]] << " " << RECOELE_PT[indexleptonfinal[3]] <<endl;
/*     if (Z1tag==2) cout << RECOELE_PT[indexleptonfinal[0]] << " " << RECOELE_PT[indexleptonfinal[1]] << " " << RECOMU_PT[indexleptonfinal[2]] << " " << RECOMU_PT[indexleptonfinal[3]] <<endl;
     
     if (Z1tag==1 && Z2tag==2){
       L11P4.SetPtEtaPhiM(RECOMU_PT[indexleptonfinal[0]], RECOMU_ETA[indexleptonfinal[0]], RECOMU_PHI[indexleptonfinal[0]], 0.105);
       L12P4.SetPtEtaPhiM(RECOMU_PT[indexleptonfinal[1]], RECOMU_ETA[indexleptonfinal[1]], RECOMU_PHI[indexleptonfinal[1]], 0.105);
       L21P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[2]], RECOELE_ETA[indexleptonfinal[2]], RECOELE_PHI[indexleptonfinal[2]], 0.000511);
       L22P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[3]], RECOELE_ETA[indexleptonfinal[3]], RECOELE_PHI[indexleptonfinal[3]], 0.000511);
     }

     if (Z1tag==2 && Z2tag==1){
       L11P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[0]], RECOELE_ETA[indexleptonfinal[0]], RECOELE_PHI[indexleptonfinal[0]], 0.000511);
       L12P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[1]], RECOELE_ETA[indexleptonfinal[1]], RECOELE_PHI[indexleptonfinal[1]], 0.000511);
       L21P4.SetPtEtaPhiM(RECOMU_PT[indexleptonfinal[2]], RECOMU_ETA[indexleptonfinal[2]], RECOMU_PHI[indexleptonfinal[2]], 0.105);
       L22P4.SetPtEtaPhiM(RECOMU_PT[indexleptonfinal[3]], RECOMU_ETA[indexleptonfinal[3]], RECOMU_PHI[indexleptonfinal[3]], 0.105);
     }
*/
     if (Z1tag==1){
       L11P4.SetPtEtaPhiM(RECOMU_PT[indexlep1Z1], RECOMU_ETA[indexlep1Z1], RECOMU_PHI[indexlep1Z1], 0.105);
       L12P4.SetPtEtaPhiM(RECOMU_PT[indexlep2Z1], RECOMU_ETA[indexlep2Z1], RECOMU_PHI[indexlep2Z1], 0.105);
     }
     if (Z1tag==2){
       L11P4.SetPtEtaPhiM(RECOELE_PT[indexlep1Z1], RECOELE_ETA[indexlep1Z1], RECOELE_PHI[indexlep1Z1], 0.000511);
       L12P4.SetPtEtaPhiM(RECOELE_PT[indexlep2Z1], RECOELE_ETA[indexlep2Z1], RECOELE_PHI[indexlep2Z1], 0.000511);
     }
/*
     //Assigning correct PID numbers
     int L11PID,L12PID,L21PID,L22PID;
     if(Z1tag==1){
       if (RECOMU_CHARGE[indexleptonfinal[0]] == 1){
         L11PID=-13;
         L12PID=+13;
       }
       else if (RECOMU_CHARGE[indexleptonfinal[0]] == -1){
         L11PID=13;
         L12PID=-13;
       }
     }
     else if(Z1tag==2){
       if (RECOELE_CHARGE[indexleptonfinal[0]] == 1){
         L11PID=-11;
         L12PID=+11;
       }
       else if (RECOELE_CHARGE[indexleptonfinal[0]] == -1){
         L11PID=11;
         L12PID=-11;
       }
     }


     if (Z1tag==1 && Z2tag==2){
       if (RECOMU_CHARGE[indexleptonfinal[0]] == 1){
	 L11PID=-13;
	 L12PID=+13;
       }
       else if (RECOMU_CHARGE[indexleptonfinal[0]] == -1){
	 L11PID=13;
	 L12PID=-13;
       }
       if (RECOELE_CHARGE[indexleptonfinal[2]] == 1){
	 L21PID=-11;
	 L22PID=11;
       }
       else if (RECOELE_CHARGE[indexleptonfinal[2]] == -1){
	 L21PID=11;
	 L22PID=-11;
       }       
     }    
     else if (Z1tag==2 && Z2tag==1){
       if (RECOELE_CHARGE[indexleptonfinal[0]] == 1){
	 L11PID=-11;
	 L12PID=+11;
       }
       else if (RECOELE_CHARGE[indexleptonfinal[0]] == -1){
	 L11PID=11;
	 L12PID=-11;
       }
       if (RECOMU_CHARGE[indexleptonfinal[2]] == 1){
	 L21PID=-13;
	 L22PID=13;
       }
       else if (RECOMU_CHARGE[indexleptonfinal[2]] == -1){
	 L21PID=13;
	 L22PID=-13;
       }
     }
     else cout << "What the hell?!?!"<<endl;
*/   
  //mark
/*        
     double massofhiggs=hP4.M();
     double massofZ1=Z1P4.M();
     double massofZ2=Z2P4.M();
     float Psig,Pbkg;
     float psPsig,psPbkg,psD,gravPsig,gravPbkg,gravD;     
     double pt4l=0.,eta4l=0.;
     cout << "Mass of Higgs passed to MELA= " << massofhiggs << endl;

     if(massofhiggs>0.){
       
       // MEM       
       vector<TLorentzVector> partP;
       partP.push_back(L11P4);
       partP.push_back(L12P4);
       partP.push_back(L21P4);
       partP.push_back(L22P4);
       
       vector<int> partId;
       partId.push_back(L11PID);
       partId.push_back(L12PID);
       partId.push_back(L21PID);
       partId.push_back(L22PID);
       
       
       vector<TLorentzVector> partPprod;
       vector<int> partIdprod;
       partPprod.push_back(L11P4);
       partPprod.push_back(L12P4);
       partPprod.push_back(L21P4);
       partPprod.push_back(L22P4);      
       if (JET1.Pt()>0.) partPprod.push_back(JET1);
       if (JET2.Pt()>0.) partPprod.push_back(JET2); // Can also use partPprod.push_back(nullFourVector) instead for integrated VBF MEs
       
       partIdprod.push_back(L11PID);
       partIdprod.push_back(L12PID);
       partIdprod.push_back(L21PID);
       partIdprod.push_back(L22PID);
       if (JET1.Pt()>0.) partIdprod.push_back(0);
       if (JET2.Pt()>0.) partIdprod.push_back(0); // For leptonic ZH in the future, this could actually be the opposite lepton flavor
       

       double p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;
       double p0minus_VAJHU, Dgg10_VAMCFM;
       double phjj_VAJHU, pvbf_VAJHU;
       
       int a=combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, partP, partId, p0plus_VAJHU); // Calculation of SM gg->H->4l JHUGen ME      
       cout << "a= "  << p0plus_VAJHU << endl;
       int b=combinedMEM.computeME(MEMNames::k0minus, MEMNames::kJHUGen, partP, partId, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME 
       int c=combinedMEM.computeME(MEMNames::kggHZZ_10, MEMNames::kMCFM, partP, partId, Dgg10_VAMCFM); // Direct calculation of Dgg (D^kin for off-shell) from MCFM MEs
       int d=combinedMEM.computeME(MEMNames::kqqZZ, MEMNames::kMCFM, partP, partId, bkg_VAMCFM); // qq->4l background calculation from MCFM
       combinedMEM.computePm4l(partP,partId, MEMNames::kNone, p0plus_m4l, bkg_m4l); // m4l probabilities for signal and background, nominal resolution
       if (njets_pass>=2){
	 int f=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partPprod, partIdprod, phjj_VAJHU); // SM gg->H+2j
	 int g=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_VAJHU);  // SM VBF->H
       }

       f_D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ); // D^kin_bkg
       f_D_bkg = p0plus_VAJHU * p0plus_m4l / ( p0plus_VAJHU * p0plus_m4l + bkg_VAMCFM * bkg_m4l ); // D^kin including superMELA
       f_D_gg = Dgg10_VAMCFM; // D_gg
       f_D_g4 = p0plus_VAJHU / ( p0plus_VAJHU + p0minus_VAJHU ); // D_0-
       if (njets_pass>=2) f_Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
       else f_Djet_VAJHU=-1;
       
       mela::computeAngles(partP[0], partId[0], partP[1], partId[1],
			   partP[2], partId[2], partP[3], partId[3],
			   angle_costhetastar,angle_costheta1,angle_costheta2,angle_Phi,angle_Phi1);
       
       
       
       TLorentzVector p4l = L11P4 + L12P4 + L21P4 + L22P4;
       pt4l = p4l.Pt();
       eta4l = p4l.Eta();
       float w; 
       //= mMEM->getMELAWeight();
       
              
       f_int_weight = w;
       f_Z1mass = massofZ1;
       f_Z2mass = massofZ2;
       f_angle_costhetastar = angle_costhetastar;
       f_angle_costheta1 = angle_costheta1;
       f_angle_costheta2 = angle_costheta2;
       f_angle_phi = angle_Phi;
       f_angle_phistar1 = angle_Phi1;
       f_pt4l = pt4l;
       f_eta4l = eta4l;
       newtree->Fill();              
       
       // Filling BNN input       
       bnn_file << newweight << " " << w << " ";
       if (Z1tag==1 && Z2tag==0){
	 printmubnn(indexleptonfinal[0]);
	 printmubnn(indexleptonfinal[1]);      
	 printelebnn(indexleptonfinal[2]);
	 printelebnn(indexleptonfinal[3]);
       }
       else if (Z1tag==0 && Z2tag==1){
	 printelebnn(indexleptonfinal[0]);
	 printelebnn(indexleptonfinal[1]);      
	 printmubnn(indexleptonfinal[2]);
	 printmubnn(indexleptonfinal[3]);
       }
       
       if(massofhiggs>100. && massofhiggs<180.){
	 hMELA_8->Fill(f_D_bkg_kin,newweight);
	 hMELA_vs_M4l_8->Fill(massofhiggs,f_D_bkg_kin,newweight);
       }
       if (massofhiggs>100. && f_D_bkg_kin > 0.1){       
	 hMELA_9->Fill(f_D_bkg_kin,newweight);
	 hMELA_vs_M4l_9->Fill(massofhiggs,f_D_bkg_kin,newweight);
	 ++N_9;
	 N_9_w=N_9_w+newweight;
       }
       if (f_D_g4 > 0.3){
	 ++N_9PS;
	 N_9PS_w=N_9PS_w+newweight;
       }
     }
     else {
       cout<<endl;
       continue;
     }
     
     // end of KD
   */ 
     //if( debug )
/*
     cout << "EVENT CANDIDATE: \n" 
	  << " N " << jentry 
	  << " RUN " << Run
	  << " EVENT " << Event	
	  << " LumiSection " << LumiSection
	  << " massZ1 " << massZ1

     output_txt << "EVENT CANDIDATE: \n" 
		<< " N " << jentry 
		<< " RUN " << Run
		<< " EVENT " << Event	
		<< " LumiSection " << LumiSection
		<< " massZ1 " << massZ1
		<< " massZ2 " << massZ2
     //<<endl;

     Char_t outformat[20000];
     float dummy=-1.;
     
     if (JET1.Pt()>0. && JET2.Pt()>0.)
       sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
		Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,JET1.Pt(),JET2.Pt(),category);
     else if (JET1.Pt()>0.)
       sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
                Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,JET1.Pt(),dummy,category);
     else sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
		Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,dummy,dummy,category);
     
     output_txt  << outformat << endl;
*/
	
     
     // fill final tree
     finaltree->Fill();

     

   } // end loop on entries

   // write on output txt file:

   
   bnn_file.close();
   output_txt.close();
   

/*
   cout << "N_0 "  << N_0  << " \n" 
	      << "N_01 " << N_01 << " \n"	
	      << "N_02 " << N_02 << " \n"	
	      << "N_1 "  << N_1  << " \n"	
	      << "N_2 "  << N_2  << " \n"	
	      << "N_3a " << N_3a << " \n"	
	      << "N_3_FSR " << N_3_FSR << " \n"	
	      << "N_3b " << N_3b << " \n"	
	      << "N_4a " << N_4a << " \n"	
	      << "N_4b " << N_4b << " \n"	
	      << "N_4c " << N_4c << " \n"	
	      << "N_4d " << N_4d << " \n"	
	      << "N_5 "  << N_5  << " \n"	
	      << "N_6 "  << N_6  << " \n"	
	      << "N_7 "  << N_7  << " \n"	
	      << "N_8 "  << N_8  << " \n"
	      << "N_8_a "<< N_8_a<< " \n"
      	      << "N_9 "  << N_9  << " \n"
      	      << "N_9_1FSR "  << N_9_1FSR  << " \n" 
	      << "N_9_2FSR "  << N_9_2FSR  << " \n" 
	      << "N_9PS "     << N_9PS << " \n"
	      << "N_9GRAV"    << N_9GRAV << "\n"
	      << "N_9a_VBF "  << N_9a_VBF << " \n"
              << "N_9b_VBF "  << N_9b_VBF << "\n"
	      << "N_VBF "     << N_VBF << " \n";
*/
/*
   nEvent_4l->GetXaxis()->SetBinLabel(1,"Init.");
   nEvent_4l->GetXaxis()->SetBinLabel(2,"MCTruth: 4mu");
   nEvent_4l->GetXaxis()->SetBinLabel(3,"MCTruth: Acc");
   nEvent_4l->GetXaxis()->SetBinLabel(4,"Init");
   nEvent_4l->GetXaxis()->SetBinLabel(5,"HLT");
   nEvent_4l->GetXaxis()->SetBinLabel(6,"Z1 lept. cuts");
   nEvent_4l->GetXaxis()->SetBinLabel(7,"Z1+#gamma");
   nEvent_4l->GetXaxis()->SetBinLabel(8,"m_{Z1}");
   nEvent_4l->GetXaxis()->SetBinLabel(9,"4#mu");
   nEvent_4l->GetXaxis()->SetBinLabel(10,"at least one Z2");
   nEvent_4l->GetXaxis()->SetBinLabel(11,"Z2 lept. cuts");
   nEvent_4l->GetXaxis()->SetBinLabel(12,"m_{Z2}");
   nEvent_4l->GetXaxis()->SetBinLabel(13,"pT cuts");
   nEvent_4l->GetXaxis()->SetBinLabel(14,"mll>4 for OS-SF");
   nEvent_4l->GetXaxis()->SetBinLabel(15,"m4l > 70");
   nEvent_4l->GetXaxis()->SetBinLabel(16,"m_{Z2} > 12");
   nEvent_4l->GetXaxis()->SetBinLabel(17,"m4l > 100");
   nEvent_4l->GetXaxis()->SetBinLabel(18,"MELA KD > 0.1");
   nEvent_4l->GetXaxis()->SetBinLabel(19,"one Z+#gamma");
   nEvent_4l->GetXaxis()->SetBinLabel(20,"two Z+#gamma");


   nEvent_4l_w->GetXaxis()->SetBinLabel(1,"Init."); 
   nEvent_4l_w->GetXaxis()->SetBinLabel(2,"MCTruth: 4mu");  
   nEvent_4l_w->GetXaxis()->SetBinLabel(3,"MCTruth: Acc");
   nEvent_4l_w->GetXaxis()->SetBinLabel(4,"Init");
   nEvent_4l_w->GetXaxis()->SetBinLabel(5,"HLT");
   nEvent_4l_w->GetXaxis()->SetBinLabel(6,"Z1 lept. cuts");
   nEvent_4l_w->GetXaxis()->SetBinLabel(7,"Z1+#gamma");
   nEvent_4l_w->GetXaxis()->SetBinLabel(8,"m_{Z1}");
   nEvent_4l_w->GetXaxis()->SetBinLabel(9,"4#mu");
   nEvent_4l_w->GetXaxis()->SetBinLabel(10,"at least one Z2");
   nEvent_4l_w->GetXaxis()->SetBinLabel(11,"Z2 lept. cuts");
   nEvent_4l_w->GetXaxis()->SetBinLabel(12,"m_{Z2}");
   nEvent_4l_w->GetXaxis()->SetBinLabel(13,"pT cuts");
   nEvent_4l_w->GetXaxis()->SetBinLabel(14,"mll>4 for OS-SF");
   nEvent_4l_w->GetXaxis()->SetBinLabel(15,"m4l > 70");
   nEvent_4l_w->GetXaxis()->SetBinLabel(16,"m_{Z2} > 12");
   nEvent_4l_w->GetXaxis()->SetBinLabel(17,"m4l > 100");
   nEvent_4l_w->GetXaxis()->SetBinLabel(18,"MELA KD > 0.1");
   nEvent_4l_w->GetXaxis()->SetBinLabel(19,"one Z+#gamma");
   nEvent_4l_w->GetXaxis()->SetBinLabel(20,"two Z+#gamma");
   
   nEvent_4l->SetBinContent(1,N_0);
   nEvent_4l->SetBinContent(2,N_01);
   nEvent_4l->SetBinContent(3,N_02);
   nEvent_4l->SetBinContent(4,N_1);
   nEvent_4l->SetBinContent(5,N_2);
   nEvent_4l->SetBinContent(6,N_3a);
   nEvent_4l->SetBinContent(7,N_3_FSR);
   nEvent_4l->SetBinContent(8,N_3b);
   nEvent_4l_w->SetBinContent(1,N_0_w);
   nEvent_4l_w->SetBinContent(2,N_01_w);
   nEvent_4l_w->SetBinContent(3,N_02_w);
   nEvent_4l_w->SetBinContent(4,N_1_w);
   nEvent_4l_w->SetBinContent(5,N_2_w);
   nEvent_4l_w->SetBinContent(6,N_3a_w);
   nEvent_4l_w->SetBinContent(7,N_3_FSR_w);
   nEvent_4l_w->SetBinContent(8,N_3b_w);
   nEvent_4l_w->SetBinContent(9,N_4a_w);
   nEvent_4l_w->SetBinContent(10,N_4b_w);
   nEvent_4l_w->SetBinContent(11,N_4c_w);
   nEvent_4l_w->SetBinContent(12,N_4d_w);
*/   
   // write on output root file:
   _filePU->Close();
   theFile->cd();
   //z1tree->Write();
   newtree->Write();
   theFile->Write();
   theFile->Close();
   skimfile->cd();
   finaltree->Write();
   skimfile->Write();
   skimfile->Close();
} // end main

double HZZ4LeptonsAnalysis::EAele(int index,bool use2011EA){

  double EffectiveArea=0.;
  if (use2011EA){
    if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.18;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.20;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.15;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.19;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.21;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.22;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.4 )                                           EffectiveArea = 0.29;
  }
  //else { // 7_4_x use eta
  // if (fabs(RECOELE_ETA[index]) >= 0.0   && fabs(RECOELE_ETA[index]) < 0.8 )   EffectiveArea = 0.1830;
  // if (fabs(RECOELE_ETA[index]) >= 0.8   && fabs(RECOELE_ETA[index]) < 1.3 )   EffectiveArea = 0.1734;
  // if (fabs(RECOELE_ETA[index]) >= 1.3   && fabs(RECOELE_ETA[index]) < 2.0 )   EffectiveArea = 0.1077;
  // if (fabs(RECOELE_ETA[index]) >= 2.0   && fabs(RECOELE_ETA[index]) < 2.2 )   EffectiveArea = 0.1565;
  // if (fabs(RECOELE_ETA[index]) >= 2.2 )                                       EffectiveArea = 0.2680;
    //}                                                                                                                                                                             
  else { // 7_6_X use eta supercluster                                                                                                                                             
    if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
  }

  return EffectiveArea;

}

double DELTAPHI( double phi1, double phi2 ){

	if( phi1 > mPI || phi1 < -mPI || phi2 > mPI || phi2 < -mPI) {
	  // cout << "Angles out of range!!! " << endl;
	  // cout << " phi1 " << phi1 << endl;
	  // cout << " phi2 " << phi2 << endl;
	  return -999;
	}
	float dp=std::abs(phi1-phi2);
	if (dp>mPI) dp-=float(2*mPI);
	return dp;
	//return  min( fabs( phi1 - phi2 ) , 2*mPI - fabs( phi1 - phi2 ) ) ;

}

double invmass (float M1, float PT1, float ETA1, float PHI1, float M2, float PT2, float ETA2, float PHI2 ){ 
 float phi1=PHI1; 
 float eta1=ETA1; 
 float pt1=PT1; 
 float m1=M1; 

 float px1=pt1*cos(phi1); 
 float py1=pt1*sin(phi1); 
 float pz1=pt1/(2.*(exp(-1*eta1))/(1.0-exp(-2.*eta1))); 

 float phi2=PHI2; 
 float eta2=ETA2; 
 float pt2=PT2; 
 float m2=M2;

 float px2=pt2*cos(phi2); 
 float py2=pt2*sin(phi2); 
 float pz2=pt2/(2.*(exp(-1*eta2))/(1.0-exp(-2.*eta2))); 

 float e1sqr=pz1*pz1+pt1*pt1+m1*m1; 
 float e2sqr=pz2*pz2+pt2*pt2+m2*m2; 
 float e1e2=sqrt(e1sqr*e2sqr); 
 float p1dotp2=px1*px2+py1*py2+pz1*pz2; 

 float m=sqrt(m1*m1+m2*m2+2.*(e1e2-p1dotp2)); 
 //cout << "Invariant mass= " << m << endl;  
 return m;
} // float invmass closed

double HZZ4LeptonsAnalysis::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){ 

  //	if(Lep.size()!=4 or pterr.size()!=4) {std::cout<<" Lepsize="<<Lep.size()<<", "<<pterr.size()<<std::endl;}
  int debug_ = 0;
  TLorentzVector compositeParticle ;
  for(unsigned int i=0; i<Lep.size(); i++){
    compositeParticle+=Lep[i];
    if(debug_) std::cout<<" in mass error :  add lep  "<<i<<endl;
  }
  double mass  =  compositeParticle.M();
  
  if(debug_) std::cout<<" in mass error :  mass "<<mass<<endl;
  double masserr = 0;
  
  for(unsigned int i=0; i<Lep.size(); i++){
    if(debug_) std::cout<<" in mass error :  varying lep "<<i<<endl;
    TLorentzVector variedLep; // = Lep[i];
    
    if(debug_) std::cout<<" in mass error : pterr = "<<pterr[i]<<endl;
    variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
    TLorentzVector compositeParticleVariation ;
    for(unsigned int j=0; j<Lep.size(); j++){
      if(i!=j)compositeParticleVariation+=Lep[j];
      else compositeParticleVariation+=variedLep;
    }
    
    masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
    if(debug_) std::cout<<" in mass error :  intermediate masserr "<<masserr<<endl;
  }
  return sqrt(masserr);
}

