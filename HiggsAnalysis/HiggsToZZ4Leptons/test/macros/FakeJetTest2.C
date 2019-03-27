#define HZZ4LeptonsAnalysis_cxx
#include "HZZ4LeptonsAnalysis_llbb.h"
#include <TH2.h>
#include <TStyle.h>
//#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TSpline.h>
#include <TRandom3.h>

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

#include "ZZMatrixElement/MELA/src/computeAngles.h"
#include "ZZMatrixElement/MELA/src/computeAngles.cc"
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "RoccoR.cc"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
using namespace std;
// using namespace RooFit;
// using namespace meMCFM;
using namespace MEMNames;
using namespace JME;

const double Zmass = 91.188; // nominal Z boson mass
const double mPI = 3.141592654; 

void HZZ4LeptonsAnalysis::Loop(Char_t *output)
{

   if (fChain == 0) return;

   TString pufile, puhist,ele_leg1,ele_leg2,id_sf,iso_sf,mu17_leg,mu8_leg,mu_track;
   string btagcal;
   string era="GH";

   if(era=="BF"){
     pufile="pileup_BCDEF.root";
     puhist="pileup_scale_BCDEF";
     ele_leg1="Leg1_BF_EGM2D.root";
     ele_leg2="Leg2_BF_EGM2D.root";
     id_sf="IDSF_BCDEF.root";
     iso_sf="ISOSF_BCDEF.root";
     mu17_leg="Mu17Leg_SF_BF.root";
     mu8_leg="Mu8Leg_SF_BF.root";
     mu_track="track_BCDEF.root";
     btagcal="CSVv2_Moriond17_B_F.csv";
     cout << "era BF" << endl;
   }
   if(era=="GH"){
     pufile="pileup.root";
     puhist="pileup_scale";
     ele_leg1="Leg1_GH_EGM2D.root";
     ele_leg2="Leg2_GH_EGM2D.root";
     id_sf="IDSF_GH.root";
     iso_sf="ISOSF_GH.root";
     mu17_leg="Mu17Leg_SF_GH.root";
     mu8_leg="Mu8Leg_SF_GH.root";
     mu_track="track_GH.root";
     btagcal="CSVv2_Moriond17_G_H.csv";
   }


   // Declare MEM class
   MEMs combinedMEM(13,125,"CTEQ6L");     
   
    // JME
   JME::JetParameters jetparameters;
   JME::JetResolution jetresolution;
   JME::JetResolutionScaleFactor jetresolution_sf;

   // BNN
   Char_t datasetChar[500],bnnOUT[500],eventsOUT[500];
  
   cout << "The output file is " << output << endl;
   TString out = output;
   TString datasetName=out.ReplaceAll(".root","");

   sprintf(datasetChar,"%s",datasetName.Data());
   sprintf(bnnOUT,"%s_bnn.txt",datasetName.Data());
   sprintf(eventsOUT,"%s_bnn.root",datasetName.Data());

   // Pileup reweighting in 80x
   TFile *_filePU;
   _filePU= TFile::Open("pileup/"+pufile);
   TH1D *puweight = (TH1D*)_filePU->Get(puhist);

   /////////////Lepton Efficiency Scale Factrons/////////////
   // Load histograms
   //
   TFile *ele_scale_factors_v3 = new TFile("SF_ELE/egammaEffi_RECO_EGM2D.root");
   TH2F *ele_scale_factors_reco = (TH2F*)gDirectory->Get("EGamma_SF2D");
   TFile *ele_scale_factors_v4 = new TFile("SF_ELE/egammaEffi_WP90_EGM2D.root");
   TH2F *ele_scale_factors_wp90 = (TH2F*)gDirectory->Get("EGamma_SF2D");

   TFile *ele_scale_factors_v1 = new TFile("SF_ELE/"+ele_leg1);
   TH2F *ele_scale_factors_leg1 = (TH2F*)gDirectory->Get("EGamma_SF2D");

   TFile *ele_scale_factors_v2 = new TFile("SF_ELE/"+ele_leg2);
   TH2F *ele_scale_factors_leg2 = (TH2F*)gDirectory->Get("EGamma_SF2D");


  TFile *mu_scale_factors3_p2 = new TFile("SF_GH/dm2/"+mu8_leg);
  TH2F *mu_scale_factors_hlt_p2 = (TH2F*)gDirectory->Get("abseta_pt_PLOT");

  TFile *mu_scale_factors1_p1 = new TFile("SF_GH/"+id_sf);
  TH2F *mu_scale_factors_id_p1 = (TH2F*)gDirectory->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

//  TFile *mu_scale_factors1_p1 = new TFile("SF_GH/IDSF_BCDEF.root");
//  TH2F *mu_scale_factors_id_p1 = (TH2F*)gDirectory->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TFile *mu_scale_factors2_p1 = new TFile("SF_GH/"+iso_sf);
  TH2F *mu_scale_factors_iso_p1 = (TH2F*)gDirectory->Get("LooseISO_LooseID_pt_eta/abseta_pt_ratio");

  TFile *mu_scale_factors3_p1 = new TFile("SF_GH/dm2/"+mu17_leg);
  TH2F *mu_scale_factors_hlt_p1 = (TH2F*)gDirectory->Get("abseta_pt_PLOT");

  TFile *mu_scale_factors4 = new TFile("SF_GH/"+mu_track); //just for GH
  TGraph *mu_scale_factors_tk = (TGraph*)gDirectory->Get("ratio_eff_eta3_dr030e030_corr");

  TFile *j_scale_factors3_p1 = new TFile("btag/iso035/fake_jet_QCD.root");
  TH2F *j_eff_p1 = (TH2F*)gDirectory->Get("Jet_9");
  TH2F *j_eff_p3 = (TH2F*)gDirectory->Get("Eta_8");
 
  


   // kfactor_ggZZ(float GENmassZZ, int finalState)     
   TString strSystTitle[9] ={
     "Nominal",
     "PDFScaleDn",
     "PDFScaleUp",
     "QCDScaleDn",
     "QCDScaleUp",
     "AsDn",
     "AsUp",
     "PDFReplicaDn",
     "PDFReplicaUp"
   };

   // Book root file (for output):
   TFile * theFile = new TFile(output,"RECREATE");

   //TString Vars("Weight:Run:Event:LumiSection:massZ1:massZ2:mass4l:Iso_max:Sip_max:MELA:FSR");
   //TNtuple * thePlots=new TNtuple("Candidates","Candidates",Vars);

    // Clone tree for Z1
   //TTree *z1tree = fChain->CloneTree(0);
   
   bool debug=false;
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
   int N_7_PFMET = 0;
    
   int N_8 = 0;
   int N_8_a = 0;
   int N_8_PFMET = 0;
   int N_9 = 0;

   int N_9_1FSR = 0;
   int N_9_2FSR = 0;

   int N_9PS = 0;
   int N_9GRAV = 0;
   
   int N_9a_VBF = 0;
   int N_9b_VBF = 0;
   int N_9_PFMET = 0;

   int N_VBF = 0;
   int N_bjets = 0;
   int N_bjets_cut = 0;
   int N_njets_cut = 0;

   int N_10 = 0;

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
   double N_7_w = 0;
   double N_7_PFMET_w = 0;


   double N_8_w = 0;
   double N_8_a_w = 0;
   double N_8_PFMET_w = 0;
   double N_9_w = 0;

   double N_9_1FSR_w = 0;
   double N_9_2FSR_w = 0;

   double N_9PS_w = 0;
   double N_9GRAV_w = 0;
   
   double N_9a_VBF_w = 0;
   double N_9b_VBF_w = 0;
   double N_9_PFMET_w = 0;
     
   double N_VBF_w = 0;
   double N_bjets_w = 0;
   double N_bjets_cut_w = 0;
   double N_njets_cut_w = 0;

   int N_10_w = 0;
   //******* SENSITIVITY ON MET
   double min_cut_PFMET = 0.; //initialize the min value of RECO_PFMET
   double max_cut_PFMET = 300.;//initialize the max value of RECO_PFMET
   vector <double> cut_PFMET;//declaration of a vector in which the cut values could be stored
   vector <double> counter_cut_PFMET;
   vector <double> counter_cut_PFMET_w;
   
   int cut_n = 61;//number of cuts on PFMET
   double step = (max_cut_PFMET - min_cut_PFMET) / (cut_n - 1);// width between two next cuts
   
   double cut_var_PFMET=0.;
   
   while(cut_var_PFMET <= max_cut_PFMET){
     cut_PFMET.push_back(cut_var_PFMET);//filling the vector with cut's value 
     cut_var_PFMET = cut_var_PFMET + step;//cuts
     
   }
   
   cout << "The size of the vector cut_PFMET is: "<<cut_PFMET.size()<<endl;
   
   //DEBUG
   // for(int i=0; i< cut_PFMET.size(); i++){
   //   cout <<"The value of the cut_PFMET vector at index "<< i <<" is "<<cut_PFMET.at(i)<<endl;
   // }
   //END DEBUG
   
   
   for(int i = 0; i < cut_PFMET.size(); i++){//initialize the counter
     counter_cut_PFMET.push_back(0);
   }
   
   for(int i=0; i< cut_PFMET.size(); i++){
     cout <<"The value of the counter_cut_PFMET vector BEFORE THE LOOP at index "<< i <<" is "<<counter_cut_PFMET.at(i)<<endl;
   }
   
   for(int i = 0; i < cut_PFMET.size(); i++){//initialize the counter
     counter_cut_PFMET_w.push_back(0);
   }
   
   //////////// END OF SENSITIVITY ON MET



   // Book Histos ***
   TH1D *nEvent_4l_w = new TH1D("nEvent_4l_w", "nEventComplete Weighted", 16, 0., 16.);
   TH1D *nEvent_4l = new TH1D("nEvent_4l", "nEventComplete", 16, 0., 16.);

   //SENSITIVITY
   
   TH1D *nEvent_CUT_w = new TH1D("nEvent_CUT_w", "nEventCUT Weightd", 61, 0., 61.);
   TH1D *nEvent_CUT = new TH1D("nEvent_CUT", "nEventCUT", 61, 0., 61.);
   
   //

   // Pileup reweighting
   TH1F *hPUvertices             = new TH1F("hPUvertices", "hPUvertices",70,0.,70.);  
   TH1F *hPUvertices_ReWeighted  = new TH1F("hPUvertices_ReWeighted", "hPUvertices_ReWeighted",70,0.,70.);  

   //step 3

   //no FSR   

   TH1D * hNjets_8 = new TH1D("hNjets_8", "Number of jets passing VBF", 10, -0.5, 9.5);
   hNjets_8->SetXTitle("# n-jets");

   TH1F * hPtJet_9 = new TH1F("hPtJet_9", "Pt of (no ID)jet after selection step 5", 300 ,  0 , 600 );
   hPtJet_9->SetXTitle("pt_jet  (GeV)");

   TH1F * hYJet_9 = new TH1F("hEtaJet_9", "Y of (no ID)jet after selection step 5", 500 , -5. , 5. );
   hYJet_9->SetXTitle("Y of jet");

   TH1F * hPtEle_8 = new TH1F("hPtEle_8", "Pt of ele after selection step 5", 300 ,  0 , 600 );

   TH1F * hEtaEle_8 = new TH1F("hEtaEle_8", "Y of ele jet after selection step 5", 500 , -5. , 5. );



   TH1F * hPtJet_8 = new TH1F("hPtJet_8", "Pt of jet after selection step 5", 300 ,  0 , 600 );
   hPtJet_8->SetXTitle("pt_jet  (GeV)");

   TH1F * hYJet_8 = new TH1F("hEtaJet_8", "Y of jet after selection step 5", 500 , -5. , 5. );
   hYJet_8->SetXTitle("Y of Jet8");

   
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

//Bottom plot
   TH1F * Mjj_6 = new TH1F("Mjj_6","invariant mass of jet pair after step 6",50,20,420);
   Mjj_6->SetXTitle("M_{jj} (GeV)");  

   TH1F * Mee_8 = new TH1F("Mee_8", "Mass jets fake ee", 200, -0.5, 499.5 );
   Mee_8->SetXTitle("Mass jets");

   TH1F * Mee_7 = new TH1F("Mee_7", "Mass ee ", 200, -0.5, 499.5 );
   Mee_7->SetXTitle("Mass jets");

   TH1F * hMZ1_5 = new TH1F("hMZ1_5", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5->SetXTitle("mass_Z1  (GeV)");

   TH1F * hN_loose_e_4 = new TH1F("hN_loose_e_4", "N_loose_e_4", 30 , 0. , 30. );
   hN_loose_e_4->SetXTitle("hN_loose_e_4");

   TH1F *mva_ele = new TH1F("mva_ele","mva_ele",200,-1.,1.);
   mva_ele->SetXTitle("mva_ele");
   // end book histo ***
      
   float newweight=1.;
   
   // New tree with clone of events passing the final selection
   // Clone tree for final events
//   TTree *finaltree = fChain->CloneTree(0);

   // loop on entries
   
   Long64_t nentries = fChain->GetEntries();

   cout << "\n****************************"  <<endl;
   cout << "Analyzing " << nentries << " entries"  <<endl;     

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

 // Initialize reduced tree variables
      //if (!(Run==1 && LumiSection==2536 && Event==486622)) continue;
      
      double mc_weight_un[9];

      if(jentry%5000 ==0) cout << "Analyzing entry: " << jentry << endl;
      

      if( RECO_NMU > 100 ) RECO_NMU = 100;
      if( RECO_NELE > 100 ) RECO_NELE = 100;
      if( RECO_NPFPHOT > 20 ) RECO_NPFPHOT = 20;
      
      bool debug=false;  //debug flag  -- default false
   
      bool zptweight=false;
      bool etaweight=false;
      bool botsf=false;
      if( datasetName.Contains("DYJetsToLL")){
        zptweight=true;
        //cout << "DY correction" << endl;
      }
      newweight=weight;
      if(jentry%5000 ==0) cout << "Starting weight= " << newweight << endl;
  
      // pileup reweighting 2012 and 2011
      if (DATA_type=="NO" && num_PU_vertices < 0) continue;                                                                                                                                              
      // pileup reweighting 2015
      hPUvertices->Fill(num_PU_vertices,weight);
   
      double pu_weight=1.;
      if (MC_type == "Spring16"||MC_type == "Summer16"){
	Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
	if(debug) cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;	
	pu_weight=double(puweight->GetBinContent(binx));
	
      }      
       
      hPUvertices_ReWeighted->Fill(num_PU_vertices,weight*pu_weight);
      if(debug) cout << "Pileup interations and weight is= " << num_PU_vertices << " " << " and weight= " << pu_weight << endl;  
      
      //if (num_PU_vertices < 0) continue;

      // Changing the weight for pileup
      newweight=weight*pu_weight;
      if(jentry%5000 ==0) cout << "Starting weight + pileup = " << newweight << endl;
           
      

      // Weight for MCNLO samples                                                                                      
      if( datasetName.Contains("amcatnlo")) {
        if(debug) cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting << endl;
        newweight=weight*pu_weight*MC_weighting;
      }


      for(int l=0; l<9; l++){
         if(MC_weighting_un[0]!=0) mc_weight_un[l]=MC_weighting_un[l];
         else mc_weight_un[l]=1;
      }
      
      float pFill[11];for(int pf=0;pf<11;pf++)pFill[11]=-999.;

      // ** Step 0:
      // simply number of entries...
      if( debug ) cout << "\n** Step 0: \nAnalyzing entry: " << jentry << " Run: " << Run << " Event: " << Event << " LumiSection: " << LumiSection << endl ;
      ++N_0 ;  // fill counter
      N_0_w=N_0_w+newweight;
      
      // ** Step 0.1:
      // number of 4L (4mu)...
      //trigger requirements (qier)
      bool is2e2mu=false;

      ++N_1 ;  // fill counter
      N_1_w=N_1_w+newweight;
      
      
      // Effective AREA
      bool tag_2011=false;
      if (DATA_type=="2010" || DATA_type=="2011" || MC_type=="Fall11"){
        tag_2011=true;
      }
      
      
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
       	
 	if(/* ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0) )
	    && RECOMU_mubesttrkType[i]!=2
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.*/
              ( RECOMU_isGlobalMu[i] || RECOMU_isTrackerMu[i] ) && RECOMU_isPFMu[i]
              && RECOMU_PT[i] > 5.
              && fabs(RECOMU_ETA[i]) < 2.4
              && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1. //loose muon (qier)
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
       	
 	if( RECOELE_PT[i] > 10. 
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
	  
	  if(/*(RECOMU_isPFMu[mu] || (RECOMU_isTrackerHighPtMu[mu] && RECOMU_PT[mu] > 200.))
	      && (RECOMU_isGlobalMu[mu] || (RECOMU_isTrackerMu[mu] && RECOMU_numberOfMatches[mu]>0))
	      && RECOMU_mubesttrkType[mu]!=2
	      && RECOMU_PT[mu] > 5. 
	      && fabs(RECOMU_ETA[mu]) < 2.4 
	      && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. 
	      && fabs(RECOMU_SIP[mu])<4. // TightID + SIP cut*/
 //              RECOMU_isMedium[mu] &&//mediumID(qier) 
              ( RECOMU_isGlobalMu[mu] || RECOMU_isTrackerMu[mu] ) && RECOMU_isPFMu[mu] 
              && RECOMU_PT[mu] > 5. 
              && fabs(RECOMU_ETA[mu]) < 2.4 
              && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. //loose muon (qier)
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
 	if( /*(RECOMU_isPFMu[i] || (RECOMU_isTrackerHighPtMu[i] && RECOMU_PT[i] > 200.))
	    && ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0))
	    && RECOMU_mubesttrkType[i]!=2	 
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.//tight muon */
 //             RECOMU_isMedium[i] &&//mediumID(qier)
              ( RECOMU_isGlobalMu[i] || RECOMU_isTrackerMu[i] ) && RECOMU_isPFMu[i]
              && RECOMU_PT[i] > 5.
              && fabs(RECOMU_ETA[i]) < 2.4
              && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1. //loose muon (qier)
        
	    ){

          Double_t Pt = RECOMU_PT[i];
          Double_t Eta = RECOMU_ETA[i];
          Int_t Q = int(RECOMU_CHARGE[i]);
          Double_t Phi = RECOMU_PHI[i];
          Double_t nl = RECOMU_mutrktrackerLayersWithMeasurement[i];
          }

	  iL[ N_good ] = i ;
	  ++N_good ;	  
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
      

    
      
      /// *** FSR
      // Photon identification & cleaning
      // ele identification is also needed
      
      //FSR photon identifications, will be used with MELA later
      int FSR_Z1_photid=-1;
      int FSR_Z2_photid=-1;
      int FSR_Z1_lepid=-1;
      int FSR_Z2_lepid=-1;

      //electrons:
      
      int Ne_good = 0 ;
      int iLe[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1}; //electrons
      int Ne_loose_4=0;


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
       	
 	if( RECOELE_PT[i] > 10. && fabs(RECOELE_ETA[i]) < 2.5 );
	  // && RECOELE_gsftrack_expected_inner_hits[i]<=1 ) /* ok */ ;
	else continue ;

        Ne_loose_4++;
        mva_ele->Fill(RECOELE_mvaNonTrigV0[i],newweight);

	bool BDT_ok = 0; // Spring16 with CMSSW_8_0_x
/*
	if( RECOELE_PT[i] > 7. &&  RECOELE_PT[i] <= 10. ){
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.211 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) < 1.479 )
						 && RECOELE_mvaNonTrigV0[i] > -0.396 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) >= 1.479 && RECOELE_mvaNonTrigV0[i] > -0.215 ) BDT_ok = 1 ;
	}
	else { 
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.870 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
						 && RECOELE_mvaNonTrigV0[i] > -0.838 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > -0.763 ) BDT_ok = 1 ;
	}	
*/
        //qier chenge to wp90
        if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaTrigV0[i] > 0.837 ) BDT_ok = 1 ;
        if( fabs(RECOELE_scl_Eta[i]) >= .8 &&fabs(RECOELE_scl_Eta[i])<1.479 && RECOELE_mvaTrigV0[i] > 0.715 ) BDT_ok = 1 ;
        if(fabs(RECOELE_scl_Eta[i])>=1.479 && RECOELE_mvaTrigV0[i] > 0.357) BDT_ok=1;
        
	if( !BDT_ok ) continue ;
	
	if( fabs(RECOELE_gsftrack_dxy[i]) < .5 
	 && fabs(RECOELE_gsftrack_dz[i])  < 1. ) /* ok */ ;
	else continue ; 
//	if (RECOELE_PFX_rho_new[iLe[i]]>=0.20) continue;
	iLe[ Ne_good ] = i ;
	++Ne_good ;

      }// end loop on electrons
      hN_loose_e_4->Fill(Ne_loose_4,newweight);

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
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 0;
	  }
	  
	}//end loop on muons  
	
	for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
	  if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[l]],2) );
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
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
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOELE_ETA[iL_loose_e[l]],2));
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    if (deltaR_ET2<min_deltaR_ET2){
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
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

      

       // Exclude that photon from the isolation cone all leptons in the event passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons
      
      if(debug) cout << "Rho for electron pileup isolation correction is= " << RHO_ele << endl;
      double EffectiveArea=-9999.;

	    
      for(int i=0.;i<Nphotons;i++) {
	if (iLp_l[i]==-1) continue;
	
	for(int e = 0; e < N_loose_e; ++e){
	  //if(!( iLp_l[i] == iL_loose_e[e] && iLp_tagEM[i] == 1 ) ) continue;
	  if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;
	  //double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[e]],2) );
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
      if( debug ) cout  << "\nStep 3: Number of good leptons: " << N_good+Ne_good << endl;
 
      if( Ne_good < 1 ) continue ; 	
      ++N_2 ;  // fill counter
      N_2_w=N_2_w+newweight;

      int Zxx_tag = 0;    // 1: Zmumu  ,  2: Zee

      int i1 = -1; //index of the first lepton (from Z1)
      int j1 = -1; //index of the second lepton (from Z1)
      int pi1 = -1; 
      int pj1 = -1;
      

      
      ++N_3a ;  // fill counter
      N_3a_w=N_3a_w+newweight;
      

      // Mass cut on Z


      ++N_3b ;  // fill counter
      N_3b_w=N_3b_w+newweight;

     


      // **** Step 4:
       // a) 4 leptons
      // b) pair #2
      // c) highest pt
      // d) mZ2 in ] 4,120 [

      int issamesign = 0;

      //if( debug ) cout  << "\nStep 4: Number of good leptons: " << N_good << endl;

      int N_Z2_pairs = 0;

      int i2 = -1; //index of the first lepton (from Z1)
      int j2 = -1; //index of the second lepton (from Z1)
      int pi2 = -1; 
      int pj2 = -1; 
      
      bool has_FSR_Z2 = 0;


      
      
      
      
      ++N_4b ;  // fill counter
      N_4b_w=N_4b_w+newweight;

      
      // **** Step 5:

      
       // Execute Efficiency Reweighting
 
      // N.B. DO NOT Update the Isolation values and correct the 4 momenta of leptons for FSR
      for(int i = 0; i < Ne_good; ++i){
	int flagFSR=0;
	int pfsr=-999;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if (iLp[p]==-1) continue;
	  if (iLp_l[p]==-1) continue;
	  
	  if(debug) cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << iLe[i] << endl;
	  if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {
	    if(debug) cout << "Electron with pT= " << RECOELE_PT[iLe[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    // RECOELE_PFX_rho_new[iLe[i]]=
	    //   (RECOELE_PFchHad[iLe[i]]+
	    //    max(0.,RECOELE_PFneuHad[iLe[i]]+
	    // 	   (RECOELE_PFphoton[iLe[i]]-RECOPFPHOT_PT[iLp[p]] )-
	    // 	   max(RHO_ele,0.0)*(EffectiveArea)))/RECOELE_PT[iLe[i]];	    
	    flagFSR=1;
	    pfsr=p;
	  }
	}
	
	if (flagFSR==1){
	  if(debug) cout << "Before correcting for FSR; electron pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	  TLorentzVector Lept,LeptCorrection;
	  Lept.SetPtEtaPhiM(RECOELE_PT[iLe[i]], RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]], 0.105);
	  LeptCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_ETA[iLp[pfsr]],RECOPFPHOT_PHI[iLp[pfsr]],0);
	  Lept+=LeptCorrection;
	  RECOELE_PT[iLe[i]]=Lept.Pt();
	  RECOELE_ETA[iLe[i]]=Lept.Eta();
	  RECOELE_PHI[iLe[i]]=Lept.Phi();
	  if(debug) cout << "After correcting for FSR; muon pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	}
      }
      
      if(debug) cout << "Kinematics of leptons corrected for FSR photons (if existing)" << endl;
      
      
      // // **** Step 6:
      //  // QCD suppression: mll>4 GeV cut on all OS-SF pairs (4/4)           
     //if( min_mass_2L <= 4 ) continue ;
     
     ++N_6 ;  // fill counter
     N_6_w=N_6_w+newweight;

     // **** Step 7:
     // mass4l > 70 
     if(Ne_good<1) continue;
     //Basic cuts to jets AND delta R section
     int njets_pass=0;
     bool lead_jet=false;
     TLorentzVector JET0,JET1,JET2,ELE0,ELE1,ELE2;
     int jet1=-999,jet2=-999,ele1=-999,ele2=-999;      
     int jetfail[100];
     float GOOD_JET_PT_MAX = 0.;
     float JET_PHI_PT_MAX = 0.;
     float max_dphi_jet_met = 0.;
     float min_dphi_jet_met = 999.;
     vector <int>  v_good_jets_index;
     ELE0.SetPtEtaPhiE(RECOELE_PT[iLe[0]],RECOELE_ETA[iLe[0]],RECOELE_PHI[iLe[0]],RECOELE_PT[iLe[0]]*TMath::CosH(RECOELE_ETA[iLe[0]]));
     if(RECOELE_PT[iLe[0]]<30) continue;
     if(!RECOELE_se_EleHLTMatch[iLe[0]]) continue;

     for(int i=0;i<100;i++) jetfail[i]=0;
     
     for(int i=0;i<RECO_PFJET_N;i++){

       if(RECO_PFJET_PT[i]<-100) continue;     

       
       double Pt=RECO_PFJET_PT[i];
       double Eta=RECO_PFJET_ETA[i];

       if(debug) cout<<i<<" Jet with pt= "<<RECO_PFJET_PT[i]<<" ETA "<<RECO_PFJET_ETA[i]<<" PUID "<<RECO_PFJET_PUID[i] << " PUID_MVA "<< RECO_PFJET_PUID_MVA[i]<<endl;
 
       float dphi_jet_met=0.;
       dphi_jet_met=RECO_PFJET_PHI[i]-RECO_PFMET_PHI;


       bool goodjet = RECO_PFJET_NHF[i] < 0.99 &&
                      RECO_PFJET_NEF[i] < 0.99 &&
                      RECO_PFJET_CHF[i] < 0.99 &&
                      RECO_PFJET_CEF[i] < 0.99 &&
                      RECO_PFJET_nconstituents[i] > 1 &&
                      RECO_PFJET_NCH[i] > 0;
       
 
       if(RECO_PFJET_PT[i]>10. && fabs(RECO_PFJET_ETA[i])<2.4 /*&& goodjet==1*/){
       
      	 //Check that jet has deltaR>0.4 away from any tight lepton corrected for FSR
	 for(int mu = 0; mu < N_good; ++mu){
//	   if (fabs(RECOMU_SIP[iL[mu]])>=4.) continue;  //qier
    	   if (RECOMU_PFX_dB_new[iL[mu]]>=0.20) continue;
	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOMU_PHI[iL[mu]]),2) + pow(RECO_PFJET_ETA[i] - RECOMU_ETA[iL[mu]],2));
	   //cout << "1st lepton muon: " << " pT=" << RECOMU_PT[iL[mu]] <<" deltaR "<< deltaR <<endl;	   
	   if (deltaR<0.4){
	     jetfail[i]=1;
     	     //cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
      	 for(int ele = 0; ele < Ne_good; ++ele){
      	//   if (fabs(RECOELE_SIP[iLe[ele]])>=4.) continue;
	   if (RECOELE_PFX_rho_new[iLe[ele]]>=0.35) continue;
      	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOELE_PHI[iLe[ele]]),2) + pow(RECO_PFJET_ETA[i] - RECOELE_ETA[iLe[ele]],2));
     	   //cout << "1st lepton electron: " << " pT=" << RECOELE_PT[iLe[ele]] <<" deltaR "<< deltaR <<endl;
	   if (deltaR<0.4){
     	     jetfail[i]=1;
     	     //cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
	 // cleaning w.r.t FSR photons attached to leptons
	 for(int j=0.;j<Nphotons;j++) {
           if (iLp_l[j]!=-1 && (iLp_tagEM[j]==0 || iLp_tagEM[j]==1) ) {
	     double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOPFPHOT_PHI[iLp[j]]),2) + pow(RECO_PFJET_ETA[i] - RECOPFPHOT_ETA[iLp[j]],2));
	     if (deltaR<0.4){
	       jetfail[i]=1;
	       //cout << " jetfail " << jetfail[i] <<endl;
	       break;
	     }
	   }
         }
	 // 
	 if (jetfail[i]==0){
	   //cout<< " PASS jet " <<i<<" PT= "<<RECO_PFJET_PT[i]<<" ETA= "<<RECO_PFJET_ETA[i]<<" PUID= "<<RECO_PFJET_PUID[i]<<endl;
	   njets_pass++;
           if(RECO_PFJET_PT[i]>50) lead_jet = true;
           hPtJet_8->Fill(RECO_PFJET_PT[i],newweight);
           hYJet_8->Fill(RECO_PFJET_ETA[i],newweight);

           int biny = j_eff_p1->GetYaxis()->FindBin(RECO_PFJET_PT[i]);
           int binx = j_eff_p1->GetXaxis()->FindBin(RECO_PFJET_ETA[i]);
           double scale_fkj=j_eff_p1->GetBinContent(binx,biny);

           hPtJet_9->Fill(RECO_PFJET_PT[i],newweight*scale_fkj);
           hYJet_9->Fill(RECO_PFJET_ETA[i],newweight*scale_fkj);
     
           JET0.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i])); 
           double mee=(ELE0+JET0).M();
           Mee_8->Fill(mee,newweight*scale_fkj);          
           hMZ1_5->Fill(mee,newweight*scale_fkj);
 
	   if (njets_pass==1){
	     jet1=i;
	     JET1.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	     GOOD_JET_PT_MAX=RECO_PFJET_PT[i];
	     JET_PHI_PT_MAX=RECO_PFJET_PHI[i];
	     //cout<<"Among the jets that pass the jet with the highet pt is the jet of index "<< i <<". It has pt "<<GOOD_JET_PT_MAX<<". The corresponding value of phi is " << JET_PHI_PT_MAX <<endl;
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
     }


     hNjets_8->Fill(njets_pass,newweight);
     hN_good_ele->Fill( Ne_good,newweight );
 

     int nele=0;
     for(int i = 1; i < Ne_good; ++i){
          if (fabs(RECOELE_PFX_rho_new[iLe[i]])>=0.35) continue; // Isolation cut
          if (RECOELE_PT[iLe[i]]<10) continue;
          nele++;

             hPtEle_8->Fill(RECOELE_PT[ iLe[i] ],newweight);
             hEtaEle_8->Fill(RECOELE_ETA[ iLe[i] ],newweight);

           if (nele==1){
             ele2=iLe[i];
             ELE2.SetPtEtaPhiE(RECOELE_PT[iLe[i]],RECOELE_ETA[iLe[i]],RECOELE_PHI[iLe[i]],RECOELE_PT[iLe[i]]*TMath::CosH(RECOELE_ETA[iLe[i]]));
           }
     }
     if(nele>=1){
       double Mee=(ELE0+ELE2).M();
       Mee_7->Fill(Mee,newweight); 
     }
 
   } // end loop on entries
   
   // write on output root file:
   _filePU->Close();
   theFile->cd();
   theFile->Write();
   theFile->Close();
  // finaltree->Write();
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


