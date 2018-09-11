#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "ZZStyle.C"
#include "TFile.h"
#include "TColor.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMath.h"
#include "TSystem.h"
#include <libgen.h>

#define nnn 1.0
using namespace std;

// Usage:
// .include /cmshome/nicola/tmp/test/Paper/last/last/CMSSW_5_3_9/src
//  gSystem->Load("libHiggsHiggs_CS_and_Width.so")
// .L PlotStack4l.C+
// PlotStack4l()

class PlotStack4l2b{
  
public: 
  PlotStack4l2b();
  void plotmZ(std::string);
  void setSamplesNamesZ();
  void printnumbers(char*, TH1F*);
  /*void createdatacards(
		       float Higgsm, float channel, float energy, 
		       float masslow, float masshigh,
		       float ggH, float qqH, float WH, float ZH, float ttH, 
		       float bkg_qqzz, float bkg_ggzz, float bkg_zjets
		       );*/
  //void getMassWindow(float Higgsm);

private:
  std::vector<string> Vdatasetnamebkg,Vdatasetnamesig,Vdatasetnamedata,Vdatasetnamebkgdata;
  std::vector<string> Vlabelbkg,Vlabelsig,Vlabeldata,Vlabelbkgdata;
  std::vector<float> Vxsectionbkg,Vxsectionsig,Vxsectiondata,Vxsectionbkgdata;
  std::vector<Color_t> Vcolorbkg, Vcolorsig/*, Vcolordata*/;

  double Nbins;
  double Xmin;
  double Xmax;
  int nRebin,nRebinZ_X;
  double Ymax;
  bool useLogY,useLogX;
  bool useDYJets,useDYJetsFromData;
  bool un_mc;
  std::string histosdir; 
  std::string inputfile;
  std::string whichchannel,whichenergy,whichsample;
  char histotitle[500];
  ofstream outputyields;
  TSystem LoadLib;
public:
  //float Higgsm, channel, energy, masslow, masshigh;
  //float ggH, qqH, WH, ZH, ttH, bkg_qqzz, bkg_ggzz, bkg_zjets;

};

PlotStack4l2b::PlotStack4l2b(){
  //TSystem LoadLib;
  //LoadLib.Load("/cmshome/nicola/slc6/MonoHiggs/Analysis13TeV/CMSSW_7_2_0/lib/slc6_amd64_gcc481/libHiggsHiggs_CS_and_Width.so");
  //getMassWindow(500.);
    
  inputfile="filelist_4mu2b_sig.txt";

  setSamplesNamesZ(); 
  cout << "\t Analysing samples for " << whichchannel << " analysis" << endl; 


  //WARNING: depending on histolabel, modify the declaration and the settings of hframe below
  //also choose a sensible value for nRebin

  //std::string histolabel = "hPUvertices";    // numPU
  //std::string histolabel = "hPUvertices_ReWeighted";    // numPY reweighted

  //std::string histolabel = "hPFMET_3"; // PFMET

  // ****** After cuts on ,Z1, mZ2 and pT >20/10: step 5 ******

//  std::string histolabel = "hMZ1_5";    // Z1 mass   
//  std::string histolabel = "hMZ2_5";    // Z2 mass 

  //std::string histolabel = "hIso_5";    // worst isolation 
  //std::string histolabel = "hSip_5";  // worst sip 
  //std::string histolabel = "hIp_5";   // worst IP

//   std::string histolabel = "Mbb_6";
//  std::string histolabel = "bdiscr_5_lead";
//  std::string histolabel = "bdiscr_5_sub";
//    std::string histolabel = "Mjj_6";

//  std::string histolabel = "hMZ1_6";    // Z1 mass   
//  std::string histolabel = "hMZ2_6";    // Z2 mass 
//  std::string histolabel = "hMZ1_7";
//  std::string histolabel = "hMZ2_7";    // Z2 mass
  // After full selection
  //std::string histolabel = "hMZ1_noFSR_8"; // Z1 mass after full selection without FSR recovery
  //std::string histolabel = "hMZ2_noFSR_8"; // Z2 mass after full selection without FSR recovery
//  std::string histolabel = "hPtZ1_5"; // Z1 pt after full selection
//    std::string histolabel = "hPtZ2_5";
//   std::string histolabel = "hPtLep4_7";
//   std::string histolabel = "hPtLep2_7"; 
// std::string histolabel = "hPtLep1_8";
//   std::string histolabel = "hPtLep2_8"; 

//   std::string histolabel = "hEtaLep1_7";
//   std::string histolabel = "hEtaLep4_7";
//   std::string histolabel = "hIsoLep1_7";
//   std::string histolabel = "hIsoLep2_7";
//    std::string histolabel = "hYZ1_5";
//    std::string histolabel = "hYZ2_5";
  //std::string histolabel = "hPtZ2_8"; // Z2 pt after full selection
  //std::string histolabel = "hYZ1_8"; // Z1 rapidity after full selection
  //std::string histolabel = "hYZ2_8"; // Z2 rapidity after full selection
  //std::string histolabel = "hIso_8"; // worst isolated lepton: isolation value after full selection
  //std::string histolabel = "hSip_8"; // worst sip lepton: sip value after full selection
  //std::string histolabel = "hMELA_8"; // MELA discriminant after full selection 
//   std::string histolabel = "hPFMET_8"; // PFMET
  //std::string histolabel = "hM4l_T_8"; // Transverse mass
  //std::string histolabel = "DPHI_8"; // DeltaPhi - 4l + MET
//   std::string histolabel = "hPtJet_7"; 
//   std::string histolabel = "hPtJet_8";
//   std::string histolabel = "hEtaJet_7";
//   std::string histolabel = "hEtaJet_8";
//  std::string histolabel = "hNjets_8";
  //std::string histolabel = "hDjj_8"; // delta eta between jets for VBF analysis
 // std::string histolabel = "hMjj_8"; // dimass between jets for VBF analysis
//  std::string histolabel = "hN_loose_e";
//  std::string histolabel = "hN_loose_mu"; 
//   std::string histolabel = "hN_good_ele";
//   std::string histolabel = "hN_good_mu";
   std::string histolabel="Mllllbb";
  

  std::string label5[11]={"hMZ1_5","hMZ2_5","hMZ1_6","hMZ2_6","hMZ1_7","hMZ2_7","Mbb_6","Mjj_6","hPtJet_8","hPtBot_8","hPFMET_8"};
//  std::string label20[9]={"hPtZ1_5","hPtZ2_5","hYZ1_5","hYZ2_5","hPtZ1_6","hPtZ2_6","hYZ1_6","hYZ2_6","hEtaJet_8"};
  std::string label20[3]={"hPtZ1_7","hPtZ2_7","ptbb_6"};
  std::string control5[5]={"hMZ1_5","hMZ2_5","hMZ1_6","hMZ2_6","Mjj_6"};
  std::string control20[8]={"hPtZ1_5","hPtZ2_5","hYZ1_5","hYZ2_5","hPtZ1_6","hPtZ2_6","hYZ1_6","hYZ2_6"};
  std::string label1[2]= {"hNjets_8","hNbjets"};
  std::string label55[12]={"hMZee_5","hMZmm_5","hMZee_6","hMZmm_6","hPtZee_5","hPtZmm_5","hYZee_5","hYZmm_5","hPtZee_6","hPtZmm_6","hYZee_6","hYZmm_6"};


  useLogY = false;
  useLogX = false;

  useDYJets=true;
  useDYJetsFromData=false;
 
  un_mc=true;
 
  nRebin=5;
  
  // Final yields
  system("mkdir plots");

  // Execute the analysis
  plotmZ(histolabel);
//  for(int i=2; i<3; i++){
//    plotmZ(label20[i]);
//  } 
 

}

void PlotStack4l2b::plotmZ(std::string histlabel){
  
  TStyle * style = getStyle("ZZ");
  //style->SetMarkerSize(0.8);
  style->cd();
  style->SetNdivisions(508, "X");
  style->SetNdivisions(508, "Y");
  style->SetMarkerSize(0.8);

//   double ytitleOffset = 1.36;
//   double xtitleOffset = 1.18;
//   double labelSize = 0.05;
//   double titleSize = 0.05;
//   double lineWidth = 2;
  
  TCanvas *c1 = new TCanvas("c1","c1",600,800);  
  c1->cd();
  c1->SetTicks(1,1);

  //  ratioPad->cd();

  // TString lumist="5.2 fb^{-1}";
//   TPaveText *ll = new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC");
//   ll->SetTextSize(0.03);
//   ll->SetTextFont(42);
//   ll->SetFillColor(0);
//   ll->SetBorderSize(0);
//   ll->SetMargin(0.01);
//   ll->SetTextAlign(12); // align left
//   TString text = "CMS Preliminary 2011";
//   ll->AddText(0.01,0.5,text);
//   text = "#sqrt{s} = 7 TeV  L = ";
//   text = text + lumist;
//   //  ll->SetTextAlign(32); // align right
//   ll->AddText(0.65, 0.6, text);
 

  TPaveText *ll = new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC");
  ll->SetTextSize(0.027);
  ll->SetTextFont(42);
  ll->SetFillColor(0);
  ll->SetBorderSize(0);
  ll->SetMargin(0.01);
  ll->SetTextAlign(12); // align left
  TString text = "CMS Preliminary";
  //TString text = "CMS";
  ll->AddText(0.01,0.5,text);
  cout << "Energy= " << whichenergy << endl;
  text = "            35.8 fb^{-1} (13TeV)" ;
  //text = "#sqrt{s} = 13 TeV, L = 14.77 fb^{-1}" ;
  ll->AddText(0.65, 0.6, text);
  //ll->Draw();



  TLegend *leg0 = new TLegend(0.6,0.40,0.8,0.90,NULL,"brNDC");
  leg0->SetTextSize(0.020);
  leg0->SetLineColor(0);
  leg0->SetLineWidth(1);
  leg0->SetFillColor(kWhite);
  leg0->SetBorderSize(0);
 
  TLegend* legend = new TLegend( 0.70, 0.76, 0.9, 0.9);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.020);


//   TLegend *leg1 = new TLegend(0.25,0.80,0.4,0.9,NULL,"brNDC");
//   leg1->SetTextSize(0.030);
//   leg1->SetLineColor(0);
//   leg1->SetLineWidth(1);
//   leg1->SetFillColor(kWhite);
//   leg1->SetBorderSize(0);
  
//   Color_t redBgColor = kGreen-5;
  Color_t ZZBgColor = kAzure-9;
//   Color_t h200Color = kOrange;
//   Color_t h350Color = kRed+1;
//   Color_t h400Color = kRed;
  



  if(useLogY) c1->SetLogy(1);
  else        c1->SetLogy(0);

  if(useLogX) c1->SetLogx(1);
  else        c1->SetLogx(0);

  
  cout << "Vdatasetnamedata " <<  Vdatasetnamedata.size() << endl;

  TFile *ff=NULL; 
 
  if(Vdatasetnamedata.size()>0)
    ff = TFile::Open(Vdatasetnamedata.at(0).c_str());  // just a random file, to determine the binning
  else if(Vdatasetnamebkg.size()>0)  
    ff = TFile::Open(Vdatasetnamebkg.at(0).c_str());  
  else if(Vdatasetnamesig.size()>0)  
    ff = TFile::Open(Vdatasetnamesig.at(0).c_str());  

  
  TH1F *hfourlepbestmass_4l_afterSel_new = (TH1F*)ff->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);

  Nbins = hfourlepbestmass_4l_afterSel_new->GetNbinsX() / nRebin;
  Xmin  = hfourlepbestmass_4l_afterSel_new->GetXaxis()->GetXmin();
  Xmax  = hfourlepbestmass_4l_afterSel_new->GetXaxis()->GetXmax() ;
  Ymax  = hfourlepbestmass_4l_afterSel_new->GetBinContent(hfourlepbestmass_4l_afterSel_new->GetMaximumBin()) * 580.;

  cout << "Ymax = " << Ymax << endl;
  TH2F *hframe=NULL,*hframe2=NULL;
 
  hframe= new TH2F("hframe","hframe",80,70.,1000.,500,5.,700.);// 4l analysis mass nrebin=10 GeV
  hframe2= new TH2F("hframe2","hframe2",80,70.,1000.,1000, 0.5, 20.);// 4l analysis mass


  if (histlabel.find("hPtLep_0")<10){
    hframe= new TH2F("hframe","hframe",80,5.,200.,500,0.01,10000000.);// pT                                                                 
    hframe2= new TH2F("hframe2","hframe2",80,5.,200.,500, 0.5, 2.);// pT                                                                                                              
  }

  if (histlabel.find("hIsoLep_0")<10){
    hframe= new TH2F("hframe","hframe",80,0.,10.,500,0.01,10000000.);// Isolation                                                                                                 
    hframe2= new TH2F("hframe2","hframe2",80,0.,10.,500, 0.5, 2.);// Isolation                                                                                                                  
  }
  
  if (histlabel.find("hTKIsoLep_0")<10){
    hframe= new TH2F("hframe","hframe",80,0.,10.,500,0.01,10000000.);// Tracker Isolation                                                                                                 
    hframe2= new TH2F("hframe2","hframe2",80,0.,10.,500, 0.5, 2.);// Tk Isolation                                                                                                                  
  }

  if (histlabel.find("hMZ1_5")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",600, 60., 120., 1000, 0.1, 1.9);// mZ1
    hframe2->SetXTitle("M_{ll_{1}} [GeV]");
  }

  if (histlabel.find("hPtZ1_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ll_{1}} [GeV]");
  }
  if (histlabel.find("hPtZ2_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ll_{2}} [GeV]");
  }

  if (histlabel.find("hPtZ1_6")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ll_{1}} [GeV]");
  }
  if (histlabel.find("hPtZ2_6")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ll_{2}} [GeV]");
  }

  if (histlabel.find("hMZee_5")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",600, 60., 120., 1000, 0.1, 1.9);// mZ1
    hframe2->SetXTitle("M_{ee} [GeV]");
  }

  if (histlabel.find("hMZmm_5")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",600, 60., 120., 1000, 0.1, 1.9);// mZ1
    hframe2->SetXTitle("M_{#mu#mu} [GeV]");
  }

  if (histlabel.find("hPtZee_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ee} [GeV]");
  }
  if (histlabel.find("hPtZmm_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{#mu#mu} [GeV]");
  }

  if (histlabel.find("hPtZee_6")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ee} [GeV]");
  }
  if (histlabel.find("hPtZmm_6")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,20000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{#mu#mu} [GeV]");
  }


  if (histlabel.find("hPtLep1_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,0.1,2000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Pt_{l_{1}} [GeV]");
  }

  if (histlabel.find("hPtLep2_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,0.1,200.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Pt_{l_{2}} [GeV]");
  }
  if (histlabel.find("hPtLep3_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,0.1,2000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Pt_{l_{3}} [GeV]");
  }

  if (histlabel.find("hPtLep4_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,0.1,400.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Pt_{l_{4}} [GeV]");
  }

  if (histlabel.find("hPtLep1_8")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hPtLep2_8")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 2.);// mZ1 
  }

  if (histlabel.find("hYZ1_5")<10 || histlabel.find("hYZee_5")<10 || histlabel.find("hYZmm_5")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,0.01,30000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Y_{ll_{1}} [GeV]");
  }

  if (histlabel.find("hYZ2_5")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,0.01,30000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Y_{ll_{2}} [GeV]");
  }

  if (histlabel.find("hYZ1_6")<10 || histlabel.find("hYZee_6")<10 || histlabel.find("hYZmm_6")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,0.01,30000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Y_{ll_{1}} [GeV]");
  }

  if (histlabel.find("hYZ2_6")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,0.01,30000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Y_{ll_{2}} [GeV]");
  }

  if (histlabel.find("hEtaLep1_7")<10 || histlabel.find("hEtaLep2_7")<10 || histlabel.find("hEtaLep3_7")<10 || histlabel.find("hEtaLep4_7")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,0.1,400.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Eta_{l_{1}} [GeV]");
  }

  if (histlabel.find("hIsoLep1_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,1.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 1., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Eta_{l_{2}} [GeV]");
  }

  if (histlabel.find("hIsoLep2_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,1.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 1., 1000, 0.5, 1.5);// mZ1 
  }
  if (histlabel.find("hMZ1_6")<10 || histlabel.find("hMZee_6")<10 ||histlabel.find("hMZmm_6")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,10000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.1, 1.9);// mZ1 
  }
  if (histlabel.find("hMZ2_6")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,10000.);// mZ2 
    hframe2= new TH2F("hframe2","hframe2",600, 60., 120., 1000, 0.1, 1.9);// mZ2 
  }

  if (histlabel.find("hMZ1_7")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,1000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("M_{ll_{1}} [GeV]");
  }
  if (histlabel.find("hMZ2_7")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,1000.);// mZ2 
    hframe2= new TH2F("hframe2","hframe2",600, 60., 120., 1000, 0.1, 1.9);// mZ2 
    hframe2->SetXTitle("M_{ll_{2}} [GeV]");
  }
  if (histlabel.find("ptbb_6")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,1000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 200., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Pt_{bb} [GeV]");
  }
  if (histlabel.find("hPtZ1_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,2000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ll_{1}} [GeV]");
  }
  if (histlabel.find("hPtZ2_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,0.01,2000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("Pt_{ll_{2}} [GeV]");
  }
  if (histlabel.find("Mbb_6")<10){
    hframe= new TH2F("hframe","hframe",190,20.,400.,500,0.001,100.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",190, 20., 400., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("M_{bb} [GeV]");
  }
  if (histlabel.find("Mllllbb")<10){
    hframe= new TH2F("hframe","hframe",150,0.,1500.,500,0.001,1.5);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",150, 0., 1500., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("M_{4l2b} [GeV]");
  }

  if (histlabel.find("Mjj_6")<10){
    hframe= new TH2F("hframe","hframe",80,20.,400.,500,0.01,100.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",80, 20., 400., 1000, 0.1, 1.9);// mZ1 
    hframe2->SetXTitle("M_{jj} [GeV]");
  }


  if (histlabel.find("bdiscr_5_lead")<10){
    hframe= new TH2F("hframe","hframe",20,0.,1.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",20, 0., 1., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("bdiscr_5_sub")<10){
    hframe= new TH2F("hframe","hframe",20,0.,1.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",20, 0., 1., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hMZ2_5")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,0.01,20000.);// mZ2 
    hframe2= new TH2F("hframe2","hframe2",600, 60., 120., 1000, 0.1, 1.9);// mZ2
    hframe2->SetXTitle("M_{ll_{2}} [GeV]"); 
  }

  if (histlabel.find("hPtJet_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,400.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",1000, 0., 400., 1000, 0.1, 1.9);// mZ1 
  }
 
  if (histlabel.find("hPtJet_8")<10){
    hframe= new TH2F("hframe","hframe",50,0.,200.,500,0.1,1000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 200., 1000, 0.1, 1.9);// mZ1
    hframe2->SetXTitle("Pt_{j} [GeV]"); 
  }

  if (histlabel.find("hPtBot_8")<10){
    hframe= new TH2F("hframe","hframe",50,0.,200.,500,0.1,1000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 200., 1000, 0.1, 1.9);// mZ1
    hframe2->SetXTitle("Pt_{b} [GeV]");
  }
 
  if (histlabel.find("hEtaJet_7")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.1, 1.9);// mZ1 
  }

  if (histlabel.find("hEtaJet_8")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,0.01,1000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.1, 1.9);// mZ1 
  }

  if (histlabel.find("hEtaBot_8")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,0.01,1000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.1, 1.9);// mZ1 
  }


  if (histlabel.find("hNjets_8")<10){
    hframe= new TH2F("hframe","hframe",10,-0.5,9.5,500,0.01,10000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",10, -0.5, 9.5, 1000, 0.1, 1.9);// mZ1
    hframe2->SetXTitle("jet multiplicity"); 
  }

  if (histlabel.find("hNbjets")<10){
    hframe= new TH2F("hframe","hframe",10,-0.5,9.5,500,0.01,10000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",10, -0.5, 9.5, 1000, 0.1, 1.9);// mZ1
    hframe2->SetXTitle("b tagging multiplicity");
  }

  if (histlabel.find("hN_loose_e")<10){
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hN_loose_mu")<10){
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hN_good_ele")<10){
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hN_good_mu")<10){
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hMZ1_8")<10 && whichchannel.find("4#mu")<20){
    hframe= new TH2F("hframe","hframe",80,40.,200.,500,0.0001,100000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 40., 160., 1000, 0.0, 2.0);// mZ1 
  }
  
  if (histlabel.find("hMZ2_8")<10 && whichchannel.find("4#mu")<20){
    hframe= new TH2F("hframe","hframe",80,40.,200.,500,0.0001,100000.);// mZ2 
    hframe2= new TH2F("hframe2","hframe2",6000, 40., 160., 1000, 0.5, 2.);// mZ2 
  }
  
  if (histlabel.find("hM4l_T_8")<10 ){
    hframe= new TH2F("hframe","hframe",80,10.,1200.,500,0.00000004,100000.);// transverse mass , M4l+MET
    hframe2= new TH2F("hframe2","hframe2",6000, 10., 1200., 1000, 0.5, 2.);// 
  }

  if (histlabel.find("DPHI_8")<10 ){
    hframe= new TH2F("hframe","hframe",80,0,3.14,500,0.00000004,1000.);// deltaphi 4l,MET
    hframe2= new TH2F("hframe2","hframe2",80,0,3.14, 1000, 0.5, 2.);// 
  }

  if (histlabel.find("hMELA_8")<10){
    hframe= new TH2F("hframe","hframe",600,0.,1.,600,0.,20.);// MELA at final stage
    hframe2= new TH2F("hframe2","hframe2",600, 0., 1., 1000, 0.5, 2.);// MELA at final stage
  }

  if (histlabel.find("hPFMET_8")<10){
    //hframe= new TH2F("hframe","hframe",1000, 0., 1000., 1000, 0.0000004, 50000.);// PFMET
    hframe= new TH2F("hframe","hframe",200,0,200., 1000, 0.001, 100000.);// PFMET
    hframe2= new TH2F("hframe2","hframe2",200, 0,200, 1000, 0.0, 2.0);// PFMET

  }

  if (histlabel.find("hPFMET_9")<10){
    //hframe= new TH2F("hframe","hframe",1000, 0., 1000., 1000, 0.0000004, 50000.);// PFMET
    hframe= new TH2F("hframe","hframe",1000, 0., 1000., 1000, 0.000001, 100.);// PFMET
    hframe2= new TH2F("hframe2","hframe2",1000, 0.,1000., 1000, 0.5, 1.5);// PFMET
  }

  if (histlabel.find("hMjj_8")<10){
    hframe= new TH2F("hframe","hframe",600,20.,500.,600,0.000004,10E4);//mass jet jet
    hframe2= new TH2F("hframe2","hframe2",6000, 20., 500., 1000, 0.5, 1.5);// mass jet jet
  }

  if (histlabel.find("hDjj_8")<10){
    hframe= new TH2F("hframe","hframe",600,0.,10.,600,0.000004,10E4);//delta eta jet jet
    hframe2= new TH2F("hframe2","hframe2",600, 0., 10., 1000, 0.5, 2.);// delta eta jet jet
  }

  //TH2F *hframe= new TH2F("hframe","hframe",6000, 0., 200., 1000, 0.004, 700000.);// ptZ

  if (nRebin==1) hframe->SetYTitle("Events/1 GeV");
  if (nRebin==2) hframe->SetYTitle("Events/2 GeV");
  if (nRebin==20) hframe->SetYTitle("Events/20 GeV");
  if (nRebin==3) hframe->SetYTitle("Events/3 GeV");
  if (nRebin==10) hframe->SetYTitle("Events/10 GeV");
  if (histlabel.find("hSip_3")<10) hframe->SetYTitle("Events/bin=0.12"); // nRebin=2 sip
  if (nRebin==4 && histlabel.find("hSip_3")<10) hframe->SetYTitle("Events/bin=0.24"); // nRebin=4 sip

  if (histlabel.find("hIso_3")<10) hframe->SetYTitle("Events/bin=0.04"); //nrebin=4 iso
  if (nRebin==5 && histlabel.find("hIso_3")<10) hframe->SetYTitle("Events/bin=0.05"); //nrebin=5 iso

  if (nRebin==2 && histlabel.find("hMjj_3")<10)  hframe->SetYTitle("Events/5 GeV"); // nRebin=2 mass jet jet
  if (nRebin==2 && histlabel.find("hDjj_3")<10)  hframe->SetYTitle("Events/bin=0.4"); // nRebin=2 deltetaeta
  if (nRebin==2 && histlabel.find("hVD_3")<10)   hframe->SetYTitle("Events/bin=0.1"); // nRebin=2 fisher
  if (nRebin==10 && histlabel.find("hMELA_8")<10) hframe->SetYTitle("Events / 0.033"); // MELA
  if (nRebin==5 && (histlabel.find("hPFMET_8")<10 ||histlabel.find("hPFMET_3")<10 )) hframe->SetYTitle("Events/5 GeV"); // PFMET
  

  
  if (histlabel.find("hM4l_T_8")<10 ){
    sprintf(histotitle,"m_{T} (%s+MET) [GeV]",whichchannel.c_str());
    hframe->SetXTitle(histotitle);
  }

  if (histlabel.find("DPHI_8")<10 ){
    sprintf(histotitle,"#Delta#phi (%s,MET) [GeV]",whichchannel.c_str());
    hframe->SetXTitle(histotitle);
    hframe->SetYTitle("Events/0.05 rad");
  }
  

  //hframe->SetXTitle("M_{Z2} [GeV]");
  if ((histlabel.find("hMZ_3")<10 || histlabel.find("hMZ1_8")<10) && whichchannel.find("4#mu")<20) 
    hframe->SetXTitle("M_{Z#rightarrow#mu#mu} [GeV]");
  if ((histlabel.find("hMZ_3")<10 || histlabel.find("hMZ1_8")<10) && whichchannel.find("4e")<20) 
    hframe->SetXTitle("M_{Z#rightarrow ee} [GeV]");
  //hframe->SetXTitle("M_{Z#rightarrow ee} (BB)  [GeV]");
  //hframe->SetXTitle("M_{Z#rightarrow ee} (EE)  [GeV]");
  if (histlabel.find("hPFMET_8")<10 || histlabel.find("hPFMET_3")<10 ) hframe2->SetXTitle("PF MET (GeV)"); // PFMET

 


  hframe->GetXaxis()->SetLabelOffset(0.005);
  hframe->GetXaxis()->SetTitleOffset(0.9);
  hframe->GetYaxis()->SetLabelOffset(0.007);

  hframe->Draw();

  TH1F *htotaldata = new TH1F("htotaldata", "htotaldata", Nbins, Xmin, Xmax);
  //TH1F *htotaldata = new TH1F("htotaldata", "htotaldata", 51, 95., 605.);
  //TH1F *htotaldata = new TH1F("htotaldata", "htotaldata", 50, 99.5, 599.5);
  htotaldata->SetMarkerColor(1);
  htotaldata->SetMarkerStyle(20);
  THStack *htotal=new THStack("Qier",""); 
  TH1F *htotalHisto = new TH1F("htotalHisto", "htotalHisto", Nbins, Xmin, Xmax);
  TH1F *htotalHistoRatio = new TH1F("htotalHistoRatio", "htotalHistoRatio", Nbins, Xmin, Xmax);

    if (histlabel.find("hPFMET_8")<10 /*|| histlabel.find("hPtJet_8")<10*/){
      double bins[16]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,100,200};
      htotaldata=(TH1F*)htotaldata->Rebin(15,histlabel.c_str(),bins);
      htotalHisto=(TH1F*)htotalHisto->Rebin(15,histlabel.c_str(),bins);
      htotalHistoRatio=(TH1F*)htotalHistoRatio->Rebin(15,histlabel.c_str(),bins);
    }
  cout << "Nbins=" << Nbins << endl;

  // data
  bool signal=false;

  for (unsigned int datasetIdData=0; datasetIdData<Vdatasetnamedata.size(); datasetIdData++){
    
    char dataset[328];
    sprintf(dataset,"%s",Vdatasetnamedata.at(datasetIdData).c_str());
    cout << "Root-ple= " << dataset << endl;
    // cout << "Counter=" << datasetIdData << " Root-ple=" << dataset << " Label=" << Vlabelbkg.at(datasetIdData) <<endl; 
    
    TFile *f1 = TFile::Open(dataset);
    hfourlepbestmass_4l_afterSel_new = (TH1F*)f1->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
    TH1 *hfourlepbestmass_4l_afterSel_new_new=NULL;
    if(histlabel.find("hPFMET_8")<10 /*|| histlabel.find("hPtJet_8")<10*/){
      double bins[16]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,100,200};
      hfourlepbestmass_4l_afterSel_new_new=(TH1F*)hfourlepbestmass_4l_afterSel_new->Rebin(15,histlabel.c_str(),bins);
    }
    else hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str());
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerColor(1);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(20);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerSize(0.95);
    // hfourlepbestmass_4l_afterSel_new_new->Draw("EPsame");
    //leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabeldata.at(datasetIdData).c_str(), "P");  
    //this part for efficiency of tight lepton calculation(qier)
/*
    TH1F *loose_e = (TH1F*)f1->Get("hN_loose_e");
    TH1F *loose_mu = (TH1F*)f1->Get("hN_loose_mu");
    TH1F *good_e = (TH1F*)f1->Get("hN_good_ele");
    TH1F *good_mu = (TH1F*)f1->Get("hN_good_mu");
    double deno1=0;
    double deno2=0;
    double nume1=0;
    double nume2=0;
    for (int nbins=2;nbins<=10; nbins++){
       deno1=deno1+(nbins-1)*loose_e->GetBinContent(nbins);
       deno2=deno2+(nbins-1)*loose_mu->GetBinContent(nbins);
       nume1=nume1+(nbins-1)*good_e->GetBinContent(nbins);
       nume2=nume2+(nbins-1)*good_mu->GetBinContent(nbins);
    }
    double eff_e = nume1/deno1;
    double eff_mu = nume2/deno2;
    cout << "eff_e = " << eff_e << "\neff_mu = " << eff_mu << endl;
*/

    //cout << "Nbins=" << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << endl;
    //cout << "htotaldata nBins = " << htotaldata->GetNbinsX() << ", hfourlepbestmass_4l_afterSel_new_new nBins = " << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << endl;
    //cout << "htotaldata lowestX = " << htotaldata->GetXaxis()->GetXmin() <<  ", htotaldata highestX = " << htotaldata->GetXaxis()->GetXmax() << ", hfourlepbestmass_4l_afterSel_new_new lowestX = " << hfourlepbestmass_4l_afterSel_new_new->GetXaxis()->GetXmin() << ", hfourlepbestmass_4l_afterSel_new_new highestX = " << hfourlepbestmass_4l_afterSel_new_new->GetXaxis()->GetXmax() << endl;
    cout << "test=" << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX()<< " test="<< htotaldata->GetNbinsX()<<endl;
    htotaldata->Add(hfourlepbestmass_4l_afterSel_new_new);
    cout << "Label= " << Vlabeldata.at(datasetIdData) <<"  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) <<endl;
    if (datasetIdData==(Vdatasetnamedata.size()-1)) {
      leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"Data 2011 + 2012", "P"); 
      legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"DATA 2016", "P"); 
      //htotaldata->Draw("EPsame");
    }
  }



  // Set Errors as in http://www-cdf.fnal.gov/physics/statistics/notes/pois_eb.txt
  //Float_t x[51],y[51],exl[51],exh[51],eyl[51],eyh[51];

  int* arraysize = new int[1];
  arraysize[0] = htotaldata->GetNbinsX();
  //std::cout << "arraysize = " << arraysize[0] << std::endl;
  Float_t x[arraysize[0]],y[arraysize[0]],exl[arraysize[0]],exh[arraysize[0]],eyl[arraysize[0]],eyh[arraysize[0]];
  delete [] arraysize;

  float totaldataentries=0.,totaldataentries100=0.;

  cout << "test bin number = " << htotaldata->GetNbinsX() << endl;
  for (int nbins=1;nbins<=htotaldata->GetNbinsX(); nbins++){
    // cout << "BinCenter=" << htotaldata->GetBinCenter(nbins) << " BinContent=" << htotaldata->GetBinContent(nbins) << " BinErrorContent=" << htotaldata->GetBinError(nbins) << endl;
    x[nbins-1]=htotaldata->GetBinCenter(nbins);
    y[nbins-1]=htotaldata->GetBinContent(nbins);
    exl[nbins-1]=0.5*htotaldata->GetBinWidth(nbins);
    exh[nbins-1]=0.5*htotaldata->GetBinWidth(nbins);
    totaldataentries=totaldataentries+htotaldata->GetBinContent(nbins);
    if (htotaldata->GetBinCenter(nbins)>100. && htotaldata->GetBinCenter(nbins)<800.) totaldataentries100=totaldataentries100+htotaldata->GetBinContent(nbins);
    if (htotaldata->GetBinContent(nbins)>0){
        eyh[nbins-1]=+0.5 + sqrt(htotaldata->GetBinContent(nbins)+0.25);   
        eyl[nbins-1]=-0.5 + sqrt(htotaldata->GetBinContent(nbins)+0.25);
    }
    else{
           x[nbins-1] = 0.;
         eyl[nbins-1] = 0.;
         eyh[nbins-1] = 0.;
    }
   cout << "x array " << nbins-1 << "= " << x[nbins-1] << endl;
   cout << "y array " << nbins-1 << "= " << y[nbins-1] << endl;
  }

  cout << "Total data= " << totaldataentries << endl;
  Nbins= htotaldata->GetNbinsX();
  TGraphAsymmErrors *gr = new TGraphAsymmErrors(Nbins,x,y,exl,exh,eyl,eyh);
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.95);
  //
 
  // Z+X from data

  // Background
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_350 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_350", "hfourlepbestmass_4l_afterSel_new_sig_350", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_400 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_400", "hfourlepbestmass_4l_afterSel_new_sig_400", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_450 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_450", "hfourlepbestmass_4l_afterSel_new_sig_450", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_500 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_500", "hfourlepbestmass_4l_afterSel_new_sig_500", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_550 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_550", "hfourlepbestmass_4l_afterSel_new_sig_550", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_600 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_600", "hfourlepbestmass_4l_afterSel_new_sig_600", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_650 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_650", "hfourlepbestmass_4l_afterSel_new_sig_650", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_700 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_700", "hfourlepbestmass_4l_afterSel_new_sig_700", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_750 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_750", "hfourlepbestmass_4l_afterSel_new_sig_750", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_sig_800 = new TH1F("hfourlepbestmass_4l_afterSel_new_sig_800", "hfourlepbestmass_4l_afterSel_new_sig_800", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_singlet = new TH1F("hfourlepbestmass_4l_afterSel_new_singlet", "hfourlepbestmass_4l_afterSel_new_singlet", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_DY      = new TH1F("hfourlepbestmass_4l_afterSel_new_DY", "hfourlepbestmass_4l_afterSel_new_DY",Nbins, Xmin, Xmax);   
  TH1F *hfourlepbestmass_4l_afterSel_new_ZZ      = new TH1F("hfourlepbestmass_4l_afterSel_new_ZZ", "hfourlepbestmass_4l_afterSel_new_ZZ",Nbins, Xmin, Xmax); 
  TH1F *hfourlepbestmass_4l_afterSel_new_DYlight = new TH1F("hfourlepbestmass_4l_afterSel_new_DYlight", "hfourlepbestmass_4l_afterSel_new_DYlight",Nbins, Xmin, Xmax);
        
  TH1F *hfourlepbestmass_4l_afterSel_new_DYbb    = new TH1F("hfourlepbestmass_4l_afterSel_new_DYbb", "hfourlepbestmass_4l_afterSel_new_DYbb",Nbins, Xmin, Xmax);        
  TH1F *hfourlepbestmass_4l_afterSel_new_DYcc    = new TH1F("hfourlepbestmass_4l_afterSel_new_DYcc", "hfourlepbestmass_4l_afterSel_new_DYcc",Nbins, Xmin, Xmax);        


  TH1F *hfourlepbestmass_4l_afterSel_new_WW    = new TH1F("hfourlepbestmass_4l_afterSel_new_WW", "hfourlepbestmass_4l_afterSel_new_WW",Nbins, Xmin, Xmax); 
  TH1F *hfourlepbestmass_4l_afterSel_new_WZ    = new TH1F("hfourlepbestmass_4l_afterSel_new_WZ", "hfourlepbestmass_4l_afterSel_new_WZ",Nbins, Xmin, Xmax);   
  
  TH1F *hfourlepbestmass_4l_afterSel_new_h    = new TH1F("hfourlepbestmass_4l_afterSel_new_h", "hfourlepbestmass_4l_afterSel_new_h",Nbins, Xmin, Xmax);                   
  TH1F *hfourlepbestmass_4l_afterSel_new_TT    = new TH1F("hfourlepbestmass_4l_afterSel_new_TT", "hfourlepbestmass_4l_afterSel_new_TT",Nbins, Xmin, Xmax); 
  TH1F *hfourlepbestmass_4l_afterSel_new_Wj    = new TH1F("hfourlepbestmass_4l_afterSel_new_Wj", "hfourlepbestmass_4l_afterSel_new_Wj",Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_TTZ    = new TH1F("hfourlepbestmass_4l_afterSel_new_TTZ", "hfourlepbestmass_4l_afterSel_new_TTZ",Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_gg    = new TH1F("hfourlepbestmass_4l_afterSel_new_gg", "hfourlepbestmass_4l_afterSel_new_gg",Nbins, Xmin, Xmax);
  if(histlabel.find("hPFMET_8")<10/*|| histlabel.find("hPtJet_8")<10*/){
      double bins[16]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,100,200};
      hfourlepbestmass_4l_afterSel_new_DY=(TH1F*)hfourlepbestmass_4l_afterSel_new_DY->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_TT=(TH1F*)hfourlepbestmass_4l_afterSel_new_TT->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_h=(TH1F*)hfourlepbestmass_4l_afterSel_new_h->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_TTZ=(TH1F*)hfourlepbestmass_4l_afterSel_new_TTZ->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_ZZ=(TH1F*)hfourlepbestmass_4l_afterSel_new_ZZ->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_350=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_350->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_400=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_400->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_450=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_450->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_500=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_500->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_550=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_550->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_600=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_600->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_650=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_650->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_700=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_700->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_750=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_750->Rebin(15,histlabel.c_str(),bins);
      hfourlepbestmass_4l_afterSel_new_sig_800=(TH1F*)hfourlepbestmass_4l_afterSel_new_sig_800->Rebin(15,histlabel.c_str(),bins);

      hfourlepbestmass_4l_afterSel_new_gg=(TH1F*)hfourlepbestmass_4l_afterSel_new_gg->Rebin(15,histlabel.c_str(),bins);
  }

  
/*
  TFile * PU= new TFile("pileup_MC_80x_271036-276811_69200.root");
  TH1F * HistoPUData= (TH1F *) PU->Get("pileup");
  TH1F * HistoPUMC= (TH1F *) PU->Get("pileup_mc");
  double nnn = HistoPUData->Integral()/HistoPUMC->Integral();
  cout << "nnn = " << nnn << endl;
*/
//signal
  for ( int datasetId=Vdatasetnamesig.size()-1; datasetId >=0; datasetId--){
    char dataset[328];
    sprintf(dataset,"%s",Vdatasetnamesig.at(datasetId).c_str());
    //cout << "Root-ple= " << dataset << "N entries= " <<  hfourlepbestmass_4l_afterSel_new->GetEntries() << endl;
    cout << "Counter=" << datasetId << " Root-ple=" << dataset << " Label=" << Vlabelsig.at(datasetId) <<endl;
 
    std::string datasetnamesig = "";
    datasetnamesig = Vdatasetnamesig.at(datasetId);

    TFile *f2 = TFile::Open(dataset);
    hfourlepbestmass_4l_afterSel_new = (TH1F*)f2->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
    TH1 *hfourlepbestmass_4l_afterSel_new_new;
      if(histlabel.find("hPFMET_8")<10/*|| histlabel.find("hPtJet_8")<10*/){
        double bins[16]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,100,200};
        hfourlepbestmass_4l_afterSel_new_new=(TH1F*)hfourlepbestmass_4l_afterSel_new->Rebin(15,histlabel.c_str(),bins);
      }
      else hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str());

      hfourlepbestmass_4l_afterSel_new_new->SetLineColor(Vcolorsig.at(datasetId));
      hfourlepbestmass_4l_afterSel_new_new->SetFillColor(Vcolorsig.at(datasetId));
      hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(24);
      hfourlepbestmass_4l_afterSel_new_new->SetLineWidth(2);
      cout << "Label= " << Vlabelsig.at(datasetId)
           << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
           << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1);
      if (hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0)
        cout << "  Error= " << sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries();
      cout << endl;


      if(datasetnamesig.find("sig_350") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.2535*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_350->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_350->SetMarkerColor(40);
        hfourlepbestmass_4l_afterSel_new_sig_350->SetLineColor(40);
        hfourlepbestmass_4l_afterSel_new_sig_350->SetLineWidth(2);
        char temp[328];
        sprintf(temp,"%s",histosdir.c_str());
          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H350", "L");
      }
      if(datasetnamesig.find("sig_400") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_400->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_400->SetMarkerColor(6);
        hfourlepbestmass_4l_afterSel_new_sig_400->SetLineColor(6);
        hfourlepbestmass_4l_afterSel_new_sig_400->SetLineWidth(3);
        cout << "fill sig histograms" << endl;

        char temp[328];
        sprintf(temp,"%s",histosdir.c_str());
          cout << "sig= " << hfourlepbestmass_4l_afterSel_new_sig_400->Integral(0,-1) << endl;
          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H400", "L");
      }

      if(datasetnamesig.find("sig_450") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.005*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_450->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_450->SetMarkerColor(41);
        hfourlepbestmass_4l_afterSel_new_sig_450->SetLineColor(41);
        hfourlepbestmass_4l_afterSel_new_sig_450->SetLineWidth(2);

          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H450", "L");
      }
      if(datasetnamesig.find("sig_500") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.019*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_500->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_500->SetMarkerColor(42);
        hfourlepbestmass_4l_afterSel_new_sig_500->SetLineColor(42);
        hfourlepbestmass_4l_afterSel_new_sig_500->SetLineWidth(2);

          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H500", "L");
      }
      if(datasetnamesig.find("sig_550") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.039*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_550->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_550->SetMarkerColor(43);
        hfourlepbestmass_4l_afterSel_new_sig_550->SetLineColor(43);
        hfourlepbestmass_4l_afterSel_new_sig_550->SetLineWidth(2);
        char temp[328];
        sprintf(temp,"%s",histosdir.c_str());
          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H550", "L");
      }
      if(datasetnamesig.find("sig_600") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.023*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_600->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_600->SetMarkerColor(7);
        hfourlepbestmass_4l_afterSel_new_sig_600->SetLineColor(7);
        hfourlepbestmass_4l_afterSel_new_sig_600->SetLineWidth(3);

        char temp[328];
        sprintf(temp,"%s",histosdir.c_str());
          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H600", "L");
      }
      if(datasetnamesig.find("sig_650") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.001*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_650->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_650->SetMarkerColor(7);
        hfourlepbestmass_4l_afterSel_new_sig_650->SetLineColor(7);
        hfourlepbestmass_4l_afterSel_new_sig_650->SetLineWidth(3);

          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H650", "L");
      }
      if(datasetnamesig.find("sig_700") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.029*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_700->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_700->SetMarkerColor(46);
        hfourlepbestmass_4l_afterSel_new_sig_700->SetLineColor(46);
        hfourlepbestmass_4l_afterSel_new_sig_700->SetLineWidth(2);
          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H700", "L");
      }
      if(datasetnamesig.find("sig_750") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.039*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_750->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_750->SetMarkerColor(47);
        hfourlepbestmass_4l_afterSel_new_sig_750->SetLineColor(47);
        hfourlepbestmass_4l_afterSel_new_sig_750->SetLineWidth(2);

          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H750", "L");
      }

      if(datasetnamesig.find("sig_800") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         signal = true;
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.024*0.00115*35812.*nnn/(99960)));
        hfourlepbestmass_4l_afterSel_new_sig_800->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_sig_800->SetMarkerColor(kYellow+1);
        hfourlepbestmass_4l_afterSel_new_sig_800->SetLineColor(kYellow+1);
        hfourlepbestmass_4l_afterSel_new_sig_800->SetLineWidth(3);

          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"H800", "L");
      }
  }

  for ( int datasetId=Vdatasetnamebkg.size()-1; datasetId >=0; datasetId--){  

    char dataset[328];
    sprintf(dataset,"%s",Vdatasetnamebkg.at(datasetId).c_str());
    //cout << "Root-ple= " << dataset << "N entries= " <<  hfourlepbestmass_4l_afterSel_new->GetEntries() << endl;
    cout << "Counter=" << datasetId << " Root-ple=" << dataset << " Label=" << Vlabelbkg.at(datasetId) <<endl; 
    
    std::string datasetnamebkg = "";
    datasetnamebkg = Vdatasetnamebkg.at(datasetId);

    TFile *f2 = TFile::Open(dataset);
    hfourlepbestmass_4l_afterSel_new = (TH1F*)f2->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
    TH1 *hfourlepbestmass_4l_afterSel_new_new;

    //this part for tight lepton efficiency calculation (qier)
/*
    TH1F *loose_e = (TH1F*)f2->Get("hN_loose_e");
    TH1F *loose_mu = (TH1F*)f2->Get("hN_loose_mu");
    TH1F *good_e = (TH1F*)f2->Get("hN_good_ele");
    TH1F *good_mu = (TH1F*)f2->Get("hN_good_mu");
    double deno1=0;
    double deno2=0;
    double nume1=0;
    double nume2=0;

    for (int nbins=2;nbins<=10; nbins++){
       deno1=deno1+(nbins-1)*loose_e->GetBinContent(nbins);
       deno2=deno2+(nbins-1)*loose_mu->GetBinContent(nbins);
       nume1=nume1+(nbins-1)*good_e->GetBinContent(nbins);
       nume2=nume2+(nbins-1)*good_mu->GetBinContent(nbins);
    }
    double eff_e = nume1/deno1;
    double eff_mu = nume2/deno2;
    cout << "eff_e = " << eff_e << "\neff_mu = " << eff_mu << endl;



    TH1 *nevent = (TH1F*)f2->Get("nEvent_4l");
    TH1 *nevent_w = (TH1F*)f2->Get("nEvent_4l_w");
*/   
    if(datasetnamebkg.find("WZ")   < 200 ||
       datasetnamebkg.find("DYJetsToLL") < 200 ||
       datasetnamebkg.find("TT_TuneCUETP8M1") < 200 ||
       datasetnamebkg.find("WWTo2L2Nu")< 200 ||
       datasetnamebkg.find("WJetsToLNu") < 200 ||
       datasetnamebkg.find("ZZTo4L") < 200 ||
       datasetnamebkg.find("TTZToLLNuNu") < 200 ||
       datasetnamebkg.find("GluGluTo4L") < 200 ||
       datasetnamebkg.find("hTo4l") < 200 ||
       datasetnamebkg.find("WZJToLLLNu") < 200
       ){
      
      if(histlabel.find("hPFMET_8")<10 /*|| histlabel.find("hPtJet_8")<10*/){
        double bins[16]={0,5,10,15,20,25,30,35,40,45,50,55,60,70,100,200};
        hfourlepbestmass_4l_afterSel_new_new=(TH1F*)hfourlepbestmass_4l_afterSel_new->Rebin(15,histlabel.c_str(),bins);
      }
      else hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str());

      hfourlepbestmass_4l_afterSel_new_new->SetLineColor(Vcolorbkg.at(datasetId)/*datasetId+2*/);
      hfourlepbestmass_4l_afterSel_new_new->SetFillColor(Vcolorbkg.at(datasetId)/*datasetId+2*/);
      hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(24);
      hfourlepbestmass_4l_afterSel_new_new->SetLineWidth(1);
      cout << "Label= " << Vlabelbkg.at(datasetId) 
	   << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
	   << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1);
      if (hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0) 
	cout << "  Error= " << sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries();
      cout << endl;

      // Higgs as background
      // DYJetsToLL check normalization
      if(datasetnamebkg.find("DYJetsToLL") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         cout << "test1" << endl;
	 hfourlepbestmass_4l_afterSel_new_new->Scale(double(5765.4*381.*nnn/(5223984*0.6678)));//*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/12138430./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
	hfourlepbestmass_4l_afterSel_new_DY->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DY->SetMarkerColor(kAzure+2);
	hfourlepbestmass_4l_afterSel_new_DY->SetFillColor(kAzure+2);                                                        
        cout << "fill DY histograms" << endl;        

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if(datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && (datasetnamebkg.find("DYJetsToLL_M-50_TuneZ2Star")<200 || datasetnamebkg.find("DYJetsToLL")<200)) {
          cout << "DY= " << hfourlepbestmass_4l_afterSel_new_DY->Integral(0,-1) << endl;
	  if (useDYJets==true) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
	}
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }


      // DYlightJetsToLL check normalization
      if(datasetnamebkg.find("DYlightJetsToLL") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){	
	hfourlepbestmass_4l_afterSel_new_DYlight->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DYlight->SetMarkerColor(kAzure+6);
	//hfourlepbestmass_4l_afterSel_new_DYlight->SetLineColor(kAzure+6);
	//hfourlepbestmass_4l_afterSel_new_DYlight->SetLineWidth(2); 
	hfourlepbestmass_4l_afterSel_new_DYlight->SetFillColor(kAzure+6);                                                        
	
	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if(datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && datasetnamebkg.find("DYlightJetsToLL_TuneZ2_M-50")<200) {
          cout << "DYlight= " << hfourlepbestmass_4l_afterSel_new_DYlight->Integral(0,-1) << endl;
	  if (useDYJets==false) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
        }
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }
      
      // DYbb
      if(datasetnamebkg.find("DYbbJetsToLL") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
	hfourlepbestmass_4l_afterSel_new_DYbb->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DYbb->SetMarkerColor(kAzure+2);
	//	hfourlepbestmass_4l_afterSel_new_DYbb->SetLineColor(kAzure+2);
	//	hfourlepbestmass_4l_afterSel_new_DYbb->SetLineWidth(2);
	hfourlepbestmass_4l_afterSel_new_DYbb->SetFillColor(kAzure+2);                                                        
	
	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if(datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && datasetnamebkg.find("DYbbJetsToLL_TuneZ2_M-50")<200) {
          cout << "DYbb= " << hfourlepbestmass_4l_afterSel_new_DYbb->Integral(0,-1) << endl;
	  if (useDYJets==false) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
        }
      }

      //DYCC
      if(datasetnamebkg.find("DYccJetsToLL") < 200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
	
	hfourlepbestmass_4l_afterSel_new_DYcc->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DYcc->SetMarkerColor(kRed+0);
	// hfourlepbestmass_4l_afterSel_new_DYcc->SetLineColor(kRed+0);
	// hfourlepbestmass_4l_afterSel_new_DYcc->SetLineWidth(2);
       	hfourlepbestmass_4l_afterSel_new_DYcc->SetFillColor(kRed+0);                                                        
	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if(datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && datasetnamebkg.find("DYccJetsToLL_M-50_TuneZ2Star")<200) {
          cout << "DYcc= " << hfourlepbestmass_4l_afterSel_new_DYcc->Integral(0,-1) << endl;
	  if (useDYJets==false) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
        }
      }

      if (useDYJetsFromData==false){
	// WW
	if(datasetnamebkg.find("WWTo2L2Nu") < 200){
	  hfourlepbestmass_4l_afterSel_new_WW->Add(hfourlepbestmass_4l_afterSel_new_new); 
	  hfourlepbestmass_4l_afterSel_new_WW->SetMarkerColor(kCyan+3); 
	  hfourlepbestmass_4l_afterSel_new_WW->SetFillColor(kCyan+3);
	  
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  if(datasetnamebkg.find(temp) < 200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && hfourlepbestmass_4l_afterSel_new_WW->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");  
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");    
	}     
	
	// WZ     
	if(datasetnamebkg.find("WZJToLLLNu") < 200){  
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(5.26*35812.*nnn/(1930828.*7.5)));
	  hfourlepbestmass_4l_afterSel_new_WZ->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_WZ->SetMarkerColor(kCyan-2);  
	  hfourlepbestmass_4l_afterSel_new_WZ->SetFillColor(kCyan-2);  
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP"); 
	}    
	
	// TTT
	if(datasetnamebkg.find("TT_TuneCUETP8M1") < 200){
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(831.76*381.*nnn/5800000.)); 
	  hfourlepbestmass_4l_afterSel_new_TT->Add(hfourlepbestmass_4l_afterSel_new_new);     
	  hfourlepbestmass_4l_afterSel_new_TT->SetMarkerColor(kTeal-6); 
	  hfourlepbestmass_4l_afterSel_new_TT->SetFillColor(kTeal-6);                     
	  
	  cout << "TT+jets= " << hfourlepbestmass_4l_afterSel_new_TT->GetEntries() << endl;
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  if(datasetnamebkg.find(temp) < 200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && hfourlepbestmass_4l_afterSel_new_TT->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
	}      

        if(datasetnamebkg.find("hTo4l") < 200){
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(0.01212*1.7*35812.*nnn/992224.)); 
          hfourlepbestmass_4l_afterSel_new_h->Add(hfourlepbestmass_4l_afterSel_new_new);
          hfourlepbestmass_4l_afterSel_new_h->SetMarkerColor(kGray+1);
          hfourlepbestmass_4l_afterSel_new_h->SetFillColor(kGray+1);

          cout << "hTo4l= " << hfourlepbestmass_4l_afterSel_new_h->Integral(0,-1) << endl;
          char temp[328];
          sprintf(temp,"%s",histosdir.c_str()); 
          if(datasetnamebkg.find(temp) < 200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && hfourlepbestmass_4l_afterSel_new_h->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
          //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
        }

        if(datasetnamebkg.find("GluGluTo4L") < 200){
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(0.01427*1.7*35812.*nnn/500000.));
          hfourlepbestmass_4l_afterSel_new_gg->Add(hfourlepbestmass_4l_afterSel_new_new);
          hfourlepbestmass_4l_afterSel_new_gg->SetMarkerColor(kAzure-1);
          hfourlepbestmass_4l_afterSel_new_gg->SetFillColor(kAzure-1);

          cout << "ggZZTo4l= " << hfourlepbestmass_4l_afterSel_new_gg->Integral(0,-1) << endl;
          char temp[328];
          sprintf(temp,"%s",histosdir.c_str());
          if(datasetnamebkg.find(temp) < 200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && hfourlepbestmass_4l_afterSel_new_gg->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
          //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
        }
	
        if(datasetnamebkg.find("TTZToLLNuNu") < 200){
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(0.2529*35812*nnn/(5934228.*0.2529)));
          hfourlepbestmass_4l_afterSel_new_TTZ->Add(hfourlepbestmass_4l_afterSel_new_new);
          hfourlepbestmass_4l_afterSel_new_TTZ->SetMarkerColor(kTeal-6);
          hfourlepbestmass_4l_afterSel_new_TTZ->SetFillColor(kTeal-6);

          cout << "TTZ " << hfourlepbestmass_4l_afterSel_new_TTZ->Integral(0,-1) << endl;
          char temp[328];
          sprintf(temp,"%s",histosdir.c_str());
          if( hfourlepbestmass_4l_afterSel_new_TTZ->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
          //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");

        }

	
	// W+jets		      					
	if(datasetnamebkg.find("_WJetsToLNu") < 200){
	  cout << "Wjets" << endl;
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(60290.*1031./871626.));
	  hfourlepbestmass_4l_afterSel_new_Wj->Add(hfourlepbestmass_4l_afterSel_new_new);        
	  hfourlepbestmass_4l_afterSel_new_Wj->SetMarkerColor(kSpring);	
	  hfourlepbestmass_4l_afterSel_new_Wj->SetFillColor(kSpring);	
	  
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  if(datasetnamebkg.find(temp) < 200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && hfourlepbestmass_4l_afterSel_new_Wj->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");	        
	}   
      } 
      
      if(datasetnamebkg.find("ZZTo4L") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         hfourlepbestmass_4l_afterSel_new_new->Scale(double(1.212*1.15*35812*nnn/(10709784*1.793)));
        hfourlepbestmass_4l_afterSel_new_ZZ->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_ZZ->SetMarkerColor(kAzure+2);
        hfourlepbestmass_4l_afterSel_new_ZZ->SetFillColor(kAzure+2);
        cout << "fill DY histograms" << endl;

        char temp[328];
        sprintf(temp,"%s",histosdir.c_str());
          cout << "ZZ= " << hfourlepbestmass_4l_afterSel_new_ZZ->Integral(0,-1) << endl;
          legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
        //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");

      }  
    }
    //
    else if (useDYJetsFromData==false){   
      // single top
      if(datasetnamebkg.find("ST_") < 200 ||  datasetnamebkg.find("Tbar_") < 200 ){
	hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_singlet->Add(hfourlepbestmass_4l_afterSel_new_new);
	// hfourlepbestmass_4l_afterSel_new_new->SetMarkerColor(datasetId+4);
	//       hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(26);
	//       hfourlepbestmass_4l_afterSel_new_new->SetLineWidth(2);
	//       hfourlepbestmass_4l_afterSel_new_new->Dra;
	hfourlepbestmass_4l_afterSel_new_singlet->SetLineColor(Vcolorbkg.at(datasetId)/*datasetId-12+5*/);
	hfourlepbestmass_4l_afterSel_new_singlet->SetFillColor(kViolet);
	//hfourlepbestmass_4l_afterSel_new_singlet->SetFillStyle(3004);  
	hfourlepbestmass_4l_afterSel_new_singlet->SetLineWidth(1);
	
	
	char temp[328];
	sprintf(temp,"%s/output_ST_",histosdir.c_str());
	
	if(datasetnamebkg.find(temp) < 200 && datasetnamebkg.find("t-channel") < 200 ){ // provided that this is the last single-top sample
	  //hfourlepbestmass_4l_afterSel_new_singlet->Draw("sameP");
	  leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_singlet,"Single Top", "F");  
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_singlet,"Single Top", "F");  
	  cout << "Label= single t     Entries= " << hfourlepbestmass_4l_afterSel_new_singlet->Integral(0,-1) <<endl;
	}  
	cout << "Label= " << Vlabelbkg.at(datasetId) << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) <<endl;
      }
    }
    
    char tempp[328];
    if (whichsample.find("8TeV")<200) sprintf(tempp,"%s/output_ZZTo2e2mu_%s",histosdir.c_str(),whichsample.c_str());
    else if (whichsample.find("7TeV")<200) sprintf(tempp,"%s/output_ZZTo2e2mu_mll4_%s",histosdir.c_str(),whichsample.c_str());    
    else if (whichsample.find("13TeV")<200) sprintf(tempp,"%s/output_ZZTo4L_%s",histosdir.c_str(),whichsample.c_str());  
    //else if (whichsample.find("13TeV")<200) sprintf(tempp,"%s/output_ZZ_TuneCUETP8M1_%s",histosdir.c_str(),whichsample.c_str());    
    
       char temppp[328];
       sprintf(temppp,"%s",histosdir.c_str());


       // other backgrounds
       if (useDYJets==true){
	 if(datasetnamebkg.find(temppp) < 200 && (
						 datasetnamebkg.find("output_DYJetsToLL") < 200 || 
						 (datasetnamebkg.find("output_DYJetsToLL_TuneZ2_M-50") < 200 && datasetnamebkg.find(whichenergy.c_str())<200)
						 ) 
	    )  {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DY);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DY);
	 }
       }
       else if (useDYJets==false){
	 if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_DYlightJetsToLL_TuneZ2_M-50") < 200 )  {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DYlight);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYlight);
	 }
	 if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_DYbbJetsToLL_TuneZ2_M-50") < 200 )  {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DYbb);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYbb);
	 }
	 if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_DYccJetsToLL_TuneZ2_M-50") < 200 )  {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DYcc);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYcc);
	 }
       }

       if(datasetnamebkg.find(temppp) < 200 && (
	  datasetnamebkg.find("output_WW") < 200 
	  )
	  )  { 
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_WW); 
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_WW); 
       }
       if(datasetnamebkg.find(temppp) < 200 && (
					       datasetnamebkg.find("output_WZJToLLLNu") < 200 
					       ) 
	  )  {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_WZ);         
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_WZ);   
       }
       if(datasetnamebkg.find(temppp) < 200 && (
					       datasetnamebkg.find("output_TT") < 200 
					       )
	  )  {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_TT);                        
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_TT);   
       }
       if(datasetnamebkg.find(temppp) < 200 && (
                                               datasetnamebkg.find("output_hTo4l") < 200 
                                               )
          )  {
         htotal->Add(hfourlepbestmass_4l_afterSel_new_h);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_h);
       }

       if(datasetnamebkg.find(temppp) < 200 && (
                                               datasetnamebkg.find("output_GluGluTo4L") < 200
                                               )
          )  {
         htotal->Add(hfourlepbestmass_4l_afterSel_new_gg);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_gg);
       }

       if(datasetnamebkg.find(temppp) < 200 && (
                                               datasetnamebkg.find("output_TTZToLLNuNu") < 200
                                               )
          )  {
         htotal->Add(hfourlepbestmass_4l_afterSel_new_TTZ);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_TTZ);
         cout << "testha" << endl;
       }
       if(datasetnamebkg.find(temppp) < 200 && ( 
          datasetnamebkg.find("output_ZZ") < 200 
          )
          )  {
         htotal->Add(hfourlepbestmass_4l_afterSel_new_ZZ);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_ZZ);
       }

       if(datasetnamebkg.find(temppp) < 200 && (
					       datasetnamebkg.find("output_WJ") < 200 
					       )
	  )  {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_Wj);   
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_Wj);
       }
       
       
       //if(datasetnamebkg.find("GluGluToZZTo4mu_BackgroundOnly") < 200 ) // provided that this is the last single-top sample      
       //htotal->Add(hfourlepbestmass_4l_afterSel_new_ggZZ);      
       
       //if(datasetnamebkg.find(temppp) < 200 datasetnamebkg.find("output_ZZTo2e2mu_8TeV") < 200 ) 
       //	htotal->Add(hfourlepbestmass_4l_afterSel_new_qqZZ); 
       ////htotal->Add(hfourlepbestmass_4l_afterSel_new_qqZZ); 
       
       if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_ST_") < 200 && datasetnamebkg.find("t-channel") < 200 ){     	 
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_singlet); 
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_singlet); 
       }     
     
     // htotal->Add(hfourlepbestmass_4l_afterSel_new_new);
     
  }

   TH1F *denom_up, *denom_dow;
   TH1F *denom_1;

   TH1F *denom_mc_up, *denom_mc_dow;
   TH1F *denom_mc_1;

   TString file_un="/eos/uscms/store/user/wangz/ntuple/unc/4mu2b/un_all.root";
   TString histname = histlabel.c_str();

   if(un_mc){
     TFile *f3 = TFile::Open(file_un);

     TH1F *hfourlepbestmass_mc_up = (TH1F*)f3->Get(histname+"_all_up_rel");
     TH1F *hfourlepbestmass_mc_dow = (TH1F*)f3->Get(histname+"_all_dow_rel");

     denom_mc_up=hfourlepbestmass_mc_up;
     denom_mc_dow=hfourlepbestmass_mc_dow;
   }


   TH1F *denom = (TH1F*)htotal->GetStack()->Last()->Clone();
   TH1F *ratiostaterr = (TH1F*)denom->Clone();
   ratiostaterr->SetStats(0);
   ratiostaterr->SetTitle("");
 //  ratiostaterr->GetYaxis()->SetTitle("Obs / Exp");
//    ratiostaterr->SetMaximum(1.2);
//    ratiostaterr->SetMinimum(0.8);
    ratiostaterr->SetMarkerSize(0);
    ratiostaterr->SetFillColor(kGray+3);
    ratiostaterr->SetFillStyle(3013);
    ratiostaterr->GetXaxis()->SetLabelSize(0.19);
    ratiostaterr->GetXaxis()->SetTitleSize(0.21);
    ratiostaterr->GetXaxis()->SetTitleOffset(1.0);
    ratiostaterr->GetXaxis()->SetLabelOffset(0.03);
    ratiostaterr->GetYaxis()->SetLabelSize(0.19);
    ratiostaterr->GetYaxis()->SetLabelOffset(0.006);
    ratiostaterr->GetYaxis()->SetTitleSize(0.21);
    ratiostaterr->GetYaxis()->SetTitleOffset(0.35);
    ratiostaterr->GetYaxis()->SetNdivisions(503);

    int nnbins=htotaldata->GetNbinsX();

    int* arraysize2 = new int[1];
    arraysize2[0] = htotaldata->GetNbinsX();
    //std::cout << "arraysize = " << arraysize[0] << std::endl;
    Float_t x2[arraysize2[0]],y2[arraysize2[0]],exl2[arraysize2[0]],exh2[arraysize2[0]],eyl2[arraysize2[0]],eyh2[arraysize2[0]];
    Float_t x3[arraysize2[0]],y3[arraysize2[0]],exl3[arraysize2[0]],exh3[arraysize2[0]],eyl3[arraysize2[0]],eyh3[arraysize2[0]],eyl4[arraysize2[0]],eyh4[arraysize2[0]];
    Float_t rx[arraysize2[0]],ry[arraysize2[0]],rexl[arraysize2[0]],rexh[arraysize2[0]],reyl[arraysize2[0]],reyh[arraysize2[0]],reyl1[arraysize2[0]],reyh1[arraysize2[0]],reyl2[arraysize2[0]],reyh2[arraysize2[0]];
    delete [] arraysize2;


  // Building the ratio
  for (int nbins=1;nbins<=htotaldata->GetNbinsX(); nbins++){
    //cout << "Total: BinCenter=" << htotalHisto->GetBinCenter(nbins) << " BinContent=" << htotalHisto->GetBinContent(nbins) << " BinErrorContent=" << htotalHisto->GetBinError(nbins) << endl;
    if (htotalHisto->GetBinContent(nbins)>0.) {
      htotalHistoRatio->SetBinContent(nbins,double(htotaldata->GetBinContent(nbins)/htotalHisto->GetBinContent(nbins)));
      htotalHistoRatio->SetBinError(nbins,double((htotaldata->GetBinError(nbins))/htotalHisto->GetBinContent(nbins)));
       
/*      htotalHistoRatio->SetBinError(nbins,double(sqrt( 
		 (1./(htotalHisto->GetBinContent(nbins)*htotalHisto->GetBinContent(nbins)))*htotaldata->GetBinContent(nbins) +
		 (htotaldata->GetBinContent(nbins)*htotaldata->GetBinContent(nbins)/pow(htotalHisto->GetBinContent(nbins),4))
		 *htotalHisto->GetBinContent(nbins) 						      
		 )));
*/		      
    }

    rx[nbins-1]=htotalHistoRatio->GetBinCenter(nbins);
    ry[nbins-1]=htotalHistoRatio->GetBinContent(nbins);
    rexl[nbins-1]=0.5*htotalHistoRatio->GetBinWidth(nbins);
    rexh[nbins-1]=0.5*htotalHistoRatio->GetBinWidth(nbins);

   ratiostaterr->SetBinContent(nbins, 1.0);
    if (denom->GetBinContent(nbins)>0.) {
     double binerror = denom->GetBinError(nbins)/denom->GetBinContent(nbins);
//     cout << "error="<< denom->GetBinError(nbins) << " content=" << denom->GetBinContent(nbins) << endl;
      ratiostaterr->SetBinError(nbins,binerror);
      reyl[nbins-1] = eyl[nbins-1]/denom->GetBinContent(nbins);
      reyh[nbins-1] = eyh[nbins-1]/denom->GetBinContent(nbins);
    }
    else {
      ratiostaterr->SetBinError(nbins,999.);
      reyl[nbins-1] = 0;
      reyh[nbins-1] = 0;
    }

    x2[nbins-1]=htotaldata->GetBinCenter(nbins);
    y2[nbins-1]=denom->GetBinContent(nbins);
    exl2[nbins-1]=0.5*htotaldata->GetBinWidth(nbins);
    exh2[nbins-1]=0.5*htotaldata->GetBinWidth(nbins);

    x3[nbins-1]=htotaldata->GetBinCenter(nbins);
    y3[nbins-1]=1.0;
    exl3[nbins-1]=0.5*htotaldata->GetBinWidth(nbins);
    exh3[nbins-1]=0.5*htotaldata->GetBinWidth(nbins);

    if (denom->GetBinContent(nbins)>0){
        eyh2[nbins-1]=denom->GetBinError(nbins);
        eyl2[nbins-1]=denom->GetBinError(nbins);
        eyh3[nbins-1]=denom->GetBinError(nbins)/denom->GetBinContent(nbins);
        eyl3[nbins-1]=denom->GetBinError(nbins)/denom->GetBinContent(nbins);

        if(un_mc){
          int nbk = denom_mc_up->GetXaxis()->FindBin(x3[nbins-1]);
          double e_mc_up=denom_mc_up->GetBinContent(nbk)-1;
          double e_mc_dow=1-denom_mc_dow->GetBinContent(nbk);
          reyl2[nbins-1]=e_mc_dow;
          reyh2[nbins-1]=e_mc_up;
          eyl4[nbins-1]=denom->GetBinContent(nbins)*e_mc_dow;
          eyh4[nbins-1]=denom->GetBinContent(nbins)*e_mc_up;
        }
   }
    else{
           x2[nbins-1] = 0.;
         eyl2[nbins-1] = 0.;
         eyh2[nbins-1] = 0.;
           x3[nbins-1] = 0.;
         eyl3[nbins-1] = 1.;
         eyh3[nbins-1] = 1.;
         reyl1[nbins-1] = 0.;
         reyh1[nbins-1] =0.;
         reyl2[nbins-1] = 0.;
         reyh2[nbins-1] = 0.;
          eyl4[nbins-1]= 0 ;
          eyh4[nbins-1]= 0 ;
    }

  }

    TGraphAsymmErrors *rgr = new TGraphAsymmErrors(Nbins,rx,ry,rexl,rexh,reyl,reyh);
    rgr->SetMarkerColor(1);
    rgr->SetMarkerStyle(20);
    rgr->SetMarkerSize(0.95);

    
    //cout << "Nbins=" << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << endl;
           
   htotal->GetHistogram()->GetXaxis()->SetLabelOffset(999);
   TH1F *temp = NULL;
   temp = (TH1F*)htotal->GetStack()->Last();
   TH1F *staterr = (TH1F*)temp->Clone("htotal");
   staterr->SetFillColor(kGray+3);
   staterr->SetLineColor(kGray+3);
   staterr->SetLineWidth(0);
   staterr->SetMarkerSize(0);
   staterr->SetFillStyle(3013);
   staterr->Draw("e2 same");

    TGraphAsymmErrors *rgr_mc = new TGraphAsymmErrors(Nbins,x3,y3,exl3,exh3,reyl2,reyh2);
    rgr_mc->SetTitle("");
    rgr_mc->SetMarkerSize(0);
    rgr_mc->SetFillColor(kAzure+1);
    rgr_mc->SetFillStyle(3021);
    rgr_mc->GetXaxis()->SetLabelSize(0.19);
    rgr_mc->GetXaxis()->SetTitleSize(0.21);
    rgr_mc->GetXaxis()->SetTitleOffset(1.0);
    rgr_mc->GetXaxis()->SetLabelOffset(0.03);
    rgr_mc->GetYaxis()->SetLabelSize(0.19);
    rgr_mc->GetYaxis()->SetLabelOffset(0.006);
    rgr_mc->GetYaxis()->SetTitleSize(0.21);
    rgr_mc->GetYaxis()->SetTitleOffset(0.35);
    rgr_mc->GetYaxis()->SetNdivisions(503);

    TGraphAsymmErrors *gr2 = new TGraphAsymmErrors(Nbins,x2,y2,exl2,exh2,eyl2,eyh2);
    gr2->SetMarkerColor(1);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(0.95);
    gr2->SetFillColor(kGray+3);
    gr2->SetLineColor(kGray+3);
    gr2->SetLineWidth(0);
    gr2->SetMarkerSize(0);
    gr2->SetFillStyle(3013);


    TGraphAsymmErrors *gr3 = new TGraphAsymmErrors(Nbins,x3,y3,exl3,exh3,eyl3,eyh3);
    gr3->SetTitle("");
    gr3->SetMarkerSize(0);
    gr3->SetFillColor(kGray+3);
    gr3->SetFillStyle(3013);
    gr3->GetXaxis()->SetLabelSize(0.19);
    gr3->GetXaxis()->SetTitleSize(0.21);
    gr3->GetXaxis()->SetTitleOffset(1.0);
    gr3->GetXaxis()->SetLabelOffset(0.03);
    gr3->GetYaxis()->SetLabelSize(0.19);
    gr3->GetYaxis()->SetLabelOffset(0.006);
    gr3->GetYaxis()->SetTitleSize(0.21);
    gr3->GetYaxis()->SetTitleOffset(0.35);
    gr3->GetYaxis()->SetNdivisions(503);
 
    TGraphAsymmErrors *gr4 = new TGraphAsymmErrors(Nbins,x2,y2,exl2,exh2,eyl4,eyh4);
    gr4->SetTitle("");
    gr4->SetMarkerSize(0);
    gr4->SetFillColor(kAzure+1);
    gr4->SetFillStyle(3021);
    gr4->GetXaxis()->SetLabelSize(0.19);
    gr4->GetXaxis()->SetTitleSize(0.21);
    gr4->GetXaxis()->SetTitleOffset(1.0);
    gr4->GetXaxis()->SetLabelOffset(0.03);
    gr4->GetYaxis()->SetLabelSize(0.19);
    gr4->GetYaxis()->SetLabelOffset(0.006);
    gr4->GetYaxis()->SetTitleSize(0.21);
    gr4->GetYaxis()->SetTitleOffset(0.35);
    gr4->GetYaxis()->SetNdivisions(503);

  htotal->Draw("hist same");
//  htotaldata->Draw("E0P0same");
  gr4->Draw("E2 same");
//  gr->Draw("E0P0Zsame");  nodata
  gr2->Draw("e2 same");
  if(signal) {
//       hfourlepbestmass_4l_afterSel_new_sig_350->Draw("hist same"); 
       hfourlepbestmass_4l_afterSel_new_sig_400->Draw("hist same");
//       hfourlepbestmass_4l_afterSel_new_sig_450->Draw("hist same");
//       hfourlepbestmass_4l_afterSel_new_sig_500->Draw("hist same");
//       hfourlepbestmass_4l_afterSel_new_sig_550->Draw("hist same");
       hfourlepbestmass_4l_afterSel_new_sig_600->Draw("hist same");
//       hfourlepbestmass_4l_afterSel_new_sig_650->Draw("hist same");
//       hfourlepbestmass_4l_afterSel_new_sig_700->Draw("hist same");
//       hfourlepbestmass_4l_afterSel_new_sig_750->Draw("hist same");
       hfourlepbestmass_4l_afterSel_new_sig_800->Draw("hist same");
  }
  legend->AddEntry(rgr_mc,"syst. unc.","F");
  legend->AddEntry(gr2,"stat. unc.","F");

  legend->Draw("same");
  ll->Draw("same");

  
  std::string saveaspdfzoom = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom_log.pdf" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom.pdf";
  std::string saveaspngzoom = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom_log.png" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom.png";
  std::string saveasepszoom = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom_log.eps" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom.eps";
//   c1->SaveAs(saveaspdfzoom.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.pdf"*/);
//   c1->SaveAs(saveaspngzoom.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.png"*/);
//   c1->SaveAs(saveasepszoom.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.eps"*/);


  // full spectrum

  gPad->RedrawAxis();

  //  c1->Update();
  double canvasratio = 0.2;
  c1->SetBottomMargin(canvasratio + (1-canvasratio)*c1->GetBottomMargin()-canvasratio*c1->GetTopMargin());
  //cout << "Canvas= " << canvasratio + (1-canvasratio)*c1->GetBottomMargin()-canvasratio*c1->GetTopMargin() << endl;
 
  // Ratio: data / total bkg 
  canvasratio = 0.16;
  TPad *ratioPad = new TPad("BottomPad","",0,0,1,1);
  ratioPad->SetTopMargin((1-canvasratio) - (1-canvasratio)*ratioPad->GetBottomMargin()+canvasratio*ratioPad->GetTopMargin());
  ratioPad->SetFillStyle(4000);
  ratioPad->SetFillColor(4000);
  ratioPad->SetFrameFillColor(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->SetFrameBorderMode(0);
  ratioPad->SetTicks(1,1);
  //ratioPad->SetLogx();
  ratioPad->Draw();
  ratioPad->cd();

  //TH2F *hframe2= new TH2F("hframe2","hframe2",6000, 0., 2.2, 1000, 0.5, 2.);// iso
 
  hframe2->GetYaxis()->SetLabelSize(0.020);
//  hframe2->GetXaxis()->SetLabelSize(0.020);
  //  hframe2->GetYaxis()->SetTitleSize(0.047);
  hframe2->SetYTitle("Data/MC");
  //  hframe2->GetYaxis()->SetRangeUser(-10,10);
  hframe2->GetYaxis()->SetNdivisions(503);
  //hframe2->GetXaxis()->SetTitleOffset(1.25);
  hframe2->Draw("");

  htotalHistoRatio->SetMarkerStyle(20);
  htotalHistoRatio->SetMarkerSize(0.95);
  htotalHistoRatio->SetMarkerColor(kBlack);
//  htotalHistoRatio->Draw("P0 E0 same");
//nodata  rgr->Draw("P0 E0 Z same");
  ratiostaterr->Draw("e2 same");
  if(un_mc) rgr_mc->Draw("E2 same");
  gr3->Draw("E2 same");

  c1->Update();

  std::string saveaspdfratio = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.pdf" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.pdf";
  std::string saveaspngratio = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.png" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.png";
  std::string saveasepsratio = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.eps" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.eps";
  std::string saveasrootratio= useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.root": "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.root";

  cout << saveasrootratio.c_str() << endl;
//  c1->SaveAs(saveaspdfratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.pdf"*/);
  c1->SaveAs(saveaspngratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.png"*/);
//  c1->SaveAs(saveasepsratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.eps"*/);
//  c1->SaveAs(saveasrootratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.root"*/);
  

  
  // Write final histogram in a file 
/*  char htotalfinal[300];
  sprintf(htotalfinal,"plots/htotalfinal_%s_%s.root",histlabel.c_str(),whichchannel.c_str());
  TFile *file1 = new TFile(htotalfinal, "RECREATE");
  file1->cd();
  htotalHisto->Write();
  htotaldata->Write();
  file1->Write();
  file1->Close();
*/
}

void PlotStack4l2b::setSamplesNamesZ()
{
  
  std::ifstream infile;
  infile.open(inputfile.c_str());
  
  if (inputfile.find("2011")<500) { 
    whichenergy="7TeV";
    whichsample="7TeV";
  }
  else if (inputfile.find("RunI_")<500){
    whichenergy="RunI_";
    whichsample="8TeV";
  }
  else if (inputfile.find("2015")<500){
    whichenergy="13TeV";
    whichsample="13TeV";
  }
  else if (inputfile.find("2016")<500){
    whichenergy="13TeV";
    whichsample="13TeV";
  }

  
  cout << "Doing plot for " << whichenergy.c_str() << "  and " << whichsample.c_str() << endl;

  whichchannel="4mu2b_sig";
  histosdir="4l2b";

  Vdatasetnamedata.clear();
  Vdatasetnamebkg.clear();
  Vdatasetnamesig.clear();

  std::string inputfilename;

  while(std::getline(infile,inputfilename)){
    cout << "Reading " << inputfilename.c_str() << endl;
    
    // DATA  
    
   
    if(inputfilename.find("_DoubleMuon")<200){ // as many times as it occurs in the input file
      Vdatasetnamedata.push_back(inputfilename);
      Vlabeldata.push_back("Double Muon - 2016");
      Vxsectiondata.push_back(1.); //pb
    }  
             
    if(inputfilename.find("_DoubleEG")<200){ // as many times as it occurs in the input file
      Vdatasetnamedata.push_back(inputfilename);
      Vlabeldata.push_back("Double EGamma - 2016");
      Vxsectiondata.push_back(1.); //pb
    }   
    
    if(inputfilename.find("_SingleElectron")<200){ // as many times as it occurs in the input file
      Vdatasetnamedata.push_back(inputfilename);
      Vlabeldata.push_back("Single Electron - 2016");
      Vxsectiondata.push_back(1.); //pb
    }    
    
    if(inputfilename.find("_MuonEG")<200){ // as many times as it occurs in the input file
      Vdatasetnamedata.push_back(inputfilename);
      Vlabeldata.push_back("Muon EGamma - 2016");
      Vxsectiondata.push_back(1.); //pb
    }    
       
    //Z+X from data
    if(inputfilename.find("Z+X")<85){
      Vdatasetnamebkgdata.push_back(inputfilename);
      Vlabelbkgdata.push_back("Z+X");
      Vxsectionbkgdata.push_back(1.); //pb
    }
    
    // background higgs
    else if(inputfilename.find("DYJetsToLL")<100 || 
	    inputfilename.find("TT_TuneCUETP8M1")<100   || 
	    inputfilename.find("WJetsToLNu")<100     ||
	    inputfilename.find("WWTo2L2Nu")<100  ||  
            inputfilename.find("ZZTo4L")<100 ||
            inputfilename.find("TTZToLLNuNu")<100 ||
            inputfilename.find("GluGluTo4L")<100 || 
            inputfilename.find("hTo4l")<100 ){ 
      //Vdatasetnamebkg.push_back(inputfilename);
      

      if(inputfilename.find("DYJetsToLL")<200){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("DY");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kAzure+2);
      }

      if(inputfilename.find("GluGluTo4L")<200 ){
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("ggZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kAzure-1);
        //cout << "VBF H" << endl;
      }


      if(inputfilename.find("hTo4l")<200 ){
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("hToZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kGray+1);
        //cout << "VBF H" << endl;
      }

      if(inputfilename.find("TT_TuneCUETP8M1")<200 ){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("TT");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-6);
	//cout << "VBF H" << endl;
      }

      if(inputfilename.find("TTZToLLNuNu")<200 ){
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("TTZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kTeal-6);
        //cout << "VBF H" << endl;
      }
      
      if( (inputfilename.find("WJetsToLNu")<200)){
	  Vdatasetnamebkg.push_back(inputfilename);
	  Vlabelbkg.push_back("W");
	  Vxsectionbkg.push_back(1.); //pb      
          Vcolorbkg.push_back(kSpring);
      }


      if(inputfilename.find("WWTo2L2Nu")<200){
	Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("WW");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kCyan+3);
        //cout << "VH" << endl;
      }

      if(inputfilename.find("ZZTo4L")<200){
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("ZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kAzure+2);
      }
 
    }
    // BACKGROUND from other sources
    else{  // provided that everything that is neither data nor signal is background
      //Vdatasetnamebkg.push_back(inputfilename);
      if(inputfilename.find("sig_350")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_350");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(40);
      }
      if(inputfilename.find("sig_400")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_400");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(6);
      }
      if(inputfilename.find("sig_450")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_450");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(41);
      }
      if(inputfilename.find("sig_500")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_500");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(42);
      }
      if(inputfilename.find("sig_550")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_550");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(43);
      }
      if(inputfilename.find("sig_600")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_600");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(7);
      }
      if(inputfilename.find("sig_650")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_650");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(45);
      }
      if(inputfilename.find("sig_700")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_700");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(46);
      }
      if(inputfilename.find("sig_750")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_750");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(47);
      }
      if(inputfilename.find("sig_800")<200){
        Vdatasetnamesig.push_back(inputfilename);
        Vlabelsig.push_back("sig_800");
        Vxsectionsig.push_back(1.); //pb
        Vcolorsig.push_back(kYellow+1);
      }

      // WZ
      if(inputfilename.find("WZJToLLLNu")<200){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("WZ");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kCyan-2);
      } 
      
      // TTTo2L2Nu2B_8TeV-powheg-pythia6 
      if(inputfilename.find("TTTo2L2Nu")<200){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("t#bar{t}");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-6);
      } 
      
      // Single Top
      // T_TuneZ2_s-channel_8TeV-madgraph
      if(inputfilename.find("T_s-channel")<200 || inputfilename.find("T_TuneZ2_s-channel")<200 ){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("T_TuneZ2_s-channel_8TeV-powheg-tauola");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      } 
      
      // T_TuneZ2_t-channel_8TeV-madgraph
      if(inputfilename.find("T_t-channel")<200 || inputfilename.find("T_TuneZ2_t-channel")<200){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("T_TuneZ2_t-channel_8TeV-powheg-tauola");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      } 
      
      // T_TuneZ2_tW-channel_8TeV-madgraph
      if(inputfilename.find("T_tW")<200 || inputfilename.find("T_TuneZ2_tW-channel")<200 ){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("T_TuneZ2_tW-channel-DR_8TeV-powheg-tauola");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      } 

      // Tbar_TuneZ2_s-channel_8TeV-madgraph
      if(inputfilename.find("Tbar_s-channel")<200 || inputfilename.find("Tbar_TuneZ2_s-channel")<200 ){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Tbar_TuneZ2_s-channel_8TeV-powheg-tauola");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      } 
      
      // Tbar_TuneZ2_t-channel_8TeV-madgraph
      if(inputfilename.find("Tbar_t-channel")<200 || inputfilename.find("Tbar_TuneZ2_t-channel")<200 ){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Tbar_TuneZ2_t-channel_8TeV-powheg-tauola");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      } 
      
      // Tbar_TuneZ2_tW-channel-DR_8TeV-madgraph
      if(inputfilename.find("Tbar_tW-channel")<200 || inputfilename.find("Tbar_TuneZ2_tW-channel")<200 ){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Tbar_TuneZ2_tW-channel-DR_8TeV-powheg-tauola");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      } 


    
    }

  }
  
  infile.close();
 
}


void PlotStack4l2b::printnumbers(char name[400], TH1F *h){

  int nSkim_MC=h->GetBinContent(1);
  int nZLQpt_MC=h->GetBinContent(2);
  int nZLQptLI_MC=h->GetBinContent(3);
  int n4LQptLI_MC=h->GetBinContent(4);
  int n4lAfterPresel_MC=h->GetBinContent(5);
  int nIso_MC=h->GetBinContent(6);
  int nIsopt_MC=h->GetBinContent(7);
  int nIP_MC=h->GetBinContent(8);
  int nKin_MC=h->GetBinContent(9);

  cout << "Sample= " <<  name << endl;
  cout << "Skim MC= " <<  nSkim_MC 
       << " \n ZLQpt_MC="         << nZLQpt_MC
       << " \n ZLQptLI_MC="       << nZLQptLI_MC 
       << " \n 4LQptLI_MC="       << n4LQptLI_MC 
       << " \n BestCand_MC="      << n4lAfterPresel_MC 
       << " \n Iso_MC="           << nIso_MC 
       << " \n Iso(pt)_MC="       << nIsopt_MC 
       << " \n IP_MC="            << nIP_MC
       << " \n Kin_MC="           << nKin_MC
       << endl;
}

/*
void PlotStack4l::getMassWindow(float Higgsm){

  std::string fileLoc = "/cmshome/nicola/tmp/test/Paper/last/last/CMSSW_5_3_9/src/Higgs/Higgs_CS_and_Width/txtFiles";
  HiggsCSandWidth myCSW = HiggsCSandWidth(fileLoc);

  double widthHVal =  myCSW.HiggsWidth(0,Higgsm);
  double windowVal = max( widthHVal, 1.0);
  double lowside = 100.0;
  double highside = 1000.0;
  
  if (Higgsm >= 275){
     lowside = 180.0;
     highside = 650.0;}
  if (Higgsm >= 350){
     lowside = 200.0;
     highside = 900.0;}
   if (Higgsm >= 500){
     lowside = 250.0;
     highside = 1000.0;}
   if (Higgsm >= 700){
     lowside = 350.0;
     highside = 1400.0;}
  
   float low_M = max( (Higgsm - 20.*windowVal), lowside);
   float high_M = min( (Higgsm + 15.*windowVal), highside);
  
   cout << "Mass window: " << low_M << " " << high_M << endl;

}
*/
