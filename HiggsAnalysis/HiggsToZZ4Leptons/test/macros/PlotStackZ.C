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

class PlotStackZ{
  
public: 
  PlotStackZ();
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

PlotStackZ::PlotStackZ(){
  //TSystem LoadLib;
  //LoadLib.Load("/cmshome/nicola/slc6/MonoHiggs/Analysis13TeV/CMSSW_7_2_0/lib/slc6_amd64_gcc481/libHiggsHiggs_CS_and_Width.so");
  //getMassWindow(500.);
    
  inputfile="filelist_2e2b_samept_input.txt";

  setSamplesNamesZ(); 
  cout << "\t Analysing samples for " << whichchannel << " analysis" << endl; 


  //WARNING: depending on histolabel, modify the declaration and the settings of hframe below
  //also choose a sensible value for nRebin

  //std::string histolabel = "hPUvertices";    // numPU
  //std::string histolabel = "hPUvertices_ReWeighted";    // numPY reweightesmC
  // Step 0
  //std::string histolabel = "hPtLep_0"; // pt not cuts (only pt>5)
  //std::string histolabel = "hIsoLep_0"; // PF isol _ delta B
  //std::string histolabel = "hSipLep_0"; // SIP
  //std::string histolabel = "hDxyLep_0"; // d_xy
  //std::string histolabel = "hDzLep_0"; // d_z
  //std::string histolabel = "hMuHitLep_0"; // n. muon hits
  //std::string histolabel = "hPxHitLep_0"; // n. pixel hits
  //std::string histolabel = "hTKLayLep_0"; // n. tracker layers
  //std::string histolabel = "hTKIsoLep_0"; // track isolation = sumpT

  // ****** Standard candle: Z1 selection: step 3 ******

  //std::string histolabel = "hMZ_3";    // Z mass 
  //std::string histolabel = "hMZBB_3";    // Z mass 
  //std::string histolabel = "hMZEE_3";    // Z mass 
  //std::string histolabel = "hPtZ_3"; // Z pt 
  //std::string histolabel = "hYZ_3"; // Z Y

  //std::string histolabel = "hPtLep_3";    // pT lepton from Z 
  //std::string histolabel = "hIsoLep_3";    // isolation of lowest pT lepton from Z 
  //std::string histolabel = "hSipLep_3";  // sip the lowest pT lepton from Z
  //std::string histolabel = "hIpLep_3";   // IP of the lowest pT lepton from Z
  //std::string histolabel = "hIpErLep_3";   // Iperror of the lowest pT lepton from Z

  //std::string histolabel = "hIso_3";    // worst isolation value of lepton not coming from Z1
  //std::string histolabel = "hSip_3";  // worst sip value of lepton not coming from Z1
  //std::string histolabel = "hIp_3";   // worst IP value of lepton not coming from Z1

  //std::string histolabel = "hMjj_3"; // mass of di-jet for VBF analysis
  //std::string histolabel = "hDjj_3"; // delta eta between jets for VBF analysis
  //std::string histolabel = "hVD_3"; // Fisher discriminant for VBF analysis

  //std::string histolabel = "hPFMET_3"; // PFMET

  // ****** After cuts on ,Z1, mZ2 and pT >20/10: step 5 ******

//  std::string histolabel = "hMZ1_5";    // Z1 mass   
  //std::string histolabel = "hMZ2_5";    // Z2 mass 
  //std::string histolabel = "hM4l_5";    // 4l mass 

  //std::string histolabel = "hIso_5";    // worst isolation 
  //std::string histolabel = "hSip_5";  // worst sip 
  //std::string histolabel = "hIp_5";   // worst IP

//   std::string histolabel = "Mbb_6";
//  std::string histolabel = "bdiscr_5_lead";
//  std::string histolabel = "bdiscr_5_sub";
//    std::string histolabel = "Mjj_6";
//    std::string histolabel = "hPtBot_8";
//    std::string histolabel = "hPtBot_7";
//   std::string histolabel = "hPtZ1_7";
//  std::string histolabel = "hMZ1_7";
  // After full selection
  //std::string histolabel = "hM4l_7"; // 4l mass after full selection but m4l > 70

  //std::string histolabel = "hM4l_8"; // 4l mass after full selection
  //std::string histolabel = "hM4l_9"; // 4l mass after full selection
  //std::string histolabel = "hM4l_8_100_800"; // 4l mass in the range [100,800] after full selection
  //std::string histolabel = "hMZ1_8"; // Z1 mass after full selection
  //std::string histolabel = "hMZ2_8"; // Z2 mass after full selection
  //std::string histolabel = "hMZ1_noFSR_8"; // Z1 mass after full selection without FSR recovery
  //std::string histolabel = "hMZ2_noFSR_8"; // Z2 mass after full selection without FSR recovery
//  std::string histolabel = "hPtZ1_5"; // Z1 pt after full selection
//   std::string histolabel = "hPtLep1_7";
//   std::string histolabel = "hPtLep2_7"; 
// std::string histolabel = "hPtLep1_8";
//   std::string histolabel = "hPtLep2_8"; 
//   std::string histolabel = "hEtaLep1_8";
//     std::string histolabel = "hEtaLep2_8";
//     std::string histolabel = "hPhiLep1_7";
//   std::string histolabel = "hPhiLep2_7";

//   std::string histolabel = "hEtaLep1_7";
//   std::string histolabel = "hEtaLep2_7";
//   std::string histolabel = "hIsoLep1_7";
//   std::string histolabel = "hIsoLep2_7";
//    std::string histolabel = "hYZ1_5";
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
//   std::string histolabel = "hPtBot_8";
//   std::string histolabel = "hEtaJet_7";
//   std::string histolabel = "hEtaJet_8";
//  std::string histolabel = "hNjets_8";
  //std::string histolabel = "hDjj_8"; // delta eta between jets for VBF analysis
 // std::string histolabel = "hMjj_8"; // dimass between jets for VBF analysis
//  std::string histolabel = "hN_loose_e";
//  std::string histolabel = "hN_loose_mu"; 
//   std::string histolabel = "hN_good_ele";
//   std::string histolabel = "hN_good_mu";
//    std::string histolabel = "hMZ1_BB_5";
// std::string histolabel = "hPtZ1_BB_5";
//  std::string histolabel = "hYZ1_BB_5";
//    std::string histolabel = "hMZ1_5_EB";
//  std::string histolabel = "hPtZ1_EB_5";
//  std::string histolabel = "hYZ1_EB_5";
//  std::string histolabel = "hMZ1_EE_5";
//  std::string histolabel = "hPtZ1_EE_5";
//  std::string histolabel = "hYZ1_EE_5";

  useLogY = true;
  useLogX = false;

  useDYJets=true;
  useDYJetsFromData=false;

  nRebin=5;
//  std::cout << "Histogram label is= " << histolabel << std::endl;
 
  std::string label5[8]={"hMZ1_5","hPtZ1_5","hYZ1_5","hEtaJet_8","hMZ1_7","hPtZ1_7","hYZ1_7","ptbb_6"};
  std::string label6[2]={"hPtJet_8","hPtBot_8"};
  std::string label2[2]={"Mjj_6","Mbb_6"};
  std::string label1[2]={"hNjets_8","hNbjets"};
  std::string label15[4]={"hPtLep1_7","hPtLep2_7","hEtaLep1_7","hEtaLep2_7"};

  un_mc=false;
 
  // Final yields
  system("mkdir plots");

  Char_t yieldsOUT[500];
  sprintf(yieldsOUT,"plots/amcatnlo/yields_%s_%s.txt",whichchannel.c_str(),whichenergy.c_str());
  
  // Execute the analysis
//  plotmZ(histolabel);
  for(int i=0; i<1; i++){
    plotmZ(label5[i]);
  }

  

}

void PlotStackZ::plotmZ(std::string histlabel){
  
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



  TLegend *leg0 = new TLegend(0.6,0.40,0.8,0.90,NULL,"brNDC");
  leg0->SetTextSize(0.020);
  leg0->SetLineColor(0);
  leg0->SetLineWidth(1);
  leg0->SetFillColor(kWhite);
  leg0->SetBorderSize(0);
 
  TLegend* legend = new TLegend(0.70, 0.80, 0.9, 0.9);
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

  //hframe= new TH2F("hframe","hframe",80,70.,810.,500,0.004,25.);// 4l analysis mass
  //hframe= new TH2F("hframe","hframe",80,70.,900.,500,0.004,25.);// 4l analysis mass // animate gif 
  //hframe2= new TH2F("hframe","hframe",80,70.,182.,1000, 0.5, 4.);// 4l analysis mass

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

  if (histlabel.find("hM4l_7")<10 && 
      (whichchannel.find("4e")<20 || whichchannel.find("4#mu")<20 || whichchannel.find("2e2#mu")<20)){
    hframe= new TH2F("hframe","hframe",80,70.,182.,500,0.004,12.);// 4l analysis mass nrebin=3 GeV
    hframe2= new TH2F("hframe2","hframe2",6000, 70., 182., 1000, 0.5, 2.);// 
  }
  if (histlabel.find("hM4l_7")<10 && whichchannel.find("4l")<20){
    hframe= new TH2F("hframe","hframe",80,70.,182.,500,0.004,15.);// 4l analysis mass nrebin=3
    hframe2= new TH2F("hframe2","hframe2",6000, 70., 182., 1000, 0.5, 2.);// 
  }

  if (histlabel.find("hM4l_8")<10 && (whichchannel.find("4l")<20)){
    hframe= new TH2F("hframe","hframe",80,70.,182.,500,0.00000004,350000000.);// 4l analysis mass nrebin=3
    hframe2= new TH2F("hframe2","hframe2",6000, 70., 182., 1000, 0.5, 2.);// 
  }


  if (histlabel.find("hM4l_9")<10 && whichchannel.find("4#mu")<20){
    hframe= new TH2F("hframe","hframe",80,70.,182.,500,0.000004,35.);// 4l analysis mass nrebin=3
    hframe2= new TH2F("hframe2","hframe2",6000, 70., 182., 1000, 0.5, 2.);// 
  }
  

  //hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,2600000.);// mZ1 ee/mumu
  //hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 ee/mumu
  
  if (histlabel.find("hMZ_3")<10 ){
    hframe= new TH2F("hframe","hframe",100,20.,120.,500,0.01,1700000000.);// mZ1 mumu
    hframe2= new TH2F("hframe2","hframe2",6000, 20., 120., 1000, 0.5, 2.);// mZ1 mumu
  }
  //hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,400000.);// mZ1 mumu  7TeV
  //hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 mumu 7TeV
  
  if (histlabel.find("hMZ_3")<10 && whichchannel.find("4e")<20){
    hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,1100000.);// mZ1 ee                                             
    hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 ee
  }

  if (histlabel.find("hMZ_3")<10 && whichchannel.find("2e2#mu")<20){
    hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,1100000.);// mZ1 ee                                             
    hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 ee
  }

  if (histlabel.find("hPtLep_3")<10){
    hframe= new TH2F("hframe","hframe",80,5.,300.,500,0.001,10000.);// pT                                                                 
    hframe2= new TH2F("hframe2","hframe2",80,5.,200.,500, 0.5, 2.);// pT                                                                   
  }
  
  if (histlabel.find("hIsoLep_3")<10){
    hframe= new TH2F("hframe","hframe",80,0.,5.,500,0.001,100000.);// Isolation                                                            
    hframe2= new TH2F("hframe2","hframe2",80,0.,5.,500, 0.5, 2.);// Isolation                                                              
  }

   if (histlabel.find("hSipLep_3")<10){
    hframe= new TH2F("hframe","hframe",600, 0., 5., 1000, 0.0004, 1000000000.);// sip
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 5., 1000, 0.5, 2.);// sip
  }
  

  //hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,1100000.);// mZ1 ee
  //hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 ee
  //hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,700000.);// mZ1 ee BB
  //hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 ee BB
  //hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,100000.);// mZ1 ee EE
  //hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 ee EE
  
  //hframe= new TH2F("hframe","hframe",80,60.,120.,500,0.4,300000.);// mZ1 ee 7TeV
  //hframe2= new TH2F("hframe2","hframe2",6000, 60., 120., 1000, 0.5, 2.);// mZ1 ee 7TeV

  if (histlabel.find("hMjj_3")<10){
    hframe= new TH2F("hframe","hframe",600,20.,500.,600,0.0004,10E7);//mass jet jet
    hframe2= new TH2F("hframe2","hframe2",6000, 20., 500., 1000, 0.5, 2.);// mass jet jet
  }

  if (histlabel.find("hDjj_3")<10){
    hframe= new TH2F("hframe","hframe",600,0.,10.,600,0.0004,10E7);//delta eta jet jet
    hframe2= new TH2F("hframe2","hframe2",600, 0., 10., 1000, 0.5, 2.);// delta eta jet jet
  }

  if (histlabel.find("hVD_3")<10){
    hframe= new TH2F("hframe","hframe",600,0.,2.,600,4E-4,5E11);// fisher
    hframe2= new TH2F("hframe2","hframe2",600, 0., 2., 1000, 0.5, 2.);// fisher
  }

  // hframe= new TH2F("hframe","hframe",Nbins*nRebin, Xmin, Xmax, Nbins*nRebin, 0.0004, Ymax);//mass

  if (histlabel.find("hSip_3")<10){
    hframe= new TH2F("hframe","hframe",600, 0., 10., 1000, 0.0004, 1000000000.);// sip
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 10., 1000, 0.5, 2.);// sip
  }
  if (histlabel.find("hIp_3")<10){
    hframe= new TH2F("hframe","hframe",600, 0., 2., 1000, 0.000004, 1000000000.);// Ip
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 2., 1000, 0.5, 2.);// Ip
  }
  //TH2F *hframe= new TH2F("hframe","hframe",600, 0., 10., 1000, 0.04, 10000000.);// chi2
  //TH2F *hframe= new TH2F("hframe","hframe",600, 0., 1., 1000, 0.004, 6.);// prob
  // TH2F *hframe= new TH2F("hframe","hframe",600, 0., 0.35, 1000, 0.04, 10000000.);// prob

  if (histlabel.find("hIso_3")<10){
    hframe= new TH2F("hframe","hframe",6000, 0., 3., 1000, 0.0004, 1000000000.);// iso
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 3., 1000, 0.5, 2.);// iso
  }

  if (histlabel.find("hPFMET_3")<10){
    hframe= new TH2F("hframe","hframe",1000, 0., 600., 1000, 0.000001, 10000.);// PFMET
    hframe2= new TH2F("hframe2","hframe2",1000, 0., 600., 1000, 0.5, 2.);// PFMET
  }
  

  if (histlabel.find("hMZ1_5")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,1.0,60000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",60, 60., 120., 10000, 0.8, 1.2);// mZ1 
    hframe2->SetXTitle("M_{ll} [GeV]");
  }

  if (histlabel.find("hPtZ1_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,1.0,80000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 200., 1000, 0.8, 1.2);// mZ1 
    hframe2->SetXTitle("Pt_{ll} [GeV]");
  }
  if (histlabel.find("ptbb_6")<10){
    hframe= new TH2F("hframe","hframe",90,20.,200.,500,0.01,8000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",90, 20., 200., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Pt_{bb} [GeV]");
  }

  if (histlabel.find("hPtLep1_7")<10){
    hframe= new TH2F("hframe","hframe",170,30.,200.,500,1.0,80000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",170, 30., 200., 500, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtLep2_7")<10){
    hframe= new TH2F("hframe","hframe",180,20.,200.,500,1.0,8000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",180, 20., 200., 500, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtLep1_8")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,30000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtLep2_8")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,30000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hYZ1_5")<10){
    hframe= new TH2F("hframe","hframe",100,-2.5,2.5,500,1.0,100000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -2.5, 2.5, 1000, 0.8, 1.2);// mZ1
    hframe2->SetXTitle("Y_{ll} [GeV]"); 
  }

  if (histlabel.find("hEtaLep1_8")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,10000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hEtaLep2_8")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,10000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hEtaLep1_7")<10){
    hframe= new TH2F("hframe","hframe",50,-2.5,2.5,500,1.0,1000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, -2.5, 2.5, 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hEtaLep2_7")<10){
    hframe= new TH2F("hframe","hframe",50,-2.5,2.5,500,1.0,1000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, -2.5, 2.5, 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPhiLep1_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,4.,500,1.0,1000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 4., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPhiLep2_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,4.,500,1.0,1000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 4., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hIsoLep1_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,1.,500,1.0,20000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 1., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hIsoLep2_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,1.,500,1.0,20000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 1., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hMZ1_7")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,1.0,200000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",60, 60., 120., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("M_{ll} [GeV]");
  }

  if (histlabel.find("hPtZ1_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,1.0,100000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 200., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Pt_{ll} [GeV]");
  }

  if (histlabel.find("hYZ1_7")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,100000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
    hframe2->SetXTitle("Y_{ll} [GeV]");
  }


  if (histlabel.find("Mbb_6")<10){
    hframe= new TH2F("hframe","hframe",80,20.,400.,500,0.1,100000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",80, 20., 400., 1000, 0.7, 1.3);// mZ1
    hframe2->SetXTitle("M_{bb} [GeV]");
  }

  if (histlabel.find("Mjj_6")<10){
    hframe= new TH2F("hframe","hframe",80,20.,400.,500,0.1,1000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",80, 20., 400., 1000, 0.8, 1.2);// mZ1 
    hframe2->SetXTitle("M_{jj} [GeV]");
  }


  if (histlabel.find("bdiscr_5_lead")<10){
    hframe= new TH2F("hframe","hframe",20,0.,1.,500,1.0,1000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",20, 0., 1., 1000, 0.7, 1.3);// mZ1 
  }

  if (histlabel.find("bdiscr_5_sub")<10){
    hframe= new TH2F("hframe","hframe",20,0.,1.,500,1.0,300000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",20, 0., 1., 1000, 0.7, 1.3);// mZ1 
  }

  if (histlabel.find("hPtJet_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,400.,500,1.0,200000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",1000, 0., 400., 1000, 0.7, 1.3);// mZ1 
  }
 
  if (histlabel.find("hPtJet_8")<10){
    hframe= new TH2F("hframe","hframe",190,20.,400.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",190, 20., 400., 1000, 0.7, 1.3);// mZ1 
   
  }

  if (histlabel.find("hPtBot_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,400.,500,1.0,200000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",1000, 0., 400., 1000, 0.7, 1.3);// mZ1 
  }

  if (histlabel.find("hPtBot_8")<10){
    hframe= new TH2F("hframe","hframe",190,20.,400.,500,1.0,2000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",190, 20., 400., 1000, 0.7, 1.3);// mZ1 
  }
 
  if (histlabel.find("hEtaJet_7")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hEtaJet_8")<10){
    hframe= new TH2F("hframe","hframe",100,-2.5,2.5,500,1.0,4000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -2.5, 2.5, 1000, 0.7, 1.3);// mZ1 
  }

  if (histlabel.find("hNjets_8")<10){
    hframe= new TH2F("hframe","hframe",8,-0.5,7.5,500,0.1,100000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",8, -0.5, 7.5, 1000, 0.7, 1.3);// mZ1
    hframe2->SetXTitle("Njets"); 
  }

  if (histlabel.find("hNbjets")<10){
    hframe= new TH2F("hframe","hframe",6,-0.5,5.5,500,0.1,100000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6, -0.5, 5.5, 1000, 0.7, 1.3);// mZ1 
    hframe2->SetXTitle("Nbjets");
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
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,2000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hN_good_mu")<10){
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,100000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.7, 1.3);// mZ1 
  }

  if (histlabel.find("hMZ1_8")<10 && whichchannel.find("4#mu")<20){
    hframe= new TH2F("hframe","hframe",80,40.,200.,500,0.0001,100000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 40., 160., 1000, 0.5, 1.5);// mZ1 
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
    hframe= new TH2F("hframe","hframe",200, 0., 200., 1000, 0.000001, 6000000.);// PFMET
    hframe2= new TH2F("hframe2","hframe2",200, 0.,200., 1000, 0.7, 1.3);// PFMET
  }

  if (histlabel.find("hPFMET_9")<10){
    //hframe= new TH2F("hframe","hframe",1000, 0., 1000., 1000, 0.0000004, 50000.);// PFMET
    hframe= new TH2F("hframe","hframe",1000, 0., 1000., 1000, 0.000001, 100.);// PFMET
    hframe2= new TH2F("hframe2","hframe2",1000, 0.,1000., 1000, 0.5, 1.5);// PFMET
  }

  if (histlabel.find("hMjj_8")<10){
    hframe= new TH2F("hframe","hframe",600,20.,500.,600,0.000004,50000);//mass jet jet
    hframe2= new TH2F("hframe2","hframe2",6000, 20., 500., 1000, 0.5, 1.5);// mass jet jet
  }

  if (histlabel.find("hDjj_8")<10){
    hframe= new TH2F("hframe","hframe",600,0.,10.,600,0.000004,10E4);//delta eta jet jet
    hframe2= new TH2F("hframe2","hframe2",600, 0., 10., 1000, 0.5, 2.);// delta eta jet jet
  }

  if (histlabel.find("hMZ1_BB_5")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,1.0,3000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",60, 60., 120., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtZ1_BB_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,1.0,2000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 200., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hYZ1_BB_5")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,1500000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hMZ1_5_EB")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,1.0,3000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",60, 60., 120., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtZ1_EB_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,1.0,2000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 200., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hYZ1_EB_5")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,1000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hMZ1_EE_5")<10){
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,1.0,600000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",60, 60., 120., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtZ1_EE_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,1.0,300000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 200., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hYZ1_EE_5")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,600000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
  }

  //TH2F *hframe= new TH2F("hframe","hframe",6000, 0., 200., 1000, 0.004, 700000.);// ptZ

  if (nRebin==1) hframe->SetYTitle("Events/8 GeV");
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
  

  if (histlabel.find("hM4l_7")<10 || histlabel.find("hM4l_8")<10  || histlabel.find("hM4l_9")<10  ){
    sprintf(histotitle,"m_{%s} [GeV]",whichchannel.c_str());
    hframe->SetXTitle(histotitle);
  }
  
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
  if (histlabel.find("hIso_3")<10 && whichchannel.find("4e")<20)   hframe->SetXTitle("electron worst rel. iso.");
  if (histlabel.find("hIso_3")<10 && whichchannel.find("4#mu")<20) hframe->SetXTitle("muon worst iso.");
  if (histlabel.find("hIso_3")<10 && whichchannel.find("2e2#mu")<20) hframe->SetXTitle("lepton worst iso.");
  //hframe->SetXTitle("lepton worst iso.");
  if (histlabel.find("hSip_3")<10 && whichchannel.find("4e")<20) hframe->SetXTitle("electron worst SIP_{3D}");
  if (histlabel.find("hSip_3")<10 && whichchannel.find("4#mu")<20) hframe->SetXTitle("muon worst sign. 3DIP");
  if (histlabel.find("hSip_3")<10 && whichchannel.find("2e2#mu")<20) hframe->SetXTitle("lepton worst sign. 3DIP");
  //hframe->SetXTitle("lepton worst sign. 3DIP");
  // hframe->SetXTitle("Prob(#chi^{2}/ndof.) of 4#mu vertex fitter (best 4#mu comb.)");
  // hframe->SetXTitle("#chi^{2}/ndof. of 4#mu vertex fitter (best 4#mu comb.)");
  if (histlabel.find("hMjj_3")<10 || histlabel.find("hMjj_8")<10) hframe->SetXTitle("di jet mass"); //Mjj
  if (histlabel.find("hDjj_3")<10 || histlabel.find("hDjj_3")<10) hframe->SetXTitle("#Delta#eta jets"); //deltaetajj
  //if (histlabel.find("hVD_3")<10) hframe->SetXTitle("Fisher discriminant"); // fisher
  if (histlabel.find("hVD_3")<10) hframe->SetXTitle("D_{jet}"); // fisher

  if (histlabel.find("hMELA_8")<10) hframe->SetXTitle("D_{bkg}^{kin}"); // MELA
  if (histlabel.find("hPFMET_8")<10 || histlabel.find("hPFMET_3")<10 ) hframe->SetXTitle("PF MET (GeV)"); // PFMET

 


  hframe->GetXaxis()->SetLabelOffset(0.007);
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

//  THStack *htotup=new THStack("up bound","");
//  THStack *htotdow=new THStack("down bound","");

  // data

  for (unsigned int datasetIdData=0; datasetIdData<Vdatasetnamedata.size(); datasetIdData++){
    
    char dataset[328];
    sprintf(dataset,"%s",Vdatasetnamedata.at(datasetIdData).c_str());
    cout << "Root-ple= " << dataset << endl;
    // cout << "Counter=" << datasetIdData << " Root-ple=" << dataset << " Label=" << Vlabelbkg.at(datasetIdData) <<endl; 
    
    TFile *f1 = TFile::Open(dataset);
    hfourlepbestmass_4l_afterSel_new = (TH1F*)f1->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
    TH1 *hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str()/*"hfourlepbestmass_4l_afterSel_new_new"*/);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerColor(1);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(20);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerSize(0.95);
    // hfourlepbestmass_4l_afterSel_new_new->Draw("EPsame");

    //leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabeldata.at(datasetIdData).c_str(), "P");  
    //this part for efficiency of tight lepton calculation(qier)

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

  for (int nbins=1;nbins<=htotaldata->GetNbinsX(); nbins++){
    // cout << "BinCenter=" << htotaldata->GetBinCenter(nbins) << " BinContent=" << htotaldata->GetBinContent(nbins) << " BinErrorContent=" << htotaldata->GetBinError(nbins) << endl;
    x[nbins-1]=htotaldata->GetBinCenter(nbins);
    y[nbins-1]=htotaldata->GetBinContent(nbins);
    exl[nbins-1]=0.;
    exh[nbins-1]=0.;
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
  }

  cout << "Total data= " << totaldataentries << endl;

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(Nbins,x,y,exl,exh,eyl,eyh);
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.95);
  //


  // Background

  TH1F *hfourlepbestmass_4l_afterSel_new_qcdDEM  = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdDEM", "hfourlepbestmass_4l_afterSel_new_qcdDEM", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcdMu   = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdMu", "hfourlepbestmass_4l_afterSel_new_qcdMu", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcdBC   = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdBC", "hfourlepbestmass_4l_afterSel_new_qcdBC", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcd     = new TH1F("hfourlepbestmass_4l_afterSel_new_qcd", "hfourlepbestmass_4l_afterSel_new_qcd", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_singlet = new TH1F("hfourlepbestmass_4l_afterSel_new_singlet", "hfourlepbestmass_4l_afterSel_new_singlet", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_ST = new TH1F("hfourlepbestmass_4l_afterSel_new_ST", "hfourlepbestmass_4l_afterSel_new_ST", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_DY      = new TH1F("hfourlepbestmass_4l_afterSel_new_DY", "hfourlepbestmass_4l_afterSel_new_DY",Nbins, Xmin, Xmax);    
  TH1F *hfourlepbestmass_4l_afterSel_new_DYlight = new TH1F("hfourlepbestmass_4l_afterSel_new_DYlight", "hfourlepbestmass_4l_afterSel_new_DYlight",Nbins, Xmin, Xmax);
        
  TH1F *hfourlepbestmass_4l_afterSel_new_DYbb    = new TH1F("hfourlepbestmass_4l_afterSel_new_DYbb", "hfourlepbestmass_4l_afterSel_new_DYbb",Nbins, Xmin, Xmax);        
  TH1F *hfourlepbestmass_4l_afterSel_new_DYcc    = new TH1F("hfourlepbestmass_4l_afterSel_new_DYcc", "hfourlepbestmass_4l_afterSel_new_DYcc",Nbins, Xmin, Xmax);        


  TH1F *hfourlepbestmass_4l_afterSel_new_WW    = new TH1F("hfourlepbestmass_4l_afterSel_new_WW", "hfourlepbestmass_4l_afterSel_new_WW",Nbins, Xmin, Xmax); 
  TH1F *hfourlepbestmass_4l_afterSel_new_WZ    = new TH1F("hfourlepbestmass_4l_afterSel_new_WZ", "hfourlepbestmass_4l_afterSel_new_WZ",Nbins, Xmin, Xmax);   
                     
  TH1F *hfourlepbestmass_4l_afterSel_new_TT    = new TH1F("hfourlepbestmass_4l_afterSel_new_TT", "hfourlepbestmass_4l_afterSel_new_TT",Nbins, Xmin, Xmax); 
  TH1F *hfourlepbestmass_4l_afterSel_new_Wj    = new TH1F("hfourlepbestmass_4l_afterSel_new_Wj", "hfourlepbestmass_4l_afterSel_new_Wj",Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_ss    = new TH1F("hfourlepbestmass_4l_afterSel_ss", "hfourlepbestmass_4l_afterSel_ss",Nbins, Xmin, Xmax);

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

//    TH1 *hfourlepbestmass_mbb_up_new, hfourlepbestmass_mbb_dow_new;
  TString label=histlabel.c_str();
  
    //this part for tight lepton efficiency calculation (qier)
    
    if(datasetnamebkg.find("WZ")   < 200 ||
       datasetnamebkg.find("DYJetsToLL") < 200 ||
       datasetnamebkg.find("TT_TuneCUETP8M1") < 200 ||
       datasetnamebkg.find("WWTo2L2Nu")< 200 ||
       datasetnamebkg.find("DYBB") < 200 ||
       datasetnamebkg.find("WJetsToLNu") < 200 ||
       datasetnamebkg.find("ST_top_t-channel") < 200 ||
       datasetnamebkg.find("ST_antitop_t-channel") < 200 ||
       datasetnamebkg.find("datadriven") < 200
       ){
      
      hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
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

      hfourlepbestmass_4l_afterSel_new_new->Sumw2();
      // Higgs as background
      // DYJetsToLL check normalization
      if(datasetnamebkg.find("DYJetsToLL") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
        hfourlepbestmass_4l_afterSel_new_new->Scale(double(6025.2*35812*nnn/(35938649*1.559*10000)));
	hfourlepbestmass_4l_afterSel_new_DY->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DY->SetMarkerColor(kAzure+2);
	hfourlepbestmass_4l_afterSel_new_DY->SetFillColor(kAzure+2);                                                        
        cout << "fill DY histograms" << endl;        

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
          cout << "DY= " << hfourlepbestmass_4l_afterSel_new_DY->Integral(0,-1) << endl;
	  if (useDYJets==true) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 

      }
      
      // DYbb
      if(datasetnamebkg.find("DYBB") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
        hfourlepbestmass_4l_afterSel_new_new->Scale(double(40.4*35812.*nnn/2554303));
	hfourlepbestmass_4l_afterSel_new_DYbb->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DYbb->SetMarkerColor(kAzure+2);
	//	hfourlepbestmass_4l_afterSel_new_DYbb->SetLineColor(kAzure+2);
	//	hfourlepbestmass_4l_afterSel_new_DYbb->SetLineWidth(2);
	hfourlepbestmass_4l_afterSel_new_DYbb->SetFillColor(kAzure+2);                                                        
	
	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if((datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200)) || datasetnamebkg.find("DYBB")<200) {
          cout << "DYbb= " << hfourlepbestmass_4l_afterSel_new_DYbb->Integral(0,-1) << endl;
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
        }
      }

      if (useDYJetsFromData==false){
	// WW
	if(datasetnamebkg.find("WWTo2L2Nu") < 200){
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(12.178*2573.*nnn/1999000.));
	  hfourlepbestmass_4l_afterSel_new_WW->Add(hfourlepbestmass_4l_afterSel_new_new); 
	  hfourlepbestmass_4l_afterSel_new_WW->SetMarkerColor(kCyan+3); 
	  hfourlepbestmass_4l_afterSel_new_WW->SetFillColor(kCyan+3);
	  
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  if(datasetnamebkg.find(temp) < 200 &&  hfourlepbestmass_4l_afterSel_new_WW->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");  
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");    
	}     
	
	// WZ     
	if(datasetnamebkg.find("WZ") < 200){ 
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(2*5.595*35812.*nnn/(20280590.*0.6643)));
	  hfourlepbestmass_4l_afterSel_new_WZ->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_WZ->SetMarkerColor(kCyan-2);  
	  hfourlepbestmass_4l_afterSel_new_WZ->SetFillColor(kCyan-2);  
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  if(datasetnamebkg.find(temp) < 200 && hfourlepbestmass_4l_afterSel_new_WZ->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP"); 
	}    
	
	// TTT
	if(datasetnamebkg.find("TT_TuneCUETP8M1") < 200){
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(87.31*35812*nnn/(5669539.)));
	  hfourlepbestmass_4l_afterSel_new_TT->Add(hfourlepbestmass_4l_afterSel_new_new);     
	  hfourlepbestmass_4l_afterSel_new_TT->SetMarkerColor(kTeal-6); 
	  hfourlepbestmass_4l_afterSel_new_TT->SetFillColor(kTeal-6);                     
	  
	  cout << "TT+jets= " << hfourlepbestmass_4l_afterSel_new_TT->GetEntries() << endl;
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
           legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
       
       }      
	
	
	// W+jets		      					
	if(datasetnamebkg.find("_WJetsToLNu") < 200){
	  cout << "Wjets" << endl;
          hfourlepbestmass_4l_afterSel_new_new->Scale(61526.7*35812/(234227690*1.53*100000));
	  hfourlepbestmass_4l_afterSel_new_Wj->Add(hfourlepbestmass_4l_afterSel_new_new);        
	  hfourlepbestmass_4l_afterSel_new_Wj->SetMarkerColor(kGray+2);	
	  hfourlepbestmass_4l_afterSel_new_Wj->SetFillColor(kGray+2);	
	  
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");	   

	}   

      //Single Top
        if(datasetnamebkg.find("ST_top_t-channel") < 200){
          cout << "single top" << endl;
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(136.1*35812./5887322.));
          hfourlepbestmass_4l_afterSel_new_ST->Add(hfourlepbestmass_4l_afterSel_new_new);
          hfourlepbestmass_4l_afterSel_new_ST->SetMarkerColor(kViolet);
          hfourlepbestmass_4l_afterSel_new_ST->SetFillColor(kViolet);

          char temp[328];
          sprintf(temp,"%s",histosdir.c_str());
          if(datasetnamebkg.find(temp) < 200 && hfourlepbestmass_4l_afterSel_new_ST->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
          //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");                
        }

      //Single antitop
        if(datasetnamebkg.find("ST_antitop_t-channel") < 200){
          cout << "single antitop" << endl;
          hfourlepbestmass_4l_afterSel_new_new->Scale(double(81.0*35812./3775674.));
          hfourlepbestmass_4l_afterSel_new_ST->Add(hfourlepbestmass_4l_afterSel_new_new);
          hfourlepbestmass_4l_afterSel_new_ST->SetMarkerColor(kViolet);
          hfourlepbestmass_4l_afterSel_new_ST->SetFillColor(kViolet);
        }


      } 

      if(datasetnamebkg.find("datadriven") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
        hfourlepbestmass_4l_afterSel_new_new->Scale(1.8);
        hfourlepbestmass_4l_afterSel_new_ss->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_ss->SetMarkerColor(kGray+2);
        hfourlepbestmass_4l_afterSel_new_ss->SetFillColor(kGray+2);

        char temp[328];
        sprintf(temp,"%s",histosdir.c_str());
        legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");

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
    
       
       char temppp[328];
       sprintf(temppp,"%s",histosdir.c_str());

       // other backgrounds
       if (useDYJets==true){ 
	 if(	 datasetnamebkg.find("output_DYJetsToLL")<200
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
	  datasetnamebkg.find("output_WW") < 200 || 
	  (datasetnamebkg.find("output_WW") < 200 && datasetnamebkg.find(whichenergy.c_str())<200)
	  )
	  )  { 
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_WW); 
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_WW); 
       }
       if(datasetnamebkg.find(temppp) < 200 && (
					       datasetnamebkg.find("output_WZ") < 200 ||
					       (datasetnamebkg.find("output_WZ") < 200 && datasetnamebkg.find(whichenergy.c_str())<200)
					       ) 
	  )  {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_WZ);         
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_WZ);   
       }
       if(   datasetnamebkg.find("output_TT") < 200  
	  )  {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_TT);                        
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_TT);   
       }
       if(datasetnamebkg.find("output_WJetsToLNu") < 200
	  )  {
         cout << "add wj" << endl;
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_Wj);   
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_Wj);
       }
       if(datasetnamebkg.find(temppp) < 200 && (
                                               datasetnamebkg.find("output_DYBB") < 200 ||  
                                               (datasetnamebkg.find("output_DYBB") < 200 && datasetnamebkg.find(whichenergy.c_str())<200)
                                               )
          )  {
         htotal->Add(hfourlepbestmass_4l_afterSel_new_DYbb);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYbb);
       }
       
       if(  datasetnamebkg.find("output_ST")<200)  {
         htotal->Add(hfourlepbestmass_4l_afterSel_new_ST);                        
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_ST);
       }
       if(  datasetnamebkg.find("datadriven")<200)  {
         htotal->Add(hfourlepbestmass_4l_afterSel_new_ss);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_ss);
       }
 
  }

   TH1F *denomup, *denomdow;
   TH1F *denomup2, *denomdow2;
   TH1F *denom1, *denom2;

   TH1F *denom_up, *denom_dow;
   TH1F *denom_1;

   TH1F *denom_mc_up, *denom_mc_dow;
   TH1F *denom_mc_1;

   TString file_un="/eos/uscms/store/user/wangz/ntuple/unc/2e2b/un_all.root";
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
    }

   //error bar set
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
      ratiostaterr->SetBinError(nbins,1.);
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
         reyl[nbins-1] = 0.;
         reyh[nbins-1] =0.;
         reyl2[nbins-1] = 0.;
         reyh2[nbins-1] = 0.;
         eyl4[nbins-1]=0.;
         eyh4[nbins-1]=0.;
    }

  }

    htotal->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    TH1F *temp = NULL;
    temp = (TH1F*)htotal->GetStack()->Last();
    TH1F *staterr = (TH1F*)temp->Clone("htotal");
    staterr->SetFillColor(kGray+3);
    staterr->SetLineColor(kGray+3);
    staterr->SetLineWidth(0);
    staterr->SetMarkerSize(0);
    staterr->SetFillStyle(3013);
//    staterr->Draw("e2 same");

    TGraphAsymmErrors *rgr = new TGraphAsymmErrors(Nbins,rx,ry,exl3,exh3,reyl,reyh);
    rgr->SetMarkerColor(1);
    rgr->SetMarkerStyle(20);
    rgr->SetMarkerSize(0.95);
 
//    cout << "bin1 l=" << reyl[0] <<" bin2 h=" << reyh[0] << endl;

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
//  htotaldata->Draw("EPsame");
  if(un_mc) gr4->Draw("E2 same");
  gr->Draw("EPsame");
  gr2->Draw("e2 same");
  
  if(un_mc) legend->AddEntry(rgr_mc,"syst. unc.","F"); 
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
  hframe2->GetXaxis()->SetLabelSize(0.020);
  //  hframe2->GetYaxis()->SetTitleSize(0.047);
  hframe2->SetYTitle("Data/MC");
  //  hframe2->GetYaxis()->SetRangeUser(-10,10);
  hframe2->GetYaxis()->SetNdivisions(503);
  //hframe2->GetXaxis()->SetTitleOffset(1.25);
  hframe2->Draw("");

  htotalHistoRatio->SetMarkerStyle(20);
  htotalHistoRatio->SetMarkerSize(0.95);
  htotalHistoRatio->SetMarkerColor(kBlack);
//  htotalHistoRatio->Draw("Psame");
//  ratiostaterr->Draw("e2 same");
  rgr->Draw("P0 E0 Z same");
  if(un_mc) rgr_mc->Draw("E2 same");
  gr3->Draw("E2 same");

  c1->Update();
  
  if (histosdir.find("4mu")<25) whichchannel="4mu";
  if (histosdir.find("2e2mu")<25) whichchannel="2e2mu";

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
  char htotalfinal[300];
  sprintf(htotalfinal,"plots/htotalfinal_%s_%s.root",histlabel.c_str(),whichchannel.c_str());
  TFile *file1 = new TFile(htotalfinal, "RECREATE");
  file1->cd();
  htotalHisto->Write();
  htotaldata->Write();
  file1->Write();
  file1->Close();

}

void PlotStackZ::setSamplesNamesZ()
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

  whichchannel="2e2b_BH";
  histosdir="llbb";

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
            inputfilename.find("DYBB")<100 ||
	    inputfilename.find("WWTo2L2Nu")<100 ||
            inputfilename.find("_ST") < 100 ||
            inputfilename.find("datadriven") < 100
           ){ 
      //Vdatasetnamebkg.push_back(inputfilename);
      

      if(inputfilename.find("DYJetsToLL")<200){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("DY");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kAzure+2);
      }

      if(inputfilename.find("DYBB")<200){
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("DYBB");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kAzure+2);
      }

      if(inputfilename.find("TT_TuneCUETP8M1")<200 ){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("TT");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-6);
	//cout << "VBF H" << endl;
      }
      
      if( (inputfilename.find("WJetsToLNu")<200)){
	  Vdatasetnamebkg.push_back(inputfilename);
	  Vlabelbkg.push_back("W");
	  Vxsectionbkg.push_back(1.); //pb      
          Vcolorbkg.push_back(kGray+2);
      }


      if(inputfilename.find("WWTo2L2Nu")<200){
	Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("WW");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kCyan+3);
        //cout << "VH" << endl;
      }

      if(inputfilename.find("_ST")<200){
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("ST");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kViolet);
        //cout << "VH" << endl;
      }

      if(inputfilename.find("datadriven")<200){
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("same sign");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kGray+2);
      }

 
    }
    // BACKGROUND from other sources
    else{  // provided that everything that is neither data nor signal is background
      //Vdatasetnamebkg.push_back(inputfilename);
      
      // WZ
      if(inputfilename.find("WZ")<200){
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


void PlotStackZ::printnumbers(char name[400], TH1F *h){

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
void PlotStack4l::createdatacards(float Higgsm, float channel, float energy, float masslow, float masshigh, float ggH, float qqH, float WH, float ZH, float ttH, float bkg_qqzz, float bkg_ggzz, float bkg_zjets){

  
  Char_t txtOUT[500];
  //sprintf(txtOUT,"datacards/%s/hzz4l_%s_txt.txt",Higgsm,datasetName.Data());
  cout << "Opening a datacard file " << txtOUT << endl;
  ofstream output_txt;
  output_txt.open(txtOUT);

  output_txt << "imax 1" << endl; 
  output_txt << "jmax 7" << endl;
  output_txt << "kmax *" << endl;
  output_txt <<"------------" << endl; 
  output_txt <<"shapes * * hzz4l_4muS_8TeV.input.root w:$PROCESS" << endl;
  output_txt <<"------------" << endl;
  output_txt <<"bin a1" << endl; 
  output_txt <<"observation 16" << endl; 
  output_txt <<"------------" << endl;
  output_txt <<"## mass window [105.0,140.0]" << endl; 
  output_txt <<"bin a1 a1 a1 a1 a1 a1 a1 a1" << endl; 
  output_txt <<"process ggH qqH WH ZH ttH bkg2d_qqzz bkg2d_ggzz bkg2d_zjets" << endl; 
  output_txt <<"process -4 -3 -2 -1 0 1 2 3" << endl; 
  output_txt <<"rate 1.0000 1.0000 1.0000 1.0000 1.0000 7.6204 0.1543 1.1796" << endl; 
  output_txt <<"------------" << endl;
  output_txt <<"lumi_8TeV lnN 1.026 1.026 1.026 1.026 1.026 1.026 1.026 -" << endl; 
  output_txt <<"pdf_gg lnN 1.0720 - - - 1.0780 - 1.0708 -" << endl;
  output_txt <<"pdf_qqbar lnN - 1.0270 1.0350 1.0350 - 1.0341 - -" << endl;
  output_txt <<"pdf_hzz4l_accept lnN 1.02 1.02 1.02 1.02 1.02 - - -" << endl; 
  output_txt <<"QCDscale_ggH lnN 1.0750 - - - - - - -" << endl; 
  output_txt <<"QCDscale_qqH lnN - 1.0020 - - - - - -" << endl; 
  output_txt <<"QCDscale_VH lnN - - 1.0040 1.0155 - - - -" << endl; 
  output_txt <<"QCDscale_ttH lnN - - - - 1.0655 - - -" << endl;   
  output_txt <<"QCDscale_ggVV lnN - - - - - - 1.2431 -" << endl;  
  output_txt <<"QCDscale_VV lnN - - - - - 1.0284 - -" << endl;  
  output_txt <<"BRhiggs_hzz4l lnN 1.02 1.02 1.02 1.02 1.02 - - -" << endl;  
  output_txt <<"CMS_eff_m lnN 1.043 1.043 1.043 1.043 1.043 1.043 1.043 -" << endl;  
  output_txt <<"CMS_hzz4mu_Zjets lnN - - - - - - - 0.6/1.4" << endl;  
  output_txt <<"CMS_zz4l_bkgMELA param 0  1  [-3,3]" << endl; 
  output_txt <<"CMS_zz4l_mean_m_sig param 0.0 1.0" << endl;  
  output_txt <<"## CMS_zz4l_mean_m_sig = 0.001" << endl;  
  output_txt <<"CMS_zz4l_sigma_m_sig param 0.0 0.2" << endl;  
  output_txt <<"CMS_zz4l_n_sig_1_8 param 0.0 0.01" << endl;  
  output_txt <<"interf_ggH param 0 1 [-1,1]" << endl;  
  
  output_txt.close();

}
*/

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
