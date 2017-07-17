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
    
  inputfile="filelist_input.txt";

  setSamplesNamesZ(); 
  cout << "\t Analysing samples for " << whichchannel << " analysis" << endl; 


  //WARNING: depending on histolabel, modify the declaration and the settings of hframe below
  //also choose a sensible value for nRebin

  //std::string histolabel = "hPUvertices";    // numPU
  //std::string histolabel = "hPUvertices_ReWeighted";    // numPY reweighted

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

  std::string histolabel = "hMZ1_5";    // Z1 mass   
  //std::string histolabel = "hMZ2_5";    // Z2 mass 
  //std::string histolabel = "hM4l_5";    // 4l mass 

  //std::string histolabel = "hIso_5";    // worst isolation 
  //std::string histolabel = "hSip_5";  // worst sip 
  //std::string histolabel = "hIp_5";   // worst IP

//   std::string histolabel = "Mbb_6";
//  std::string histolabel = "bdiscr_5_lead";
//  std::string histolabel = "bdiscr_5_sub";
//    std::string histolabel = "Mjj_6";

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

//   std::string histolabel = "hEtaLep1_7";
//   std::string histolabel = "hEtaLep2_7";
//   std::string histolabel = "hIsoLep1_7";
//   std::string histolabel = "hIsoLep2_7";
 //   std::string histolabel = "hYZ1_5";
  //std::string histolabel = "hPtZ2_8"; // Z2 pt after full selection
  //std::string histolabel = "hYZ1_8"; // Z1 rapidity after full selection
  //std::string histolabel = "hYZ2_8"; // Z2 rapidity after full selection
  //std::string histolabel = "hIso_8"; // worst isolated lepton: isolation value after full selection
  //std::string histolabel = "hSip_8"; // worst sip lepton: sip value after full selection
  //std::string histolabel = "hMELA_8"; // MELA discriminant after full selection 
 //  std::string histolabel = "hPFMET_8"; // PFMET
  //std::string histolabel = "hM4l_T_8"; // Transverse mass
  //std::string histolabel = "DPHI_8"; // DeltaPhi - 4l + MET
//   std::string histolabel = "hPtJet_7"; 
//   std::string histolabel = "hPtJet_8";
//   std::string histolabel = "hEtaJet_7";
//   std::string histolabel = "hEtaJet_8";
 // std::string histolabel = "hNjets_8";
  //std::string histolabel = "hDjj_8"; // delta eta between jets for VBF analysis
 // std::string histolabel = "hMjj_8"; // dimass between jets for VBF analysis
//  std::string histolabel = "hN_loose_e";
//  std::string histolabel = "hN_loose_mu"; 
//   std::string histolabel = "hN_good_ele";
//   std::string histolabel = "hN_good_mu";

  useLogY = false;
  useLogX = false;

  useDYJets=true;
  useDYJetsFromData=false;
  
  nRebin=2;
  std::cout << "Histogram label is= " << histolabel << std::endl;
  
  // Final yields
  system("mkdir plots");

  Char_t yieldsOUT[500];
  sprintf(yieldsOUT,"plots/amcatnlo/yields_%s_%s.txt",whichchannel.c_str(),whichenergy.c_str());
  if (histolabel.find("hM4l_9")<10 ) {
    cout << "Opening a file for final numbers= " << yieldsOUT << endl;
    outputyields.open(yieldsOUT);
  }
  
  // Execute the analysis
  plotmZ(histolabel);
  
  // close file for final yields
  if (histolabel.find("hM4l_9")<10 ) outputyields.close();

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
  if (whichenergy.find("RunI")<50) {
    text = "#sqrt{s} = 7 TeV, L = 5.05 fb^{-1} ; #sqrt{s} = 8 TeV, L = 19.71 fb^{-1}" ;
    ll->AddText(0.3, 0.6, text);
  }
  else if (whichenergy.find("7TeV")<50) {
    text = "#sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" ;
    ll->AddText(0.65, 0.6, text);
  }
  else if (whichenergy.find("8TeV")<50) {
    text = "#sqrt{s} = 8 TeV, L = 19.71 fb^{-1}" ;
    ll->AddText(0.65, 0.6, text);
  }
  else if (whichenergy.find("13TeV")<50) {
    text = "#sqrt{s} = 13 TeV, L = 36.46 fb^{-1}" ;
    //text = "#sqrt{s} = 13 TeV, L = 14.77 fb^{-1}" ;
    ll->AddText(0.65, 0.6, text);
  }
  //ll->Draw();



  TLegend *leg0 = new TLegend(0.6,0.40,0.8,0.90,NULL,"brNDC");
  leg0->SetTextSize(0.020);
  leg0->SetLineColor(0);
  leg0->SetLineWidth(1);
  leg0->SetFillColor(kWhite);
  leg0->SetBorderSize(0);
 
  TLegend* legend = new TLegend( 0.70, 0.76, 0.86, 0.86);
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
    hframe= new TH2F("hframe","hframe",60,60.,120.,500,1.0,200000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",60, 60., 120., 1000, 0.7, 1.3);// mZ1 
  }

  if (histlabel.find("hPtZ1_5")<10){
    hframe= new TH2F("hframe","hframe",100,0.,200.,500,1.0,60000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 200., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtLep1_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,200000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",200, 0., 200., 500, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtLep2_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,200000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",200, 0., 200., 500, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hPtLep1_8")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hPtLep2_8")<10){
    hframe= new TH2F("hframe","hframe",200,0.,200.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 0., 200., 1000, 0.5, 2.);// mZ1 
  }

  if (histlabel.find("hYZ1_5")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,30000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hEtaLep1_7")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,50000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hEtaLep2_7")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,50000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hIsoLep1_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,1.,500,1.0,500000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 1., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hIsoLep2_7")<10){
    hframe= new TH2F("hframe","hframe",100,0.,1.,500,1.0,500000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, 0., 1., 1000, 0.8, 1.2);// mZ1 
  }

  if (histlabel.find("hMZ1_7")<10){
    hframe= new TH2F("hframe","hframe",80,40.,200.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",6000, 40., 160., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("Mbb_6")<10){
    hframe= new TH2F("hframe","hframe",80,20.,400.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",80, 20., 400., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("Mjj_6")<10){
    hframe= new TH2F("hframe","hframe",80,20.,400.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",80, 20., 400., 1000, 0.5, 1.5);// mZ1 
  }


  if (histlabel.find("bdiscr_5_lead")<10){
    hframe= new TH2F("hframe","hframe",20,0.,1.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",20, 0., 1., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("bdiscr_5_sub")<10){
    hframe= new TH2F("hframe","hframe",20,0.,1.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",20, 0., 1., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hMZ2_5")<10 && whichchannel.find("4#mu")<20){
    hframe= new TH2F("hframe","hframe",80,40.,200.,500,0.0000000001,1000000.);// mZ2 
    hframe2= new TH2F("hframe2","hframe2",6000, 40., 200., 1000, 0.5, 2.);// mZ2 
  }

  if (histlabel.find("hPtJet_7")<10){
    hframe= new TH2F("hframe","hframe",200,0.,400.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",1000, 0., 400., 1000, 0.5, 1.5);// mZ1 
  }
 
  if (histlabel.find("hPtJet_8")<10){
    hframe= new TH2F("hframe","hframe",200,0.,400.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",1000, 0., 400., 1000, 0.5, 1.5);// mZ1 
  }
 
  if (histlabel.find("hEtaJet_7")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hEtaJet_8")<10){
    hframe= new TH2F("hframe","hframe",100,-5.,5.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",100, -5., 5., 1000, 0.5, 1.5);// mZ1 
  }

  if (histlabel.find("hNjets_8")<10){
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.5, 1.5);// mZ1 
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
    hframe= new TH2F("hframe","hframe",10,0.,10.,500,1.0,10000000.);// mZ1 
    hframe2= new TH2F("hframe2","hframe2",50, 0., 10., 1000, 0.5, 1.5);// mZ1 
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
    hframe= new TH2F("hframe","hframe",300, 0., 300., 1000, 0.000001, 100000.);// PFMET
    hframe2= new TH2F("hframe2","hframe2",300, 0.,300., 1000, 0.5, 1.5);// PFMET
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
  //hframe->SetXTitle("M_{Z#rightarrow ee} (BB)  [GeV]");
  //hframe->SetXTitle("M_{Z#rightarrow ee} (EE)  [GeV]");
  // hframe->SetXTitle("Sign. 3DIP");  
  // hframe->SetXTitle("R^{iso}_{12} [GeV/c^{2}]");
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


    //cout << "Nbins=" << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << endl;
    //cout << "htotaldata nBins = " << htotaldata->GetNbinsX() << ", hfourlepbestmass_4l_afterSel_new_new nBins = " << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << endl;
    //cout << "htotaldata lowestX = " << htotaldata->GetXaxis()->GetXmin() <<  ", htotaldata highestX = " << htotaldata->GetXaxis()->GetXmax() << ", hfourlepbestmass_4l_afterSel_new_new lowestX = " << hfourlepbestmass_4l_afterSel_new_new->GetXaxis()->GetXmin() << ", hfourlepbestmass_4l_afterSel_new_new highestX = " << hfourlepbestmass_4l_afterSel_new_new->GetXaxis()->GetXmax() << endl;

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
  if (histlabel.find("hM4l_9")<10 ) outputyields << "Data " << totaldataentries << " +/- 0" << endl;

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(Nbins,x,y,exl,exh,eyl,eyh);
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.95);
  //
 
  // Z+X from data
  if (useDYJets==false && useDYJetsFromData==true){
    if (histlabel.find("hM4l_7")<10 || histlabel.find("hM4l_8")<10){
      for (unsigned int datasetIdbkgData=0; datasetIdbkgData<Vdatasetnamebkgdata.size(); datasetIdbkgData++){
	char dataset[328];
	sprintf(dataset,"%s",Vdatasetnamebkgdata.at(datasetIdbkgData).c_str());
	cout << "Root-ple= " << dataset << endl;
	TFile *f3 = TFile::Open(dataset);
	TH1F *hfourlepbestmass_4l_afterSel_orig = new TH1F("hfourlepbestmass_4l_afterSel_orig", "Mass of four leptons after fullselection", 2460, 0.,1230. );
        hfourlepbestmass_4l_afterSel_orig = (TH1F*)f3->Get("h_3P1F_2P2F");

        hfourlepbestmass_4l_afterSel_new = new TH1F("hfourlepbestmass_4l_afterSel_new", "Mass of four leptons after fullselection", 2400, 4.5,1204.5 );

        int mbins=1;
        for (int nbins=1;nbins<=hfourlepbestmass_4l_afterSel_orig->GetNbinsX(); nbins++){
          if (hfourlepbestmass_4l_afterSel_orig->GetBinCenter(nbins)>4.5 && hfourlepbestmass_4l_afterSel_orig->GetBinCenter(nbins)<1204.5){
            hfourlepbestmass_4l_afterSel_new->SetBinContent(mbins,double(hfourlepbestmass_4l_afterSel_orig->GetBinContent(nbins)));
            mbins++;
          }
        }

	TH1 *hfourlepbestmass_4l_afterSel_new_new;
	
	nRebinZ_X=nRebin*2;
	hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebinZ_X, "h_3P1F_2P2F");      
	hfourlepbestmass_4l_afterSel_new_new->SetLineColor(kCyan-2);
	hfourlepbestmass_4l_afterSel_new_new->SetFillColor(kCyan-2);
	hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(24);
	hfourlepbestmass_4l_afterSel_new_new->SetLineWidth(1);
	
	char temp[328];
        sprintf(temp,"%s",histosdir.c_str());

	if ( 
	    histlabel.find("hM4l_7")<10 && Vdatasetnamebkgdata.at(datasetIdbkgData).find("m4l_gt_70")<85 ) {
	  cout << "Adding Z+X for m4l > 70. GeV" << endl;
	  cout << "N bins Z+X= " << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX()<< endl;
	  cout << "Z+X entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << endl;
	  htotal->Add(hfourlepbestmass_4l_afterSel_new_new);
	  htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_new);
	  if ( 
	      Vdatasetnamebkgdata.at(datasetIdbkgData).find(temp) <200 && 
	      (Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichenergy) < 200 || Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichsample) < 200) &&	      
	      Vdatasetnamebkgdata.at(datasetIdbkgData).find("2mu2e")>85 ) {
	    cout << "Adding legend for Z+X" << endl;
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkgdata.at(datasetIdbkgData).c_str(), "F"); }
	}
	else if ( 
		 Vdatasetnamebkgdata.at(datasetIdbkgData).find(temp) <200 && 
		 (Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichenergy) < 200 || Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichsample) < 200) &&
		 histlabel.find("hM4l_8")<10 && !(Vdatasetnamebkgdata.at(datasetIdbkgData).find("m4l_gt_70")<85) ) {
	  cout << "Adding Z+X for m4l > 100. GeV" << endl;
	  cout << "Z+X entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << endl;
	  htotal->Add(hfourlepbestmass_4l_afterSel_new_new);
	  htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_new);	  
	  if (
	      Vdatasetnamebkgdata.at(datasetIdbkgData).find(temp) <200 && 
	      (Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichenergy) < 200 || Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichsample) < 200) &&
	      Vdatasetnamebkgdata.at(datasetIdbkgData).find("2mu2e")>85 ){
	    cout << "Adding legend for Z+X" << endl;
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkgdata.at(datasetIdbkgData).c_str(), "F");}
	}
      }
    }
  }

  //cout << "Total Z+X is= " << htotal->GetHistogram()->GetEntries() << endl;
  cout << "Total Z+X is= " << htotalHisto->Integral(0,-1)  << endl;
  //outputyields << "Z+X "   << htotal->GetHistogram()->GetEntries() << " +/- 0"<< endl;
  if (histlabel.find("hM4l_9")<10 ) outputyields << "Z+X " << htotalHisto->Integral(0,-1) << " +/- 0" << endl;


  // Background

  TH1F *hfourlepbestmass_4l_afterSel_new_qcdDEM  = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdDEM", "hfourlepbestmass_4l_afterSel_new_qcdDEM", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcdMu   = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdMu", "hfourlepbestmass_4l_afterSel_new_qcdMu", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcdBC   = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdBC", "hfourlepbestmass_4l_afterSel_new_qcdBC", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcd     = new TH1F("hfourlepbestmass_4l_afterSel_new_qcd", "hfourlepbestmass_4l_afterSel_new_qcd", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_singlet = new TH1F("hfourlepbestmass_4l_afterSel_new_singlet", "hfourlepbestmass_4l_afterSel_new_singlet", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_DY      = new TH1F("hfourlepbestmass_4l_afterSel_new_DY", "hfourlepbestmass_4l_afterSel_new_DY",Nbins, Xmin, Xmax);    
  TH1F *hfourlepbestmass_4l_afterSel_new_DYlight = new TH1F("hfourlepbestmass_4l_afterSel_new_DYlight", "hfourlepbestmass_4l_afterSel_new_DYlight",Nbins, Xmin, Xmax);
        
  TH1F *hfourlepbestmass_4l_afterSel_new_DYbb    = new TH1F("hfourlepbestmass_4l_afterSel_new_DYbb", "hfourlepbestmass_4l_afterSel_new_DYbb",Nbins, Xmin, Xmax);        
  TH1F *hfourlepbestmass_4l_afterSel_new_DYcc    = new TH1F("hfourlepbestmass_4l_afterSel_new_DYcc", "hfourlepbestmass_4l_afterSel_new_DYcc",Nbins, Xmin, Xmax);        


  TH1F *hfourlepbestmass_4l_afterSel_new_WW    = new TH1F("hfourlepbestmass_4l_afterSel_new_WW", "hfourlepbestmass_4l_afterSel_new_WW",Nbins, Xmin, Xmax); 
  TH1F *hfourlepbestmass_4l_afterSel_new_WZ    = new TH1F("hfourlepbestmass_4l_afterSel_new_WZ", "hfourlepbestmass_4l_afterSel_new_WZ",Nbins, Xmin, Xmax);   
                     
  TH1F *hfourlepbestmass_4l_afterSel_new_TT    = new TH1F("hfourlepbestmass_4l_afterSel_new_TT", "hfourlepbestmass_4l_afterSel_new_TT",Nbins, Xmin, Xmax); 
  TH1F *hfourlepbestmass_4l_afterSel_new_Wj    = new TH1F("hfourlepbestmass_4l_afterSel_new_Wj", "hfourlepbestmass_4l_afterSel_new_Wj",Nbins, Xmin, Xmax);

/*
  TFile * PU= new TFile("pileup_MC_80x_271036-276811_69200.root");
  TH1F * HistoPUData= (TH1F *) PU->Get("pileup");
  TH1F * HistoPUMC= (TH1F *) PU->Get("pileup_mc");
  double nnn = HistoPUData->Integral()/HistoPUMC->Integral();
  cout << "nnn = " << nnn << endl;
*/

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
    
    if(datasetnamebkg.find("WZ")   < 200 ||
       datasetnamebkg.find("DYJetsToLL") < 200 ||
       datasetnamebkg.find("TT_TuneCUETP8M1") < 200 ||
       datasetnamebkg.find("WWTo2L2Nu")< 200 ||
       datasetnamebkg.find("WJetsToLNu") < 200  
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

      // Higgs as background
      // DYJetsToLL check normalization
      if(datasetnamebkg.find("DYJetsToLL") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries()>0 ){
         cout << "test1" << endl;
	 hfourlepbestmass_4l_afterSel_new_new->Scale(double(5765.4*1152.*nnn/(3416613*0.6643)));//*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/12138430./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
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
	if(datasetnamebkg.find("WZ") < 200){  
	  hfourlepbestmass_4l_afterSel_new_WZ->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_WZ->SetMarkerColor(kCyan-2);  
	  hfourlepbestmass_4l_afterSel_new_WZ->SetFillColor(kCyan-2);  
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  if(datasetnamebkg.find(temp) < 200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && hfourlepbestmass_4l_afterSel_new_WZ->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
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
      
      
    }
    //
    else if (useDYJetsFromData==false){      
      if(datasetnamebkg.find("_MuPt5Enriched") < 200 ){
	hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcdMu->Add(hfourlepbestmass_4l_afterSel_new_new);      
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001); 
	hfourlepbestmass_4l_afterSel_new_qcdMu->SetLineColor(1);
	hfourlepbestmass_4l_afterSel_new_qcdMu->SetFillColor(kTeal-8);
	hfourlepbestmass_4l_afterSel_new_qcdMu->SetLineWidth(1);
	
	char temp[328];
	sprintf(temp,"%s/output_QCD_Pt-15to20_MuPt5Enriched",histosdir.c_str());
	
	if(datasetnamebkg.find(temp) < 200 ){ // provided that this is the last single-top sample
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcdMu,"QCD MuPt5", "F");  
	  cout << "Label= QCD MuPt5    Entries= " << hfourlepbestmass_4l_afterSel_new_qcdMu->Integral(0,-1) <<endl;
	}  
	cout << "Label= " << Vlabelbkg.at(datasetId) << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) <<endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");  
      } 
      else if(datasetnamebkg.find("_doubleEMEnriched") < 200 ){
	hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcdDEM->Add(hfourlepbestmass_4l_afterSel_new_new);      
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001); 
	hfourlepbestmass_4l_afterSel_new_qcdDEM->SetLineColor(1);
	hfourlepbestmass_4l_afterSel_new_qcdDEM->SetFillColor(kTeal+8);
	hfourlepbestmass_4l_afterSel_new_qcdDEM->SetLineWidth(1);
	
	char temp[328];
	sprintf(temp,"%s/output_QCD_Pt-80_doubleEMEnriched",histosdir.c_str());
	
	if(datasetnamebkg.find(temp) < 80 ){ 
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcdDEM,"QCD doubleEM", "F");  
	  cout << "Label= QCD double EM    Entries= " << hfourlepbestmass_4l_afterSel_new_qcdMu->Integral(0,-1) <<endl;
	}  
	cout << "Label= " << Vlabelbkg.at(datasetId) << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) <<endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");  
      }  
      else if(datasetnamebkg.find("_BCtoE") < 200 ){
	hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcdBC->Add(hfourlepbestmass_4l_afterSel_new_new);      
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001); 
	hfourlepbestmass_4l_afterSel_new_qcdBC->SetLineColor(1);
	hfourlepbestmass_4l_afterSel_new_qcdBC->SetFillColor(kTeal-2);
	hfourlepbestmass_4l_afterSel_new_qcdBC->SetLineWidth(1);
	
	char temp[328];
	sprintf(temp,"%s/output_QCD_Pt-20to30_BCtoE",histosdir.c_str());
	
	if(datasetnamebkg.find(temp) < 200){ // provided that this is the last single-top sample
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcdBC,"QCD BCtoE", "F");  
	  cout << "Label= QCD BCtoE    Entries= " << hfourlepbestmass_4l_afterSel_new_qcdBC->Integral(0,-1) <<endl;
	}  
	cout << "Label= " << Vlabelbkg.at(datasetId) << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) <<endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");  
      }
       else if(datasetnamebkg.find("QCD_Pt") ){
	
	hfourlepbestmass_4l_afterSel_new_new=hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcd->Add(hfourlepbestmass_4l_afterSel_new_new);      
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001); 
	hfourlepbestmass_4l_afterSel_new_qcd->SetLineColor(kTeal-2);
	hfourlepbestmass_4l_afterSel_new_qcd->SetFillColor(kTeal-2);
	hfourlepbestmass_4l_afterSel_new_qcd->SetLineWidth(1);

	char temp[328];
	//sprintf(temp,"%s/output_QCD_Pt_10to15",histosdir.c_str());
        sprintf(temp,"%s/output_QCD_Pt_1000to1400",histosdir.c_str());
	
	cout << "alpha" << temp << datasetnamebkg.find(temp) << endl;
	//sprintf(temp,"%s",histosdir.c_str());
	//if(datasetnamebkg.find(temp) < 200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && hfourlepbestmass_4l_afterSel_new_Wj->GetEntries()>0. ) legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 

	if(datasetnamebkg.find(temp) < 200 && hfourlepbestmass_4l_afterSel_new_qcd->GetEntries()>0.){ // provided that this is the last single-top sample
	  leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_qcd,"QCD", "F");  
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcd,"QCD", "F");  
	  cout << "Label= QCD Entries= " << hfourlepbestmass_4l_afterSel_new_qcd->Integral(0,-1) <<endl;
	}  
	cout << "Label= " << Vlabelbkg.at(datasetId) << "  Entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) <<endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F"); 
      }       
      // single top
      //else if (datasetId>=12 && datasetId<=14){
      else if(datasetnamebkg.find("ST_") < 200 ||  datasetnamebkg.find("Tbar_") < 200 ){
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
    
    if(datasetnamebkg.find(tempp) < 500) {
    }
    else {
       
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
       if(datasetnamebkg.find(temppp) < 200 && (
					       datasetnamebkg.find("output_TT") < 200 || 
					       (datasetnamebkg.find("output_TT") < 200 && datasetnamebkg.find(whichenergy.c_str())<200)
					       )
	  )  {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_TT);                        
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_TT);   
       }
       if(datasetnamebkg.find(temppp) < 200 && (
					       datasetnamebkg.find("output_WJ") < 200 ||
					       (datasetnamebkg.find("output_WJ") < 200 && datasetnamebkg.find(whichenergy.c_str())<200)
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
       
       if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt-15to20_MuPt5Enriched") < 200 ){ 
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_qcdMu); 
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcdMu); 
       }
       
       if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt-40_doubleEMEnriched") < 200 ){ 
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_qcdDEM); 
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcdDEM); 
       }
       
       if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt_20to30_BCtoE") < 200 ){     
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_qcdBC);   
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcdBC);
       }   

       if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt_1000to1400") < 200 ){
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_qcd);   
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcd);
       }
       
       if(datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_ST_") < 200 && datasetnamebkg.find("t-channel") < 200 ){     	 
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_singlet); 
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_singlet); 
       }     
     }
     
     // htotal->Add(hfourlepbestmass_4l_afterSel_new_new);
     
  }


  // Building the ratio
  for (int nbins=1;nbins<=htotaldata->GetNbinsX(); nbins++){
    //cout << "Total: BinCenter=" << htotalHisto->GetBinCenter(nbins) << " BinContent=" << htotalHisto->GetBinContent(nbins) << " BinErrorContent=" << htotalHisto->GetBinError(nbins) << endl;
    if (htotalHisto->GetBinContent(nbins)>0.) {
      htotalHistoRatio->SetBinContent(nbins,double(htotaldata->GetBinContent(nbins)/htotalHisto->GetBinContent(nbins)));
      //htotalHistoRatio->SetBinError(nbins,double(sqrt(htotaldata->GetBinContent(nbins))/htotalHisto->GetBinContent(nbins)));
      htotalHistoRatio->SetBinError(nbins,double(sqrt( 
		 (1./(htotalHisto->GetBinContent(nbins)*htotalHisto->GetBinContent(nbins)))*htotaldata->GetBinContent(nbins) +
		 (htotaldata->GetBinContent(nbins)*htotaldata->GetBinContent(nbins)/pow(htotalHisto->GetBinContent(nbins),4))
		 *htotalHisto->GetBinContent(nbins) 						      
		 )));
		      
    }
  }


    
    //cout << "Nbins=" << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << endl;
           

 

  htotal->Draw("hist same");
//  htotaldata->Draw("EPsame");
  gr->Draw("EPsame");
  

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
  double canvasratio = 0.3;
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
  htotalHistoRatio->Draw("Psame");


  c1->Update();
  
  if (histosdir.find("4mu")<25) whichchannel="4mu";
  if (histosdir.find("2e2mu")<25) whichchannel="2e2mu";

  std::string saveaspdfratio = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.pdf" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.pdf";
  std::string saveaspngratio = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.png" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.png";
  std::string saveasepsratio = useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.eps" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.eps";
  std::string saveasrootratio= useLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.root": "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.root";

  cout << saveasrootratio.c_str() << endl;
  c1->SaveAs(saveaspdfratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.pdf"*/);
  c1->SaveAs(saveaspngratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.png"*/);
  c1->SaveAs(saveasepsratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.eps"*/);
  c1->SaveAs(saveasrootratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.root"*/);
  

  
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

  whichchannel="2mu2b";
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
	    inputfilename.find("WWTo2L2Nu")<100){ 
      //Vdatasetnamebkg.push_back(inputfilename);
      

      if(inputfilename.find("DYJetsToLL")<200){
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("DY");
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
          Vcolorbkg.push_back(kSpring);
      }


      if(inputfilename.find("WWTo2L2Nu")<200){
	Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("WW");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kCyan+3);
        //cout << "VH" << endl;
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
