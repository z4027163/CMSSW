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

#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>

using namespace std;

void draw_genZ(){

  TH2F *hframe=NULL,*hframe2=NULL;
//  float bins[12]={0,5,8,10,15,20,30,40,50,65,85,150};
  float bins[13]={0,2,5,8,10,15,20,30,40,50,65,85,150};
  TH1F *h1 = new TH1F("genZpt_mu_2017","genZpt_mu_2017",12,bins);
  TH1F *h2 = new TH1F("genZpt_e_2017","genZpt_e_2017",12,bins);
  TH1F *h3 = new TH1F("genZpt_mu_2016","genZpt_e_2016",12,bins);
  TH1F *h4 = new TH1F("genZpt_e_2016","genZpt_e_2016",12,bins);

  TFile *f1 = new TFile("/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_ext_v1/190430_164515/1.root");
  TTree *t1 = (TTree*)f1->Get("HZZ4LeptonsAnalysis");

  int entries1 = t1->GetEntries();

  Float_t         MC_Z_PT1[2][5];
  Float_t         MC_Z_PDGID1[2][5];
  Float_t         MC_weight1;

  t1->SetBranchAddress("MC_weighting", &MC_weight1);
  t1->SetBranchAddress("MC_Z_PT", MC_Z_PT1);
  t1->SetBranchAddress("MC_Z_PDGID", MC_Z_PDGID1);
  
//  entries1=100;

  double scale1=0, scale2=0,scale3=0,scale4=0;
  for(int i=0; i < entries1 ; i++){
    t1->GetEntry(i);
    if(i%10000==0) cout << "entry " << i << endl;
    if(MC_Z_PT1[0][0]>=0 && abs(MC_Z_PDGID1[0][1])==13) {h1->Fill(MC_Z_PT1[0][0],MC_weight1);scale1+=MC_weight1;}
    if(MC_Z_PT1[0][0]>=0 && abs(MC_Z_PDGID1[0][1])==11) {h2->Fill(MC_Z_PT1[0][0],MC_weight1);scale2+=MC_weight1;}
  }

  TFile *f2 = new TFile("/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v8/180913_180905/345.root");

  TTree *t2 = (TTree*)f2->Get("HZZ4LeptonsAnalysis");

  Float_t         MC_Z_PT2[2][5];
  Float_t         MC_Z_PDGID2[2][5];
  Float_t         MC_weight2;

  int entries2 = t2->GetEntries();
  t2->SetBranchAddress("MC_weighting", &MC_weight2);
  t2->SetBranchAddress("MC_Z_PT", MC_Z_PT2);
  t2->SetBranchAddress("MC_Z_PDGID", MC_Z_PDGID2);

//  entries2=100;
  for(int i=0; i < entries2 ; i++){
    t2->GetEntry(i);
    if(i%10000==0) cout << "entry " << i << endl;
    if(MC_Z_PT2[0][0]>=0 && abs(MC_Z_PDGID2[0][1])==13) {h3->Fill(MC_Z_PT2[0][0],MC_weight2);scale3+=MC_weight2;}
    if(MC_Z_PT2[0][0]>=0 && abs(MC_Z_PDGID2[0][1])==11) {h4->Fill(MC_Z_PT2[0][0],MC_weight2);scale4+=MC_weight2;}
  }

  h1->Scale(1/scale1);
  h2->Scale(1/scale2);
  h3->Scale(1/scale3);
  h4->Scale(1/scale4);

  h1->SetLineColor(1);
  h1->SetLineWidth(2);
  h2->SetLineColor(2);
  h2->SetLineWidth(2);
  h3->SetLineColor(3);
  h3->SetLineWidth(2);
  h4->SetLineColor(4);
  h4->SetLineWidth(2);

  TLegend *legend = new TLegend(0.6,0.70,0.85,0.87);
  legend->AddEntry(h1,"2017 Zmm","lep");
  legend->AddEntry(h2,"2017 Zee","lep");
  legend->AddEntry(h3,"2016 Zmm","lep");
  legend->AddEntry(h4,"2016 Zee","lep");
  legend->SetBorderSize(0);

  TH1F *hd1 = NULL, *hn1=NULL;
  hn1 = (TH1F*)h1->Clone();
  hd1 = (TH1F*)h3->Clone();
  hn1->Divide(hd1);

  hn1->SetStats(kFALSE);
  hd1->SetStats(kFALSE); 

  TH1F *hd2 = NULL, *hn2=NULL;
  hn2 = (TH1F*)h2->Clone();
  hd2 = (TH1F*)h4->Clone();
  hn2->Divide(hd2);

  hn2->SetStats(kFALSE);
  hd2->SetStats(kFALSE);

  h1->SetStats(kFALSE);
  h2->SetStats(kFALSE);
  h3->SetStats(kFALSE);
  h4->SetStats(kFALSE);

  TCanvas *c1 = new TCanvas("c1","c1",600,800);
  c1->cd();
  c1->SetLogy(1);
  c1->SetTicks(1,1);
 
  h1->Draw("HIST");
  h2->Draw("HIST same");
  h3->Draw("HIST same");
  h4->Draw("HIST same");

  legend->Draw();

  double canvasratio = 0.3;
  c1->SetBottomMargin(canvasratio + (1-canvasratio)*c1->GetBottomMargin()-canvasratio*c1->GetTopMargin());

  canvasratio = 0.16;
  TPad *ratioPad = new TPad("BottomPad","",0,0,1,1);
  ratioPad->SetTopMargin((1-canvasratio) - (1-canvasratio)*ratioPad->GetBottomMargin()+canvasratio*ratioPad->GetTopMargin());
  ratioPad->SetFillStyle(4000);
  ratioPad->SetFillColor(4000);
  ratioPad->SetFrameFillColor(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->SetFrameBorderMode(0);
  ratioPad->SetTicks(1,1);
  ratioPad->Draw();
  ratioPad->cd();

  hn1->SetMarkerStyle(20);
  hn1->SetMarkerSize(0.8);
  hn1->SetMarkerColor(1);
  hn1->SetMinimum(0);
  hn1->SetMaximum(2);
  hn1->SetTitle("2017 vs 2016");
  hn1->Draw("Psame");

  hn2->SetMarkerStyle(20);
  hn2->SetMarkerSize(0.8);
  hn2->SetMarkerColor(2);
  hn2->Draw("Psame");

  c1->Update();
  c1->SaveAs("genZPt.png");

  TFile * theFile = new TFile("genZpt_2017_2016.root","RECREATE");
  theFile->cd();
  h1->Write(0,TObject::kOverwrite);
  h2->Write(0,TObject::kOverwrite);
  h3->Write(0,TObject::kOverwrite);
  h4->Write(0,TObject::kOverwrite);
  hn1->Write("sf_Zmm",TObject::kOverwrite);
  hn2->Write("sf_zee",TObject::kOverwrite);
  theFile->Close();
}
  
