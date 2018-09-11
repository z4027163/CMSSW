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

void checkweight(){

   TString infilename="/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_dm_mcweight_noMuCal_v4/171128_184521/012.root";
   TFile *f1 = new TFile(infilename);
   TTree *t1 = (TTree*)f1->Get("HZZ4LeptonsAnalysis");

   float mc_weight;
   float temp=0;

   t1->SetBranchAddress("MC_weighting", &mc_weight);

   int entries = t1->GetEntries();
   for(int i=0; i < entries; i++){
       if(i%1000==0) cout << i << endl;
       t1->GetEntry(i);
       if(i>0){ if(abs(mc_weight)!=abs(temp)) cout << "different " << mc_weight << " and " << temp << endl;} 
       temp=mc_weight;
   }
} 
