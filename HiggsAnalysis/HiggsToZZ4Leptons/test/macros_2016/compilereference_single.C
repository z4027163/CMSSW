#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <TDCacheFile.h>

#endif

using namespace std;

int main(int argc, char ** argv){
  

  string dataconf=argv[1];
  cout << dataconf << " dataconfiguration" <<endl;

  string mcconf=argv[2];
  cout << mcconf << " mc configuration" <<endl;
  
  Char_t nome[300];

//  sprintf(nome,"/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v2/170508_174010/DYJetsToLL_merge.root");
//  sprintf(nome, "/eos/uscms/store/user/wangz/TT_TuneCUETP8M1_13TeV-powheg-pythia8-evtgen/TT_TuneCUETP8M1_13TeV-powheg-pythia8-evtgen_v2_part1/170508_173643/TT_TuneCUETP8M1_merge.root");
  sprintf(nome,"/eos/uscms/store/user/wangz/DoubleMuon/DoubleMuon_Run2016D-23Sep2016-v2_part1/170512_000801/DoubleMuon.root");
    

  cout << "test01" << endl;
  TFile *file3;
  file3 = TFile::Open(nome);
  
  cout << "test02" << endl;
  TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
  cout << "test03" << endl;
  HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf);
  //HZZ4LeptonsAnalysis make3(tree3);

  sprintf(nome,"output_DoubleMuon.root");
  make3.Loop(nome);

  cout << "Create file with name: " << nome << endl;
  delete tree3;
  file3 -> Close();

  return 0; 

}

