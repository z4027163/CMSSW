#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_llbb.h"
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

//    sprintf(nome,"/eos/uscms/store/user/wangz/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_regression_v7/170712_071915/DYJetsToLL.root");
  sprintf(nome,"/eos/uscms/store/user/wangz/DoubleEG/DoubleEG_Run2016E-03Feb2017-v2_part2/170714_174048/DoubleEG_E_part2.root");
//    sprintf(nome,"/eos/uscms/store/user/wangz/DoubleMuon/DoubleMuon_Run2016E-03Feb2017-v2/170714_163811/DoubleMuon_E_par2.root");
//    sprintf(nome,"/eos/uscms/store/user/wangz/DoubleEG/DoubleEG_Run2016D-03Feb2017-v2/170710_194246/DoubleEG.root");
//  sprintf(nome,"/eos/uscms/store/user/wangz/DoubleMuon/DoubleMuon_Run2016D-03Feb2017-v1/170712_072212/DoubleMuon.root");  
//  sprintf(nome,"/eos/uscms/store/user/wangz/DoubleMuon/DoubleMuon_Run2016C-03Feb2017-v5/170712_211425/678.root");

  cout << "test01" << endl;
  TFile *file3;
  file3 = TFile::Open(nome);
  
  cout << "test02" << endl;
  TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
  cout << "test03" << endl;
  HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf);
  //HZZ4LeptonsAnalysis make3(tree3);

  sprintf(nome,"output_DoubleEG_E_part2.root");

  make3.Loop(nome);

  cout << "Create file with name: " << nome << endl;
  delete tree3;
  file3 -> Close();

  return 0; 

}

