#include "TRandom3.h"
#include "TMath.h"
#include "RoccoR.cc"
#include <iostream>

void random_test(){
    
    printf("test1\n");   
    RoccoR  rc("/uscms/home/zwang4/nobackup/WORKSPCACE/ntuple/CMSSW_8_0_24/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros/roccor/rcdata.2016.v3");
    int charge = -1;
    double pt = 50.0;
    double eta = 1.3;
    double phi = 1.0;
    printf("test2\n");
    double dataSF = rc.kScaleDT(charge, pt, eta, phi);
    std::cout << "SF=" << dataSF << endl;
}
