# HiggsAnalysis-HiggsToZZ4LeptonsRun2
Package for  H->ZZ->4l anaysis for Run2

git clone https://github.com/cms-analysis/HiggsAnalysis-HiggsToZZ4LeptonsRun2 HiggsAnalysis/HiggsToZZ4Leptons -b origin-CMSSW_7_6_3
git clone -b 74x-root6 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit 
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement -b V00-02-01-patch1
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V2 KaMuCa 
cp /afs/cern.ch/user/n/ndefilip/public/KalmanMuonCalibrator.cc $CMSSW_BASE/src/KaMuCa/Calibration/src/.
git clone https://github.com/tocheng/KinZfitter.git
git cms-merge-topic -u matteosan1:smearer_76X

cp -r $CMSSW_RELEASE_BASE/src/RecoEgamma/ElectronIdentification .

cd  RecoEgamma/ElectronIdentification/data/Spring15
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml .
cp $CMSSW_RELEASE_BASE/external/${SCRAM_ARCH}/data/RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml .
