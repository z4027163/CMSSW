#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>

// #include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/GsfElectronRegressionEnergyProducer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/PhiPlusPlusMCGenProducer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/SimpleCounter.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/GenAnalyzer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/ZbbGenAnalyzer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/ZccGenAnalyzer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/ZjetsGenAnalyzer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsRunEventFilter.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsHLTAnalysis.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsHLTAnalysisFilter.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsHLTInfo.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronSelector.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronOrdering.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsPFtoRECOMuon.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsPFfsrPhoton.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonCalibrator.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonSelector.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronIsolationProducerEgamma.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronIsolationTest.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonIsolationProducerMu.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonPreselection.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonPreselectionFilter.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsVtxProducer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsIpToVtxProducer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsTipLipToVtxProducer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonOfflineSelection.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonOfflineSelectionFilter.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsBestCandidate.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCP.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMCGenFilter.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonRootTree.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsRootTree.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsConstraintFitProducer.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/ConvValueMapProd.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCandViewCleaner.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsPFJetSelector.h>
#include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsZZMassErrors.h>
// #include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/GsfCalibratedElectronProducer.h>
// #include <HiggsAnalysis/HiggsToZZ4Leptons/plugins/RegressionElectronProducer.h>

DEFINE_FWK_MODULE(PhiPlusPlusMCGenProducer);
DEFINE_FWK_MODULE(SimpleCounter);
DEFINE_FWK_MODULE(HZZ4LeptonsCandViewCleaner);
DEFINE_FWK_MODULE(ConvValueMapProd);
DEFINE_FWK_MODULE(GenAnalyzer);
DEFINE_FWK_MODULE(ZbbGenAnalyzer);
DEFINE_FWK_MODULE(ZccGenAnalyzer);
DEFINE_FWK_MODULE(ZjetsGenAnalyzer);
DEFINE_FWK_MODULE(HZZ4LeptonsRunEventFilter);
DEFINE_FWK_MODULE(HZZ4LeptonsHLTAnalysis);
DEFINE_FWK_MODULE(HZZ4LeptonsHLTAnalysisFilter);
DEFINE_FWK_MODULE(HZZ4LeptonsHLTInfo);
DEFINE_FWK_MODULE(HZZ4LeptonsElectronSelector);
DEFINE_FWK_MODULE(HZZ4LeptonsElectronOrdering);
DEFINE_FWK_MODULE(HZZ4LeptonsPFtoRECOMuon);
DEFINE_FWK_MODULE(HZZ4LeptonsPFfsrPhoton);
DEFINE_FWK_MODULE(HZZ4LeptonsMuonCalibrator);
DEFINE_FWK_MODULE(HZZ4LeptonsMuonSelector);
DEFINE_FWK_MODULE(HZZ4LeptonsElectronIsolationProducerEgamma);
DEFINE_FWK_MODULE(HZZ4LeptonsElectronIsolationTest);
DEFINE_FWK_MODULE(HZZ4LeptonsMuonIsolationProducerMu);
DEFINE_FWK_MODULE(HZZ4LeptonsCommonPreselection);
DEFINE_FWK_MODULE(HZZ4LeptonsCommonPreselectionFilter);
DEFINE_FWK_MODULE(HZZ4LeptonsVtxProducer);
DEFINE_FWK_MODULE(HZZ4LeptonsIpToVtxProducer);
DEFINE_FWK_MODULE(HZZ4LeptonsTipLipToVtxProducer);
DEFINE_FWK_MODULE(HZZ4LeptonsCommonOfflineSelection);
DEFINE_FWK_MODULE(HZZ4LeptonsCommonOfflineSelectionFilter);
DEFINE_FWK_MODULE(HZZ4LeptonsBestCandidate);
DEFINE_FWK_MODULE(HZZ4LeptonsCP);
DEFINE_FWK_MODULE(HZZ4LeptonsMCGenFilter);
DEFINE_FWK_MODULE(HZZ4LeptonsCommonRootTree);
DEFINE_FWK_MODULE(HZZ4LeptonsRootTree);
DEFINE_FWK_MODULE(HZZ4LeptonsConstraintFitProducer);
DEFINE_FWK_MODULE(HZZ4LeptonsPFJetSelector);
// DEFINE_FWK_MODULE(HZZ4LeptonsZZMassErrors);
// DEFINE_FWK_MODULE(CalibratedElectronProducer);
// DEFINE_FWK_MODULE(RegressionElectronProducer);
// DEFINE_FWK_MODULE(GsfElectronRegressionEnergyProducer);
