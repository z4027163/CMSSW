#ifndef HZZ4LeptonsCommonRootTree_h
#define HZZ4LeptonsCommonRootTree_h

/** \class  HZZ4LeptonsCommonRootTree
 *
 *  Root Tree for H->ZZ->4l analysis.
 *
 *  Author: N. De Filippis - Politecnico and INFN Bari
 *          
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/DataKeyTags.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/FileBlock.h"

// MCTruth
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

// Data format
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"    
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
//#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
//#include "DataFormats/METReco/interface/PFMET.h"
//#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/GenJet.h"


// Trigger
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

// Maps
#include "DataFormats/Common/interface/ValueMap.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

// user include files
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

// Transient tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//PU Jet ID
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"

//Full Error
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

#include <TMatrixD.h>

//chisquare
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
class MultiTrajectoryStateMode ;
class EgammaTowerIsolation ;


// Class to create TTree variables
#include <TFile.h> 
#include <TTree.h> 

//DEBUG REMOVE ME!!!
#include <exception>

// Namespaces
using namespace reco;
using namespace std;
using namespace pat;




class HZZ4LeptonsCommonRootTree : public edm::EDAnalyzer {
  
 public:
  
  HZZ4LeptonsCommonRootTree(const edm::ParameterSet& pset);
  
  ~HZZ4LeptonsCommonRootTree();
  
  void analyze(const edm::Event& e, const edm::EventSetup& c);
  void beginJob();
  void endJob();
  
  void respondToOpenInputFile(edm::FileBlock const& fb) {
    inputfileName = fb.fileName();
    cout << "Input Filename is=" << inputfileName.c_str() << endl;
  }
	
  void ReadParameters(const edm::ParameterSet& pset){ 
    std::cout << "Reading parameters from cfg" << std::endl;
    
    typedef std::vector<edm::InputTag> vtag;
    // Get the various input parameters
    decaychannel              = pset.getParameter<std::string>("decaychannel");
    rootFileName              = pset.getUntrackedParameter<std::string>("rootFileName");
    useRECOformat             = pset.getUntrackedParameter<bool>("useRECOformat");

    // search param in cfg
    module_to_search=pset.getUntrackedParameter<vector<std::string> >("module_to_search");
    par_to_search= pset.getUntrackedParameter<std::string>("par_to_search");
    

    // Get PU simulation info
    fillPUinfo                = pset.getUntrackedParameter<bool>("fillPUinfo");
    PileupSrc_                = consumes<std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("PileupSrc"));

    // Generator
    generator_                = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("Generator"));
    lheEventProductToken_     = consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("lheEventProduct"));
    // Get HLT flags
    fillHLTinfo               = pset.getUntrackedParameter<bool>("fillHLTinfo");
    HLTInfoFired              = pset.getParameter<edm::InputTag>("HLTInfoFired");
    HLTAnalysisinst           = pset.getParameter<string>("HLTAnalysisinst");
    flagHLTnames              = pset.getParameter<vtag>("flagHLTnames");
    HLTFilter_                 = pset.getParameter<std::vector<std::string>>("HLTFilter");
    fillLHEinfo               = pset.getUntrackedParameter<bool>("fillLHEinfo");
    filljec               = pset.getUntrackedParameter<bool>("filljec");

    // Get HLT matching
//AOD    triggerEvent              = consumes<trigger::TriggerEvent >(pset.getParameter<edm::InputTag>("triggerEvent"));
    triggerObjects_              = consumes<pat::TriggerObjectStandAloneCollection>(pset.getParameter<edm::InputTag>("triggerobjects"));

    triggerBits_              = consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("triggerbits"));
    triggerPrescales_         = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescales"));
    triggerPrescalesL1min_         = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescalesl1min"));
    triggerPrescalesL1max_         = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescalesl1max"));
   
    triggerFilter             = pset.getParameter<std::string>("triggerFilter");
    triggerMatchObject        = consumes<edm::Association<std::vector<pat::TriggerObjectStandAlone> > >(pset.getParameter<edm::InputTag>("triggerMatchObject"));
    triggerMatchObjectEle     = pset.getParameter<edm::InputTag>("triggerMatchObjectEle");
    triggerHLTcollection      = pset.getParameter<std::string>("triggerHLTcollection");

    // Get flags
    flaginst                  = pset.getParameter<std::string>("flaginst");
    flagtags                  = pset.getParameter<std::vector<std::string> >("flagtags");
    // MCtruth tags
    fillMCTruth               = pset.getUntrackedParameter<bool>("fillMCTruth");
    MCcollName                = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("MCcollName"));
    genParticles_             = consumes<vector<reco::GenParticle> >(pset.getParameter<edm::InputTag>("genParticles"));
    fourgenleptons_           = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("fourgenleptons"));
    digenZ_                   = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("digenZ"));

    
    // RECO tags
    RECOcollNameBest2e2mu     = pset.getParameter<vtag>("RECOcollNameBest2e2mu");
    RECOcollNameBest4mu       = pset.getParameter<vtag>("RECOcollNameBest4mu");
    RECOcollNameBest4e        = pset.getParameter<vtag>("RECOcollNameBest4e");
 
    // RECO additional tags
    useAdditionalRECO         = pset.getUntrackedParameter<bool>("useAdditionalRECO");
    RECOcollNameZ             = pset.getParameter<vtag>("RECOcollNameZ");
    RECOcollNameZss           = pset.getParameter<vtag>("RECOcollNameZss");
    RECOcollNameDiLep_         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameDiLep"));

    RECOcollNameDiLep         = pset.getParameter<edm::InputTag>("RECOcollNameDiLep");
    RECOcollNameMMMM_         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameMMMM"));
    RECOcollNameEEEE          = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameEEEE")); 
    RECOcollNameEEMM          = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameEEMM"));
    RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
    RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
    RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
    RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
    RECOcollNameLLLL         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameLLLL"));   

//    RECOcollNameMMMM          = pset.getParameter<vtag>("RECOcollNameMMMM");
//    RECOcollNameEEEE          = pset.getParameter<vtag>("RECOcollNameEEEE");
//    RECOcollNameEEMM          = pset.getParameter<vtag>("RECOcollNameEEMM");
    RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
    RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
    RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
    RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
//    RECOcollNameLLLL          = pset.getParameter<edm::InputTag>("RECOcollNameLLLL");

    // electrons and muons tags
    use2011EA                 = pset.getUntrackedParameter<bool>("use2011EA");

 
    muonTag_                  = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("MuonsLabel"));
    muonCorrPtErrorMapTag_    = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsCorrPtErrorMapLabel"));
    
    muonPFTag_                  = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("PFMuonsLabel"));
/*
    muonPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueChargedAll"));
    muonPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueCharged"));
    muonPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueNeutral"));
    muonPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueGamma"));
    muonPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValuePU"));
 */   
    
    clusterCollectionTag_    = pset.getParameter<edm::InputTag>("SuperClustersLabel");
    gsftrackCollection_      = pset.getParameter<edm::InputTag>("GsfTracksElectronsLabel");

 electronEgmTag_          = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("ElectronsEgmLabel"));
    electronEgmTkMapTag_     = pset.getParameter<edm::InputTag>("ElectronsEgmTkMapLabel");
    electronEgmEcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmEcalMapLabel");
    electronEgmHcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmHcalMapLabel");
    
    mvaElectronTag_          = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("mvaElectronTag"));
//    mvaTrigV0MapTag_         = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaTrigV0MapTag"));
//    mvaNonTrigV0MapTag_      = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaNonTrigV0MapTag"));
/*    
    electronPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueChargedAll"));
    electronPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueCharged"));
    electronPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueNeutral"));
    electronPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueGamma"));
    electronPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValuePU"));
  */  
//    eleRegressionEnergyErrorTag_= pset.getParameter<edm::InputTag>("eleRegressionEnergyErrorLabel");
//    eleRegressionEnergyTag_     = pset.getParameter<edm::InputTag>("eleRegressionEnergyLabel");
    
    
    // PF photons
//    pfphotonsTag_                 = consumes<edm::View<reco::PFCandidate>>(pset.getParameter<edm::InputTag>("PFPhotonsLabel"));
    pfTag_                 = consumes<edm::View<pat::PackedCandidate>>(pset.getParameter<edm::InputTag>("pfCands"));
/*    photonPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueChargedAll"));
    photonPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueCharged"));
    photonPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueNeutral"));
    photonPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueGamma"));
    photonPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValuePU"));
*/   

    // vertexing 
    // 3D w.r.t primary vertex DA
    muonTag_Vert           = pset.getParameter<edm::InputTag>("MuonsLabelVert");


    muonMapTag_Vert        = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVert"));
    muonMapTag_VertValue   = consumes<edm::ValueMap<float> >( pset.getParameter<edm::InputTag>("MuonsMapLabelVertValue"));
    muonMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertError"));
/*
    // KF
    muonMapTag_VertKF        = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertKF"));
    muonMapTag_VertValueKF   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertValueKF"));
    muonMapTag_VertErrorKF   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertErrorKF"));
    
    // STIP SLIP
    muonSTIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVert"));
    muonSLIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVert"));
    
    muonSTIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertValue"));
    muonSLIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertValue"));
    muonSTIPMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertError"));
    muonSLIPMapTag_VertError   =consumes<edm::ValueMap<float> >( pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertError"));
*/    
    // 3D w.r.t primary vertex DA
    
/*    
    // geom. discr.
    StandardFitVertexDiLep   = pset.getParameter<edm::InputTag>("StandardFitVertexDiLep");
*/
    //electronID
    eleIDTag_                 = pset.getParameter<vtag>("eleIDLabel");

    // Other objets
    photonsTag_              = consumes<edm::View<pat::Photon> >(pset.getParameter<edm::InputTag>("PhotonsLabel"));
//    tracksTag_               = consumes<vector<reco::Track> >(pset.getParameter<edm::InputTag>("TracksLabel"));
    jetsTag_                 = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsLabel"));
 //   jetsDataTag_             = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsDataLabel"));
    jetsMVATag_              = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsMVALabel"));
//    PuJetMvaMCfullDiscr_     = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaMCfullDiscrLabel"));
//    PuJetMvaMCfullId_        = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaMCfullIdLabel"));
//    PuJetMvaDatafullDiscr_   = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaDatafullDiscrLabel"));
//    PuJetMvaDatafullId_      = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaDatafullIdLabel"));
 
    rhojetsTag_              = consumes<double>(pset.getParameter<edm::InputTag>("RhoJetsLabel"));
    verticesTag_             = consumes<vector<reco::Vertex> >(pset.getParameter<edm::InputTag>("VerticesLabel"));

    // GenJet 
    genjetTag_              = consumes<vector<reco::GenJet> >(pset.getParameter<edm::InputTag>("GenJetLabel")); // GenJet

/*    
    // MET reco
    genmetTag_              = consumes<vector<reco::GenMET> >(pset.getParameter<edm::InputTag>("GenMETLabel")); // GenMET
    trackermetTag_          = pset.getParameter<edm::InputTag>("TrackerMETLabel"); // TrackerMET
    // Calo MET
    calometTag_             = pset.getParameter<edm::InputTag>("CaloMETLabel");
    calometnohfTag_         = pset.getParameter<edm::InputTag>("CaloMET_NoHFLabel");
    useAdditionalMET_       = pset.getUntrackedParameter<bool>("useAdditionalMET");
    calomethoTag_           = pset.getParameter<edm::InputTag>("CaloMET_HOLabel");
    calometoptTag_          = pset.getParameter<edm::InputTag>("CaloMET_OptLabel");
    calometoptnohfTag_      = pset.getParameter<edm::InputTag>("CaloMET_OptNoHFLabel");
    calometoptnohfhoTag_    = pset.getParameter<edm::InputTag>("CaloMET_OptNoHFHOLabel");
    calometopthoTag_        = pset.getParameter<edm::InputTag>("CaloMET_OptHOLabel");
    calometnohfhoTag_       = pset.getParameter<edm::InputTag>("CaloMET_NoHFHOLabel");
*/    
    // PF MET
    pfmetTag_               = consumes<vector<pat::MET> >(pset.getParameter<edm::InputTag>("PfMETLabel"));
/*
    // htMET
    htmetic5Tag_            = pset.getParameter<edm::InputTag>("HtMET_IC5Label");
    htmetkt4Tag_            = pset.getParameter<edm::InputTag>("HtMET_KT4Label");
    htmetkt6Tag_            = pset.getParameter<edm::InputTag>("HtMET_KT6Label");
    htmetsc5Tag_            = pset.getParameter<edm::InputTag>("HtMET_SC5Label");
    htmetsc7Tag_            = pset.getParameter<edm::InputTag>("HtMET_SC7Label");
    // JES correction MET
    jescormetic5Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorIC5CaloJetLabel");
    jescormetkt4Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorKT4CaloJetLabel");
    jescormetkt6Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorKT6CaloJetLabel");
    jescormetsc5Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorSC5CaloJetLabel");
    jescormetsc7Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorSC7CaloJetLabel");
    // Type I Muon Correction MET
    cormetMuTag_            = pset.getParameter<edm::InputTag>("CorMETGlobalMuLabel");
*/    
    // btagging
//   bDiscriminators_         =pset.getParameter<std::vector<std::string> >("bDiscriminators");

   tCHighEff_bTag_         =pset.getParameter<std::string>("tCHighEff_bTagLabel");
   tCHighPur_bTag_         =pset.getParameter<std::string>("tCHighPur_bTagLabel");
   cSV_bTag_         =pset.getParameter<std::string>("cSV_bTagLabel");
//    tCHighPur_bTag_         = consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(pset.getParameter<edm::InputTag>("tCHighPur_bTagLabel"));
//    cSV_bTag_               = consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(pset.getParameter<edm::InputTag>("cSV_bTagLabel"));

    // Conversion finder
//    ConvMapDistTag_       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ConvMapDist"));
//    ConvMapDcotTag_       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ConvMapDcot"));
    
    // Matching

    goodElectronMCMatch_  = consumes<edm::Association<std::vector<reco::GenParticle> >>(pset.getParameter<edm::InputTag>("goodElectronMCMatch"));
    myElectrons_          = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myElectrons"));
    goodMuonMCMatch_      = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodMuonMCMatch"));
    myMuons_              = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myMuons"));
    goodGammaMCMatch_     = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodGammaMCMatch"));
    myGammas_             = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myGammas"));


   //mc 

   prunedGenToken_  =consumes<edm::View<reco::GenParticle> >(pset.getParameter<edm::InputTag>("pruned"));
   packedGenToken_  =consumes<edm::View<pat::PackedGenParticle> >(pset.getParameter<edm::InputTag>("packed"));

   goodZtoMuMuMCMatch_    = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodZtoMuMuMCMatch"));
 
   goodZtoEEMCMatch_    = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodZtoEEMCMatch"));

   goodHiggsTozzToEEMMMCMatch_ = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodHiggsTozzToEEMMMCMatch"));

   goodHiggsTozzToMMMMMCMatch_ = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodHiggsTozzToMMMMMCMatch"));
  
   goodHiggsTozzToEEEEMCMatch_ = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodHiggsTozzToEEEEMCMatch"));

    // Beam Spot
    offlineBeamSpot_      = consumes<reco::BeamSpot>(pset.getParameter<edm::InputTag>("offlineBeamSpot"));
  }
  
  
  void DefineBranches(TTree *Tree_){
    // Run event
    Tree_->Branch("Run",&irun,"irun/i");
    Tree_->Branch("Event",&ievt,"ievt/i");
    Tree_->Branch("LumiSection",&ils,"ils/i");
    Tree_->Branch("Avginstlumi",&Avginstlumi,"Avginstlumi/F");
    
    // MC
    Tree_->Branch("num_PU_vertices",&num_PU_vertices,"num_PU_vertices/I");
    Tree_->Branch("PU_BunchCrossing",&PU_BunchCrossing,"PU_BunchCrossing/I");
    Tree_->Branch("MC_weighting",&MC_weighting,"MC_weighting/F");
    Tree_->Branch("MC_weighting_un",&MC_weighting_un,"MC_weighting_un[9]/F");
    Tree_->Branch("PDF_weighting_un",&PDF_weighting_un,"PDF_weighting_un/F");    
    // HLT 
    Tree_->Branch("RECO_nMuHLTMatch",&RECO_nMuHLTMatch,"RECO_nMuHLTMatch/I");
    Tree_->Branch("RECOMU_PT_MuHLTMatch",RECOMU_PT_MuHLTMatch,"RECOMU_PT_MuHLTMatch[50]/F");
    Tree_->Branch("RECOMU_ETA_MuHLTMatch",RECOMU_ETA_MuHLTMatch,"RECOMU_ETA_MuHLTMatch[50]/F");
    Tree_->Branch("RECOMU_PHI_MuHLTMatch",RECOMU_PHI_MuHLTMatch,"RECOMU_PHI_MuHLTMatch[50]/F");
    Tree_->Branch("RECOMU_sm_MuHLTMatch",RECOMU_sm_MuHLTMatch,"RECOMU_sm_MuHLTMatch[50]/b");
    Tree_->Branch("RECOMU_dm_MuHLTMatch",RECOMU_dm_MuHLTMatch,"RECOMU_dm_MuHLTMatch[50]/I");
   
    Tree_->Branch("dm_trig",&dm_trig,"dm_trig/b");
    Tree_->Branch("sm_trig",&sm_trig,"sm_trig/b");
    Tree_->Branch("de_trig",&de_trig,"de_trig/b");
    Tree_->Branch("se_trig",&se_trig,"se_trig/b");
    Tree_->Branch("tri_trig",&tri_trig,"tri_trig/b");
    Tree_->Branch("jet_trig",&jet_trig,"jet_trig/I");

    Tree_->Branch("RECO_nEleHLTMatch",&RECO_nEleHLTMatch,"RECO_nEleHLTMatch/I");
    Tree_->Branch("RECOELE_PT_EleHLTMatch",RECOELE_PT_EleHLTMatch,"RECOELE_PT_EleHLTMatch[50]/F");
    Tree_->Branch("RECOELE_ETA_EleHLTMatch",RECOELE_ETA_EleHLTMatch,"RECOELE_ETA_EleHLTMatch[50]/F");
    Tree_->Branch("RECOELE_PHI_EleHLTMatch",RECOELE_PHI_EleHLTMatch,"RECOELE_PHI_EleHLTMatch[50]/F");
    Tree_->Branch("RECOELE_se_MuHLTMatch",RECOELE_se_EleHLTMatch,"RECOELE_se_EleHLTMatch[50]/b");
    Tree_->Branch("RECOELE_de_MuHLTMatch",RECOELE_de_EleHLTMatch,"RECOELE_de_EleHLTMatch[50]/I");
    Tree_->Branch("RECOJET_JetHLTMatch",RECOJET_JetHLTMatch,"RECOJET_JetHLTMatch[50]/I");

//    Tree_->Branch("HLTPathsFired",HLTPathsFired,"HLTPathsFired/C");

   
    // MC block 
    Tree_->Branch("MC_E", MC_E, "MC_E[7]/F"); 
    Tree_->Branch("MC_PT", MC_PT, "MC_PT[7]/F"); 
    Tree_->Branch("MC_ETA", MC_ETA, "MC_ETA[7]/F"); 
    Tree_->Branch("MC_THETA", MC_THETA, "MC_THETA[7]/F");
    Tree_->Branch("MC_PHI", MC_PHI, "MC_PHI[7]/F");
    Tree_->Branch("MC_MASS", MC_MASS, "MC_MASS[7]/F");
    Tree_->Branch("MC_PDGID", &MC_PDGID, "MC_PDGID[7]/F");
   

    //  Leptons stable ordered in PT
    Tree_->Branch("MC_LEPT_PT", MC_LEPT_PT, "MC_LEPT_PT[4]/F"); 
    Tree_->Branch("MC_LEPT_ETA", MC_LEPT_ETA, "MC_LEPT_ETA[4]/F"); 
    Tree_->Branch("MC_LEPT_PHI", MC_LEPT_PHI, "MC_LEPT_PHI[4]/F"); 
    Tree_->Branch("MC_LEPT_THETA", MC_LEPT_THETA, "MC_LEPT_THETA[4]/F"); 
    Tree_->Branch("MC_LEPT_PDGID", MC_LEPT_PDGID, "MC_LEPT_PDGID[4]/F");

    // MC Z1, Z2 (first index) and daughter leptons and photon (second index, order: L1, L2, P1, P2 with PT ordering of the leptons and P1 associated to L1)
    Tree_->Branch("MC_Z_PT", MC_Z_PT, "MC_Z_PT[2][5]/F");
    Tree_->Branch("MC_Z_ETA", MC_Z_ETA, "MC_Z_ETA[2][5]/F");
    Tree_->Branch("MC_Z_PHI", MC_Z_PHI, "MC_Z_PHI[2][5]/F");
    Tree_->Branch("MC_Z_THETA", MC_Z_THETA, "MC_Z_THETA[2][5]/F");
    Tree_->Branch("MC_Z_MASS", MC_Z_MASS, "MC_Z_MASS[2][5]/F");
    Tree_->Branch("MC_Z_PDGID", MC_Z_PDGID, "MC_Z_PDGID[2][5]/F");

    // GenJet
    Tree_->Branch( "MC_GENJET_PT",  MC_GENJET_PT,  "MC_GENJET_PT[50]/F");
    Tree_->Branch( "MC_GENJET_ETA", MC_GENJET_ETA, "MC_GENJET_ETA[50]/F");
    Tree_->Branch( "MC_GENJET_PHI", MC_GENJET_PHI, "MC_GENJET_PHI[50]/F");

    // GenMET  
    Tree_->Branch("MC_GENMET", &genmet, "MC_GENMET/F");
    
  
    // Electron block
    Tree_->Branch("RECOELE_E", RECOELE_E, "RECOELE_E[50]/F"); 
    Tree_->Branch("RECOELE_PT",RECOELE_PT,"RECOELE_PT[50]/F");
    Tree_->Branch("RECOELE_PTError",RECOELE_PTError,"RECOELE_PTError[50]/F");
    Tree_->Branch("RECOELE_P",RECOELE_P,"RECOELE_P[50]/F");
    Tree_->Branch("RECOELE_ETA",RECOELE_ETA,"RECOELE_ETA[50]/F"); 
    Tree_->Branch("RECOELE_THETA",RECOELE_THETA,"RECOELE_THETA[50]/F"); 
    Tree_->Branch("RECOELE_PHI",RECOELE_PHI,"RECOELE_PHI[50]/F"); 
    Tree_->Branch("RECOELE_MASS",RECOELE_MASS,"RECOELE_MASS[50]/F"); 
    Tree_->Branch("RECOELE_CHARGE",RECOELE_CHARGE,"RECOELE_CHARGE[50]/F");  
    Tree_->Branch("RECOELE_ID",RECOELE_ID,"RECOELE_ID[50]/F");
    Tree_->Branch("RECOELE_PT_uncorr",RECOELE_PT_uncorr,"RECOELE_PT_uncorr[50]/F");
 
    // Core attributes
    Tree_->Branch("RECOELE_isEcalDriven",   RECOELE_isEcalDriven,   "RECOELE_isEcalDriven[50]/b");   
    Tree_->Branch("RECOELE_isTrackerDriven",RECOELE_isTrackerDriven,"RECOELE_isTrackerDriven[50]/b");  
    Tree_->Branch("RECOELE_gsftrack_NPixHits",RECOELE_gsftrack_NPixHits,"RECOELE_gsftrack_NPixHits[50]/f");
    Tree_->Branch("RECOELE_gsftrack_NStripHits",RECOELE_gsftrack_NStripHits,"RECOELE_gsftrack_NStripHits[50]/f");
    Tree_->Branch("RECOELE_gsftrack_chi2",  RECOELE_gsftrack_chi2, "RECOELE_gsftrack_chi2[50]/F");
    Tree_->Branch("RECOELE_gsftrack_dxyB",  RECOELE_gsftrack_dxyB,"RECOELE_gsftrack_dxyB[50]/F");
    Tree_->Branch("RECOELE_gsftrack_dxy",   RECOELE_gsftrack_dxy,"RECOELE_gsftrack_dxy[50]/F");
    Tree_->Branch("RECOELE_gsftrack_dxyError", RECOELE_gsftrack_dxyError,"RECOELE_gsftrack_dxyError[50]/F");
    Tree_->Branch("RECOELE_gsftrack_dzB",   RECOELE_gsftrack_dzB,"RECOELE_gsftrack_dzB[50]/F");
    Tree_->Branch("RECOELE_gsftrack_dz",    RECOELE_gsftrack_dz,"RECOELE_gsftrack_dz[50]/F");
    Tree_->Branch("RECOELE_gsftrack_dzError", RECOELE_gsftrack_dzError,"RECOELE_gsftrack_dzError[50]/F");
    Tree_->Branch("RECOELE_scl_E",  RECOELE_scl_E,"RECOELE_scl_E[50]/F");
    Tree_->Branch("RECOELE_scl_Et", RECOELE_scl_Et,"RECOELE_scl_Et[50]/F");
    Tree_->Branch("RECOELE_scl_Eta",RECOELE_scl_Eta,"RECOELE_scl_Eta[50]/F");
    Tree_->Branch("RECOELE_scl_Phi",RECOELE_scl_Phi,"RECOELE_scl_Phi[50]/F");    
    Tree_->Branch("RECOELE_ecalEnergy",RECOELE_ecalEnergy,"RECOELE_ecalEnergy[50]/F");

    // Track-Cluster Matching    
/*
    Tree_->Branch("RECOELE_ep",          RECOELE_ep,"RECOELE_ep[50]/F");
    Tree_->Branch("RECOELE_eSeedp",      RECOELE_eSeedp,"RECOELE_eSeedp[50]/F");
    Tree_->Branch("RECOELE_eSeedpout",   RECOELE_eSeedpout,"RECOELE_eSeedpout[50]/F");
    Tree_->Branch("RECOELE_eElepout",    RECOELE_eElepout,"RECOELE_eElepout[50]/F");
    Tree_->Branch("RECOELE_deltaEtaIn",  RECOELE_deltaEtaIn,"RECOELE_deltaEtaIn[50]/F");
    Tree_->Branch("RECOELE_deltaEtaSeed",RECOELE_deltaEtaSeed,"RECOELE_deltaEtaSeed[50]/F");
    Tree_->Branch("RECOELE_deltaEtaEle", RECOELE_deltaEtaEle,"RECOELE_deltaEtaEle[50]/F");
    Tree_->Branch("RECOELE_deltaPhiIn",  RECOELE_deltaPhiIn,"RECOELE_deltaPhiIn[50]/F");
    Tree_->Branch("RECOELE_deltaPhiSeed",RECOELE_deltaPhiSeed,"RECOELE_deltaPhiSeed[50]/F");
    Tree_->Branch("RECOELE_deltaPhiEle", RECOELE_deltaPhiEle,"RECOELE_deltaPhiEle[50]/F");
    // Fiducial flags
    Tree_->Branch("RECOELE_isbarrel",    RECOELE_isbarrel,   "RECOELE_isbarrel[50]/I");   
    Tree_->Branch("RECOELE_isendcap",    RECOELE_isendcap,   "RECOELE_isendcap[50]/I");   
    Tree_->Branch("RECOELE_isGap",       RECOELE_isGap,      "RECOELE_isGap[50]/I");
    Tree_->Branch("RECOELE_isEBetaGap",  RECOELE_isEBetaGap, "RECOELE_isEBetaGap[50]/I");   
    Tree_->Branch("RECOELE_isEBphiGap",  RECOELE_isEBphiGap, "RECOELE_isEBphiGap[50]/I");   
    Tree_->Branch("RECOELE_isEEdeeGap",  RECOELE_isEEdeeGap, "RECOELE_isEEdeeGap[50]/I");   
    Tree_->Branch("RECOELE_isEEringGap", RECOELE_isEEringGap,"RECOELE_isEEringGap[50]/I");   
*/
    // Shower shape
    Tree_->Branch("RECOELE_sigmaIetaIeta", RECOELE_sigmaIetaIeta, "RECOELE_sigmaIetaIeta[50]/F");   
    Tree_->Branch("RECOELE_sigmaEtaEta",   RECOELE_sigmaEtaEta,   "RECOELE_sigmaEtaEta[50]/F");   
    Tree_->Branch("RECOELE_e15",           RECOELE_e15,           "RECOELE_e15[50]/F");   
    Tree_->Branch("RECOELE_e25max",        RECOELE_e25max,        "RECOELE_e25max[50]/F");   
    Tree_->Branch("RECOELE_e55",           RECOELE_e55,           "RECOELE_e55[50]/F");   
    Tree_->Branch("RECOELE_he",            RECOELE_he,            "RECOELE_he[50]/F");   
    Tree_->Branch("RECOELE_r9",            RECOELE_r9,            "RECOELE_r9[50]/F");   
    // Particle flow
    Tree_->Branch("RECOELE_mva", RECOELE_mva,"RECOELE_mva[50]/F");   
    // Brem & Classifaction
    Tree_->Branch("RECOELE_fbrem",  RECOELE_fbrem,  "RECOELE_fbrem[50]/F");   
    Tree_->Branch("RECOELE_nbrems", RECOELE_nbrems, "RECOELE_nbrems[50]/I");   
    //  golden/bigbrem/(narrow)/showering/crack
    Tree_->Branch("RECOELE_Class",  RECOELE_Class,  "RECOELE_Class[50]/I");  
    //fBrem addition
    Tree_->Branch("RECOELE_fbrem_mode", &RECOELE_fbrem_mode,"RECOELE_fbrem_mode[50]/D");
    Tree_->Branch("RECOELE_fbrem_mean", &RECOELE_fbrem_mean,"RECOELE_fbrem_mean[50]/D");
    
    // Isolation
    Tree_->Branch("RECOELE_EGMTRACKISO",RECOELE_EGMTRACKISO,"RECOELE_EGMTRACKISO[50]/F");  
    Tree_->Branch("RECOELE_EGMHCALISO",RECOELE_EGMHCALISO,"RECOELE_EGMHCALISO[50]/F");  
    Tree_->Branch("RECOELE_EGMECALISO",RECOELE_EGMECALISO,"RECOELE_EGMECALISO[50]/F"); 
    Tree_->Branch("RECOELE_EGMX",RECOELE_EGMX,"RECOELE_EGMX[50]/F"); 

    // PF isolation
    Tree_->Branch("RECOELE_PFchAllPart",  RECOELE_PFchAllPart,  "RECOELE_PFchAllPart[50]/D");
    Tree_->Branch("RECOELE_PFchHad",         RECOELE_PFchHad,  "RECOELE_PFchHad[50]/D");
    Tree_->Branch("RECOELE_PFneuHad",        RECOELE_PFneuHad, "RECOELE_PFneuHad[50]/D");
    Tree_->Branch("RECOELE_PFphoton",        RECOELE_PFphoton, "RECOELE_PFphoton[50]/D");
    Tree_->Branch("RECOELE_PFPUchAllPart",   RECOELE_PFPUchAllPart,"RECOELE_PFPUchAllPart[50]/D");
    Tree_->Branch("RECOELE_PFX_dB",          RECOELE_PFX_dB,   "RECOELE_PFX_dB[50]/D");
    Tree_->Branch("RECOELE_PFX_rho",         RECOELE_PFX_rho,  "RECOELE_PFX_rho[50]/D");

    // Electron Regression
    Tree_->Branch("RECOELE_regEnergy",RECOELE_regEnergy,"RECOELE_regEnergy[50]/D");
    Tree_->Branch("RECOELE_regEnergyError",RECOELE_regEnergyError,"RECOELE_regEnergyError[50]/D");
    
    // Vertexing DA and KF
    Tree_->Branch("RECOELE_SIP",RECOELE_SIP,"RECOELE_SIP[50]/F"); 
    Tree_->Branch("RECOELE_IP",RECOELE_IP,"RECOELE_IP[50]/F"); 
    Tree_->Branch("RECOELE_IPERROR",RECOELE_IPERROR,"RECOELE_IPERROR[50]/F"); 
/*
    Tree_->Branch("RECOELE_STIP",RECOELE_STIP,"RECOELE_STIP[50]/F"); 
    Tree_->Branch("RECOELE_SLIP",RECOELE_SLIP,"RECOELE_SLIP[50]/F"); 
    Tree_->Branch("RECOELE_TIP",RECOELE_TIP,"RECOELE_TIP[50]/F"); 
    Tree_->Branch("RECOELE_LIP",RECOELE_LIP,"RECOELE_LIP[50]/F"); 
    Tree_->Branch("RECOELE_TIPERROR",RECOELE_TIPERROR,"RECOELE_TIPERROR[50]/F"); 
    Tree_->Branch("RECOELE_LIPERROR",RECOELE_LIPERROR,"RECOELE_LIPERROR[50]/F"); 
*/
    Tree_->Branch("RECOELE_sclRawE",ele_sclRawE,"RECOELE_sclRawE[50]/D") ;
    Tree_->Branch("RECOELE_sclX",ele_sclX,"RECOELE_sclX[50]/D") ; 
    Tree_->Branch("RECOELE_sclY",ele_sclY,"RECOELE_sclY[50]/D") ; 
    Tree_->Branch("RECOELE_sclZ",ele_sclZ,"RECOELE_sclZ[50]/D") ;
/*
    Tree_->Branch("RECOELE_seedSubdet1",ele_seedSubdet1,"RECOELE_seedSubdet1[50]/D") ;
    Tree_->Branch("RECOELE_seedDphi1",ele_seedDphi1,"RECOELE_seedDphi[50]/D") ; 
    Tree_->Branch("RECOELE_seedDrz1",ele_seedDrz1,"RECOELE_seedDrz1[50]/D") ;
    Tree_->Branch("RECOELE_seedSubdet2",ele_seedSubdet2,"RECOELE_seedSubdet2[50]/D") ;
    Tree_->Branch("RECOELE_seedDphi2",ele_seedDphi2,"RECOELE_seedDphi2[50]/D") ; 
    Tree_->Branch("RECOELE_seedDrz2",ele_seedDrz2,"RECOELE_seedDrz2[50]/D") ;
    Tree_->Branch("RECOELE_eidVeryLoose",ele_eidVeryLoose,"RECOELE_eidVeryLoose[50]/D") ; 
    Tree_->Branch("RECOELE_eidLoose",ele_eidLoose,"RECOELE_eidLoose[50]/D") ; 
    Tree_->Branch("RECOELE_eidMedium",ele_eidMedium,"RECOELE_eidMedium[50]/D") ; 
    Tree_->Branch("RECOELE_eidTight",ele_eidTight,"RECOELE_eidTight[50]/D") ; 
    Tree_->Branch("RECOELE_eidHZZVeryLoose",ele_eidHZZVeryLoose,"RECOELE_eidHZZVeryLoose[50]/D") ; 
    Tree_->Branch("RECOELE_eidHZZLoose",ele_eidHZZLoose,"RECOELE_eidHZZLoose[50]/D") ; 
    Tree_->Branch("RECOELE_eidHZZMedium",ele_eidHZZMedium,"RECOELE_eidHZZMedium[50]/D") ; 
    Tree_->Branch("RECOELE_eidHZZTight",ele_eidHZZTight,"RECOELE_eidHZZTight[50]/D") ; 
*/
    Tree_->Branch("RECOELE_mvaTrigV0",RECOELE_mvaTrigV0,"RECOELE_mvaTrigV0[50]/D") ;     
    Tree_->Branch("RECOELE_mvaNonTrigV0",RECOELE_mvaNonTrigV0,"RECOELE_mvaNonTrigV0[50]/D") ; 
    Tree_->Branch("RECOELE_COV",RECOELE_COV,"RECOELE_COV[50][3][3]/D"); 


    Tree_->Branch("RECOELE_ecalTrkEnergyPreCorr",RECOELE_ecalTrkEnergyPreCorr,"RECOELE_ecalTrkEnergyPreCorr[50]/F");
    Tree_->Branch("RECOELE_ecalTrkEnergyErrPreCorr",RECOELE_ecalTrkEnergyErrPreCorr,"RECOELE_ecalTrkEnergyErrPreCorr[50]/F");
    Tree_->Branch("RECOELE_ecalTrkEnergyErrPostCorr",RECOELE_ecalTrkEnergyErrPostCorr,"RECOELE_ecalTrkEnergyErrPostCorr[50]/F");
    Tree_->Branch("RECOELE_energyScaleValue",RECOELE_energyScaleValue,"RECOELE_energyScaleValue[50]/F");       
    Tree_->Branch("RECOELE_energySigmaValue",RECOELE_energySigmaValue, "RECOELE_energySigmaValue[50]/F");
    Tree_->Branch("RECOELE_energyScaleUp",RECOELE_energyScaleUp,"RECOELE_energyScaleUp[50]/F");     
    Tree_->Branch("RECOELE_energyScaleDown",RECOELE_energyScaleDown,"RECOELE_energyScaleDown[50]/F");       

    Tree_->Branch("RECOELE_energyScaleStatUp",RECOELE_energyScaleStatUp,"RECOELE_energyScaleStatUp[50]/F");       
    Tree_->Branch("RECOELE_energyScaleStatDown",RECOELE_energyScaleStatDown,"RECOELE_energyScaleStatDown[50]/F");        
    Tree_->Branch("RECOELE_energyScaleSystUp",RECOELE_energyScaleSystUp,"RECOELE_energyScaleSystUp[50]/F");        
    Tree_->Branch("RECOELE_energyScaleSystDown",RECOELE_energyScaleSystDown,"RECOELE_energyScaleSystDown[50]/F");        
/*    Tree_->Branch("RECOELE_energyScaleGainUp",RECOELE_energyScaleGainUp,"RECOELE_energyScaleGainUp[50]/F");        
    Tree_->Branch("RECOELE_energyScaleGainDown",RECOELE_energyScaleGainDown,"RECOELE_energyScaleGainDown[50]/F");      
    Tree_->Branch("RECOELE_energyScaleEtUp",RECOELE_energyScaleEtUp,"RECOELE_energyScaleEtUp[50]/F");       
    Tree_->Branch("RECOELE_energyScaleEtDown",RECOELE_energyScaleEtDown,"RECOELE_energyScaleEtDown[50]/F");       
*/
    Tree_->Branch("RECOELE_energySigmaUp",RECOELE_energySigmaUp,"RECOELE_energySigmaUp[50]/F");         
    Tree_->Branch("RECOELE_energySigmaDown",RECOELE_energySigmaDown,"RECOELE_energySigmaDown[50]/F");       
/*
    Tree_->Branch("RECOELE_energySigmaPhiUp",RECOELE_energySigmaPhiUp,"RECOELE_energySigmaPhiUp[50]/F");        
    Tree_->Branch("RECOELE_energySigmaPhiDown",RECOELE_energySigmaPhiDown,"RECOELE_energySigmaPhiDown[50]/F");     
    Tree_->Branch("RECOELE_energySigmaRhoUp",RECOELE_energySigmaRhoUp,"RECOELE_energySigmaRhoUp[50]/F");        
    Tree_->Branch("RECOELE_energySigmaRhoDown",RECOELE_energySigmaRhoDown,"RECOELE_energySigmaRhoDown[50]/F"); 
*/
    Tree_->Branch("RECOELE_TLE_ParentSC_X",RECOELE_TLE_ParentSC_X,"RECOELE_TLE_ParentSC_X[20]/F");
    Tree_->Branch("RECOELE_TLE_ParentSC_Y",RECOELE_TLE_ParentSC_Y,"RECOELE_TLE_ParentSC_Y[20]/F");
    Tree_->Branch("RECOELE_TLE_ParentSC_Z",RECOELE_TLE_ParentSC_Z,"RECOELE_TLE_ParentSC_X[20]/F");

    // Muon block
    Tree_->Branch("RECOMU_isPFMu", RECOMU_isPFMu, "RECOMU_isPFMu[50]/b");
    Tree_->Branch("RECOMU_isGlobalMu", RECOMU_isGlobalMu, "RECOMU_isGlobalMu[50]/b");
    Tree_->Branch("RECOMU_isStandAloneMu", RECOMU_isStandAloneMu, "RECOMU_isStandAloneMu[50]/b");
    Tree_->Branch("RECOMU_isTrackerMu", RECOMU_isTrackerMu, "RECOMU_isTrackerMu[50]/b");
    Tree_->Branch("RECOMU_isCaloMu", RECOMU_isCaloMu, "RECOMU_isCaloMu[50]/b");
    Tree_->Branch("RECOMU_isTrackerHighPtMu", RECOMU_isTrackerHighPtMu, "RECOMU_isTrackerHighPtMu[50]/b");
    Tree_->Branch("RECOMU_isME0Muon", RECOMU_isME0Muon, "RECOMU_isME0Muon[50]/b");
    Tree_->Branch("RECOMU_E",RECOMU_E,"RECOMU_E[50]/F"); 
    Tree_->Branch("RECOMU_PT",RECOMU_PT,"RECOMU_PT[50]/F"); 
    Tree_->Branch("RECOMU_P",RECOMU_P,"RECOMU_P[50]/F"); 
    Tree_->Branch("RECOMU_ETA",RECOMU_ETA,"RECOMU_ETA[50]/F"); 
    Tree_->Branch("RECOMU_THETA",RECOMU_THETA,"RECOMU_THETA[50]/F"); 
    Tree_->Branch("RECOMU_PHI",RECOMU_PHI,"RECOMU_PHI[50]/F"); 
    Tree_->Branch("RECOMU_MASS",RECOMU_MASS,"RECOMU_MASS[50]/F"); 
    Tree_->Branch("RECOMU_CHARGE",RECOMU_CHARGE,"RECOMU_CHARGE[50]/F");

    Tree_->Branch("RECOMU_COV",RECOMU_COV,"RECOMU_COV[50][3][3]/D"); 

    Tree_->Branch("RECOMU_TRACKISO",RECOMU_TRACKISO,"RECOMU_TRACKISO[50]/F");  
    Tree_->Branch("RECOMU_TRACKISO_SUMPT",RECOMU_TRACKISO_SUMPT,"RECOMU_TRACKISO_SUMPT[50]/F");  
    Tree_->Branch("RECOMU_HCALISO",RECOMU_HCALISO,"RECOMU_HCALISO[50]/F");  
    Tree_->Branch("RECOMU_ECALISO",RECOMU_ECALISO,"RECOMU_ECALISO[50]/F"); 
    Tree_->Branch("RECOMU_X",RECOMU_X,"RECOMU_X[50]/F");

    Tree_->Branch("RECOMU_PFchHad",  RECOMU_PFchHad,  "RECOMU_PFchHad[50]/D");
    Tree_->Branch("RECOMU_PFneuHad", RECOMU_PFneuHad, "RECOMU_PFneuHad[50]/D");
    Tree_->Branch("RECOMU_PFphoton", RECOMU_PFphoton, "RECOMU_PFphoton[50]/D");
    Tree_->Branch("RECOMU_PFPUchAllPart",RECOMU_PFPUchAllPart,"RECOMU_PFPUchAllPart[50]/D");
    Tree_->Branch("RECOMU_PFX_dB",   RECOMU_PFX_dB,   "RECOMU_PFX_dB[50]/D");
    Tree_->Branch("RECOMU_PFX_rho",  RECOMU_PFX_rho,  "RECOMU_PFX_rho[50]/D");


    // photon
    Tree_->Branch("RECOPFPHOT_PFchHad",  RECOPFPHOT_PFchHad,  "RECOPFPHOT_PFchHad[20]/D");
    Tree_->Branch("RECOPFPHOT_PFneuHad", RECOPFPHOT_PFneuHad, "RECOPFPHOT_PFneuHad[20]/D");
    Tree_->Branch("RECOPFPHOT_PFphoton", RECOPFPHOT_PFphoton, "RECOPFPHOT_PFphoton[20]/D");
    Tree_->Branch("RECOPFPHOT_PFPUchAllPart",RECOPFPHOT_PFPUchAllPart,"RECOPFPHOT_PFPUchAllPart[20]/D");
    Tree_->Branch("RECOPFPHOT_PFX_rho",  RECOPFPHOT_PFX_rho,  "RECOPFPHOT_PFX_rho[20]/D");

    // vertexing DA and KF
/*
    Tree_->Branch("RECOMU_SIP",RECOMU_SIP,"RECOMU_SIP[50]/F"); 
    Tree_->Branch("RECOMU_IP",RECOMU_IP,"RECOMU_IP[50]/F"); 
    Tree_->Branch("RECOMU_IPERROR",RECOMU_IPERROR,"RECOMU_IPERROR[50]/F"); 

    Tree_->Branch("RECOMU_SIP_KF",RECOMU_SIP_KF,"RECOMU_SIP_KF[50]/F"); 
    Tree_->Branch("RECOMU_IP_KF",RECOMU_IP_KF,"RECOMU_IP_KF[50]/F"); 
    Tree_->Branch("RECOMU_IPERROR_KF",RECOMU_IPERROR_KF,"RECOMU_IPERROR_KF[50]/F"); 

    // GD vertex
    Tree_->Branch("RECOMU_SIP_GD",RECOMU_SIP_GD,"RECOMU_SIP_GD[50]/F"); //2e2mu
    // Std vertex
    Tree_->Branch("RECOMU_SIP_Std",RECOMU_SIP_Std,"RECOMU_SIP_Std[50]/F"); //2e2mu
    // Kin vertex
    Tree_->Branch("RECOMU_SIP_Kin",RECOMU_SIP_Kin,"RECOMU_SIP_Kin[50]/F"); //2e2mu




    Tree_->Branch("RECOMU_STIP",RECOMU_STIP,"RECOMU_STIP[50]/F"); 
    Tree_->Branch("RECOMU_SLIP",RECOMU_SLIP,"RECOMU_SLIP[50]/F"); 
    Tree_->Branch("RECOMU_TIP",RECOMU_TIP,"RECOMU_TIP[50]/F"); 
    Tree_->Branch("RECOMU_LIP",RECOMU_LIP,"RECOMU_LIP[50]/F"); 
    Tree_->Branch("RECOMU_TIPERROR",RECOMU_TIPERROR,"RECOMU_TIPERROR[50]/F"); 
    Tree_->Branch("RECOMU_LIPERROR",RECOMU_LIPERROR,"RECOMU_LIPERROR[50]/F"); 
    
 

    Tree_->Branch("RECOMU_caloCompatibility",RECOMU_caloCompatibility,"RECOMU_caloCompatibility[50]/F");
    Tree_->Branch("RECOMU_segmentCompatibility",RECOMU_segmentCompatibility,"RECOMU_segmentCompatibility[50]/F"); 
    Tree_->Branch("RECOMU_numberOfMatches",RECOMU_numberOfMatches,"RECOMU_numberOfMatches[50]/i");
    Tree_->Branch("RECOMU_numberOfMatchedStations",RECOMU_numberOfMatchedStations,"RECOMU_numberOfMatchedStations[50]/i");
    Tree_->Branch("RECOMU_glbmuPromptTight",RECOMU_glbmuPromptTight,"RECOMU_glbmuPromptTight[50]/b");
    Tree_->Branch("RECOMU_chi2LocalPosition",RECOMU_chi2LocalPosition,"RECOMU_chi2LocalPosition[50]/F");
    Tree_->Branch("RECOMU_trkKink",RECOMU_trkKink,"RECOMU_trkKink[50]/F");
*/
    Tree_->Branch("RECOMU_isMedium",RECOMU_isMedium,"RECOMU_isMedium[50]/b");
 
    // track variables from muons:
/*
    Tree_->Branch( "RECOMU_trkmuArbitration", RECOMU_trkmuArbitration, "RECOMU_trkmuArbitration[50]/b");
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityLoose", RECOMU_trkmu2DCompatibilityLoose, "RECOMU_trkmu2DCompatibilityLoose[50]/b");
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityTight", RECOMU_trkmu2DCompatibilityTight, "RECOMU_trkmu2DCompatibilityTight[50]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationLoose", RECOMU_trkmuOneStationLoose, "RECOMU_trkmuOneStationLoose[50]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationTight", RECOMU_trkmuOneStationTight, "RECOMU_trkmuOneStationTight[50]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationLoose", RECOMU_trkmuLastStationLoose, "RECOMU_trkmuLastStationLoose[50]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationTight", RECOMU_trkmuLastStationTight, "RECOMU_trkmuLastStationTight[50]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationAngLoose", RECOMU_trkmuOneStationAngLoose, "RECOMU_trkmuOneStationAngLoose[50]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationAngTight", RECOMU_trkmuOneStationAngTight, "RECOMU_trkmuOneStationAngTight[50]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationAngLoose", RECOMU_trkmuLastStationAngLoose, "RECOMU_trkmuLastStationAngLoose[50]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationAngTight", RECOMU_trkmuLastStationAngTight, "RECOMU_trkmuLastStationAngTight[50]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtLoose",RECOMU_trkmuLastStationOptimizedLowPtLoose , "RECOMU_trkmuLastStationOptimizedLowPtLoose[50]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtTight",RECOMU_trkmuLastStationOptimizedLowPtTight , "RECOMU_trkmuLastStationOptimizedLowPtTight[50]/b");
*/
    Tree_->Branch( "RECOMU_mutrkPT", RECOMU_mutrkPT, "RECOMU_mutrkPT[50]/F");
    Tree_->Branch( "RECOMU_mutrkPTError", RECOMU_mutrkPTError, "RECOMU_mutrkPTError[50]/F");
    Tree_->Branch( "RECOMU_mutrkDxy", RECOMU_mutrkDxy, "RECOMU_mutrkDxy[50]/F");
    Tree_->Branch( "RECOMU_mutrkDxyError", RECOMU_mutrkDxyError, "RECOMU_mutrkDxyError[50]/F");
    Tree_->Branch( "RECOMU_mutrkDxyB", RECOMU_mutrkDxyB, "RECOMU_mutrkDxyB[50]/F");
    Tree_->Branch( "RECOMU_mutrkDz", RECOMU_mutrkDz, "RECOMU_mutrkDz[50]/F");
    Tree_->Branch( "RECOMU_mutrkDzError", RECOMU_mutrkDzError, "RECOMU_mutrkDzError[50]/F");
    Tree_->Branch( "RECOMU_mutrkDzB", RECOMU_mutrkDzB, "RECOMU_mutrkDzB[50]/F");
    Tree_->Branch( "RECOMU_mutrkChi2PerNdof", RECOMU_mutrkChi2PerNdof, "RECOMU_mutrkChi2PerNdof[50]/F");
    Tree_->Branch( "RECOMU_mutrkCharge", RECOMU_mutrkCharge, "RECOMU_mutrkCharge[50]/F");
    Tree_->Branch( "RECOMU_mutrkNHits", RECOMU_mutrkNHits, "RECOMU_mutrkNHits[50]/F");
    Tree_->Branch( "RECOMU_mutrkNStripHits", RECOMU_mutrkNStripHits, "RECOMU_mutrkNStripHits[50]/F");
    Tree_->Branch( "RECOMU_mutrkNPixHits", RECOMU_mutrkNPixHits, "RECOMU_mutrkNPixHits[50]/F");
    Tree_->Branch( "RECOMU_mutrkNMuonHits", RECOMU_mutrkNMuonHits, "RECOMU_mutrkNMuonHits[50]/F");
    Tree_->Branch( "RECOMU_mutrktrackerLayersWithMeasurement",RECOMU_mutrktrackerLayersWithMeasurement,"RECOMU_mutrktrackerLayersWithMeasurement[50]/F");
/*    
    Tree_->Branch( "RECOMU_muInnertrkDxy", RECOMU_muInnertrkDxy, "RECOMU_muInnertrkDxy[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkDxyError", RECOMU_muInnertrkDxyError, "RECOMU_muInnertrkDxyError[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkDxyB", RECOMU_muInnertrkDxyB, "RECOMU_muInnertrkDxyB[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkDz", RECOMU_muInnertrkDz, "RECOMU_muInnertrkDz[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkDzError", RECOMU_muInnertrkDzError, "RECOMU_muInnertrkDzError[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkDzB", RECOMU_muInnertrkDzB, "RECOMU_muInnertrkDzB[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkChi2PerNdof", RECOMU_muInnertrkChi2PerNdof, "RECOMU_muInnertrkChi2PerNdof[50]/F");
    Tree_->Branch( "RECOMU_muInnertrktrackerLayersWithMeasurement",RECOMU_muInnertrktrackerLayersWithMeasurement,"RECOMU_muInnertrktrackerLayersWithMeasurement[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkPT", RECOMU_muInnertrkPT, "RECOMU_muInnertrkPT[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkPTError", RECOMU_muInnertrkPTError, "RECOMU_muInnertrkPTError[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkCharge", RECOMU_muInnertrkCharge, "RECOMU_muInnertrkCharge[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkNHits", RECOMU_muInnertrkNHits, "RECOMU_muInnertrkNHits[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkNStripHits", RECOMU_muInnertrkNStripHits, "RECOMU_muInnertrkNStripHits[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkNPixHits", RECOMU_muInnertrkNPixHits, "RECOMU_muInnertrkNPixHits[50]/F");
    Tree_->Branch( "RECOMU_muInnertrkvalidFraction", RECOMU_muInnertrkvalidFraction, "RECOMU_muInnertrkvalidFraction[50]/F");
*/
    // best tracks for 13 TeV analysis
    Tree_->Branch( "RECOMU_mubesttrkType", RECOMU_mubesttrkType, "RECOMU_mubesttrkType[50]/I");
    Tree_->Branch( "RECOMU_mubesttrkDxy", RECOMU_mubesttrkDxy, "RECOMU_mubesttrkDxy[50]/F");
    Tree_->Branch( "RECOMU_mubesttrkDxyError", RECOMU_mubesttrkDxyError, "RECOMU_mubesttrkDxyError[50]/F");
    Tree_->Branch( "RECOMU_mubesttrkDz", RECOMU_mubesttrkDz, "RECOMU_mubesttrkDz[50]/F");
    Tree_->Branch( "RECOMU_mubesttrkDzError", RECOMU_mubesttrkDzError, "RECOMU_mubesttrkDzError[50]/F");
    Tree_->Branch( "RECOMU_mubesttrkPTError", RECOMU_mubesttrkPTError, "RECOMU_mubesttrkPTError[50]/F");

/*
    // Geom. Discri.
    Tree_->Branch("ftsigma",        &ftsigma,        "ftsigma[50]/D");
    Tree_->Branch("gdX",            &gdX,            "gdX[50]/D");
    Tree_->Branch("gdY",            &gdY,            "gdY[50]/D");
    Tree_->Branch("gdZ",            &gdZ,            "gdZ[50]/D");
    Tree_->Branch("ftsigmalag",     &ftsigmalag,     "ftsigmalag[50]/D");
    Tree_->Branch("gdlagX",         &gdlagX,         "gdlagX[50]/D");
    Tree_->Branch("gdlagY",         &gdlagY,         "gdlagY[50]/D");
    Tree_->Branch("gdlagZ",         &gdlagZ,         "gdlagZ[50]/D");
    Tree_->Branch("gdlagProb",      &gdlagProb,      "gdlagProb[50]/D");
    Tree_->Branch("gdlagNdof",      &gdlagNdof,      "gdlagNdof[50]/D");
    Tree_->Branch("ftsigmaMMMM",    &ftsigmaMMMM,    "ftsigmaMMMM[50]/D");
    Tree_->Branch("gdXMMMM",        &gdXMMMM,        "gdXMMMM[50]/D");
    Tree_->Branch("gdYMMMM",        &gdYMMMM,        "gdYMMMM[50]/D");
    Tree_->Branch("gdZMMMM",        &gdZMMMM,        "gdZMMMM[50]/D");
    Tree_->Branch("ftsigmalagMMMM", &ftsigmalagMMMM, "ftsigmalagMMMM[50]/D");
    Tree_->Branch("gdlagXMMMM",     &gdlagXMMMM,     "gdlagXMMMM[50]/D");
    Tree_->Branch("gdlagYMMMM",     &gdlagYMMMM,     "gdlagYMMMM[50]/D");
    Tree_->Branch("gdlagZMMMM",     &gdlagZMMMM,     "gdlagZMMMM[50]/D");
    Tree_->Branch("gdlagProbMMMM",  &gdlagProbMMMM,  "gdlagProbMMMM[50]/D");
    Tree_->Branch("gdlagNdofMMMM",  &gdlagNdofMMMM,  "gdlagNdofMMMM[50]/D");
    Tree_->Branch("ftsigmaEEEE",    &ftsigmaEEEE,    "ftsigmaEEEE[50]/D");
    Tree_->Branch("gdXEEEE",        &gdXEEEE,        "gdXEEEE[50]/D");
    Tree_->Branch("gdYEEEE",        &gdYEEEE,        "gdYEEEE[50]/D");
    Tree_->Branch("gdZEEEE",        &gdZEEEE,        "gdZEEEE[50]/D");
    Tree_->Branch("ftsigmalagEEEE", &ftsigmalagEEEE, "ftsigmalagEEEE[50]/D");
    Tree_->Branch("gdlagXEEEE",     &gdlagXEEEE,     "gdlagXEEEE[50]/D");
    Tree_->Branch("gdlagYEEEE",     &gdlagYEEEE,     "gdlagYEEEE[50]/D");
    Tree_->Branch("gdlagZEEEE",     &gdlagZEEEE,     "gdlagZEEEE[50]/D");
    Tree_->Branch("gdlagProbEEEE",  &gdlagProbEEEE,  "gdlagProbEEEE[50]/D");
    Tree_->Branch("gdlagNdofEEEE",  &gdlagNdofEEEE,  "gdlagNdofEEEE[50]/D");
    
    // ConstraintFit 4l
    Tree_->Branch("StdFitVertexX",        StdFitVertexX,        "StdFitVertexX[50]/D");
    Tree_->Branch("StdFitVertexY",        StdFitVertexY,        "StdFitVertexY[50]/D");
    Tree_->Branch("StdFitVertexZ",        StdFitVertexZ,        "StdFitVertexZ[50]/D");
    Tree_->Branch("StdFitVertexChi2r",    StdFitVertexChi2r,    "StdFitVertexChi2r[50]/D");
    Tree_->Branch("StdFitVertexProb",     StdFitVertexProb,     "StdFitVertexProb[50]/D");
    Tree_->Branch("StdFitVertexTrack_PT", StdFitVertexTrack_PT, "StdFitVertexTrack_PT[4][50]/F");
    Tree_->Branch("StdFitVertexTrack_ETA",StdFitVertexTrack_ETA,"StdFitVertexTrack_ETA[4][50]/F");
    Tree_->Branch("StdFitVertexTrack_PHI",StdFitVertexTrack_PHI,"StdFitVertexTrack_PHI[4][50]/F");
    Tree_->Branch("KinFitVertexX",        KinFitVertexX,        "KinFitVertexX[50]/D");
    Tree_->Branch("KinFitVertexY",        KinFitVertexY,        "KinFitVertexY[50]/D");
    Tree_->Branch("KinFitVertexZ",        KinFitVertexZ,        "KinFitVertexZ[50]/D");
    Tree_->Branch("KinFitVertexChi2r",    KinFitVertexChi2r,    "KinFitVertexChi2r[50]/D");
    Tree_->Branch("KinFitVertexProb",     KinFitVertexProb,     "KinFitVertexProb[50]/D");

    Tree_->Branch("StdFitVertexXMMMM",        StdFitVertexXMMMM,        "StdFitVertexXMMMM[50]/D");
    Tree_->Branch("StdFitVertexYMMMM",        StdFitVertexYMMMM,        "StdFitVertexYMMMM[50]/D");
    Tree_->Branch("StdFitVertexZMMMM",        StdFitVertexZMMMM,        "StdFitVertexZMMMM[50]/D");
    Tree_->Branch("StdFitVertexChi2rMMMM",    StdFitVertexChi2rMMMM,    "StdFitVertexChi2rMMMM[50]/D");
    Tree_->Branch("StdFitVertexProbMMMM",     StdFitVertexProbMMMM,     "StdFitVertexProbMMMM[50]/D");
    Tree_->Branch("StdFitVertexTrackMMMM_PT", StdFitVertexTrackMMMM_PT, "StdFitVertexTrackMMMM_PT[4][50]/F");
    Tree_->Branch("StdFitVertexTrackMMMM_ETA",StdFitVertexTrackMMMM_ETA,"StdFitVertexTrackMMMM_ETA[4][50]/F");
    Tree_->Branch("StdFitVertexTrackMMMM_PHI",StdFitVertexTrackMMMM_PHI,"StdFitVertexTrackMMMM_PHI[4][50]/F");
    Tree_->Branch("KinFitVertexXMMMM",        KinFitVertexXMMMM,        "KinFitVertexXMMMM[50]/D");
    Tree_->Branch("KinFitVertexYMMMM",        KinFitVertexYMMMM,        "KinFitVertexYMMMM[50]/D");
    Tree_->Branch("KinFitVertexZMMMM",        KinFitVertexZMMMM,        "KinFitVertexZMMMM[50]/D");
    Tree_->Branch("KinFitVertexChi2rMMMM",    KinFitVertexChi2rMMMM,    "KinFitVertexChi2rMMMM[50]/D");
    Tree_->Branch("KinFitVertexProbMMMM",     KinFitVertexProbMMMM,     "KinFitVertexProbMMMM[50]/D");
    

    Tree_->Branch("StdFitVertexXEEEE",        StdFitVertexXEEEE,        "StdFitVertexXEEEE[50]/D");
    Tree_->Branch("StdFitVertexYEEEE",        StdFitVertexYEEEE,        "StdFitVertexYEEEE[50]/D");
    Tree_->Branch("StdFitVertexZEEEE",        StdFitVertexZEEEE,        "StdFitVertexZEEEE[50]/D");
    Tree_->Branch("StdFitVertexChi2rEEEE",    StdFitVertexChi2rEEEE,    "StdFitVertexChi2rEEEE[50]/D");
    Tree_->Branch("StdFitVertexProbEEEE",     StdFitVertexProbEEEE,     "StdFitVertexProbEEEE[50]/D");
    Tree_->Branch("StdFitVertexTrackEEEE_PT", StdFitVertexTrackEEEE_PT, "StdFitVertexTrackEEEE_PT[4][50]/F");
    Tree_->Branch("StdFitVertexTrackEEEE_ETA",StdFitVertexTrackEEEE_ETA,"StdFitVertexTrackEEEE_ETA[4][50]/F");
    Tree_->Branch("StdFitVertexTrackEEEE_PHI",StdFitVertexTrackEEEE_PHI,"StdFitVertexTrackEEEE_PHI[4][50]/F");
    Tree_->Branch("KinFitVertexXEEEE",        KinFitVertexXEEEE,        "KinFitVertexXEEEE[50]/D");
    Tree_->Branch("KinFitVertexYEEEE",        KinFitVertexYEEEE,        "KinFitVertexYEEEE[50]/D");
    Tree_->Branch("KinFitVertexZEEEE",        KinFitVertexZEEEE,        "KinFitVertexZEEEE[50]/D");
    Tree_->Branch("KinFitVertexChi2rEEEE",    KinFitVertexChi2rEEEE,    "KinFitVertexChi2rEEEE[50]/D");
    Tree_->Branch("KinFitVertexProbEEEE",     KinFitVertexProbEEEE,     "KinFitVertexProbEEEE[50]/D");

    // constrintFit 3l
    Tree_->Branch("StdFitVertexChi2rMMM",     StdFitVertexChi2rMMM,    "StdFitVertexChi2rMMM[50]/D");
    Tree_->Branch("StdFitVertexProbMMM",      StdFitVertexProbMMM,     "StdFitVertexProbMMM[50]/D");
    Tree_->Branch("StdFitVertexChi2rMME",     StdFitVertexChi2rMME,    "StdFitVertexChi2rMME[50]/D");
    Tree_->Branch("StdFitVertexProbMME",      StdFitVertexProbMME,     "StdFitVertexProbMME[50]/D");
    Tree_->Branch("StdFitVertexChi2rEEE",     StdFitVertexChi2rEEE,    "StdFitVertexChi2rEEE[50]/D");
    Tree_->Branch("StdFitVertexProbEEE",      StdFitVertexProbEEE,     "StdFitVertexProbEEE[50]/D");
    Tree_->Branch("StdFitVertexChi2rMEE",     StdFitVertexChi2rMEE,    "StdFitVertexChi2rMEE[50]/D");
    Tree_->Branch("StdFitVertexProbMEE",      StdFitVertexProbMEE,     "StdFitVertexProbMEE[50]/D");


     // constrintFit Dileptons
    Tree_->Branch("StdFitVertexChi2rDiLep",   StdFitVertexChi2rDiLep,    "StdFitVertexChi2rDiLep[40]/D");
    Tree_->Branch("StdFitVertexProbDiLep",    StdFitVertexProbDiLep,     "StdFitVertexProbDiLep[40]/D");

    // Conversions
    Tree_->Branch("ConvMapDist",              ConvMapDist,              "ConvMapDist[50]/F");
    Tree_->Branch("ConvMapDcot",              ConvMapDcot,              "ConvMapDcot[50]/F");
*/


    //MatchingMC:
    //Muons:
    Tree_->Branch("RECOMU_MatchingMCTruth", RECOMU_MatchingMCTruth, "RECOMU_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOMU_MatchingMCpT", RECOMU_MatchingMCpT, "RECOMU_MatchingMCpT[50]/F");
    Tree_->Branch("RECOMU_MatchingMCEta", RECOMU_MatchingMCEta, "RECOMU_MatchingMCEta[50]/F");
    Tree_->Branch("RECOMU_MatchingMCPhi", RECOMU_MatchingMCPhi, "RECOMU_MatchingMCPhi[50]/F");

    //Electrons:
    Tree_->Branch("RECOELE_MatchingMCTruth", RECOELE_MatchingMCTruth, "RECOELE_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOELE_MatchingMCpT", RECOELE_MatchingMCpT, "RECOELE_MatchingMCpT[50]/F");
    Tree_->Branch("RECOELE_MatchingMCEta", RECOELE_MatchingMCEta, "RECOELE_MatchingMCEta[50]/F");
    Tree_->Branch("RECOELE_MatchingMCPhi", RECOELE_MatchingMCPhi, "RECOELE_MatchingMCPhi[50]/F");

    //Bottom
    Tree_->Branch("RECOBOT_MatchingMCTruth", RECOBOT_MatchingMCTruth, "RECOBOT_MatchingMCTruth[50]/I");
    Tree_->Branch("RECOBOT_MatchingMCpT", RECOBOT_MatchingMCpT, "RECOBOT_MatchingMCpT[50]/F");
    Tree_->Branch("RECOBOT_MatchingMCEta", RECOBOT_MatchingMCEta, "RECOBOT_MatchingMCEta[50]/F");
    Tree_->Branch("RECOBOT_MatchingMCPhi", RECOBOT_MatchingMCPhi, "RECOBOT_MatchingMCPhi[50]/F");

    //Gamma:
    Tree_->Branch("RECOPHOT_MatchingMCTruth", RECOPHOT_MatchingMCTruth, "RECOPHOT_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOPHOT_MatchingMCpT", RECOPHOT_MatchingMCpT, "RECOPHOT_MatchingMCpT[50]/F");
    Tree_->Branch("RECOPHOT_MatchingMCEta", RECOPHOT_MatchingMCEta, "RECOPHOT_MatchingMCEta[50]/F");
    Tree_->Branch("RECOPHOT_MatchingMCPhi", RECOPHOT_MatchingMCPhi, "RECOPHOT_MatchingMCPhi[50]/F");

/*
    //ZtoMuMu:
    Tree_->Branch("RECOzMuMu_MatchingMCTruth", RECOzMuMu_MatchingMCTruth, "RECOzMuMu_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOzMuMu_MatchingMCpT", RECOzMuMu_MatchingMCpT, "RECOzMuMu_MatchingMCpT[50]/F");
    Tree_->Branch("RECOzMuMu_MatchingMCmass", RECOzMuMu_MatchingMCmass, "RECOzMuMu_MatchingMCmass[50]/F");
    Tree_->Branch("RECOzMuMu_MatchingMCEta", RECOzMuMu_MatchingMCEta, "RECOzMuMu_MatchingMCEta[50]/F");
    Tree_->Branch("RECOzMuMu_MatchingMCPhi", RECOzMuMu_MatchingMCPhi, "RECOzMuMu_MatchingMCPhi[50]/F");

    //ZtoEE:
    Tree_->Branch("RECOzEE_MatchingMCTruth", RECOzEE_MatchingMCTruth, "RECOzEE_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOzEE_MatchingMCpT", RECOzEE_MatchingMCpT, "RECOzEE_MatchingMCpT[50]/F");
    Tree_->Branch("RECOzEE_MatchingMCmass", RECOzEE_MatchingMCmass, "RECOzEE_MatchingMCmass[50]/F");
    Tree_->Branch("RECOzEE_MatchingMCEta", RECOzEE_MatchingMCEta, "RECOzEE_MatchingMCEta[50]/F");
    Tree_->Branch("RECOzEE_MatchingMCPhi", RECOzEE_MatchingMCPhi, "RECOzEE_MatchingMCPhi[50]/F");

    //HtoZtoEEEE:
    Tree_->Branch("RECOHzzEEEE_MatchingMCTruth", RECOHzzEEEE_MatchingMCTruth, "RECOHzzEEEE_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOHzzEEEE_MatchingMCpT", RECOHzzEEEE_MatchingMCpT, "RECOHzzEEEE_MatchingMCpT[50]/F");
    Tree_->Branch("RECOHzzEEEE_MatchingMCmass", RECOHzzEEEE_MatchingMCmass, "RECOHzzEEEE_MatchingMCmass[50]/F");
    Tree_->Branch("RECOHzzEEEE_MatchingMCEta", RECOHzzEEEE_MatchingMCEta, "RECOHzzEEEE_MatchingMCEta[50]/F");
    Tree_->Branch("RECOHzzEEEE_MatchingMCPhi", RECOHzzEEEE_MatchingMCPhi, "RECOHzzEEEE_MatchingMCPhi[50]/F");

    //HtoZtoEEMM:
    Tree_->Branch("RECOHzzEEMM_MatchingMCTruth", RECOHzzEEMM_MatchingMCTruth, "RECOHzzEEMM_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOHzzEEMM_MatchingMCpT", RECOHzzEEMM_MatchingMCpT, "RECOHzzEEMM_MatchingMCpT[50]/F");
    Tree_->Branch("RECOHzzEEMM_MatchingMCmass", RECOHzzEEMM_MatchingMCmass, "RECOHzzEEMM_MatchingMCmass[50]/F");
    Tree_->Branch("RECOHzzEEMM_MatchingMCEta", RECOHzzEEMM_MatchingMCEta, "RECOHzzEEMM_MatchingMCEta[50]/F");
    Tree_->Branch("RECOHzzEEMM_MatchingMCPhi", RECOHzzEEMM_MatchingMCPhi, "RECOHzzEEMM_MatchingMCPhi[50]/F");

    //HtoZtoMMMM:
    Tree_->Branch("RECOHzzMMMM_MatchingMCTruth", RECOHzzMMMM_MatchingMCTruth, "RECOHzzMMMM_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOHzzMMMM_MatchingMCpT", RECOHzzMMMM_MatchingMCpT, "RECOHzzMMMM_MatchingMCpT[50]/F");
    Tree_->Branch("RECOHzzMMMM_MatchingMCmass", RECOHzzMMMM_MatchingMCmass, "RECOHzzMMMM_MatchingMCmass[50]/F");
    Tree_->Branch("RECOHzzMMMM_MatchingMCEta", RECOHzzMMMM_MatchingMCEta, "RECOHzzMMMM_MatchingMCEta[50]/F");
    Tree_->Branch("RECOHzzMMMM_MatchingMCPhi", RECOHzzMMMM_MatchingMCPhi, "RECOHzzMMMM_MatchingMCPhi[50]/F");
*/


    //Global Event 
    Tree_->Branch( "RECO_NMU", &RECO_NMU, "RECO_NMU/I"); 
    Tree_->Branch( "RECO_NELE", &RECO_NELE, "RECO_NELE/I"); 
    
    // Tracks
    Tree_->Branch( "RECO_NTRACK", &RECO_NTRACK, "RECO_NTRACK/I");
    Tree_->Branch( "RECO_TRACK_PT", &RECO_TRACK_PT, "RECO_TRACK_PT[100]/F");
    Tree_->Branch( "RECO_TRACK_ETA", &RECO_TRACK_ETA, "RECO_TRACK_ETA[100]/F");
    Tree_->Branch( "RECO_TRACK_PHI", &RECO_TRACK_PHI, "RECO_TRACK_PHI[100]/F");
    Tree_->Branch( "RECO_TRACK_CHI2", &RECO_TRACK_CHI2, "RECO_TRACK_CHI2[100]/F");
    Tree_->Branch( "RECO_TRACK_CHI2RED", &RECO_TRACK_CHI2RED, "RECO_TRACK_CHI2RED[100]/F");
    Tree_->Branch( "RECO_TRACK_CHI2PROB", &RECO_TRACK_CHI2PROB, "RECO_TRACK_CHI2PROB[100]/F");
    Tree_->Branch( "RECO_TRACK_NHITS", &RECO_TRACK_NHITS, "RECO_TRACK_NHITS[100]/I");
    Tree_->Branch( "RECO_TRACK_DXY", &RECO_TRACK_DXY, "RECO_TRACK_DXY[100]/F");
    Tree_->Branch( "RECO_TRACK_DXYERR", &RECO_TRACK_DXYERR, "RECO_TRACK_DXYERR[100]/F");
    Tree_->Branch( "RECO_TRACK_DZ", &RECO_TRACK_DZ, "RECO_TRACK_DZ[100]/F");
    Tree_->Branch( "RECO_TRACK_DZERR", &RECO_TRACK_DZERR, "RECO_TRACK_DZERR[100]/F");
    
    // Photons
    Tree_->Branch("RECO_NPHOT", &RECO_NPHOT, "RECO_NPHOT/I");
    Tree_->Branch("RECOPHOT_PT",RECOPHOT_PT,"RECOPHOT_PT[20]/F"); 
    Tree_->Branch("RECOPHOT_ETA",RECOPHOT_ETA,"RECOPHOT_ETA[20]/F"); 
    Tree_->Branch("RECOPHOT_PHI",RECOPHOT_PHI,"RECOPHOT_PHI[20]/F"); 
    Tree_->Branch("RECOPHOT_THETA",RECOPHOT_THETA,"RECOPHOT_THETA[20]/F"); 
    Tree_->Branch("RECOPHOT_TLE_ParentSC_X",RECOPHOT_TLE_ParentSC_X,"RECOPHOT_TLE_ParentSC_X[20]/F");
    Tree_->Branch("RECOPHOT_TLE_ParentSC_Y",RECOPHOT_TLE_ParentSC_Y,"RECOPHOT_TLE_ParentSC_Y[20]/F");
    Tree_->Branch("RECOPHOT_TLE_ParentSC_Z",RECOPHOT_TLE_ParentSC_Z,"RECOPHOT_TLE_ParentSC_X[20]/F");

    Tree_->Branch("RECO_NPFPHOT", &RECO_NPFPHOT, "RECO_NPFPHOT/I");
    Tree_->Branch("RECOPFPHOT_PT",RECOPFPHOT_PT,"RECOPFPHOT_PT[20]/F"); 
    Tree_->Branch("RECOPFPHOT_PTError",RECOPFPHOT_PTError,"RECOPFPHOT_PTError[20]/F");  
    Tree_->Branch("RECOPFPHOT_ETA",RECOPFPHOT_ETA,"RECOPFPHOT_ETA[20]/F"); 
    Tree_->Branch("RECOPFPHOT_PHI",RECOPFPHOT_PHI,"RECOPFPHOT_PHI[20]/F"); 
    Tree_->Branch("RECOPFPHOT_THETA",RECOPFPHOT_THETA,"RECOPFPHOT_THETA[20]/F"); 
    
    
    //Beam Spot position
    Tree_->Branch("BeamSpot_X",&BeamSpot_X,"BeamSpot_X/D");
    Tree_->Branch("BeamSpot_Y",&BeamSpot_Y,"BeamSpot_Y/D");
    Tree_->Branch("BeamSpot_Z",&BeamSpot_Z,"BeamSpot_Z/D");
    // Vertices
    Tree_->Branch( "RECO_NVTX", &RECO_NVTX, "RECO_NVTX/I");
    Tree_->Branch( "RECO_VERTEX_x", RECO_VERTEX_x, "RECO_VERTEX_x[15]/F");
    Tree_->Branch( "RECO_VERTEX_y", RECO_VERTEX_y, "RECO_VERTEX_y[15]/F");
    Tree_->Branch( "RECO_VERTEX_z", RECO_VERTEX_z, "RECO_VERTEX_z[15]/F");
    Tree_->Branch( "RECO_VERTEX_ndof", RECO_VERTEX_ndof, "RECO_VERTEX_ndof[15]/F");
    Tree_->Branch( "RECO_VERTEX_chi2", RECO_VERTEX_chi2, "RECO_VERTEX_chi2[15]/F");
    Tree_->Branch( "RECO_VERTEX_ntracks", RECO_VERTEX_ntracks, "RECO_VERTEX_ntracks[15]/I");
    Tree_->Branch( "RECO_VERTEXPROB", RECO_VERTEXPROB, "RECO_VERTEXPROB[15]/F");
    Tree_->Branch( "RECO_VERTEX_isValid", RECO_VERTEX_isValid, "RECO_VERTEX_isValid[15]/b");
    Tree_->Branch( "RECO_VERTEX_TRACK_PT",RECO_VERTEX_TRACK_PT,"RECO_VERTEX_TRACK_PT[15][50]/F");
    
    // PFJets
    Tree_->Branch( "RECO_PFJET_N",   &RECO_PFJET_N,   "RECO_PFJET_N/I");
    Tree_->Branch( "RECO_PFJET_CHARGE",  RECO_PFJET_CHARGE,  "RECO_PFJET_CHARGE[50]/I");
    Tree_->Branch( "RECO_PFJET_ET",  RECO_PFJET_ET,  "RECO_PFJET_ET[50]/F");
    Tree_->Branch( "RECO_PFJET_PT",  RECO_PFJET_PT,  "RECO_PFJET_PT[50]/F");
    Tree_->Branch( "RECO_PFJET_PT_UP",  RECO_PFJET_PT_UP,  "RECO_PFJET_PT_UP[50]/F");
    Tree_->Branch( "RECO_PFJET_PT_DOW",  RECO_PFJET_PT_DOW,  "RECO_PFJET_PT_DOW[50]/F");
    Tree_->Branch( "RECO_PFJET_ETA", RECO_PFJET_ETA, "RECO_PFJET_ETA[50]/F");
    Tree_->Branch( "RECO_PFJET_PHI", RECO_PFJET_PHI, "RECO_PFJET_PHI[50]/F");
    Tree_->Branch( "RECO_PFJET_PUID", RECO_PFJET_PUID, "RECO_PFJET_PUID[50]/I");
    Tree_->Branch( "RECO_PFJET_PUID_MVA", RECO_PFJET_PUID_MVA, "RECO_PFJET_PUID_MVA[50]/F");
    //jet ID
    Tree_->Branch( "RECO_PFJET_nconstituents",  RECO_PFJET_nconstituents,  "RECO_PFJET_nconstituents[50]/I");
    Tree_->Branch( "RECO_PFJET_NCH",  RECO_PFJET_NCH,  "RECO_PFJET_NCH[50]/I");
    Tree_->Branch( "RECO_PFJET_NHF",  RECO_PFJET_NHF,  "RECO_PFJET_NHF[50]/F");
    Tree_->Branch( "RECO_PFJET_NEF",  RECO_PFJET_NEF,  "RECO_PFJET_NEF[50]/F");
    Tree_->Branch( "RECO_PFJET_CHF", RECO_PFJET_CHF, "RECO_PFJET_CHF[50]/F");
    Tree_->Branch( "RECO_PFJET_CEF", RECO_PFJET_CEF, "RECO_PFJET_CEF[50]/F");
    Tree_->Branch( "RECO_PFJET_MUF", RECO_PFJET_MUF, "RECO_PFJET_MUF[50]/F");
    Tree_->Branch( "RHO_ele", &RHO_ele, "RHO_ele/D");
    Tree_->Branch( "RHO_mu", &RHO_mu, "RHO_mu/D");
    
    //CaloMET
    Tree_->Branch( "RECO_CALOMET",          &calomet,          "RECO_CALOMET/F");
 /*    Tree_->Branch( "RECO_CALOMETHO",        &calometho,        "RECO_CALOMETHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHFHO",    &calometnohfho,    "RECO_CALOMETNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHF",      &calometnohf,      "RECO_CALOMETNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTHO",     &calometoptho,     "RECO_CALOMETOPTHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHFHO", &calometoptnohfho, "RECO_CALOMETOPTNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHF",   &calometoptnohf,   "RECO_CALOMETOPTNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPT",       &calometopt,       "RECO_CALOMETOPT/F"); */
    //Particle Flow MET
    Tree_->Branch( "RECO_PFMET", &pfmet, "RECO_PFMET/F");
    Tree_->Branch( "RECO_PFMET_X", &pfmet_x, "RECO_PFMET_X/F");
    Tree_->Branch( "RECO_PFMET_Y", &pfmet_y, "RECO_PFMET_Y/F");
    Tree_->Branch( "RECO_PFMET_PHI", &pfmet_phi, "RECO_PFMET_PHI/F");
    Tree_->Branch( "RECO_PFMET_THETA", &pfmet_theta, "RECO_PFMET_THETA/F");
    //Track Corrected MET
    Tree_->Branch( "RECO_TCMET", &tcmet, "RECO_TCMET/F");
    //Type I correction MET
    Tree_->Branch( "RECO_CORMETMUONS",  &cormetmuons,  "RECO_CORMETMUONS/F");
   

    // Btagging jets and discriminators
    Tree_->Branch("tCHighEff_BTagJet_PT",tCHighEff_BTagJet_PT,"tCHighEff_BTagJet_PT[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_PT", tCHighPur_BTagJet_PT,"tCHighPur_BTagJet_PT[50]/F");
    Tree_->Branch("cSV_BTagJet_PT",cSV_BTagJet_PT,"cSV_BTagJet_PT[50]/F");
    Tree_->Branch("tCHighEff_BTagJet_ETA",tCHighEff_BTagJet_ETA,"tCHighEff_BTagJet_ETA[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_ETA", tCHighPur_BTagJet_ETA,"tCHighPur_BTagJet_ETA[50]/F");
    Tree_->Branch("cSV_BTagJet_ETA",cSV_BTagJet_ETA,"cSV_BTagJet_ETA[50]/F");
    Tree_->Branch("tCHighEff_BTagJet_PHI",tCHighEff_BTagJet_PHI,"tCHighEff_BTagJet_PHI[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_PHI", tCHighPur_BTagJet_PHI,"tCHighPur_BTagJet_PHI[50]/F");
    Tree_->Branch("cSV_BTagJet_PHI",cSV_BTagJet_PHI,"cSV_BTagJet_PHI[50]/F");
    Tree_->Branch("cSV_BTagJet_ET",cSV_BTagJet_ET,"cSV_BTagJet_ET[50]/F");
    Tree_->Branch("tCHighEff_BTagJet_DISCR",tCHighEff_BTagJet_DISCR,"tCHighEff_BTagJet_DISCR[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_DISCR", tCHighPur_BTagJet_DISCR,"tCHighPur_BTagJet_DISCR[50]/F");
    Tree_->Branch("cSV_BTagJet_DISCR",cSV_BTagJet_DISCR,"cSV_BTagJet_DISCR[50]/F");     
  }
  
  
  void Initialize(){
    
    irun=-999,ievt=-999,ils=-999;
    Avginstlumi=-999.;
    RHO=-999.,RHO_ele=-999.,RHO_mu=-999.;

    // PU
    num_PU_vertices=-999;
    PU_BunchCrossing=-999;

    
    // HLT flags   

    RECO_nMuHLTMatch=0;
    RECO_nEleHLTMatch=0;
    
    //for (int ii=0;ii<200;ii++){
    //  HLTPathsFired[ii]="";
    //}
/*
    leptonscands2e2mu_= new (CandidateCollection);
    leptonscands2e2murf_= new (CandidateCollection);
    leptonscands4mu_= new (CandidateCollection);
    leptonscands4murf_= new (CandidateCollection);
    leptonscands4e_= new (CandidateCollection);
    leptonscands4erf_= new (CandidateCollection);

    leptonscands_Z0= new (CandidateCollection);
    leptonscands_Z1= new (CandidateCollection);
    leptonscands_Zss0= new (CandidateCollection);
    leptonscands_Zss1= new (CandidateCollection);
    leptonscands_Zcross= new (CandidateCollection);
    leptonscands_DiLep= new (CandidateCollection);
    leptonscands_MMMM= new (CandidateCollection);
    leptonscands_EEEE= new (CandidateCollection);
    leptonscands_EEMM= new (CandidateCollection);
    leptonscands_LLL0= new (CandidateCollection);
    leptonscands_LLL1= new (CandidateCollection);
    leptonscands_LLL2= new (CandidateCollection);
    leptonscands_LLL3= new (CandidateCollection);
    leptonscands_LLLLss0= new (CandidateCollection);
    leptonscands_LLLLss1= new (CandidateCollection);
    leptonscands_LLLLss2= new (CandidateCollection);
    leptonscands_LLLl0= new (CandidateCollection);
    leptonscands_LLLl1= new (CandidateCollection);
    leptonscands_LLLL= new (CandidateCollection);
*/
    for  (int i=0; i<4;i++){ 
      MC_LEPT_PT[i]=-999.;
      MC_LEPT_ETA[i]=-999.;
      MC_LEPT_PHI[i]=-999.;
      MC_LEPT_THETA[i]=-999.;
      MC_LEPT_PDGID[i]=-999.;
    }

    for(int i=0; i<2; ++i) {
      for(int j=0; j<5; ++j) {
	MC_Z_PT[i][j]=-999;
	MC_Z_ETA[i][j]=-999;
	MC_Z_PHI[i][j]=-999;
	MC_Z_THETA[i][j]=-999;
	MC_Z_MASS[i][j]=-999;
	MC_Z_PDGID[i][j]=-999;
      }
    }

    for (int i=0; i<7;i++){
      MC_E[i]=-999.;
      MC_PT[i]=-999.;
      MC_ETA[i]=-999.;
      MC_THETA[i]=-999.;
      MC_PHI[i]=-999.;
      MC_MASS[i]=-999.;
      MC_PDGID[i]=-999.;
      


  

    }

    genmet=-999.,calomet=-999.;  
    RECO_NMU=0,RECO_NELE=0;
    RECO_NTRACK=0;
    
    
    RECO_NPHOT=0,RECO_NJET=0,RECO_NVTX=0;
    RECO_NPFPHOT=0;

    for (int i=0; i<20;i++){
      RECOPHOT_PT[i]=-999.;
      RECOPHOT_ETA[i]=-999.;
      RECOPHOT_PHI[i]=-999.;
      RECOPHOT_THETA[i]=-999.;
      RECOPHOT_TLE_ParentSC_X[i]=-999.;
      RECOPHOT_TLE_ParentSC_Y[i]=-999.;
      RECOPHOT_TLE_ParentSC_Z[i]=-999.;
      RECOPFPHOT_PT[i]=-999.;
      RECOPFPHOT_PTError[i]=-999.;
      RECOPFPHOT_ETA[i]=-999.;
      RECOPFPHOT_PHI[i]=-999.;
      RECOPFPHOT_THETA[i]=-999.;

      RECOPFPHOT_PFX_rho[i]=-999.;

      RECOPFPHOT_PFchAllPart[i]=-999.;
      RECOPFPHOT_PFchHad[i]=-999.;
      RECOPFPHOT_PFneuHad[i]=999.;
      RECOPFPHOT_PFphoton[i]=-999.;
      RECOPFPHOT_PFPUchAllPart[i]=-999.;
      RECOPFPHOT_PFX_rho[i]=-999.;


    }
		
    for(int ivtx=0;ivtx<15;ivtx++) {
      RECO_VERTEX_x[ivtx] = -999;
      RECO_VERTEX_y[ivtx] = -999;
      RECO_VERTEX_z[ivtx] = -999;
      RECO_VERTEXPROB[ivtx]=-999.;
      RECO_VERTEX_ndof[ivtx] = -999;
      RECO_VERTEX_chi2[ivtx] = -999;
      RECO_VERTEX_ntracks[ivtx]=-999;
      for (int i=0;i<50;i++){
	RECO_VERTEX_TRACK_PT[ivtx][i]=-999;
      }
    }
    
    RECO_PFJET_N = 0;
    for (int ijets=0;ijets<50;ijets++) {
      RECO_PFJET_CHARGE[ijets]   = -999.; 
      RECO_PFJET_ET[ijets]   = -999.; 
      RECO_PFJET_PT[ijets]  = -999.;
      RECO_PFJET_PT_UP[ijets] = -999.;
      RECO_PFJET_PT_DOW[ijets]= -999; 
      RECO_PFJET_ETA[ijets] = -999.; 
      RECO_PFJET_PHI[ijets] = -999.;
      RECO_PFJET_PUID[ijets] = -999;
      RECO_PFJET_PUID_MVA[ijets] = -999.;
      RECO_PFJET_NHF[ijets] = -999.;
      RECO_PFJET_NEF[ijets]     = -999.;
      RECO_PFJET_CHF[ijets]     = -999.;
      RECO_PFJET_CEF[ijets]    = -999.;
      RECO_PFJET_nconstituents[ijets]    = -999;
      RECO_PFJET_NCH[ijets] = -999;
      RECO_PFJET_MUF[ijets] = -999;
    }
    
/*     calomet=-999.; */
/*     calometopt=-999.; */
/*     calometoptnohf=-999.; */
/*     calometoptnohfho=-999.; */
/*     calometoptho=-999.; */
/*     calometnohf=-999.; */
/*     calometnohfho=-999.; */
/*     calometho=-999.; */
    pfmet=-999.;
    pfmet_x=-999.;
    pfmet_y=-999.;
    pfmet_phi=-999.;
    pfmet_theta=-999.;
    tcmet=-999.;
    cormetmuons=-999.;
    
    BeamSpot_X=-999.;
    BeamSpot_Y=-999.;
    BeamSpot_Z=-999.;
    
    
    for (int i=0; i<50;i++){
      MC_GENJET_PT[i]=-999;
      MC_GENJET_ETA[i]=-999;
      MC_GENJET_PHI[i]=-999;
/*
      ftsigma[i]=-999.;
      ftsigmalag[i]=-999.;
      gdX[i]=-999.;
      gdY[i]=-999.;
      gdZ[i]=-999.;
      gdlagX[i]=-999.;
      gdlagY[i]=-999.;
      gdlagZ[i]=-999.;
      gdlagProb[i]=-999.;
      gdlagNdof[i]=-999.;
      ftsigmaMMMM[i]=-999.;
      ftsigmalagMMMM[i]=-999.;
      gdXMMMM[i]=-999.;
      gdYMMMM[i]=-999.;
      gdZMMMM[i]=-999.;
      gdlagXMMMM[i]=-999.;
      gdlagYMMMM[i]=-999.;
      gdlagZMMMM[i]=-999.;
      gdlagProbMMMM[i]=-999.;
      gdlagNdofMMMM[i]=-999.;
      ftsigmaEEEE[i]=-999.;
      ftsigmalagEEEE[i]=-999.;
      gdXEEEE[i]=-999.;
      gdYEEEE[i]=-999.;
      gdZEEEE[i]=-999.;
      gdlagXEEEE[i]=-999.;
      gdlagYEEEE[i]=-999.;
      gdlagZEEEE[i]=-999.;
      gdlagProbEEEE[i]=-999.;
      gdlagNdofEEEE[i]=-999.;

      StdFitVertexX[i]=-999.;
      StdFitVertexY[i]=-999.;
      StdFitVertexZ[i]=-999.;
      StdFitVertexChi2r[i]=-999.;
      StdFitVertexProb[i]=-999.;
      KinFitVertexX[i]=-999.;
      KinFitVertexY[i]=-999.;
      KinFitVertexZ[i]=-999.;
      KinFitVertexChi2r[i]=-999.;
      KinFitVertexProb[i]=-999.;

      StdFitVertexXMMMM[i]=-999.;
      StdFitVertexYMMMM[i]=-999.;
      StdFitVertexZMMMM[i]=-999.;
      StdFitVertexChi2rMMMM[i]=-999.;
      StdFitVertexProbMMMM[i]=-999.;
      KinFitVertexXMMMM[i]=-999.;
      KinFitVertexYMMMM[i]=-999.;
      KinFitVertexZMMMM[i]=-999.;
      KinFitVertexChi2rMMMM[i]=-999.;
      KinFitVertexProbMMMM[i]=-999.;

      StdFitVertexXEEEE[i]=-999.;
      StdFitVertexYEEEE[i]=-999.;
      StdFitVertexZEEEE[i]=-999.;
      StdFitVertexChi2rEEEE[i]=-999.;
      StdFitVertexProbEEEE[i]=-999.;
      KinFitVertexXEEEE[i]=-999.;
      KinFitVertexYEEEE[i]=-999.;
      KinFitVertexZEEEE[i]=-999.;
      KinFitVertexChi2rEEEE[i]=-999.;
      KinFitVertexProbEEEE[i]=-999.;

      StdFitVertexChi2rMMM[i]=-999.;
      StdFitVertexProbMMM[i]=-999.;
      StdFitVertexChi2rMME[i]=-999.;
      StdFitVertexProbMME[i]=-999.;
      StdFitVertexChi2rEEE[i]=-999.;
      StdFitVertexProbEEE[i]=-999.;
      StdFitVertexChi2rMEE[i]=-999.;
      StdFitVertexProbMEE[i]=-999.;

      RECO_MMMM_MASS_REFIT[i]=-999.;
      RECO_EEMM_MASS_REFIT[i]=-999.;
      RECO_EEEE_MASS_REFIT[i]=-999.;

      for (int k=0;k<4;k++){
	StdFitVertexTrack_PT[k][i]=-999.;
	StdFitVertexTrack_ETA[k][i]=-999.;
	StdFitVertexTrack_PHI[k][i]=-999.;
	
	StdFitVertexTrackMMMM_PT[k][i]=-999.;
	StdFitVertexTrackMMMM_ETA[k][i]=-999.;
	StdFitVertexTrackMMMM_PHI[k][i]=-999.;
	
	StdFitVertexTrackEEEE_PT[k][i]=-999.;
	StdFitVertexTrackEEEE_ETA[k][i]=-999.;
	StdFitVertexTrackEEEE_PHI[k][i]=-999.;
      }
*/
    }
    
    dm_trig=false;
    sm_trig=false;
    de_trig=false;
    se_trig=false;
    tri_trig=false;
    jet_trig=0;   
 
    for (int i=0; i<50;i++){
      RECOELE_E[i]     = -999.;
      RECOELE_PT[i]=-999.;
      RECOELE_PTError[i]=-999.;
      RECOELE_P[i]=-999.;
      RECOELE_PHI[i]=-999.;
      RECOELE_ETA[i]=-999.;
      RECOELE_THETA[i]=-999.;
      RECOELE_MASS[i]=-999.;
      RECOELE_CHARGE[i]=-999.;

      RECOELE_ID[i]=-999.;
      RECOELE_PT_uncorr[i]=-999.;
  
      //

      RECOELE_ecalTrkEnergyPreCorr[i] =-999.;
      RECOELE_ecalTrkEnergyErrPreCorr[i]=-999.;
      RECOELE_ecalTrkEnergyErrPostCorr[i]=-999.;
      RECOELE_energyScaleValue[i]=-999.;       
      RECOELE_energySigmaValue[i]=-999.;
      RECOELE_energyScaleUp[i]=-999.;    
      RECOELE_energyScaleDown[i]=-999.;       
      RECOELE_energyScaleStatUp[i]=-999.;       
      RECOELE_energyScaleStatDown[i]=-999.;        
      RECOELE_energyScaleSystUp[i]=-999.;        
      RECOELE_energyScaleSystDown[i]=-999.;        
      RECOELE_energyScaleGainUp[i]=-999.;        
      RECOELE_energyScaleGainDown[i]=-999.;      
      RECOELE_energyScaleEtUp[i]=-999.;       
      RECOELE_energyScaleEtDown[i]=-999.;       
      RECOELE_energySigmaUp[i]=-999.;         
      RECOELE_energySigmaDown[i]=-999.;       
      RECOELE_energySigmaPhiUp[i]=-999.;        
      RECOELE_energySigmaPhiDown[i]=-999.;     
      RECOELE_energySigmaRhoUp[i]=-999.;        
      RECOELE_energySigmaRhoDown[i]=-999.;  
      
      // Core attributes
      RECOELE_isEcalDriven[i]    = false;
      RECOELE_isTrackerDriven[i] = false;
      RECOELE_gsftrack_NPixHits[i]=-999.;
      RECOELE_gsftrack_NStripHits[i]=-999.;
      RECOELE_gsftrack_chi2[i]   = -999;
      RECOELE_gsftrack_dxyB[i]   = -999;
      RECOELE_gsftrack_dxy[i]    = -999;
      RECOELE_gsftrack_dxyError[i]    = -999;
      RECOELE_gsftrack_dzB[i]    = -999;
      RECOELE_gsftrack_dz[i]     = -999;
      RECOELE_gsftrack_dzError[i]   = -999;
      RECOELE_scl_E[i]   = -999;
      RECOELE_scl_Et[i]  = -999;
      RECOELE_scl_Eta[i] = -999;
      RECOELE_scl_Phi[i] = -999;
      RECOELE_ecalEnergy[i] = -999;

      // Track-Cluster matching attributes
      RECOELE_ep[i]             = -999.;
      RECOELE_eSeedp[i]         = -999.;
      RECOELE_eSeedpout[i]      = -999.;
      RECOELE_eElepout[i]       = -999.;
      //
      RECOELE_deltaEtaIn[i]     = -999.;
      RECOELE_deltaEtaSeed[i]   = -999.;
      RECOELE_deltaEtaEle[i]    = -999.;
      RECOELE_deltaPhiIn[i]     = -999.;
      RECOELE_deltaPhiSeed[i]   = -999.;
      RECOELE_deltaPhiEle[i]    = -999.;
      // Fiducial flags 
      RECOELE_isbarrel[i]    = 0;
      RECOELE_isendcap[i]    = 0;
      RECOELE_isGap[i]  = 0;
      RECOELE_isEBetaGap[i]  = 0;
      RECOELE_isEBphiGap[i]  = 0;
      RECOELE_isEEdeeGap[i]  = 0;
      RECOELE_isEEringGap[i] = 0;
      // Shower shape
      RECOELE_sigmaIetaIeta[i] = -999.;
      RECOELE_sigmaEtaEta[i]   = -999.;
      RECOELE_e15[i]           = -999.;
      RECOELE_e25max[i]        = -999.;
      RECOELE_e55[i]           = -999.;
      RECOELE_he[i]            = -999.;
      RECOELE_r9[i]            = -999.;
      // Particle flow
      RECOELE_mva[i] = -999.; 
      // Brem & Classifaction
      RECOELE_fbrem[i]  = -999.; 
      RECOELE_nbrems[i] = -999;
      RECOELE_Class[i]  = -999;
      RECOELE_fbrem_mode[i]=-999.;
      RECOELE_fbrem_mean[i]=-999.;
      // isolation
      RECOELE_EGMTRACKISO[i]=-999.; 
      RECOELE_EGMHCALISO[i]=-999.;
      RECOELE_EGMECALISO[i]=-999.;
      RECOELE_EGMX[i]=-999.;

      // PF isolation
      RECOELE_PFchAllPart[i]=-999.;
      RECOELE_PFchHad[i]=-999.;
      RECOELE_PFneuHad[i]=999.;
      RECOELE_PFphoton[i]=-999.;
      RECOELE_PFPUchAllPart[i]=-999.;
      RECOELE_PFX_dB[i]=-999.;
      RECOELE_PFX_rho[i]=-999.;

      // Electron regression
      RECOELE_regEnergy[i]=-999.;
      RECOELE_regEnergyError[i]=-999.;

      // TLE electrons
      RECOELE_TLE_ParentSC_X[i]=-999.;
      RECOELE_TLE_ParentSC_Y[i]=-999.;
      RECOELE_TLE_ParentSC_Z[i]=-999.;

      // IP
      RECOELE_IP[i]=-9999.;
      RECOELE_SIP[i]=-9999.;
      RECOELE_IPERROR[i]=-9999.;
      RECOELE_STIP[i]=-999.;
      RECOELE_TIP[i]=-999.;
      RECOELE_TIPERROR[i]=-999.;
      RECOELE_SLIP[i]=-999.;
      RECOELE_LIP[i]=-999.;
      RECOELE_LIPERROR[i]=-999.;
      ele_sclRawE[i]=-999. ;
      ele_sclX[i]=-999.; 
      ele_sclY[i]=-999.; 
      ele_sclZ[i]=-999.;
      ele_seedSubdet1[i]=-999.;
      ele_seedDphi1[i]=-999.; 
      ele_seedDrz1[i]=-999.;
      ele_seedSubdet2[i]=-999.;
      ele_seedDphi2[i]=-999.; 
      ele_seedDrz2[i]=-999.;
      ele_eidVeryLoose[i]=-999.; 
      ele_eidLoose[i]=-999.; 
      ele_eidMedium[i]=-999.; 
      ele_eidTight[i]=-999. ;
      ele_eidHZZVeryLoose[i]=-999.; 
      ele_eidHZZLoose[i]=-999.; 
      ele_eidHZZMedium[i]=-999.; 
      ele_eidHZZTight[i]=-999. ;
      RECOELE_mvaTrigV0[i]=-999.;
      RECOELE_mvaNonTrigV0[i]=-999.;
    
      // Conversion
      ConvMapDist[i]=-999.;
      ConvMapDcot[i]=-999.;
      
      // Muon block
      RECOMU_PT_MuHLTMatch[i]=-999.;
      RECOMU_ETA_MuHLTMatch[i]=-999.;
      RECOMU_PHI_MuHLTMatch[i]=-999.;
      RECOELE_PT_EleHLTMatch[i]=-999.;
      RECOELE_ETA_EleHLTMatch[i]=-999.;
      RECOELE_PHI_EleHLTMatch[i]=-999.;
      RECOMU_sm_MuHLTMatch[i]=false;
      RECOMU_dm_MuHLTMatch[i]=-999;
      RECOELE_se_EleHLTMatch[i]=false;
      RECOELE_de_EleHLTMatch[i]=-999;
      RECOJET_JetHLTMatch[i]=0;   

      RECOMU_isPFMu[i]=false;
      RECOMU_isGlobalMu[i]=false;
      RECOMU_isStandAloneMu[i]=false;
      RECOMU_isTrackerMu[i]=false;
      RECOMU_isCaloMu[i]=false;
      RECOMU_isTrackerHighPtMu[i]=false;   
      RECOMU_isME0Muon[i]=false;

      RECOMU_E[i]=-999.;
      RECOMU_PT[i]=-999.;
      RECOMU_P[i]=-999.;
      RECOMU_PHI[i]=-999.;
      RECOMU_ETA[i]=-999.;
      RECOMU_THETA[i]=-999.;
      RECOMU_MASS[i]=-999.;      
      RECOMU_CHARGE[i]=-999.;
      
      for (int j=0; j<3; j++){
	for (int k=0;k<3;k++){
	  RECOMU_COV[i][j][k]=-999.;
	  RECOELE_COV[i][j][k]=-999.;
	}
      }

      RECOMU_TRACKISO[i]=-999.; 
      RECOMU_TRACKISO_SUMPT[i]=-999.; 
      RECOMU_ECALISO[i]=-999.;
      RECOMU_HCALISO[i]=-999.;
      RECOMU_X[i]=-999.;

      RECOMU_PFchHad[i]=-999.;
      RECOMU_PFneuHad[i]=999.;
      RECOMU_PFphoton[i]=-999.;
      RECOMU_PFX_dB[i]=-999.;
      RECOMU_PFX_rho[i]=-999.;

/*
      RECOMU_IP[i]=-9999.;
      RECOMU_SIP[i]=-9999.;
      RECOMU_IPERROR[i]=-9999.;
      RECOMU_IP_KF[i]=-999.;
      RECOMU_SIP_KF[i]=-999.;
      RECOMU_IPERROR_KF[i]=-999.;

      RECOMU_SIP_GD[i]=-999.;
      RECOMU_SIP_GDMMMM[i]=-999.;
      RECOMU_SIP_Std[i]=-999.;
      RECOMU_SIP_StdMMMM[i]=-999.;
      RECOMU_SIP_Kin[i]=-999.;
      RECOMU_SIP_KinMMMM[i]=-999.;
      RECOMU_STIP[i]=-999.;
      RECOMU_TIP[i]=-999.;
      RECOMU_TIPERROR[i]=-999.;
      RECOMU_SLIP[i]=-999.;
      RECOMU_LIP[i]=-999.;
      RECOMU_LIPERROR[i]=-999.;
 */     
 

      RECOMU_numberOfMatches[i]=-999;
      RECOMU_numberOfMatchedStations[i]=-999.;
      RECOMU_caloCompatibility[i]=-999.;
      RECOMU_segmentCompatibility[i]=-999.;
      RECOMU_glbmuPromptTight[i]=false;
      
      RECOMU_mutrkPT[i]=-999.;
      RECOMU_mutrkPTError[i]=-999.;
      RECOMU_mutrkDxy[i]=-999.;
      RECOMU_mutrkDxyError[i]=-999.;
      RECOMU_mutrkDxyB[i]=-999.;
      RECOMU_mutrkDz[i]=-999.;
      RECOMU_mutrkDzError[i]=-999.;
      RECOMU_mutrkDzB[i]=-999.;
      RECOMU_mutrkChi2PerNdof[i]=-999.;
      RECOMU_mutrktrackerLayersWithMeasurement[i]=-999.;

      RECOMU_muInnertrkDxy[i]=-999.;
      RECOMU_muInnertrkDxyError[i]=-999.;
      RECOMU_muInnertrkDxyB[i]=-999.;
      RECOMU_muInnertrkDz[i]=-999.;
      RECOMU_muInnertrkDzError[i]=-999.;
      RECOMU_muInnertrkDzB[i]=-999.;
      RECOMU_muInnertrkChi2PerNdof[i]=-999.;
      RECOMU_muInnertrktrackerLayersWithMeasurement[i]=-999.;
      RECOMU_muInnertrkCharge[i]=-999.;
      RECOMU_muInnertrkNHits[i]=-999.;
      RECOMU_muInnertrkNPixHits[i]=-999.;
      RECOMU_muInnertrkNStripHits[i]=-999.;
      RECOMU_muInnertrkPT[i]=-999.;
      RECOMU_muInnertrkPTError[i]=-999.;
      RECOMU_muInnertrkvalidFraction[i]=-999; 
     
      RECOMU_mubesttrkType[i]=-999;
      RECOMU_mubesttrkDxy[i]=-999.;
      RECOMU_mubesttrkDxyB[i]=-999.;
      RECOMU_mubesttrkDxyError[i]=-999.;
      RECOMU_mubesttrkDz[i]=-999.;
      RECOMU_mubesttrkDzB[i]=-999.;
      RECOMU_mubesttrkDzError[i]=-999.;
      RECOMU_mubesttrkPTError[i]=-999.;

      RECOMU_mutrkCharge[i]=-999.;
      RECOMU_mutrkNHits[i]=-999.;
      RECOMU_mutrkNPixHits[i]=-999.;
      RECOMU_mutrkNMuonHits[i]=-999.;
      RECOMU_mutrkNStripHits[i]=-999.;
      RECOMU_trkmuArbitration[i]=false;
      RECOMU_trkmu2DCompatibilityLoose[i]=false;
      RECOMU_trkmu2DCompatibilityTight[i]=false;
      RECOMU_trkmuOneStationLoose[i]=false;
      RECOMU_trkmuOneStationTight[i]=false;
      RECOMU_trkmuLastStationLoose[i]=false;
      RECOMU_trkmuLastStationTight[i]=false;
      RECOMU_trkmuOneStationAngLoose[i]=false;
      RECOMU_trkmuOneStationAngTight[i]=false;
      RECOMU_trkmuLastStationAngLoose[i]=false;
      RECOMU_trkmuLastStationAngTight[i]=false;
      RECOMU_trkmuLastStationOptimizedLowPtLoose[i]=false;
      RECOMU_trkmuLastStationOptimizedLowPtTight[i]=false;

    }
    

    for (int i=0; i<50;i++){      
      //Matching 

      // Muons
      RECOMU_MatchingMCTruth[i]=false;
      RECOMU_MatchingMCpT[i]=-999.;
      RECOMU_MatchingMCEta[i]=-999.;
      RECOMU_MatchingMCPhi[i]=-999.;
      
      // Electrons
      RECOELE_MatchingMCTruth[i]=false;
      RECOELE_MatchingMCpT[i]=-999.;
      RECOELE_MatchingMCEta[i]=-999.;
      RECOELE_MatchingMCPhi[i]=-999.;
    
      // Bottoms
      RECOBOT_MatchingMCTruth[i]=0;
      RECOBOT_MatchingMCpT[i]=-999.;
      RECOBOT_MatchingMCEta[i]=-999.;
      RECOBOT_MatchingMCPhi[i]=-999.;
      }
    for (int i=0; i<50;i++){
      //Gamma:
      RECOPHOT_MatchingMCTruth[i]=false;
      RECOPHOT_MatchingMCpT[i]=-999.;
      RECOPHOT_MatchingMCEta[i]=-999.;
      RECOPHOT_MatchingMCPhi[i]=-999.;

/*      
      //zToMuMu:
      RECOzMuMu_MatchingMCTruth[i]=false;
      RECOzMuMu_MatchingMCpT[i]=-999.;
      RECOzMuMu_MatchingMCmass[i]=-999.;
      RECOzMuMu_MatchingMCEta[i]=-999.;
      RECOzMuMu_MatchingMCPhi[i]=-999.;
      
      //zToEE:
      RECOzEE_MatchingMCTruth[i]=false;
      RECOzEE_MatchingMCpT[i]=-999.;
      RECOzEE_MatchingMCmass[i]=-999.;*/
    }

    for (int i=0; i<100;i++){
      RECO_TRACK_PT[i]=-999.;
      RECO_TRACK_ETA[i]=-999.;
      RECO_TRACK_PHI[i]=-999.;
      RECO_TRACK_CHI2[i]=-999.;
      RECO_TRACK_CHI2RED[i]=-999.;
      RECO_TRACK_CHI2PROB[i]=-999.;
      RECO_TRACK_NHITS[i]=0;
      RECO_TRACK_DXY[i]=-999.;
      RECO_TRACK_DXYERR[i]=-999.;
      RECO_TRACK_DZ[i]=-999.;
      RECO_TRACK_DZERR[i]=-999.;   			       					
    }

    for (int i=0; i<50;i++){
      tCHighEff_BTagJet_PT[i]=-999.;
//      tCHighPur_BTagJet_PT[i]=-999.;
      cSV_BTagJet_PT[i]=-999.;
      tCHighEff_BTagJet_ETA[i]=-999.;
 //     tCHighPur_BTagJet_ETA[i]=-999.;
      cSV_BTagJet_ETA[i]=-999.;
      tCHighEff_BTagJet_PHI[i]=-999.;
 //     tCHighPur_BTagJet_PHI[i]=-999.;
      cSV_BTagJet_PHI[i]=-999.;
      tCHighEff_BTagJet_DISCR[i]=-999.;
 //     tCHighPur_BTagJet_DISCR[i]=-999.;
      cSV_BTagJet_DISCR[i]=-999.;
      cSV_BTagJet_ET[i]=-999.;
    }
    
  }
  
  void fillPU(const edm::Event& iEvent){
      edm::Handle<vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByToken(PileupSrc_, PupInfo);

      if(!PupInfo.isValid()) return;
      for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand ) { 
	std::cout << " Pileup Information: bunchXing, nvtx: " << cand->getBunchCrossing() << " " << cand->getPU_NumInteractions() << std::endl;
	if (cand->getBunchCrossing() == 0) num_PU_vertices=cand->getTrueNumInteractions();;
	// num_PU_vertices=cand->getPU_NumInteractions(); in-time,out-of-time pileup
	PU_BunchCrossing=cand->getBunchCrossing();
      }	
  }

  void EventsMCReWeighting(const edm::Event& iEvent){

    // get the weight                                                                                                                                                     
    MC_weighting=0.;
    for(int j=0; j<9; j++) MC_weighting_un[j]=0;
    PDF_weighting_un=0;
    float EventWeight = 1.0;
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    iEvent.getByToken(generator_, gen_ev_info);


    edm::Handle<LHEEventProduct> lheInfo;
    iEvent.getByToken(lheEventProductToken_, lheInfo);

    if(!gen_ev_info.isValid()) return;
    EventWeight = gen_ev_info->weight();
    std::cout<<"mc_weight = "<< gen_ev_info->weight() <<std::endl;
                                                                                                                                                                        
    float mc_weight = EventWeight;
   // ( EventWeight > 0 ) ? 1 : -1;
    //std::cout<<"mc_weight = "<< mc_weight <<std::endl;                                                                                                                  
    MC_weighting=mc_weight;
    if(fillLHEinfo){
      double nomlheweight;
      if (lheInfo.isValid()) nomlheweight = lheInfo->weights()[0].wgt;
      if (lheInfo.isValid()) {
        for (size_t i=0; i<9; ++i) {
            MC_weighting_un[i]=lheInfo->weights()[i].wgt/nomlheweight;
         }

      //pdf weight

      double NNPDF3wgtAvg = 0.0;

      double NNPDF3wgtRMS = 0.0;
      double NNPDF3wgt = 0.0;
      double NNPDF3wgt_frac = 0.0;
      double centralWgt=nomlheweight;

      unsigned int PDFstart = 9;
      unsigned int PDFend = 109;

      for (unsigned int i_lhePDF = PDFstart; i_lhePDF < PDFend; ++i_lhePDF){
         NNPDF3wgt = lheInfo->weights()[i_lhePDF].wgt;
         NNPDF3wgt_frac = NNPDF3wgt/(centralWgt);
         NNPDF3wgtAvg += NNPDF3wgt_frac;
      }

      NNPDF3wgtAvg = NNPDF3wgtAvg/(PDFend - PDFstart);
      for (unsigned int i_lhePDF = PDFstart; i_lhePDF < PDFend; ++i_lhePDF){
         NNPDF3wgt = lheInfo->weights()[i_lhePDF].wgt;
         NNPDF3wgt_frac = NNPDF3wgt/(centralWgt);
         NNPDF3wgtRMS += (NNPDF3wgt_frac - NNPDF3wgtAvg)*(NNPDF3wgt_frac - NNPDF3wgtAvg);
      }

      NNPDF3wgtRMS = sqrt(NNPDF3wgtRMS/(PDFend - PDFstart - 1));
      PDF_weighting_un = NNPDF3wgtRMS;
      cout << "pdf_un=" <<PDF_weighting_un << endl; 
     } 
    } 
  }


  void fillHLTFired(const edm::Event& iEvent){
    edm::Handle<vector<std::string> > HLTfired_;
    iEvent.getByLabel(HLTInfoFired,HLTfired_);

    vector<string> HLTimported;
    string tmpstring="";

    for (vector<string>::const_iterator cand=HLTfired_->begin(); cand!=HLTfired_->end(); ++cand){
      unsigned int i=cand-HLTfired_->begin();
      HLTimported.push_back(cand->c_str());
      string newstr=HLTimported.at(i) + ":" + tmpstring;
      tmpstring=newstr;
    }

    cout << "HLTFiredString= " << tmpstring.c_str() << endl;
    if (!tmpstring.empty()) sprintf(HLTPathsFired,tmpstring.c_str());
        

  }


  void triggermatching(const edm::Event& iEvent){

    cout << "Start Trigger matching for muon" << endl;
    // check HLTrigger/Configuration/python/HLT_GRun_cff.py

//AOD    edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales,triggerPrescalesl1min,triggerPrescalesl1max;
    iEvent.getByToken(triggerObjects_, triggerObjects );
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);
    iEvent.getByToken(triggerPrescalesL1min_, triggerPrescalesl1min);
    iEvent.getByToken(triggerPrescalesL1max_, triggerPrescalesl1max);
//    const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
    size_t nMuHLT =0, nEleHLT=0;


    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);    


    std::vector<pat::TriggerObjectStandAlone>  HLTMuMatched_sm, HLTMuMatched_dm,HLTEleMatched_se,HLTEleMatched_de,HLTMu17,HLTMu8,HLTEleLeg1,HLTEleLeg2;
    std::vector<string> HLTMuMatchedNames_sm,HLTMuMatchedNames_dm,HLTEleMatchedNames_se,HLTEleMatchedNames_de;
    std::vector<pat::TriggerObjectStandAlone> HLTJetMatched_40, HLTJetMatched_60, HLTJetMatched_80, HLTJetMatched_140, HLTJetMatched_200, HLTJetMatched_260, HLTJetMatched_320;

/*  
    std::cout << "\n == TRIGGER PATHS= " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        std::cout << "Trigger " << names.triggerName(i) <<
                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
                ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << "\n" 
                << triggerPrescalesl1min->getPrescaleForIndex(i) <<"l1min " 
                << triggerPrescalesl1max->getPrescaleForIndex(i) <<"l1max "
                << std::endl;
    }
*/
//MiniAOD

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
       obj.unpackPathNames(names);
//       std::cout << "\t   Collection: " << obj.collection() << std::endl;
//       for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] << endl;
       obj.unpackFilterLabels(iEvent,*triggerBits);
//        std::cout << std::endl;
        // Print associated trigger filters
//        std::cout << "\t   Filters:    ";
        bool mu17_pass = false;
        bool mu8_pass = false;
        bool ele_leg1 = false;
        bool ele_leg2 = false;
        int jet_pass =0;
        for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
             TString flt = obj.filterLabels()[h].c_str();
//             cout <<"test flt=" << flt << endl;
             if(flt=="hltL3fL1DoubleMu155fFiltered17") mu17_pass=true;
             if(flt=="hltL3fL1DoubleMu155fPreFiltered8") mu8_pass=true;
             if(flt=="hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter") ele_leg1=true;
             if(flt=="hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter") ele_leg2=true;
             if(flt=="hltSinglePFJet40"&& jet_pass<1) {jet_pass=1;}
             if(flt=="hltSinglePFJet60"&& jet_pass<2) {jet_pass=2;}
             if(flt=="hltSinglePFJet80"&& jet_pass<3) {jet_pass=3;}
             if(flt=="hltSinglePFJet140"&& jet_pass<4) {jet_pass=4;}
             if(flt=="hltSinglePFJet200"&& jet_pass<5) {jet_pass=5;}
             if(flt=="hltSinglePFJet260"&& jet_pass<6) {jet_pass=6;}
             if(flt=="hltSinglePFJet320"&& jet_pass<7) {jet_pass=7;}       

        }

        std::vector<string> pathNamesAll = obj.pathNames(false);
        std::vector<string> pathNamesLast = obj.pathNames(true);
/*
        std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
            bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
            bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
            bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
            std::cout << "   " << pathNamesAll[h];
            if (isBoth) std::cout << "(L,3)";
            if (isL3 && !isBoth) std::cout << "(*,3)";
            if (isLF && !isBoth) std::cout << "(L,*)";
            if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
        }
*/
        bool hlt_pass_sm = false;
        bool hlt_pass_dm = false;
        bool hlt_pass_se = false;
        bool hlt_pass_de = false;
        bool hlt_pass_tri = false;

        for(unsigned h = 0; h < pathNamesLast.size(); h++)
         for(unsigned h2 = 0; h2 < HLTFilter_.size(); h2++){
           TString hlt = pathNamesLast[h].c_str();
           if(hlt.Contains(HLTFilter_[h2].c_str())){ 
            if(HLTFilter_[h2]=="HLT_IsoMu24_v"||HLTFilter_[h2]=="HLT_IsoTkMu24_v") hlt_pass_sm=true;
            if(HLTFilter_[h2]=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v" || HLTFilter_[h2]=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") hlt_pass_dm=true;
            if(HLTFilter_[h2]=="HLT_Ele25_eta2p1_WPTight_Gsf_v"||HLTFilter_[h2]=="HLT_Ele27_WPTight_Gsf_v") hlt_pass_se=true;
            if(HLTFilter_[h2]=="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") hlt_pass_de=true;
            if(HLTFilter_[h2]=="HLT_TripleMu_12_10_5_v"||HLTFilter_[h2]=="HLT_TripleMu_10_5_5_DZ_v"||HLTFilter_[h2]=="HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"||HLTFilter_[h2]=="HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"||HLTFilter_[h2]=="HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") hlt_pass_tri=true;
            if(HLTFilter_[h2]=="HLT_PFJet40_v") jet_trig=1;
            if(HLTFilter_[h2]=="HLT_PFJet60_v") jet_trig=2;
            if(HLTFilter_[h2]=="HLT_PFJet80_v") jet_trig=3;
            if(HLTFilter_[h2]=="HLT_PFJet140_v") jet_trig=4;
            if(HLTFilter_[h2]=="HLT_PFJet200_v") jet_trig=5;
            if(HLTFilter_[h2]=="HLT_PFJet260_v") jet_trig=6;
            if(HLTFilter_[h2]=="HLT_PFJet320_v") jet_trig=7;
            
            cout << "Matching " << HLTFilter_[h2].c_str()  << endl;
           }
           }
        
       if(hlt_pass_tri) tri_trig=true;
       if(hlt_pass_sm)  {HLTMuMatched_sm.push_back(obj); HLTMuMatchedNames_sm.push_back("IsoMu24");sm_trig=true;}
       if(hlt_pass_dm)  {HLTMuMatched_dm.push_back(obj); HLTMuMatchedNames_dm.push_back("DoubleMu");dm_trig=true;}
       if(hlt_pass_se)  {HLTEleMatched_se.push_back(obj); HLTEleMatchedNames_se.push_back("SingleEG");se_trig=true;}
       if(hlt_pass_de)  {HLTEleMatched_de.push_back(obj); HLTEleMatchedNames_de.push_back("DoubleEG");de_trig=true;}
       if(mu17_pass)     HLTMu17.push_back(obj);
       if(mu8_pass)      HLTMu8.push_back(obj);
       if(ele_leg1)      HLTEleLeg1.push_back(obj);
       if(ele_leg2)      HLTEleLeg2.push_back(obj);
       if(hlt_pass_sm||hlt_pass_dm) nMuHLT++;
       if(hlt_pass_se||hlt_pass_de) nEleHLT++;

       if(jet_pass==7) HLTJetMatched_320.push_back(obj);
       if(jet_pass==6) HLTJetMatched_260.push_back(obj);
       if(jet_pass==5) HLTJetMatched_200.push_back(obj);
       if(jet_pass==4) HLTJetMatched_140.push_back(obj);
       if(jet_pass==3) HLTJetMatched_80.push_back(obj);
       if(jet_pass==2) HLTJetMatched_60.push_back(obj);
       if(jet_pass==1) HLTJetMatched_40.push_back(obj);

      }
    


    edm::Handle<edm::View<pat::Muon> > MuCandidates;
    iEvent.getByToken(muonTag_, MuCandidates);
    float maxDeltaR_=0.15;
    float maxDPtRel_=1.0;
    int nMuHLTMatch=0;
    float minDR=0.5;
    int Ni=-1;
    for (edm::View<pat::Muon>::const_iterator iCand = MuCandidates->begin(); iCand != MuCandidates->end(); ++iCand){
      unsigned int i=iCand-MuCandidates->begin();

      cout << "Muon with pt= " << iCand->pt() << ": check trigger matching" << endl;
      if (IsMuMatchedToHLTMu(*iCand,  HLTMuMatched_sm , maxDeltaR_, maxDPtRel_)==true){
        RECOMU_sm_MuHLTMatch[i]=true;
        RECOMU_PT_MuHLTMatch[i] =iCand->pt();
        RECOMU_ETA_MuHLTMatch[i]=iCand->eta();
        RECOMU_PHI_MuHLTMatch[i]=iCand->phi();
      }
      if (IsMuMatchedToHLTMu(*iCand,  HLTMu8 , maxDeltaR_, maxDPtRel_)==true){
        RECOMU_PT_MuHLTMatch[i] =iCand->pt();
        RECOMU_ETA_MuHLTMatch[i]=iCand->eta();
        RECOMU_PHI_MuHLTMatch[i]=iCand->phi();
        RECOMU_dm_MuHLTMatch[i]=1;
      }
      if (IsMuMatchedToHLTMu(*iCand,  HLTMu17 , maxDeltaR_, maxDPtRel_)==true){
        RECOMU_dm_MuHLTMatch[i]=2;
        RECOMU_PT_MuHLTMatch[i] =iCand->pt();
        RECOMU_ETA_MuHLTMatch[i]=iCand->eta();
        RECOMU_PHI_MuHLTMatch[i]=iCand->phi();
      }
      
      if(IsMuMatchedToHLTMu(*iCand,  HLTMuMatched_sm , maxDeltaR_, maxDPtRel_)==true
       ||IsMuMatchedToHLTMu(*iCand,  HLTMuMatched_dm , maxDeltaR_, maxDPtRel_)==true){
        nMuHLTMatch++;
        cout << "Muon HLT Matched with pT= " << iCand->pt() << endl;
      }

    }
//    cout << "DR = " << minDR << " index= " << Ni << endl;

    cout << "N. Muons HLT Matched= " << nMuHLTMatch << endl;
    RECO_nMuHLTMatch    = nMuHLTMatch;


    cout << "Start Trigger matching for electron" << endl;

    int nEleHLTMatch=0;
    edm::Handle<edm::View<pat::Electron> > EleCandidates;
    iEvent.getByToken(electronEgmTag_, EleCandidates);

    minDR=0.5;
    Ni=-1;

    for (edm::View<pat::Electron>::const_iterator iCand = EleCandidates->begin(); iCand != EleCandidates->end(); ++iCand){

      unsigned int i=iCand-EleCandidates->begin();
      cout << "Electron with pt= " << iCand->pt() << ": check trigger matching" << endl;
      for(size_t k=0; k < HLTEleMatched_se.size(); k++)
            if(deltaR(HLTEleMatched_se[k],*iCand)<minDR){minDR=deltaR(HLTEleMatched_se[k],*iCand); Ni=i;}

      if (IsEleMatchedToHLTEle(*iCand,  HLTEleMatched_se ,  maxDeltaR_, maxDPtRel_)==true){
            RECOELE_se_EleHLTMatch[i]=true;
            RECOELE_PT_EleHLTMatch[i]=iCand->pt();
            RECOELE_ETA_EleHLTMatch[i]=iCand->eta();
            RECOELE_PHI_EleHLTMatch[i]=iCand->phi();
      }
      if (IsEleMatchedToHLTEle(*iCand,  HLTEleLeg2 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOELE_de_EleHLTMatch[i]=1;
            RECOELE_PT_EleHLTMatch[i]=iCand->pt();
            RECOELE_ETA_EleHLTMatch[i]=iCand->eta();
            RECOELE_PHI_EleHLTMatch[i]=iCand->phi();
      }
     if (IsEleMatchedToHLTEle(*iCand,  HLTEleLeg1 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOELE_de_EleHLTMatch[i]=2;
            RECOELE_PT_EleHLTMatch[i]=iCand->pt();
            RECOELE_ETA_EleHLTMatch[i]=iCand->eta();
            RECOELE_PHI_EleHLTMatch[i]=iCand->phi();
      }
     if (IsEleMatchedToHLTEle(*iCand,  HLTEleMatched_se , maxDeltaR_, maxDPtRel_)==true
       ||IsEleMatchedToHLTEle(*iCand,  HLTEleMatched_de , maxDeltaR_, maxDPtRel_)==true){
            cout << "Electron HLT Matched with pT= " << iCand->pt() << endl;
            nEleHLTMatch++;
      }

    cout << "DR = " << minDR << " index= " << Ni << endl;
    for(int i=0; i < 50; i++)
      if(i!=Ni) RECOELE_se_EleHLTMatch[i] = false;

    }

    cout << "N. Electrons HLT Matched= " << nEleHLTMatch << endl;

    RECO_nEleHLTMatch = nEleHLTMatch; 

    edm::Handle<pat::JetCollection> pfjets;
    iEvent.getByToken(jetsTag_, pfjets);

    for (pat::JetCollection::const_iterator iCand = pfjets->begin(); iCand != pfjets->end(); ++iCand){

      unsigned int i=iCand-pfjets->begin();

      if (IsJetMatchedToHLTJet(*iCand,  HLTJetMatched_320 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOJET_JetHLTMatch[i]=7;
      }
      else if (IsJetMatchedToHLTJet(*iCand,  HLTJetMatched_260 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOJET_JetHLTMatch[i]=6;
      }
      else if (IsJetMatchedToHLTJet(*iCand,  HLTJetMatched_200 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOJET_JetHLTMatch[i]=5;
      }
      else if (IsJetMatchedToHLTJet(*iCand,  HLTJetMatched_140 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOJET_JetHLTMatch[i]=4;
      }
      else if (IsJetMatchedToHLTJet(*iCand,  HLTJetMatched_80 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOJET_JetHLTMatch[i]=3;
      }
      else if (IsJetMatchedToHLTJet(*iCand,  HLTJetMatched_60 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOJET_JetHLTMatch[i]=2;
      }
      else if (IsJetMatchedToHLTJet(*iCand,  HLTJetMatched_40 ,  maxDeltaR_, maxDPtRel_)==true){
            RECOJET_JetHLTMatch[i]=1;
      }
    }


  }

  bool IsMuMatchedToHLTMu ( const pat::Muon &mu, std::vector<pat::TriggerObjectStandAlone> HLTMu , double DR, double DPtRel ) {
    size_t dim =  HLTMu.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      //cout << "HLT mu filter is= " << HLTMuNames[k].c_str() << " Delta R= " << deltaR(HLTMu[k], mu) << " Delta pT= " << fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt() << endl;
      if (  (deltaR(HLTMu[k], mu) < DR)   && (fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt()<DPtRel)){ 
	cout << "HLT mu filter pt is= " << mu.pt() << " Delta R= " << deltaR(HLTMu[k], mu) << " Delta pT= " << fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt() << endl;
	nPass++ ;
      }
    }
    return (nPass>0);
  }

  bool IsEleMatchedToHLTEle ( const pat::Electron &ele, std::vector<pat::TriggerObjectStandAlone> HLTEle , double DR, double DPtRel ) {
    size_t dim =  HLTEle.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      //cout << "HLT ele filter is= " << HLTEleNames[k].c_str() << " Delta R= " << deltaR(HLTEle[k], ele) << " Delta pT= " << fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt() << endl;
      if (  (deltaR(HLTEle[k], ele) < DR)   && (fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt()<DPtRel)){ 
	cout << "HLT ele filter is= " << " Delta R= " << deltaR(HLTEle[k], ele) << " Delta pT= " << fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt() << endl;
	nPass++ ;
      }
    }
    return (nPass>0);
  }

  bool IsJetMatchedToHLTJet ( const pat::Jet &jet, std::vector<pat::TriggerObjectStandAlone> HLTJet , double DR, double DPtRel ) {
    size_t dim =  HLTJet.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      //cout << "HLT ele filter is= " << HLTEleNames[k].c_str() << " Delta R= " << deltaR(HLTEle[k], ele) << " Delta pT= " << fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt() << endl;
      if (  (deltaR(HLTJet[k], jet) < DR)   && (fabs(HLTJet[k].pt() - jet.pt())/ HLTJet[k].pt()<DPtRel)){
        nPass++ ;
      }
    }
    return (nPass>0);
  }

  bool isTrackerHighPtMu (const reco::Muon &mu, math::XYZPoint primaryVertex){
    return( mu.numberOfMatchedStations() > 1 &&
            (mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt()) < 0.3
            && std::abs(mu.muonBestTrack()->dxy(primaryVertex)) < 0.2
            && std::abs(mu.muonBestTrack()->dz(primaryVertex)) < 0.5
            && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
            && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
            );
  }

  
  // PDT
  std::string getParticleName(int id) const{
    const ParticleData * pd = pdt_->particle( id );
    if (!pd) {
      std::ostringstream ss;
      ss << "P" << id;
      return ss.str();
    } else
      return pd->name();
  }

  
  void fillgenjets(const edm::Event& iEvent){
    edm::Handle<reco::GenJetCollection> genjetHandle;
    iEvent.getByToken(genjetTag_,genjetHandle);
    int i=0;
    int tem=0;
    if (genjetHandle.isValid()) {
    for ( GenJetCollection::const_iterator igen=genjetHandle->begin(); igen!=genjetHandle->end(); igen++) {
	std::exception_ptr eptr=nullptr;
      if (i>49) break;
      MC_GENJET_PT[i]=igen->pt();
      MC_GENJET_ETA[i]=igen->eta();
      MC_GENJET_PHI[i]=igen->phi();     
      if(abs(MC_GENJET_ETA[i])<4.7 && MC_GENJET_PT[i]>30) {tem++;cout << "hellgen" << endl;}
//      cout <<"genpt=" << MC_GENJET_PT[i] << " geneta=" << MC_GENJET_ETA[i] << " phi=" << MC_GENJET_PHI[i] << endl;
//      cout <<"daughters=" << igen->numberOfDaughters() << endl;
//      cout <<"mothers=" << igen->numberOfMothers() << endl;
      i++;
	cout << i << endl;
        bool matchlep=false;
        for(unsigned l=0; l < igen->numberOfDaughters(); l++){
          // Select status 3 electron, muon or tau
          cout << "l=" << l << endl;
		try {
		  eptr=nullptr;
		  auto tmp=igen->daughter(l);
//	          cout << "daughter ID=" << tmp->pdgId() << endl;
		} catch(const std::exception &e) { eptr=std::current_exception();  cout << e.what(); 
		} 
		if(eptr) { break; }
//          cout << "daughter ID=" << (igen->daughter(l))->pdgId() << endl;
          if(abs((igen->daughter(l))->pdgId())==11 || abs((igen->daughter(l))->pdgId())==13){
      cout << "status=" << (igen->daughter(l))->status() << endl;
      double deltaphi=abs((igen->daughter(l))->phi()-MC_GENJET_PHI[i-1]);
      if(deltaphi>3.1416) deltaphi-=2*3.1416;
      double deltaR = sqrt( pow( deltaphi,2) + pow((igen->daughter(l))->eta() - MC_GENJET_ETA[i-1],2) ); 
//      cout << "deltaphi=" << deltaphi << " deltaEta="<<abs((igen->daughter(l))->eta() - MC_GENJET_ETA[i]) << endl; 
      cout << "deltaR=" << deltaR << endl;
      if(deltaR<0.3) matchlep=true;
      }     
         }
      if(matchlep) MC_GENJET_PT[i-1]=-MC_GENJET_PT[i-1];
    }}
/*
   if (genjetHandle.isValid()) {
   for(unsigned i=0; i < genjetHandle->size(); i++){
     if(i>99) break;
      const reco::GenJet& igen = (*genjetHandle) [i];
      MC_GENJET_PT[i]=igen.pt();
      MC_GENJET_ETA[i]=igen.eta();
      MC_GENJET_PHI[i]=igen.phi();
      if(abs(MC_GENJET_ETA[i])<4.7 && MC_GENJET_PT[i]>30) {tem++;cout << "hellgen" << endl;}
      cout <<"genpt=" << MC_GENJET_PT[i] << " geneta=" << MC_GENJET_ETA[i] << " phi=" << MC_GENJET_PHI[i] << endl;
      cout <<"daughters=" << igen.numberOfDaughters() << endl;
      cout <<"mothers=" << igen.numberOfMothers() << endl;

   }
  }*/

   cout << "genjet=" << tem << endl;
  }

  // GenParticles  
  void fillgenparticles(const edm::Event& iEvent, const edm::EventSetup &es){    
    
    // get gen particle candidates 
    edm::Handle<reco::GenParticleCollection> genCandidates;           
    iEvent.getByToken(genParticles_, genCandidates);

    es.getData( pdt_ );
/*
    vector<float> leptonpt,leptoneta,leptonphi,leptontheta,leptonpdgid;
    
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      // lepton stable      
   
      if ( ( abs(mcIter->pdgId())==11 || abs(mcIter->pdgId())==13 || abs(mcIter->pdgId())==15  ) && mcIter->status()==1 ){	
	
	std::cout << " \n Found lepton with Id= " << mcIter->pdgId() << "  in the final state with mother= " ; 
	leptonpt.push_back(mcIter->pt());
	leptoneta.push_back(mcIter->eta());
	leptonphi.push_back(mcIter->phi());
	leptontheta.push_back(mcIter->theta());
        leptonpdgid.push_back(float(mcIter->pdgId()));

	if (mcIter->numberOfMothers()>0){ 
	  std::cout << getParticleName(int(mcIter->mother(0)->pdgId())) << " << ";
	  if (mcIter->mother(0)->status()==3) continue;
	  if (mcIter->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->status()!=3 ){
	    std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->pdgId())) << " << ";
	    if (mcIter->mother(0)->mother(0)->status()==3) continue;
	    if (mcIter->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->status()!=3 ){
	      std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
	      if (mcIter->mother(0)->mother(0)->mother(0)->status()==3) continue;
	      if (mcIter->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		  std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && 
mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		    std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << " << std::endl;
		    if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		  }
		}
	      }
	    }
	  }
	}

    
      }
    }

    sort(leptonpt.rbegin(),leptonpt.rend());
    
    if (leptonpt.size()>0) { 
      cout << "\n The 4 highest pt leptons in MCtruth are= ";
      for (unsigned int i=0;i<leptonpt.size();i++){
	if (i>3) continue;
	cout << leptonpt.at(i) << " ";
	MC_LEPT_PT[i]=leptonpt.at(i);
	MC_LEPT_ETA[i]=leptoneta.at(i);
	MC_LEPT_PHI[i]=leptonphi.at(i);
	MC_LEPT_THETA[i]=leptontheta.at(i);
        MC_LEPT_PDGID[i]=leptonpdgid.at(i);
      }
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      cout << endl;
    }
*/
    // Check if Z bosons are generated (specially for Z + jet MC)
    bool thereisZ=false;
    int i=0;
    double phi0=-100,eta0=-100;
    int nlep=0;
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      if ( abs(mcIter->pdgId())==25){
        phi0=mcIter->eta();   eta0=mcIter->phi();
       }
      if((abs(mcIter->pdgId())==11 || abs(mcIter->pdgId())==13)&&mcIter->status()==1){
       MC_LEPT_PT[nlep]=mcIter->pt();
       MC_LEPT_ETA[nlep]=mcIter->eta();
       MC_LEPT_PHI[nlep]=mcIter->phi();
    }
    }
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      // Select the Z decay: status 3
      if ( abs(mcIter->pdgId())==23/*&&mcIter->status()==22*/) {
	int j=0;
//        std::cout << "\n Found generated Z with PDGId = " << mcIter->pdgId() << " and status = "<< mcIter->status() << " and " << mcIter->numberOfDaughters() << " daughts |";
//        std::cout << " pt = "<<std::setw(12)<<mcIter->pt()<<" GeV/c\n";
        if(/*mcIter->status()==62*/mcIter->numberOfDaughters()==2){
//	std::cout << "   ===> Filled Z Bo at ["<<i<<"]["<<j<<"]\n";
	MC_Z_PT[i][0]=mcIter->pt();       MC_Z_ETA[i][0]=mcIter->eta();   MC_Z_PHI[i][0]=mcIter->phi();
        MC_Z_THETA[i][0]=mcIter->theta(); MC_Z_MASS[i][0]=mcIter->mass(); MC_Z_PDGID[i][0]=mcIter->pdgId();  
        thereisZ=true;
	// Check daughters of Z decay: status 3
	reco::GenParticle::daughters dst3 = mcIter->daughterRefVector();
	for (reco::GenParticle::daughters::const_iterator it_dst3 = dst3.begin(), est3 = dst3.end(); it_dst3 != est3; ++it_dst3) {
	  // Select status 3 electron, muon or tau
	  if (((abs((**it_dst3).pdgId()) == 11) || (abs((**it_dst3).pdgId()) == 13) || (abs((**it_dst3).pdgId()) == 15))) {
	    // daughters of the status 3 electron, muon or tau
		  ++j;
//		  std::cout<<"   ===> Filled E/MU at ["<<i<<"]["<<j<<"]";
		  MC_Z_PT[i][j]=(**it_dst3).pt();          MC_Z_ETA[i][j]=(**it_dst3).eta();      MC_Z_PHI[i][j]=(**it_dst3).phi();
		  MC_Z_THETA[i][j]=(**it_dst3).theta();    MC_Z_MASS[i][j]=(**it_dst3).mass();    MC_Z_PDGID[i][j]=(**it_dst3).pdgId();
		}
		  else {
//		    std::cout<<"   ===> NOT Filled, pt = "<<MC_Z_PT[i][j]<<" GeV/c";
		  }
                }
             i++;
            }
		else {} // other particle than e or mu or photon
        }
      //signal only
      if ( abs(mcIter->pdgId())==35&&mcIter->numberOfDaughters()==2){
       reco::GenParticle::daughters dst4 = mcIter->daughterRefVector();
       for (reco::GenParticle::daughters::const_iterator it_dst3 = dst4.begin(), est3 = dst4.end(); it_dst3 != est3; ++it_dst3) {
       cout << "daughter id=" << (**it_dst3).pdgId() << endl;
       if (abs((**it_dst3).pdgId()) == 23){
                  MC_Z_PT[0][3]=(**it_dst3).pt();          MC_Z_ETA[0][3]=(**it_dst3).eta();      MC_Z_PHI[0][3]=(**it_dst3).phi();
                  MC_Z_MASS[0][3]=(**it_dst3).mass();    MC_Z_PDGID[0][3]=(**it_dst3).pdgId();
        if(phi0!=-100&&eta0!=-100){
             double DELTAPHI;
             if(abs(phi0-MC_Z_PHI[0][3])<3.14159) DELTAPHI=abs(phi0-MC_Z_PHI[0][3]);
             else DELTAPHI=abs(phi0-MC_Z_PHI[0][3])-2*3.14159; 
             MC_Z_THETA[0][3]=sqrt( pow( DELTAPHI,2) + pow(eta0-MC_Z_ETA[0][3],2) );
        }
//        std::cout << "H->ZA Z pt = "<<(**it_dst3).pt()<<" GeV deltaR= "<< MC_Z_THETA[0][3] << endl;
       }}
      }
      if ( abs(mcIter->pdgId())==36&&mcIter->numberOfDaughters()==2){
       reco::GenParticle::daughters dst4 = mcIter->daughterRefVector();
       for (reco::GenParticle::daughters::const_iterator it_dst3 = dst4.begin(), est3 = dst4.end(); it_dst3 != est3; ++it_dst3) {
       if (abs((**it_dst3).pdgId()) == 23){
                  MC_Z_PT[0][4]=(**it_dst3).pt();          MC_Z_ETA[0][4]=(**it_dst3).eta();      MC_Z_PHI[0][4]=(**it_dst3).phi();
                  MC_Z_MASS[0][4]=(**it_dst3).mass();    MC_Z_PDGID[0][4]=(**it_dst3).pdgId();
        if(phi0!=-100&&eta0!=-100){
             double DELTAPHI;
             if(abs(phi0-MC_Z_PHI[0][4])<3.14159) DELTAPHI=abs(phi0-MC_Z_PHI[0][4]);
             else DELTAPHI=abs(phi0-MC_Z_PHI[0][4])-2*3.14159;
             MC_Z_THETA[0][4]=sqrt( pow( DELTAPHI,2) + pow(eta0-MC_Z_ETA[0][4],2) );
        }

//       std::cout << "A->ZA Z pt = "<<(**it_dst3).pt()<<" GeV deltaR= "<< MC_Z_THETA[0][4] << endl;
       }}
      }
    }
//    std::cout<<"Z number =" << i << "\n"<<std::endl;

//test part
   if(!thereisZ){
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      if ( abs(mcIter->pdgId())==1&&mcIter->status()==21) {
//        cout <<"fake Z status=" <<  mcIter->status() << endl;
//        cout <<"daughter N=" << mcIter->numberOfDaughters()<<endl;
//        std::cout << "   ===> Filled Z Bo at ["<<i<<"]["<<j<<"]";
        int j=0;
        MC_Z_PT[i][0]=mcIter->pt();       MC_Z_ETA[i][0]=mcIter->eta();   MC_Z_PHI[i][0]=mcIter->phi();
        MC_Z_THETA[i][0]=mcIter->theta(); MC_Z_MASS[i][0]=mcIter->mass(); MC_Z_PDGID[i][0]=mcIter->pdgId();
        reco::GenParticle::daughters dst3 = mcIter->daughterRefVector();
        for (reco::GenParticle::daughters::const_iterator it_dst3 = dst3.begin(), est3 = dst3.end(); it_dst3 != est3; ++it_dst3) {
          if (((abs((**it_dst3).pdgId()) == 11) || (abs((**it_dst3).pdgId()) == 13) || (abs((**it_dst3).pdgId()) == 15))) {
 //           cout << "daughter status=" << (**it_dst3).status() << endl;
//            std::cout<<"\n |--> Z daughter Particle  "<<std::setw(18)<<"| id = "<<std::setw(5)<<(**it_dst3).pdgId()<<" | st = "<<std::setw(5)<<(**it_dst3).status()<<" | pt =\n ";
                  ++j;
//                std::cout<<"   ===> Filled E/MU at ["<<i<<"]["<<j<<"]";
                  MC_Z_PT[i][j]=(**it_dst3).pt();          MC_Z_ETA[i][j]=(**it_dst3).eta();      MC_Z_PHI[i][j]=(**it_dst3).phi();
                  MC_Z_THETA[i][j]=(**it_dst3).theta();    MC_Z_MASS[i][j]=(**it_dst3).mass();    MC_Z_PDGID[i][j]=(**it_dst3).pdgId();
                }

      }
    }
   }
   }
    i=0; 
    
  }
  

  // MC Higgs 
  void fillmc(const edm::Event& iEvent){
/*

        // Pruned particles are the one containing "important" stuff
        edm::Handle<edm::View<reco::GenParticle> > pruned;
        iEvent.getByToken(prunedGenToken_,pruned);

        // Packed particles are all the status 1, so usable to remake jets
        // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
        edm::Handle<edm::View<pat::PackedGenParticle> > packed;
        iEvent.getByToken(packedGenToken_,packed);

        int index=0;
        for(size_t i=0; i<pruned->size();i++){
           if(abs((*pruned)[i].pdgId()) == 11 || abs((*pruned)[i].pdgId()) == 13){
               const Candidate * mclep = &(*pruned)[i];
               std::cout << "PdgID: " << mclep->pdgId() << " pt " << mclep->pt() << " eta: " << mclep->eta() << " phi: " << mclep->phi() << std::endl;
//if we did not return yet, then particle and ancestor are not relatives
                 MC_LEPT_PT[index]=mclep.pt();
                 MC_LEPT_ETA[index]=mclep.eta();
                 MC_LEPT_PHI[index]=mclep.phi();
                 MC_LEPT_PDGID[index]=mclep.pdgId();
                 index++;
                 if(index>=4) break;
                }
           }*/


/*
    edm::Handle<edm::View<Candidate> > Candidates;
    iEvent.getByToken(MCcollName, Candidates);
    cout << "running fillmc" << endl;
    cout << "size= " << Candidates->end()-Candidates->begin() << endl;
    for( edm::View<Candidate>::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
      cout << "Filling MC variables" << endl;
      const reco::Candidate& theParticle = (*cand);
      SetMCValues(theParticle,0);   
      int i=0,l=0;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	if (cand->daughter(j)->pdgId()==23){	  
	  const reco::Candidate& theParticle = (*cand->daughter(j));
	  SetMCValues(theParticle,j+1);
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	    if ( abs(cand->daughter(j)->daughter(k)->pdgId())==13 || 
		 abs(cand->daughter(j)->daughter(k)->pdgId())==15 || 
		 abs(cand->daughter(j)->daughter(k)->pdgId())==11){
	      l=i+j+k+3; 	      
	      const reco::Candidate& theParticle = (*cand->daughter(j)->daughter(k));
	      SetMCValues(theParticle,l);
	    }
	  }
	  i++;
	}
      }
    }    
    //fillMCCP(iEvent);*/
  }

  
  
  struct SortCandByDecreasingPt {
    bool operator()( const reco::Candidate &c1, const reco::Candidate &c2) const {
      return c1.pt() > c2.pt();
    }
  };
  

 /* 
  void fillAdditionalRECO(const edm::Event& iEvent){
   
    //Matching ZtoMuMu:
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZMuMu;
    iEvent.getByToken(goodZtoMuMuMCMatch_, GenParticlesMatchZMuMu);
    if (GenParticlesMatchZMuMu.isValid() ){
      cout << "The matched map collection has size= " <<  GenParticlesMatchZMuMu->size() << endl;
    }
    
    //Matching ZtoEE:
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZEE;
    iEvent.getByToken(goodZtoEEMCMatch_, GenParticlesMatchZEE);
    if (GenParticlesMatchZEE.isValid() ){
      cout << "The matched map collection has size= " <<  GenParticlesMatchZEE->size() << endl;
    }

    // di-leptons OS
    leptonscands_Z0->clear();
    leptonscands_Z1->clear();
   
    
    cout << "RECOcollNameZ size " << RECOcollNameZ.size() << endl;
    cout << "RECOcollNameZ string" << RECOcollNameZ.at(0) << endl;
    for (unsigned int i=0; i<RECOcollNameZ.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZ;
      if(i==0) iEvent.getByToken(zToMuMu, CandidatesZ); 
      else if(i==1) iEvent.getByToken(zToEE, CandidatesZ);
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZ->begin();cand != CandidatesZ->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_ZMM_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_ZMM_PT[0][kk]=cand->p4().pt();
	if (i==0) RECO_ZMM_ETA[0][kk]=cand->p4().eta();
	if (i==0) RECO_ZMM_PHI[0][kk]=cand->p4().phi();
	if (i==1) RECO_ZEE_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_ZEE_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_ZEE_ETA[0][kk]=cand->p4().eta();
	if (i==1) RECO_ZEE_PHI[0][kk]=cand->p4().phi();
	cout << "di-lepton candidate of type=" << RECOcollNameZ.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEE_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEE_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEE_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z1->push_back( cand->daughter(j)->clone());
	  }
	}
	
	if (i==0 && fillMCTruth==true){
	  // Matching ZtoMuMu
	  edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	  edm::Ref<std::vector<reco::GenParticle> > genrefZMuMu = (*GenParticlesMatchZMuMu)[Ref]; 
	  if (!genrefZMuMu.isNull()){
	    cout << "Gen Z with pT= " << genrefZMuMu->p4().pt() << " and mass="<< genrefZMuMu->p4().mass() << endl;	  
	    RECOzMuMu_MatchingMCTruth[kk]= true;
	    RECOzMuMu_MatchingMCpT[kk]= genrefZMuMu->p4().pt();
	    RECOzMuMu_MatchingMCmass[kk]= genrefZMuMu->p4().mass();
	    RECOzMuMu_MatchingMCEta[kk]= genrefZMuMu->p4().eta();
	    RECOzMuMu_MatchingMCPhi[kk]= genrefZMuMu->p4().phi();
	  }
	}
	if (i==1 && fillMCTruth==true){
	  // Matching ZtoEE
	  edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	  edm::Ref<std::vector<reco::GenParticle> > genrefZEE = (*GenParticlesMatchZEE)[Ref]; 
	  if (!genrefZEE.isNull()){
	    cout << "Gen Z with pT= " << genrefZEE->p4().pt() << " and mass="<< genrefZEE->p4().mass() << endl;	  
	    RECOzEE_MatchingMCTruth[kk]= true;
	    RECOzEE_MatchingMCpT[kk]= genrefZEE->p4().pt();
	    RECOzEE_MatchingMCmass[kk]= genrefZEE->p4().mass();
	    RECOzEE_MatchingMCEta[kk]= genrefZEE->p4().eta();
	    RECOzEE_MatchingMCPhi[kk]= genrefZEE->p4().phi();
	  }
	}


	kk++;
      }
    }

    // di-leptons SS and cross-leptons
    leptonscands_Zss0->clear();
    leptonscands_Zss1->clear();
    leptonscands_Zcross->clear();


    cout << "RECOcollNameZss size " << RECOcollNameZss.size() << endl; 
    for (unsigned int i=0; i<RECOcollNameZss.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZss;
      if(i==0)   iEvent.getByToken(zToMuMussmerge, CandidatesZss);
      if(i==1)   iEvent.getByToken(zToEEssmerge, CandidatesZss);
      if(i==2)   iEvent.getByToken(zToCrossLeptons, CandidatesZss); 
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZss->begin();cand != CandidatesZss->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0)  RECO_ZMMss_MASS[kk]=cand->p4().mass();
	if (i==0)  RECO_ZMMss_PT[0][kk]=cand->p4().pt();
	if (i==0)  RECO_ZMMss_ETA[0][kk]=cand->p4().eta();
	if (i==0)  RECO_ZMMss_PHI[0][kk]=cand->p4().phi();

	if (i==1)  RECO_ZEEss_MASS[kk]=cand->p4().mass();
	if (i==1)  RECO_ZEEss_PT[0][kk]=cand->p4().pt();
	if (i==1)  RECO_ZEEss_ETA[0][kk]=cand->p4().eta();
	if (i==1)  RECO_ZEEss_PHI[0][kk]=cand->p4().phi();

	if (i==2)  RECO_ZEM_MASS[kk]=cand->p4().mass();
	if (i==2)  RECO_ZEM_PT[0][kk]=cand->p4().pt();
	if (i==2)  RECO_ZEM_ETA[0][kk]=cand->p4().eta();
	if (i==2)  RECO_ZEM_PHI[0][kk]=cand->p4().phi();

	cout << "di-lepton candidate of type=" << RECOcollNameZss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMMss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMMss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMMss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEEss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEEss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEEss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss1->push_back( cand->daughter(j)->clone());
	  }
	  if (i==2) {
	    RECO_ZEM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zcross->push_back( cand->daughter(j)->clone());	  
	  }
	}
	kk++;
      }
    }


    // di-leptons ALL
    leptonscands_DiLep->clear();

    edm::Handle<edm::View<Candidate> > CandidatesDiLep;
    iEvent.getByToken(RECOcollNameDiLep_, CandidatesDiLep); 
    int kkk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesDiLep->begin();cand != CandidatesDiLep->end(); ++ cand ) { 
      if (kkk>49) continue;
      RECO_DiLep_MASS[kkk]=cand->p4().mass();
      RECO_DiLep_PT[0][kkk]=cand->p4().pt();
      RECO_DiLep_ETA[0][kkk]=cand->p4().eta();
      RECO_DiLep_PHI[0][kkk]=cand->p4().phi();
      
      cout << "di-lepton candidate of type=" << RECOcollNameDiLep.label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl;
	RECO_DiLep_PT[j+1][kkk]=cand->daughter(j)->p4().pt();
	RECO_DiLep_ETA[j+1][kkk]=cand->daughter(j)->p4().eta();
	RECO_DiLep_PHI[j+1][kkk]=cand->daughter(j)->p4().phi();
	leptonscands_DiLep->push_back( cand->daughter(j)->clone());
      }
      kkk++;
    }



    // MuMuMuMu
    int i=1; 
    leptonscands_MMMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMMM;
    iEvent.getByToken(RECOcollNameMMMM_, CandidatesMMMM);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHMMMM;
    iEvent.getByToken(goodHiggsTozzToMMMMMCMatch_, GenParticlesMatchHMMMM);
    if (GenParticlesMatchHMMMM.isValid()){
      cout << "The matched map collection has size= " <<  GenParticlesMatchHMMMM->size() << endl;
    }

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) { 
      if (kk>99) break;
      RECO_MMMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMMM_PT[i-1][kk]=cand->p4().pt();
      RECO_MMMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMMM_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  leptonscands_MMMM->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	} 
	l++;
      }


      // Matching HiggsToMMMM
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesMMMM,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzMMMM = (*GenParticlesMatchHMMMM)[Ref]; 
	if (!genrefHzzMMMM.isNull()){
	  cout << "Gen Z with pT= " << genrefHzzMMMM->p4().pt() << " and mass="<< genrefHzzMMMM->p4().mass() << endl;	  
	  RECOHzzMMMM_MatchingMCTruth[kk]= true;
	  RECOHzzMMMM_MatchingMCpT[kk]= genrefHzzMMMM->p4().pt();
	  RECOHzzMMMM_MatchingMCmass[kk]= genrefHzzMMMM->p4().mass();
	  RECOHzzMMMM_MatchingMCEta[kk]= genrefHzzMMMM->p4().eta();
	  RECOHzzMMMM_MatchingMCPhi[kk]= genrefHzzMMMM->p4().phi();
	}
      }

      kk++;
    }

    
    // EEEE
    i=1;
    leptonscands_EEEE->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEEE;
    iEvent.getByToken(RECOcollNameEEEE, CandidatesEEEE);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEEE;
    iEvent.getByToken(goodHiggsTozzToEEEEMCMatch_, GenParticlesMatchHEEEE);
    if (GenParticlesMatchHEEEE.isValid()){
      cout << "The matched map collection has size= " <<  GenParticlesMatchHEEEE->size() << endl;
    }

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) {
      if (kk>99) break;
      RECO_EEEE_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEEE_PT[i-1][kk]=cand->p4().pt();
      RECO_EEEE_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEEE_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEEE_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEEE_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEEE_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEEE_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEEE_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEEE_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEEE_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEEE_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  leptonscands_EEEE->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }

      // Matching HiggsToEEEE
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesEEEE,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzEEEE = (*GenParticlesMatchHEEEE)[Ref]; 
	if (!genrefHzzEEEE.isNull()){
	  cout << "Gen Z with pT= " << genrefHzzEEEE->p4().pt() << " and mass="<< genrefHzzEEEE->p4().mass() << endl;	  
	  RECOHzzEEEE_MatchingMCTruth[kk]= true;
	  RECOHzzEEEE_MatchingMCpT[kk]= genrefHzzEEEE->p4().pt();
	  RECOHzzEEEE_MatchingMCmass[kk]= genrefHzzEEEE->p4().mass();
	  RECOHzzEEEE_MatchingMCEta[kk]= genrefHzzEEEE->p4().eta();
	  RECOHzzEEEE_MatchingMCPhi[kk]= genrefHzzEEEE->p4().phi();
	}
      }


      kk++;
    }


    // EEMM
    i=1;
    leptonscands_EEMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEMM;
    iEvent.getByToken(RECOcollNameEEMM, CandidatesEEMM);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEMM;
    iEvent.getByToken(goodHiggsTozzToEEMMMCMatch_, GenParticlesMatchHEEMM);
    if (GenParticlesMatchHEEMM.isValid()){
      cout << "The matched map collection has size= " <<  GenParticlesMatchHEEMM->size() << endl;
    }

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) {
      if (kk>99) break;
      RECO_EEMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEMM_PT[i-1][kk]=cand->p4().pt();
      RECO_EEMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEMM_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();

	  leptonscands_EEMM->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }

      // Matching HiggsToEEMM
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesEEMM,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzEEMM = (*GenParticlesMatchHEEMM)[Ref]; 
	if (!genrefHzzEEMM.isNull()){
	  cout << "Gen Z with pT= " << genrefHzzEEMM->p4().pt() << " and mass="<< genrefHzzEEMM->p4().mass() << endl;	  
	  RECOHzzEEMM_MatchingMCTruth[kk]= true;
	  RECOHzzEEMM_MatchingMCpT[kk]= genrefHzzEEMM->p4().pt();
	  RECOHzzEEMM_MatchingMCmass[kk]= genrefHzzEEMM->p4().mass();
	  RECOHzzEEMM_MatchingMCEta[kk]= genrefHzzEEMM->p4().eta();
	  RECOHzzEEMM_MatchingMCPhi[kk]= genrefHzzEEMM->p4().phi();
	}
      }

      kk++;
    }
    
   
    // tri-leptons
    leptonscands_LLL0->clear();
    leptonscands_LLL1->clear();
    leptonscands_LLL2->clear();
    leptonscands_LLL3->clear();

    for (unsigned int i=0; i<RECOcollNameLLL.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesLLL;
      if(i==0) iEvent.getByToken(triLeptonsMuMuMu, CandidatesLLL); 
      if(i==1) iEvent.getByToken(triLeptonsMuMuE, CandidatesLLL);
      if(i==2) iEvent.getByToken(triLeptonsMuEE, CandidatesLLL);
      if(i==3) iEvent.getByToken(triLeptonsEEE, CandidatesLLL);
      int k=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLL->begin();cand != CandidatesLLL->end(); ++ cand ) { 
	if (k>49) break;
	if (i==0) RECO_LLL0_MASS[k]=cand->p4().mass();
	if (i==0) RECO_LLL0_PT[0][k]=cand->p4().pt();
	if (i==1) RECO_LLL1_MASS[k]=cand->p4().mass();
	if (i==1) RECO_LLL1_PT[0][k]=cand->p4().pt();
	if (i==2) RECO_LLL2_MASS[k]=cand->p4().mass();
	if (i==2) RECO_LLL2_PT[0][k]=cand->p4().pt();
	if (i==3) RECO_LLL3_MASS[k]=cand->p4().mass();
	if (i==4) RECO_LLL3_PT[0][k]=cand->p4().pt();
	cout << "Tri-lepton candidate of type=" << RECOcollNameLLL.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  if (i==0) {
	    RECO_LLL0_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_LLL1_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL1->push_back( cand->daughter(j)->clone());
	  }
	  if (i==2) {
	    RECO_LLL2_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL2->push_back( cand->daughter(j)->clone()); 
	  }
	  if (i==3) {
	    RECO_LLL3_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL3->push_back( cand->daughter(j)->clone());
	  }
	}
	k++;
      }
    }
    

    // 4-leptons SS
    leptonscands_LLLLss0->clear();
    leptonscands_LLLLss1->clear();
    leptonscands_LLLLss2->clear();
    
    for (unsigned int i=0; i<RECOcollNameLLLLss.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesLLLLss;
      if(i==0) iEvent.getByToken(quadLeptons4Mu, CandidatesLLLLss);
      if(i==1) iEvent.getByToken(quadLeptons2Mu2E, CandidatesLLLLss);
      if(i==2) iEvent.getByToken(quadLeptons4E, CandidatesLLLLss);
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLLLss->begin();cand != CandidatesLLLLss->end(); ++ cand ) { 
	if (kk>49) break;
	if (i==0) RECO_LLLL0ss_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_LLLL0ss_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_LLLL1ss_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_LLLL1ss_PT[0][kk]=cand->p4().pt();
	if (i==2) RECO_LLLL2ss_MASS[kk]=cand->p4().mass();
	if (i==2) RECO_LLLL2ss_PT[0][kk]=cand->p4().pt();
	cout << "4-lepton ss candidate of type=" << RECOcollNameLLLLss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	    if (i==0) {
	      if (j==0) RECO_LLLL0ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL0ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss0->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==1) {
	      if (j==0) RECO_LLLL1ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL1ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss1->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==2) {
	      if (j==0) RECO_LLLL2ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL2ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss2->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	  }
	}
	kk++;
      }      
    }


    // 4-leptons 3l+l
    leptonscands_LLLl0->clear();
    leptonscands_LLLl1->clear(); 
    
    for (unsigned int i=0; i<RECOcollNameLLLl.size(); i++){ 
      edm::Handle<edm::View<Candidate> > CandidatesLLLl; 
      if(i==0) iEvent.getByToken(quadLeptons3Mu1E, CandidatesLLLl);
      if(i==1) iEvent.getByToken(quadLeptons3E1Mu, CandidatesLLLl);  
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLLl->begin();cand != CandidatesLLLl->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_LLLl0_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_LLLl0_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_LLLl1_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_LLLl1_PT[0][kk]=cand->p4().pt();
	cout << "4-lepton from 3l+l candidate of type=" << RECOcollNameLLLl.at(i).label()
	     << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  //cout << "j= " << j << endl;
	  if ( i==0 && j==1) {
	    RECO_LLLl0_PT[j+3][kk]=cand->daughter(j)->p4().pt();
	    leptonscands_LLLl0->push_back( cand->daughter(j)->clone());
	  }
	  if ( i==1 && j==1) {
	    RECO_LLLl1_PT[j+3][kk]=cand->daughter(j)->p4().pt(); 
	    leptonscands_LLLl1->push_back( cand->daughter(j)->clone());
	  }
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) { 
	    //cout << "k= " << k << endl;
	    if (i==0) {
	      if (j==0) RECO_LLLl0_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLl0->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==1) {
	      if (j==0) RECO_LLLl1_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLl1->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    
	  } 
	} 
	kk++;
      }
    }


    // LLLL merged (no flavour, no charge)

    i=1;
    leptonscands_LLLL->clear();
    edm::Handle<edm::View<Candidate> > CandidatesLLLL;
    iEvent.getByToken(RECOcollNameLLLL, CandidatesLLLL);
    kk=0;    
    for( edm::View<Candidate>::const_iterator cand = CandidatesLLLL->begin();cand != CandidatesLLLL->end(); ++ cand ) {
      if (kk>49) break;
      cout << "4lepton (any flavour and charge) candidate of type= allLLLL"  << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      
      RECO_LLLL_MASS[0][kk]=cand->p4().mass();
      RECO_LLLL_PT[0][kk]=cand->p4().pt();
      RECO_LLLL_ETA[0][kk]=cand->p4().eta();
      RECO_LLLL_PHI[0][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_LLLL_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_LLLL_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_LLLL_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_LLLL_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_LLLL_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_LLLL_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_LLLL_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_LLLL_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  
	  leptonscands_LLLL->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }
      kk++;
    }
    
  }
*/

  void SetMCValues(const reco::Candidate& cand, int nMC){
    MC_E[nMC]     = cand.p4().energy();
    MC_PT[nMC]    = cand.p4().pt();
    MC_ETA[nMC]   = cand.p4().eta();
    MC_THETA[nMC] = cand.p4().theta();
    MC_PHI[nMC]   = cand.p4().phi();
    MC_MASS[nMC]  = cand.p4().mass();
    MC_PDGID[nMC] = cand.pdgId();
  }

 
 
  
  bool match(double mass, double pt, int charge, const reco::CandidateCollection *c1Coll){
    
    bool found=false;
    for( reco::CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++ pp ) {
      if ((fabs(pp->p4().mass()-mass)  <0.001 ) &&
	  (fabs(pp->p4().pt()  -pt)    <0.001 ) &&
	  (fabs(pp->charge()   -charge)<0.001 )  ){
	found=true;
	//std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
      }
    }
    return found;
  }
  
  bool matchParticle(double mass, double pt, int charge, const reco::Candidate *c1){
    
    bool found=false;
    if ((fabs(c1->p4().mass()-mass)  <0.001 ) &&
	(fabs(c1->p4().pt()  -pt)    <0.001 ) &&
	(fabs(c1->charge()   -charge)<0.001 )  ){
      found=true;
      //std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
    }    
    return found;
  }
  
  
  void fillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

    // Supercluster collection
    edm::Handle<reco::SuperClusterCollection> clusters;
    iEvent.getByLabel(clusterCollectionTag_,clusters);
    

    // EG isolation
    edm::Handle<edm::View<pat::Electron> > EleRefs;
    iEvent.getByToken(electronEgmTag_, EleRefs);   

    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (PV.size() != 0);
    //@// cout << "pvfound=" << pvfound << endl;
    if(pvfound){       
      for(reco::VertexCollection::const_iterator it=PV.begin() ; it!=PV.end() ; ++it){
	if(!it->isFake() && it->ndof() > 4 && fabs(it->position().z()) <= 24 && fabs(it->position().rho()) <= 2){ 	  	  
	  pVertex = math::XYZPoint(it->position().x(), it->position().y(), it->position().z());
	  //@//
	  /*cout << "P vertex position used for computing dxy and dz for electron and muons is (x,y,z)= " << 
	    pVertex.x() << " " <<
	    pVertex.y() << " " <<
	    pVertex.z() << endl;*/
	  break;
	}
      }
    } 
    
    int index=0;
    RECO_NELE=EleRefs->size();

    // Matching
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchEle;
    iEvent.getByToken(goodElectronMCMatch_, GenParticlesMatchEle);
    edm::Handle<reco::CandidateCollection > CollEle;
    iEvent.getByToken(myElectrons_, CollEle);
    
     for (edm::View<pat::Electron>::const_iterator cand = EleRefs->begin(); cand != EleRefs->end(); ++cand) {
      
     if(index>49) break;
     edm::Ref<edm::View<pat::Electron> > eletrackref(EleRefs,index);

      float corrEt = cand->et() * cand->userFloat("ecalTrkEnergyPostCorr")/cand->energy();
      auto corrP4  = cand->p4() * cand->userFloat("ecalTrkEnergyPostCorr")/cand->energy();
      cout<<"corrEt= "<<corrEt<<"corrP4 = "<<corrP4<<endl;
     
//      cout<<"TESTTTT Kinematics corrPT = "<<corrP4.pt()<<" corrE = "<<corrP4.energy()<<"corrEta = "<<corrP4.eta()<<"corr phi = "<<corrP4.phi()<<"corrTheta = "<<corrP4.theta()<<"corrcharge = "<<cand->charge()<<endl;

//      cout<<"TESTTTT Kinematics BEFORE corr. PT = "<< cand->p4().pt()<<" E = "<<cand->p4().energy()<<"Eta = "<<cand->p4().eta()<<" phi = "<<cand->p4().phi()<<"Theta = "<<cand->p4().theta()<<"charge = "<<cand->charge()<<endl;
 

      if (corrP4.pt() <7 )continue;
      if (corrP4.eta() >2.5)continue;
     
      RECOELE_E[index]       = corrP4.energy();
      RECOELE_PT[index]      = corrP4.pt();   
      RECOELE_P[index]=sqrt(corrP4.px()*corrP4.px()+corrP4.py()*corrP4.py()+corrP4.pz()*corrP4.pz());
      RECOELE_ETA[index]     = corrP4.eta();
      RECOELE_THETA[index]   = corrP4.theta();
      RECOELE_PHI[index]     = corrP4.phi();
      RECOELE_MASS[index]    = corrP4.mass();
      RECOELE_CHARGE[index]  = cand->charge();
      RECOELE_PT_uncorr[index]   = cand->pt();

	// 
      
      // cout<<"ecalTrkEnergyErrPostCorr = "<<cand->userFloat("ecalTrkEnergyErrPostCorr")<<endl;
      // cout<<"ecalTrkEnergyPreCorr"<<cand->userFloat("ecalTrkEnergyPreCorr")<<" ecalTrkEnergyErrPreCorr = "<<cand->userFloat("ecalTrkEnergyErrPreCorr");
   
      RECOELE_ecalTrkEnergyPreCorr[index]  = cand->userFloat("ecalTrkEnergyPreCorr");
      RECOELE_ecalTrkEnergyErrPreCorr[index]  = cand->userFloat("ecalTrkEnergyErrPreCorr");
      RECOELE_ecalTrkEnergyErrPostCorr[index]  = cand->userFloat("ecalTrkEnergyErrPostCorr");
      RECOELE_energyScaleValue[index]          = cand->userFloat("energyScaleValue");
      RECOELE_energySigmaValue[index]          = cand->userFloat("energySigmaValue");
      RECOELE_energyScaleUp[index]          = cand->userFloat("energyScaleUp");
      RECOELE_energyScaleDown[index]          = cand->userFloat("energyScaleDown");
      RECOELE_energyScaleStatUp[index]          = cand->userFloat("energyScaleStatUp");
      RECOELE_energyScaleStatDown[index]          = cand->userFloat("energyScaleStatDown");
      RECOELE_energyScaleSystUp[index]          = cand->userFloat("energyScaleSystUp");
      RECOELE_energyScaleSystDown[index]          = cand->userFloat("energyScaleSystDown");
      RECOELE_energyScaleGainUp[index]          = cand->userFloat("energyScaleGainUp");
      RECOELE_energyScaleGainDown[index]          = cand->userFloat("energyScaleGainDown");
      RECOELE_energyScaleEtUp[index]          = cand->userFloat("energyScaleEtUp");
      RECOELE_energyScaleEtDown[index]          = cand->userFloat("energyScaleEtDown");
      RECOELE_energySigmaUp[index]               = cand->userFloat("energySigmaUp");
      RECOELE_energySigmaDown[index]          = cand->userFloat("energySigmaDown");
      RECOELE_energySigmaPhiUp[index]          = cand->userFloat("energySigmaPhiUp");
      RECOELE_energySigmaPhiDown[index]          = cand->userFloat("energySigmaPhiDown");
      RECOELE_energySigmaRhoUp[index]          = cand->userFloat("energySigmaRhoUp");
      RECOELE_energySigmaRhoDown[index]          = cand->userFloat("energySigmaRhoDown");

      cout <<"ele test scaleup=" << RECOELE_energyScaleUp[index] <<" scale dow=" << RECOELE_energyScaleDown[index] << " sigmaUp=" << RECOELE_energySigmaUp[index] << "sigmaDow" << RECOELE_energySigmaDown[index] << " center_corr=" << RECOELE_E[index] << endl;
     // Global variables
      RECOELE_isEcalDriven[index]    = cand->ecalDrivenSeed();
      RECOELE_isTrackerDriven[index] = cand->trackerDrivenSeed();

     //@//
       std::cout << "\n Electron in the event: "
	 << "  isEcalDriven="    << RECOELE_isEcalDriven[index]  
	 << "  isTrackerDriven=" << RECOELE_isTrackerDriven[index] 
<< std::endl;

     cout << cand->electronID("mvaEleID-Fall17-iso-V2-wpHZZ") << " " << std::abs(cand->dB(pat::Electron::PV3D))/cand->edB(pat::Electron::PV3D) << endl;
     cout<<">>>>>MVA value " <<cand->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values")<<endl;

     RECOELE_ID[index] = cand->electronID("mvaEleID-Fall17-iso-V2-wpHZZ");
     RECOELE_mvaNonTrigV0[index] =  cand->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values");
     RECOELE_mvaTrigV0[index] =  cand->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values");
     cout<<"RECOELE_ID[index] = "<<RECOELE_ID[index]<<" RECOELE_mvaNonTrigV0[index] = "<< RECOELE_mvaNonTrigV0[index]<<endl;

     //@@@///  

     RECOELE_SIP[index]= std::abs(cand->dB(pat::Electron::PV3D))/cand->edB(pat::Electron::PV3D);
     RECOELE_IP[index]= std::abs(cand->dB(pat::Electron::PV3D));
     RECOELE_IPERROR[index]= cand->edB(pat::Electron::PV3D);

     // cout<<"ID = "<<RECOELE_ID[index]<<" SIP = "<<RECOELE_SIP[index]<<endl;
      
      // SuperCluster
      reco::SuperClusterRef sclRef = cand->superCluster();
      math::XYZPoint sclPos = cand->superClusterPosition();

      //if(!cand->ecalDrivenSeed() && cand->trackerDrivenSeed())
      //clRef = cand->pflowSuperCluster();
      double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
      double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
      ele_sclRawE[index]   = sclRef->rawEnergy() ;
      RECOELE_scl_E[index]   = sclRef->energy() ;
      RECOELE_scl_Et[index]  = sclRef->energy()*(Rt/R) ;
      RECOELE_scl_Eta[index] = sclRef->eta() ;
      RECOELE_scl_Phi[index] = sclRef->phi() ;
      ele_sclX[index]  = sclPos.X();
      ele_sclY[index] =  sclPos.Y();
      ele_sclZ[index] = sclPos.Z();

      TMatrixDSym bigCov;
      double dp = 0.;
      if (cand->ecalDriven()) {
	dp = cand->p4Error(pat::Electron::P4_COMBINATION);	
	RECOELE_PTError[index]=dp;
      }
      
      else {
	// Parametrization from Claude Charlot, 
	// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CJLST/ZZAnalysis/AnalysisStep/src/ZZMassErrors.cc?revision=1.2&view=markup
#if CMSSW_VERSION<500
	double ecalEnergy = cand->ecalEnergy() ;
#else
	double ecalEnergy = cand->correctedEcalEnergy() ;
#endif
	double err2 = 0.0;
	if (cand->isEB()) {
	  err2 += (5.24e-02*5.24e-02)/ecalEnergy;  
	  err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
	  err2 += 1.00e-02*1.00e-02;
	} else if (cand->isEE()) {
	  err2 += (1.46e-01*1.46e-01)/ecalEnergy;  
	  err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
	  err2 += 1.94e-03*1.94e-03;
	}
	dp = ecalEnergy * sqrt(err2);
	RECOELE_PTError[index]=dp;
      }
      // In order to produce a 3x3 matrix, we need a jacobian from (p) to (px,py,pz), i.e.
      //            [ Px/P  ]
      //  C_(3x3) = [ Py/P  ] * sigma^2(P) * [ Px/P Py/P Pz/P  ]
      //            [ Pz/P  ]
      AlgebraicMatrix31 ptop3;
      ptop3(0,0) = cand->px()/cand->p();
      ptop3(1,0) = cand->py()/cand->p();
      ptop3(2,0) = cand->pz()/cand->p();
      AlgebraicSymMatrix33 mat = ROOT::Math::Similarity(ptop3, AlgebraicSymMatrix11(dp*dp) );
      for (int i = 0; i < 3; ++i) { 
	for (int j = 0; j < 3; ++j) {
	  RECOELE_COV[index][i][j]=mat[i][j];
	}
      }
      
      cout<<"#########Electron step 2 "<<endl;
      cout<<"RECOELE_PTError[index]= "<<RECOELE_PTError[index]<<endl;
      cout<<"cand ecal energy = "<<cand->correctedEcalEnergy()<<endl;

      RECOELE_EGMECALISO[index]=(eletrackref->dr03EcalRecHitSumEt())/eletrackref->pt();
      RECOELE_EGMHCALISO[index]=(eletrackref->dr03HcalTowerSumEt())/eletrackref->pt();
      RECOELE_EGMX[index]=RECOELE_EGMTRACKISO[index]+RECOELE_EGMECALISO[index]+RECOELE_EGMHCALISO[index];
   
      //@//   
      /* std::cout << "--EG isolation: electron"
		<< "  X="     << RECOELE_EGMX[index]
		<< "  Track=" << RECOELE_EGMTRACKISO[index]
		<< "  Ecal= " << RECOELE_EGMECALISO[index]
		<< "  Hcal=" << RECOELE_EGMHCALISO[index]
		<< std::endl;*/


      // PF isolation  (R unknown, TODO)
      RECOELE_PFchHad[index]  = (cand->pfIsolationVariables().sumChargedHadronPt);
      RECOELE_PFneuHad[index] = (cand->pfIsolationVariables().sumNeutralHadronEt);
      RECOELE_PFphoton[index] = (cand->pfIsolationVariables().sumPhotonEt);
      RECOELE_PFPUchAllPart[index]= (cand->pfIsolationVariables().sumPUPt);
      RECOELE_PFchAllPart[index]= (cand->pfIsolationVariables().sumChargedParticlePt);

      //RECOELE_PFchAllPart[index]      =  (*isoPFChargedAllelemap)[eletrackref];
      //RECOELE_PFchHad[index]          =  (*isoPFChargedelemap)[eletrackref];
      //RECOELE_PFneuHad[index]         =  (*isoPFNeutralelemap)[eletrackref];
      //RECOELE_PFphoton[index]         =  (*isoPFGammaelemap)[eletrackref];
     // RECOELE_PFPUchAllPart[index]    =  (*isoPFPUelemap)[eletrackref];

      RECOELE_PFX_dB[index]           =  (RECOELE_PFchHad[index] + max(RECOELE_PFphoton[index]+RECOELE_PFneuHad[index]-0.5*RECOELE_PFPUchAllPart[index],0.))/eletrackref->pt();

      float EffectiveArea=0.;
      if (use2011EA){
	if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.18;
	if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.20;
	if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.15;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.19;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.21;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.22;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.4 )                                       EffectiveArea = 0.29;
      }
      else { 
	// 7_4_X use eta 
	//if (fabs(RECOELE_ETA[index]) >= 0.0   && fabs(RECOELE_ETA[index]) < 0.8 )   EffectiveArea = 0.1830;
	//if (fabs(RECOELE_ETA[index]) >= 0.8   && fabs(RECOELE_ETA[index]) < 1.3 )   EffectiveArea = 0.1734;
	//if (fabs(RECOELE_ETA[index]) >= 1.3   && fabs(RECOELE_ETA[index]) < 2.0 )   EffectiveArea = 0.1077;
	//if (fabs(RECOELE_ETA[index]) >= 2.0   && fabs(RECOELE_ETA[index]) < 2.2 )   EffectiveArea = 0.1565;
	//if (fabs(RECOELE_ETA[index]) >= 2.2 )                                       EffectiveArea = 0.2680;

	// 7_6_X use eta supercluster   
//	if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
//        if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
//        if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
        // 8_0_X  
        if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
        if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
        if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
      }

      cout<<"#########Electron step 3 "<<endl;

      //RECOELE_PFX_rho[index]=(RECOELE_PFchHad[index]+
      //		      max( (RECOELE_PFneuHad[index]+RECOELE_PFphoton[index]-max(RHO_ele,0.0)*EffectiveArea),0.0) )/cand->p4().pt(); 

      RECOELE_PFX_rho[index]= (cand->pfIsolationVariables().sumChargedHadronPt+
			       std::max(
					cand->pfIsolationVariables().sumPhotonEt+
					cand->pfIsolationVariables().sumNeutralHadronEt-
					max(RHO_ele,0.0)*EffectiveArea,
					0.0)
			       )/cand->p4().pt();

      RECOELE_PFchHad[index]  = (cand->pfIsolationVariables().sumChargedHadronPt);
      RECOELE_PFneuHad[index] = (cand->pfIsolationVariables().sumNeutralHadronEt); 
      RECOELE_PFphoton[index] = (cand->pfIsolationVariables().sumPhotonEt); 
      //RECOELE_PFsumPUPt[index]= (cand->pfIsolationVariables().sumPUPt); 
      RECOELE_PFPUchAllPart[index]=(cand->pfIsolationVariables().sumPUPt);   




      cout<<"#########Electron step 4 "<<endl;

      // GsfTrack
      RECOELE_gsftrack_NPixHits[index]  = cand->gsfTrack()->hitPattern().numberOfValidPixelHits();
      RECOELE_gsftrack_NStripHits[index]= cand->gsfTrack()->hitPattern().numberOfValidStripHits();
      RECOELE_gsftrack_chi2[index]      = cand->gsfTrack()->normalizedChi2();
      RECOELE_gsftrack_dxyB[index]      = cand->gsfTrack()->dxy(bs.position()) ;
      RECOELE_gsftrack_dxy[index]       = cand->gsfTrack()->dxy(pVertex);
      RECOELE_gsftrack_dxyError[index]  = cand->gsfTrack()->dxyError();
      RECOELE_gsftrack_dzB[index]       = cand->gsfTrack()->dz(bs.position());
      RECOELE_gsftrack_dz[index]        = cand->gsfTrack()->dz(pVertex);
      RECOELE_gsftrack_dzError[index]   = cand->gsfTrack()->dzError();

      //Conversion variables

      // Track-Cluster matching attributes
      RECOELE_ep[index]             = cand->eSuperClusterOverP(); 
      RECOELE_eSeedp[index]         = cand->eSeedClusterOverP();     
      RECOELE_eSeedpout[index]      = cand->eSeedClusterOverPout(); 
      RECOELE_eElepout[index]       = cand->eEleClusterOverPout();        
      
      RECOELE_deltaEtaIn [index]  = cand->deltaEtaSuperClusterTrackAtVtx(); 
      RECOELE_deltaEtaSeed[index] = cand->deltaEtaSeedClusterTrackAtCalo(); 
      RECOELE_deltaEtaEle[index]  = cand->deltaEtaEleClusterTrackAtCalo();  
      RECOELE_deltaPhiIn[index]   = cand->deltaPhiSuperClusterTrackAtVtx();
      RECOELE_deltaPhiSeed[index] = cand->deltaPhiSeedClusterTrackAtCalo(); 
      RECOELE_deltaPhiEle[index]  = cand->deltaPhiEleClusterTrackAtCalo() ;  

      cout<<"#########Electron step 5 "<<endl;

      //@// 
     /* std::cout << "--track-cluster matching: " 
	<< "  eSC_p="        << RECOELE_ep[index] 
	<< "  eSeed_p="      << RECOELE_eSeedp[index] 
	<< "  eSeed_pout="   << RECOELE_eSeedpout[index] 
	<< "  eEle_pout="    << RECOELE_eElepout[index] 
	<< "  deltaEtaIn="   << RECOELE_deltaEtaIn[index] 
        << "  deltaEtaSeed=" << RECOELE_deltaEtaSeed[index] 
        << "  deltaEtaEle="  << RECOELE_deltaEtaEle[index] 
        << "  deltaPhiIn="   << RECOELE_deltaPhiIn[index] 
        << "  deltaPhiSeed=" << RECOELE_deltaPhiSeed[index] 
	<< "  deltaPhiEle="  << RECOELE_deltaPhiEle[index] 
	<< std::endl;*/
      
       // Fiducial flags
      if (cand->isEB()) RECOELE_isbarrel[index] = 1 ; 
      else  RECOELE_isbarrel[index] = 0 ;
      if (cand->isEE()) RECOELE_isendcap[index] = 1 ; 
      else  RECOELE_isendcap[index] = 0 ;
      if (cand->isEBEtaGap())  RECOELE_isEBetaGap[index]  = 1 ;  
      if (cand->isEBPhiGap())  RECOELE_isEBphiGap[index]  = 1 ;  
      if (cand->isEEDeeGap())  RECOELE_isEEdeeGap[index]  = 1 ;  
      if (cand->isEERingGap()) RECOELE_isEEringGap[index] = 1 ;
      if (cand->isGap())       RECOELE_isGap[index] = 1 ;

      //@//
      /*  std::cout << "--fiducial flags: " 
	<< "  isEB="        << RECOELE_isbarrel[index] 
	<< "  isEBetagap="  << RECOELE_isEBetaGap[index] 
	<< "  isEBphigap="  << RECOELE_isEBphiGap[index] 
	<< "  isEE="        << RECOELE_isendcap[index] 
        << "  isEEdeegap="  << RECOELE_isEEdeeGap[index] 
	<< "  isEEringgap=" << RECOELE_isEEringGap[index] 
        << "  isGap="       << RECOELE_isGap[index]
	<< std::endl;*/
      

      // Shower shape
      RECOELE_sigmaIetaIeta[index] = cand->sigmaIetaIeta() ; 
      RECOELE_sigmaEtaEta[index]   = cand->sigmaEtaEta() ;
      RECOELE_e15[index]           = cand->e1x5() ;
      RECOELE_e25max[index]        = cand->e2x5Max() ;
      RECOELE_e55[index]           = cand->e5x5() ;
      RECOELE_he[index]            = cand->hadronicOverEm() ;
      // RECOELE_r9[index]            = cand->r9() ;

      //@//
      /*  std::cout << "--shower shape:" 
	<< "  sigmaIetaIeta=" << RECOELE_sigmaIetaIeta[index] 
	<< "  sigmaEtaEta=" << RECOELE_sigmaEtaEta[index]  
	<< "  e15=" << RECOELE_e15[index] 
	<< "  e25max=" << RECOELE_e25max[index] 
	<< "  e55=" << RECOELE_e55[index]  
	<< "  he=" << RECOELE_he[index] 
	<< "  r9=" << RECOELE_r9[index] 
	<< std::endl;*/
      
      // also variables on hcal
      // Isolation
      // Particle flow
      //RECOELE_mva[index] = cand->mva();
      //@//  std::cout << "--PF: mva = " << RECOELE_mva[index] << std::endl;

      // Brem & Classifaction
      RECOELE_fbrem[index]         = cand->fbrem() ;
      RECOELE_nbrems[index]        = cand->numberOfBrems();
      RECOELE_Class[index]         = cand->classification() ;
      if (useRECOformat) RECOELE_fbrem_mean[index]=1. - cand->gsfTrack()->outerMomentum().R()/cand->gsfTrack()->innerMomentum().R();
      RECOELE_fbrem_mode[index]=cand->fbrem();

      //@//
      /* std::cout << "--brem & classification: fbrem/nbrems/Class/fbrem_mean/fbrem_mode = "
		<< "  fbrem="      << RECOELE_fbrem[index]
		<< "  nbrems="     << RECOELE_nbrems[index]
        	<< "  class="      << RECOELE_Class[index]
		<< "  fbrem_mean=" << RECOELE_fbrem_mean[index]
		<< "  fbrem_mode=" << RECOELE_fbrem_mode[index]
		<< std::endl;*/

      // Corrections
      /*        RECOELE_isEcalEnergyCorrected[index] = cand->isEcalEnergyCorrected(); */
              RECOELE_ecalEnergy[index]            = cand->ecalEnergy(); 
      /*        RECOELE_ecalEnergyError[index]       = cand->ecalEnergyError(); */
      /*        RECOELE_isMomentumCorrected[index]   = cand->isMomentumCorrected(); */
      /*        RECOELE_trackMomentumError[index]    = cand->trackMomentumError(); */
      /*        RECOELE_electronMomentumError[index] = cand->electronMomentumError(); */
      
	      cout<<"#########Electron step 6 "<<endl;
    
      // Seed Collection
/*
      if (useRECOformat) {
	edm::RefToBase<TrajectorySeed> seed = cand->gsfTrack()->extra()->seedRef();
	reco::ElectronSeedRef MyS = seed.castTo<reco::ElectronSeedRef>();
	ele_seedSubdet2[index] = int(MyS->subDet2());
	if(fabs(MyS->dPhi2()) < 100.) ele_seedDphi2[index] = double(MyS->dPhi2());
	if(fabs(MyS->dRz2()) < 100.)  ele_seedDrz2[index]  = double(MyS->dRz2());
	
	ele_seedSubdet1[index] = int(MyS->subDet1());
	if(fabs(MyS->dPhi1()) < 100.) ele_seedDphi1[index] = double(MyS->dPhi1());
	if(fabs(MyS->dRz1()) < 100.)  ele_seedDrz1[index]  = double(MyS->dRz1());
      }
*/
      // Matching
      if (fillMCTruth==true){
      	    edm::Ref<std::vector<reco::GenParticle> > genrefEle = (*GenParticlesMatchEle)[eletrackref];
      	    if (!genrefEle.isNull()){
      	      cout << "GenElectron with pT= " << genrefEle->p4().pt() << " and mass="<< genrefEle->p4().mass()<< endl;
      	      RECOELE_MatchingMCTruth[index]= true;
      	      RECOELE_MatchingMCpT[index]= genrefEle->p4().pt();
      	      RECOELE_MatchingMCEta[index]= genrefEle->p4().eta();
      	      RECOELE_MatchingMCPhi[index]= genrefEle->p4().phi();
      	    }
      	    else {
      	      cout << "There is no a reference to a genElectron" << endl;
      	    }
         }//end match

      cout<<"RECOELE_MatchingMCTruth[index] = "<<RECOELE_MatchingMCTruth[index]<<endl;
      cout<<" RECOELE_MatchingMCpT[index] = "<< RECOELE_MatchingMCpT[index]<<endl;
      cout<<"RECOELE_MatchingMCEta[index]= "<<RECOELE_MatchingMCEta[index]<<endl;
      cout<<" RECOELE_MatchingMCPhi[index]= "<< RECOELE_MatchingMCPhi[index]<<endl;
      ///////

      index ++;
    } 
  }
  

  void fillMuons(const edm::Event& iEvent,const edm::EventSetup& iSetup){

    // Get the B-field
    //edm::ESHandle<MagneticField> magfield_;
    iSetup.get<IdealMagneticFieldRecord>().get( magfield_ );        


    // Muons
    edm::Handle<edm::View<pat::Muon> > MuCandidates;
    iEvent.getByToken(muonTag_, MuCandidates);

    edm::Handle<edm::ValueMap<float> > corrpterrormumap;
    iEvent.getByToken(muonCorrPtErrorMapTag_,corrpterrormumap);
    
    // Particle Flow Isolation
/*
    edm::Handle<edm::ValueMap<double> > isoPFChargedAllmumap;
    iEvent.getByToken(muonPFIsoValueChargedAllTag_, isoPFChargedAllmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFChargedmumap;
    iEvent.getByToken(muonPFIsoValueChargedTag_, isoPFChargedmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFNeutralmumap;
    iEvent.getByToken(muonPFIsoValueNeutralTag_, isoPFNeutralmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFGammamumap;
    iEvent.getByToken(muonPFIsoValueGammaTag_, isoPFGammamumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFPUmumap;
    iEvent.getByToken(muonPFIsoValuePUTag_, isoPFPUmumap);
*/

    // Vertexing
    // 3D
    edm::Handle<edm::View<pat::Muon> > VertMuCandidates;
    iEvent.getByLabel(muonTag_Vert, VertMuCandidates);
    
    // 3D w.r.t primary vertex DA
/*
    edm::Handle<edm::ValueMap<float> > vertexmumap;
    iEvent.getByToken(muonMapTag_Vert, vertexmumap);
    
    edm::Handle<edm::ValueMap<float> > vertexmumapvalue;
    iEvent.getByToken(muonMapTag_VertValue, vertexmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > vertexmumaperror;
    iEvent.getByToken(muonMapTag_VertError, vertexmumaperror);

    // 3D w.r.t primary vertex KF
    edm::Handle<edm::ValueMap<float> > vertexmumapKF;
//    iEvent.getByToken(muonMapTag_VertKF, vertexmumapKF);
    
    edm::Handle<edm::ValueMap<float> > vertexmumapvalueKF;
//    iEvent.getByToken(muonMapTag_VertValueKF, vertexmumapvalueKF);
    
    edm::Handle<edm::ValueMap<float> > vertexmumaperrorKF;
//    iEvent.getByToken(muonMapTag_VertErrorKF, vertexmumaperrorKF);


    // STIP SLIP
    edm::Handle<edm::ValueMap<float> > stipmumap;
//    iEvent.getByToken(muonSTIPMapTag_Vert, stipmumap);
    
    edm::Handle<edm::ValueMap<float> > slipmumap;
//    iEvent.getByToken(muonSLIPMapTag_Vert, slipmumap);
    
    edm::Handle<edm::ValueMap<float> > stipmumapvalue;
//    iEvent.getByToken(muonSTIPMapTag_VertValue, stipmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > slipmumapvalue;
//    iEvent.getByToken(muonSLIPMapTag_VertValue, slipmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > stipmumaperror;
//    iEvent.getByToken(muonSTIPMapTag_VertError, stipmumaperror);
    
    edm::Handle<edm::ValueMap<float> > slipmumaperror;
//    iEvent.getByToken(muonSLIPMapTag_VertError, slipmumaperror);
*/        
    // primary vertex
 
    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (PV.size() != 0);
    cout << "pvfound=" << pvfound << endl;
    if(pvfound){       
      for(reco::VertexCollection::const_iterator it=PV.begin() ; it!=PV.end() ; ++it){
	if(!it->isFake() && it->ndof() > 4 && fabs(it->position().z()) <= 24 && fabs(it->position().rho()) <= 2){ 	  	  
	  pVertex = math::XYZPoint(it->position().x(), it->position().y(), it->position().z());
	  cout << "P vertex position used for computing dxy and dz for electron and muons is (x,y,z)= " << 
	    pVertex.x() << " " <<
	    pVertex.y() << " " <<
	    pVertex.z() << endl;
	  break;
	}
      }
    }
    

    int indexbis=0;
    RECO_NMU=MuCandidates->size();
    cout << "Event has muon numbers of " << RECO_NMU << endl;

/*
    // Matching
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchMu;
    iEvent.getByToken(goodMuonMCMatch_, GenParticlesMatchMu);
    edm::Handle<reco::CandidateCollection > CollMu;
    iEvent.getByToken(myMuons_, CollMu);

    if (GenParticlesMatchMu.isValid()){
      cout << endl<< "Muons:"<<endl<<"The reco collection to be matched has size= " <<  CollMu->size() << endl;
      cout << "The matched map collection has size= " <<  GenParticlesMatchMu->size() << endl;
    }
*/
    //


    for (edm::View<pat::Muon>::const_iterator cand = MuCandidates->begin(); cand != MuCandidates->end(); ++cand) {
            
      if(indexbis>99) break;
      
      edm::Ref<edm::View<pat::Muon> > mutrackref(MuCandidates,indexbis); 
      edm::Ref<edm::View<pat::Muon> > mutrackrefv(VertMuCandidates,indexbis); 
      RECOMU_isPFMu[indexbis]=cand->isPFMuon();

      RECOMU_isGlobalMu[indexbis]=cand->isGlobalMuon();
      RECOMU_isStandAloneMu[indexbis]=cand->isStandAloneMuon();
      RECOMU_isTrackerMu[indexbis]=cand->isTrackerMuon();
      RECOMU_isCaloMu[indexbis]=cand->isCaloMuon();
      
     if(cand->globalTrack().isAvailable()) RECOMU_isTrackerHighPtMu[indexbis]=isTrackerHighPtMu(*cand,pVertex);   

      std::cout << "\n Muon in the event: "
	        <<   "  isPF=" << RECOMU_isPFMu[indexbis]
		<<   "  isGB=" << RECOMU_isGlobalMu[indexbis]
		<< "   isSTA=" << RECOMU_isStandAloneMu[indexbis]
		<<   "  isTM=" << RECOMU_isTrackerMu[indexbis]
		<< "  isCalo=" << RECOMU_isCaloMu[indexbis]
		<< std::endl;

      // Kinematic of muon
      RECOMU_E[indexbis]=cand->p4().energy();
      RECOMU_PT[indexbis]=cand->p4().pt();
      RECOMU_P[indexbis]=sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz());
      RECOMU_ETA[indexbis]=cand->p4().eta();
      RECOMU_THETA[indexbis]=cand->p4().theta();
      RECOMU_PHI[indexbis]=cand->p4().phi();
      RECOMU_MASS[indexbis]=cand->p4().mass();
      RECOMU_CHARGE[indexbis]=cand->charge();

      std::cout << "--kinematic:"
		<< "  pT="     << RECOMU_PT[indexbis]
	        << "  E="      << RECOMU_E[indexbis]
		<< "  p="      << RECOMU_P[indexbis]
		<< "  eta="    << RECOMU_ETA[indexbis]
		<< "  theta="  << RECOMU_THETA[indexbis]
		<< "  phi="    << RECOMU_PHI[indexbis]
		<< "  mass="   << RECOMU_MASS[indexbis]
		<< "  charge=" << RECOMU_CHARGE[indexbis]
		<< std::endl;
	
      // Covariance matrix
      
      if(cand->innerTrack().isAvailable()){

	GlobalTrajectoryParameters gp(GlobalPoint(cand->innerTrack()->vx(), cand->innerTrack()->vy(),  cand->innerTrack()->vz()),
				      GlobalVector(cand->innerTrack()->px(),cand->innerTrack()->py(),cand->innerTrack()->pz()),
				      cand->innerTrack()->charge(),
				      magfield_.product());
	JacobianCurvilinearToCartesian curv2cart(gp);
	CartesianTrajectoryError cartErr= ROOT::Math::Similarity(curv2cart.jacobian(), cand->innerTrack()->covariance());
	const AlgebraicSymMatrix66 mat = cartErr.matrix();
	for (int i = 0; i < 3; ++i) { 
	  for (int j = 0; j < 3; ++j) { 
	    //bigCov(i,j) = mat[i+3][j+3]; 
	    //	  cout << "index= " << indexbis << "i= " << i << "j= " << j <<endl;
	    RECOMU_COV[indexbis][i][j]=mat[i+3][j+3]; 
	  }  
	}  
	
      }

      // Isolation
      //      RECOMU_TRACKISO[indexbis]=(*isoTkmumap)[mutrackref]/cand->p4().pt();
      RECOMU_TRACKISO_SUMPT[indexbis]=(cand->isolationR03().sumPt)/cand->p4().pt();
      //RECOMU_ECALISO[indexbis]=(*isoEcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_HCALISO[indexbis]=(*isoHcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_X[indexbis]=[(*isoTkmumap)[mutrackref]+(*isoEcalmumap)[mutrackref]+(*isoHcalmumap)[mutrackref]]/cand->p4().pt();
      // temporary solution for reducedrechit problem
      RECOMU_ECALISO[indexbis]=(cand->isolationR03().emEt)/cand->p4().pt();
      RECOMU_HCALISO[indexbis]=(cand->isolationR03().hadEt)/cand->p4().pt();
      RECOMU_X[indexbis]=RECOMU_TRACKISO[indexbis]+RECOMU_ECALISO[indexbis]+RECOMU_HCALISO[indexbis];

      RECOMU_PFchHad[indexbis]  = (cand->pfIsolationR03().sumChargedHadronPt);
      RECOMU_PFneuHad[indexbis] = (cand->pfIsolationR03().sumNeutralHadronEt);
      RECOMU_PFphoton[indexbis] = (cand->pfIsolationR03().sumPhotonEt);
      RECOMU_PFPUchAllPart[indexbis]= (cand->pfIsolationR03().sumPUPt);

//      RECOMU_PFchHad[indexbis]          =  (*isoPFChargedmumap)[mutrackref];
//      RECOMU_PFneuHad[indexbis]         =  (*isoPFNeutralmumap)[mutrackref];
//      RECOMU_PFphoton[indexbis]         =  (*isoPFGammamumap)[mutrackref];
//      RECOMU_PFPUchAllPart[indexbis]    =  (*isoPFPUmumap)[mutrackref];

      RECOMU_PFX_dB[indexbis]=(RECOMU_PFchHad[indexbis]+max(0.,RECOMU_PFneuHad[indexbis]+RECOMU_PFphoton[indexbis]-0.5*RECOMU_PFPUchAllPart[indexbis]))/cand->p4().pt();

      float EffectiveArea=0.;
      if (use2011EA){
	if (fabs(RECOMU_ETA[indexbis]) >= 0.0 && fabs(RECOMU_ETA[indexbis]) < 1.0 ) EffectiveArea = 0.132;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.0 && fabs(RECOMU_ETA[indexbis]) < 1.5 ) EffectiveArea = 0.120;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.5 && fabs(RECOMU_ETA[indexbis]) < 2.0 ) EffectiveArea = 0.114;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.0 && fabs(RECOMU_ETA[indexbis]) < 2.2 ) EffectiveArea = 0.139;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.2 && fabs(RECOMU_ETA[indexbis]) < 2.3 ) EffectiveArea = 0.168;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.3 )                                     EffectiveArea = 0.189;
      }
      else {       
	if (fabs(RECOMU_ETA[indexbis]) >= 0.0 && fabs(RECOMU_ETA[indexbis]) < 1.0 ) EffectiveArea = 0.674;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.0 && fabs(RECOMU_ETA[indexbis]) < 1.5 ) EffectiveArea = 0.565;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.5 && fabs(RECOMU_ETA[indexbis]) < 2.0 ) EffectiveArea = 0.442;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.0 && fabs(RECOMU_ETA[indexbis]) < 2.2 ) EffectiveArea = 0.515;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.2 && fabs(RECOMU_ETA[indexbis]) < 2.3 ) EffectiveArea = 0.821;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.3 )                                     EffectiveArea = 0.660;
      }
       
      RECOMU_PFX_rho[indexbis]=(RECOMU_PFchHad[indexbis]+max( (RECOMU_PFneuHad[indexbis]+RECOMU_PFphoton[indexbis]-max(RHO_mu,0.0)*(EffectiveArea)),0.0) 
)/double(cand->p4().pt());


      std::cout << "--isolation: muon"
		<< "  X="       << RECOMU_X[indexbis]
		<< "  Tracker=" << RECOMU_TRACKISO[indexbis] <<  "  no veto=" << RECOMU_TRACKISO_SUMPT[indexbis]
		<< "  Ecal="    << RECOMU_ECALISO[indexbis]
		<< "  Hcal="    << RECOMU_HCALISO[indexbis]
	        << endl;

      std::cout << "--PF isolation"
	        << "  chHad="   << RECOMU_PFchHad[indexbis] 
	        << "  neuHad="  << RECOMU_PFneuHad[indexbis]
	        << "  photon="  << RECOMU_PFphoton[indexbis] 
	        << "  sumPUPt=" << RECOMU_PFPUchAllPart[indexbis]
                << "  PFX_dB="  << RECOMU_PFX_dB[indexbis] 
                << "  PFX_rho=" << RECOMU_PFX_rho[indexbis]
		<< std::endl;

      // Vertexing
/*
      RECOMU_SIP[indexbis]=(*vertexmumap)[mutrackrefv];
      RECOMU_IP[indexbis]=(*vertexmumapvalue)[mutrackrefv];
      RECOMU_IPERROR[indexbis]=(*vertexmumaperror)[mutrackrefv];
*/
      
      // Other properties
      RECOMU_numberOfMatches[indexbis]=cand->numberOfMatches();
      RECOMU_numberOfMatchedStations[indexbis]=cand->numberOfMatchedStations();
      RECOMU_caloCompatibility[indexbis]=cand->caloCompatibility();
      RECOMU_segmentCompatibility[indexbis]=(muon::segmentCompatibility( (*cand)));
      RECOMU_glbmuPromptTight[indexbis]=(muon::isGoodMuon( (*cand),muon::GlobalMuonPromptTight));
//      RECOMU_chi2LocalPosition[indexbis]=cand.combinedQuality().chi2LocalPosition;
//      RECOMU_trkKink[indexbis]=cand->combinedQuality()->trkKink;
      RECOMU_isMedium[indexbis] = muon::isMediumMuon(*cand);

      std::cout	<< "--other properties:"
		<< "  n.matches="   << RECOMU_numberOfMatches[indexbis]
		<< "  caloComp="    << RECOMU_caloCompatibility[indexbis]
		<< "  segmentComp=" << RECOMU_segmentCompatibility[indexbis]
	        << "  glbmuPromptTight=" << RECOMU_glbmuPromptTight[indexbis]
                << "  isMediumMuon= " << RECOMU_isMedium[indexbis]
//                << "  chi2LocalPosition=" << RECOMU_chi2LocalPosition[indexbis]
//                << "  trkKink=" << RECOMU_trkKink[indexbis]
		<< std::endl;

  
      // Track properties
      if(cand->muonBestTrack().isAvailable()){
	RECOMU_mubesttrkType[indexbis]=cand->muonBestTrackType();
	RECOMU_mubesttrkDxy[indexbis]=cand->muonBestTrack()->dxy(pVertex);
	RECOMU_mubesttrkDxyB[indexbis]=cand->muonBestTrack()->dxy(bs.position());
	RECOMU_mubesttrkDxyError[indexbis]=cand->muonBestTrack()->dxyError();
	RECOMU_mubesttrkDz[indexbis]=cand->muonBestTrack()->dz(pVertex);
	RECOMU_mubesttrkDzB[indexbis]=cand->muonBestTrack()->dz(bs.position());
	RECOMU_mubesttrkDzError[indexbis]=cand->muonBestTrack()->dzError();
	//RECOMU_mubesttrkPTError[indexbis]=(*corrpterrormumap)[mutrackref];;
      }

      if(cand->globalTrack().isAvailable()){
	RECOMU_mutrkPT[indexbis]=cand->globalTrack()->pt();
	RECOMU_mutrkPTError[indexbis]=cand->globalTrack()->ptError();
	//RECOMU_mutrkPTError[indexbis]=cand->bestTrack()->ptError(); // bestTrack
	RECOMU_mutrkDxy[indexbis]=cand->globalTrack()->dxy(pVertex);
	RECOMU_mutrkDxyError[indexbis]=cand->globalTrack()->dxyError();
	RECOMU_mutrkDxyB[indexbis]=cand->globalTrack()->dxy(bs.position()) ;
	RECOMU_mutrkDz[indexbis]=cand->globalTrack()->dz(pVertex);
	RECOMU_mutrkDzError[indexbis]=cand->globalTrack()->dzError();
	RECOMU_mutrkDzB[indexbis]=cand->globalTrack()->dz(bs.position());
	RECOMU_mutrkChi2PerNdof[indexbis]=cand->globalTrack()->normalizedChi2();
	RECOMU_mutrkCharge[indexbis]=cand->globalTrack()->charge();
 	RECOMU_mutrkNHits[indexbis]=cand->globalTrack()->numberOfValidHits(); 
	RECOMU_mutrkNPixHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_mutrkNStripHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidStripHits();
	RECOMU_mutrkNMuonHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidMuonHits(); 
	RECOMU_mutrktrackerLayersWithMeasurement[indexbis]=cand->globalTrack()->hitPattern().trackerLayersWithMeasurement(); 

	RECOMU_muInnertrkDxy[indexbis]=cand->innerTrack()->dxy(pVertex);
	RECOMU_muInnertrkDxyError[indexbis]=cand->innerTrack()->dxyError();
	RECOMU_muInnertrkDxyB[indexbis]=cand->innerTrack()->dxy(bs.position()) ;
	RECOMU_muInnertrkDz[indexbis]=cand->innerTrack()->dz(pVertex);
	RECOMU_muInnertrkDzError[indexbis]=cand->innerTrack()->dzError();
	RECOMU_muInnertrkDzB[indexbis]=cand->innerTrack()->dz(bs.position());
	RECOMU_muInnertrkChi2PerNdof[indexbis]=cand->innerTrack()->normalizedChi2();
	RECOMU_muInnertrktrackerLayersWithMeasurement[indexbis]=cand->innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
	RECOMU_muInnertrkPT[indexbis]=cand->innerTrack()->pt();	
	//RECOMU_muInnertrkPTError[indexbis]=cand->innerTrack()->ptError();
	RECOMU_muInnertrkPTError[indexbis]=cand->bestTrack()->ptError(); // Besttrack

	RECOMU_muInnertrkCharge[indexbis]=cand->innerTrack()->charge();
 	RECOMU_muInnertrkNHits[indexbis]=cand->innerTrack()->numberOfValidHits(); 
	RECOMU_muInnertrkNPixHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_muInnertrkNStripHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidStripHits();
        RECOMU_muInnertrkvalidFraction[indexbis]=cand->innerTrack()->validFraction();
      }
      else if(cand->innerTrack().isAvailable()){
	RECOMU_muInnertrkDxy[indexbis]=cand->innerTrack()->dxy(pVertex);
	RECOMU_muInnertrkDxyError[indexbis]=cand->innerTrack()->dxyError();
	RECOMU_muInnertrkDxyB[indexbis]=cand->innerTrack()->dxy(bs.position()) ;
	RECOMU_muInnertrkDz[indexbis]=cand->innerTrack()->dz(pVertex);
	RECOMU_muInnertrkDzError[indexbis]=cand->innerTrack()->dzError();
	RECOMU_muInnertrkDzB[indexbis]=cand->innerTrack()->dz(bs.position());
	RECOMU_muInnertrkChi2PerNdof[indexbis]=cand->innerTrack()->normalizedChi2();
	RECOMU_muInnertrktrackerLayersWithMeasurement[indexbis]=cand->innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
	RECOMU_muInnertrkPT[indexbis]=cand->innerTrack()->pt();	
	//RECOMU_muInnertrkPTError[indexbis]=cand->innerTrack()->ptError();
	RECOMU_muInnertrkPTError[indexbis]=cand->bestTrack()->ptError(); // Besttrack
	RECOMU_muInnertrkCharge[indexbis]=cand->innerTrack()->charge();
 	RECOMU_muInnertrkNHits[indexbis]=cand->innerTrack()->numberOfValidHits(); 
	RECOMU_muInnertrkNPixHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_muInnertrkNStripHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidStripHits();
        RECOMU_muInnertrkvalidFraction[indexbis]=cand->innerTrack()->validFraction();
      }

      if(cand->globalTrack().isAvailable() || cand->innerTrack().isAvailable() ){
	std::cout << "--muon track properties: "
	          << "  pt="        << RECOMU_mutrkPT[indexbis] 
	          << "  bestTrackType= " << RECOMU_mubesttrkType[indexbis]
		  << "  dxy="       << RECOMU_mubesttrkDxy[indexbis]
		  << "  dxyError="  << RECOMU_mubesttrkDxyError[indexbis]
		  << "  dxyB="      << RECOMU_mubesttrkDxyB[indexbis] 
		  << "  dz="        << RECOMU_mubesttrkDz[indexbis] 
		  << "  dzError="   << RECOMU_mubesttrkDzError[indexbis]
		  << "  dzB="       << RECOMU_mubesttrkDzB[indexbis]
		  << "  PtError="   << RECOMU_mubesttrkPTError[indexbis]
		  << "  chi2_nodf=" << RECOMU_mutrkChi2PerNdof[indexbis]
		  << "  charge="    << RECOMU_mutrkCharge[indexbis]
		  << "  nhits="     << RECOMU_mutrkNHits[indexbis] 
		  << "  nPixhits="   << RECOMU_mutrkNPixHits[indexbis]
	          << "  nStriphits=" << RECOMU_mutrkNStripHits[indexbis]
	          << "  nMuonhits="  << RECOMU_mutrkNMuonHits[indexbis]
		  << std::endl;
	
	
	// Tracker muon properties
	// RECOMU_trkmuArbitration[indexbis]=(muon::isGoodMuon( (*cand),muon::TrackerMuonArbitrated));
	RECOMU_trkmuArbitration[indexbis]=(muon::segmentCompatibility((*cand),pat::Muon::SegmentAndTrackArbitration));
	RECOMU_trkmu2DCompatibilityLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityLoose));
	RECOMU_trkmu2DCompatibilityTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityTight));
	RECOMU_trkmuOneStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationLoose));
	RECOMU_trkmuOneStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationTight));
	RECOMU_trkmuLastStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationLoose));
	RECOMU_trkmuLastStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationTight));
	RECOMU_trkmuOneStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngLoose));
	RECOMU_trkmuOneStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngTight));
	RECOMU_trkmuLastStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngLoose));
	RECOMU_trkmuLastStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngTight));
	RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtLoose));
	RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtTight));
	
	std::cout << "--tracker muon properties:"
		  << "  arbitration="              << RECOMU_trkmuArbitration[indexbis]
		  << "  2DCompLoose="              << RECOMU_trkmu2DCompatibilityLoose[indexbis]
		  << "  2DCompTight="              << RECOMU_trkmu2DCompatibilityTight[indexbis]
		  << "  1StationLoose="            << RECOMU_trkmuOneStationLoose[indexbis]
		  << "  1StationTight="            << RECOMU_trkmuOneStationTight[indexbis]
		  << "  LastStationLoose="         << RECOMU_trkmuLastStationLoose[indexbis]
		  << "  LastStationTight="         << RECOMU_trkmuLastStationTight[indexbis]
		  << "  1StationAngLoose="         << RECOMU_trkmuOneStationAngLoose[indexbis]
		  << "  1StationAngTight="         << RECOMU_trkmuOneStationAngTight[indexbis]
		  << "  LastStationAngLoose="      << RECOMU_trkmuLastStationAngLoose[indexbis]
		  << "  LastStationAngTight="      << RECOMU_trkmuLastStationAngTight[indexbis]
		  << "  LastStationOptLowptLoose=" << RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]
		  << "  LastStationOptLowptTight=" << RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]
		  << std::endl;
      }
  
      if (fillMCTruth==true){
        int k=-1;
          //cout << "Reco Electron with pT= " << hIter->pt() << " " << RECOELE_PT[index] << " " << fabs(hIter->pt()-RECOELE_PT[index]) << " and mass="<< hIter->mass()<< endl;
        edm::Handle<edm::View<reco::GenParticle> > pruned;
        iEvent.getByToken(prunedGenToken_,pruned);
             
        double dRmin=10;
        for(size_t j=0; j<pruned->size();j++){
           if(abs((*pruned)[j].pdgId()) == 13&&(*pruned)[j].isPromptFinalState()){
             const Candidate * genmu = &(*pruned)[j];
             double phi1 = genmu->p4().phi();
             double phi2 = cand->p4().phi();
             double eta1 = genmu->p4().eta();
             double eta2 = cand->p4().eta();
   //          double pt1 = genmu->p4().pt();
   //          double pt2 = cand->p4().pt();
             double DELTAPHI;
             if(abs(phi1-phi2)<3.14159) DELTAPHI=abs(phi1-phi2);
             else DELTAPHI=abs(phi1-phi2)-2*3.14159;
             double deltaR = sqrt( pow( DELTAPHI,2) + pow(eta1-eta2,2) );
             if(deltaR<dRmin) {k=j;dRmin=deltaR;}
           }
         }
        if(k>=0&&dRmin<0.15){
              RECOMU_MatchingMCTruth[indexbis]= true;
              RECOMU_MatchingMCpT[indexbis]= (*pruned)[k].p4().pt();
              RECOMU_MatchingMCEta[indexbis]= (*pruned)[k].p4().eta();
              RECOMU_MatchingMCPhi[indexbis]= (*pruned)[k].p4().phi();
        }
      }
 
      indexbis++;
    }
  }

  void fillPhotons(const edm::Event& iEvent){
    // Photons
    //edm::Handle<edm::View<reco::Candidate> > photons;
    edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByToken(photonsTag_, photons);
    RECO_NPHOT=photons->size();
/*
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchPhot;
    iEvent.getByToken(goodGammaMCMatch_, GenParticlesMatchPhot);
    edm::Handle<reco::CandidateCollection > CollPhot;
    iEvent.getByToken(myGammas_, CollPhot);
    bool ismyGammas=false;
    if (CollPhot.isValid()) ismyGammas=true;
*/    
    int iphot=0;
    for (edm::View<pat::Photon>::const_iterator cand = photons->begin(); cand != photons->end(); ++cand) {
      if (iphot>19) break;
      RECOPHOT_PT[iphot]=cand->pt();
      RECOPHOT_ETA[iphot]=cand->eta();
      RECOPHOT_PHI[iphot]=cand->phi();
      RECOPHOT_THETA[iphot]=cand->theta();
/*
      if (!cand->parentSuperCluster().isNull()){
        RECOPHOT_TLE_ParentSC_X[iphot]=cand->parentSuperCluster()->x();
        RECOPHOT_TLE_ParentSC_Y[iphot]=cand->parentSuperCluster()->y();
        RECOPHOT_TLE_ParentSC_Z[iphot]=cand->parentSuperCluster()->z();
      }

      cout << "Parent SC x,y,x= " <<RECOPHOT_TLE_ParentSC_X[iphot] << " " << RECOPHOT_TLE_ParentSC_Y[iphot] << " " << RECOPHOT_TLE_ParentSC_Z[iphot] << endl;
*/

      cout << "Reco Photon with pT= " << RECOPHOT_PT[iphot] << " eta= " << RECOPHOT_ETA[iphot] << " phi= " << RECOPHOT_PHI[iphot] << endl;
/*
      // Matching
      int i=0;
      if (ismyGammas){
	for ( reco::CandidateCollection::const_iterator hIter=CollPhot->begin(); hIter!= CollPhot->end(); ++hIter ){
	  //cout << "Reco Photon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl;
	  if (fabs(hIter->pt()-RECOPHOT_PT[iphot])<0.01){
	    i=hIter-(CollPhot->begin());
	    CandidateRef Ref( CollPhot, i );
	    edm::Ref<std::vector<reco::GenParticle> > genrefPhot = (*GenParticlesMatchPhot)[Ref];
	    if (!genrefPhot.isNull()){
	      cout << "GenMuon with pT= " << genrefPhot->p4().pt() << " and mass="<< genrefPhot->p4().mass()<< endl;
	      RECOPHOT_MatchingMCTruth[i]= true;
	      RECOPHOT_MatchingMCpT[i]= genrefPhot->p4().pt();
	      RECOPHOT_MatchingMCEta[i]= genrefPhot->p4().eta();
	      RECOPHOT_MatchingMCPhi[i]= genrefPhot->p4().phi();	    
	    } 
	  }   
	}
      }  */ 
      iphot++;
    }

    // PF ISR photon

    // Photons
//    edm::Handle<edm::View<pat::PackedCandidate> > pfphotons;
//    iEvent.getByToken(pfphotonsTag_, pfphotons);

/*
    edm::Handle<edm::ValueMap<double> > isoPFChargedAllphmap; 
    iEvent.getByToken(photonPFIsoValueChargedAllTag_, isoPFChargedAllphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFChargedphmap; 
    iEvent.getByToken(photonPFIsoValueChargedTag_, isoPFChargedphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFNeutralphmap; 
    iEvent.getByToken(photonPFIsoValueNeutralTag_, isoPFNeutralphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFGammaphmap; 
    iEvent.getByToken(photonPFIsoValueGammaTag_, isoPFGammaphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFPUphmap; 
    iEvent.getByToken(photonPFIsoValuePUTag_, isoPFPUphmap); 
*/

    iphot=0;
    RECO_NPFPHOT=photons->size();
    cout << "There are " << photons->size() << " photons" << endl;
    for (edm::View<pat::Photon>::const_iterator cand = photons->begin(); cand != photons->end(); ++cand) {

      if (iphot>19) break;

      edm::Ref<edm::View<pat::Photon> > phtrackref(photons,iphot); 

      RECOPFPHOT_PT[iphot]=cand->pt();

      // error on pt of photon taken from #include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"
      double C;
      double S;
      double N;
      if(TMath::Abs(cand->eta())<1.48){C=0.35/100; S=5.51/100; N=98./1000.;}
      else{C=0; S=12.8/100; N=440./1000.;} 
      double  perr = TMath::Sqrt(C*C*cand->p4().e()*cand->p4().e() + S*S*cand->p4().e() + N*N);
      double pterr = perr*cand->pt()/cand->p(); 
      RECOPFPHOT_PTError[iphot]=float(pterr);
      
      cout << "PF Photon : pT= " << RECOPFPHOT_PT[iphot] << " pTerr= " << RECOPFPHOT_PTError[iphot]<< endl;

      RECOPFPHOT_ETA[iphot]=cand->eta();
      RECOPFPHOT_PHI[iphot]=cand->phi();
      RECOPFPHOT_THETA[iphot]=cand->theta();
/*      
      RECOPFPHOT_PFchAllPart[iphot]      =  (*isoPFChargedAllphmap)[phtrackref]; 
      RECOPFPHOT_PFchHad[iphot]          =  (*isoPFChargedphmap)[phtrackref]; 
      RECOPFPHOT_PFneuHad[iphot]         =  (*isoPFNeutralphmap)[phtrackref]; 
      RECOPFPHOT_PFphoton[iphot]         =  (*isoPFGammaphmap)[phtrackref]; 
      RECOPFPHOT_PFPUchAllPart[iphot]    =  (*isoPFPUphmap)[phtrackref]; 
*/
//TODO      RECOPFPHOT_PFchAllPart[iphot]      =  (cand->PflowIsolationVariables().sumChargedParticlePt); 
      RECOPFPHOT_PFchHad[iphot]          =  (cand->chargedHadronIso());
      RECOPFPHOT_PFneuHad[iphot]         =  (cand->neutralHadronIso()); 
      RECOPFPHOT_PFphoton[iphot]         =  (cand->photonIso()); 
      RECOPFPHOT_PFPUchAllPart[iphot]    =  (cand->puChargedHadronIso());


      
      
      RECOPFPHOT_PFX_rho[iphot]=( 
 				 RECOPFPHOT_PFchHad[iphot]+ 
 				 RECOPFPHOT_PFneuHad[iphot]+ 
 				 RECOPFPHOT_PFphoton[iphot]+ 
 				 RECOPFPHOT_PFPUchAllPart[iphot] 
				  )/double(cand->p4().pt()); 
      
      iphot++;
      
    }
    
  }


void fillTracks(const edm::Event& iEvent){
    // Tracks
  using namespace edm; using namespace std; using namespace reco;
  edm::Handle<edm::View<pat::PackedCandidate>> cands;
  iEvent.getByToken(pfTag_,cands);
 
  int countk=0;
  RECO_NTRACK=0;
  
  for(unsigned int i=0;i<cands->size();i++){
     if (countk>199) break;
     const pat::PackedCandidate & c = (*cands)[i];
     if(!(c.charge() != 0 && c.numberOfHits()> 0)) continue;
     RECO_NTRACK++;
    RECO_TRACK_PT[countk]=c.pseudoTrack().pt();
    RECO_TRACK_ETA[countk]=c.pseudoTrack().eta();
    RECO_TRACK_PHI[countk]=c.pseudoTrack().phi();
    RECO_TRACK_CHI2[countk]=c.pseudoTrack().chi2();
    RECO_TRACK_CHI2RED[countk]=c.pseudoTrack().normalizedChi2();
    //RECO_TRACK_CHI2PROB=TMath::Prob(i->chi2(),i->ndof());
    RECO_TRACK_CHI2PROB[countk]=ChiSquaredProbability(c.pseudoTrack().chi2(),c.pseudoTrack().ndof());
    RECO_TRACK_NHITS[countk]=c.pseudoTrack().numberOfValidHits();
    RECO_TRACK_DXY[countk]=c.pseudoTrack().dxy();
    RECO_TRACK_DXYERR[countk]=c.pseudoTrack().dxyError();
    RECO_TRACK_DZ[countk]=c.pseudoTrack().dz();
    RECO_TRACK_DZERR[countk]=c.pseudoTrack().dzError();
    countk++;
  }
 cout << "Number of Tracks in the event= " << RECO_NTRACK << endl;
}
        

  
 /* 
  double SetMET(const edm::Event& iEvent, edm::InputTag myTag_ ){
    double met=-999.;
    edm::Handle<reco::METCollection> metHandle;
    iEvent.getByLabel(myTag_,metHandle);
    for ( pat::METCollection::const_iterator iMet=metHandle->begin(); iMet!=metHandle->end(); iMet++) {
      met = iMet->pt();
    }
    return met;
  }
  
  double SetCaloMET(const edm::Event& iEvent, edm::InputTag myTag_ ){
    double met=-999.;
    edm::Handle<reco::CaloMETCollection> metHandle;
    iEvent.getByLabel(myTag_,metHandle);
    for ( CaloMETCollection::const_iterator iMet=metHandle->begin(); iMet!=metHandle->end(); iMet++) {
      met = iMet->pt();
    }
    return met;
  }
  
*/
  void fillMET(const edm::Event& iEvent){
    //TCMET
    //tcmet = SetMET(iEvent,trackermetTag_);
    //CALOMET
    //calomet = SetCaloMET(iEvent,calometTag_);
/*     calometnohf = SetCaloMET(iEvent,calometnohfTag_); */

/*     if (useAdditionalMET_==true){ */
/*       calometho = SetCaloMET(iEvent,calomethoTag_); */
/*       calometopt = SetCaloMET(iEvent,calometoptTag_); */
/*       calometoptnohf =  SetCaloMET(iEvent,calometoptnohfTag_); */
/*       calometoptnohfho =  SetCaloMET(iEvent,calometoptnohfhoTag_); */
/*       calometoptho = SetCaloMET(iEvent,calomethoTag_); */
/*       calometnohfho = SetCaloMET(iEvent,calometnohfhoTag_); */
/*     } */

    //Type I Muon Correction MET
    //cormetmuons = SetCaloMET(iEvent,cormetMuTag_);
    
    //PFMET
    edm::Handle<pat::METCollection> pfmetHandle;
    iEvent.getByToken(pfmetTag_,pfmetHandle);
    for ( pat::METCollection::const_iterator i=pfmetHandle->begin(); i!=pfmetHandle->end(); i++) {
      cormetmuons = i->pt();
      pfmet       = i->uncorPt();     
/*
      pfmet_x     = i->uncorPx();
      pfmet_y     = i->uncorPy();
      pfmet_phi   = i->uncorPhi();
*/
      pfmet_x     = i->px();
      pfmet_y     = i->py();
      pfmet_phi   = i->phi();
      pfmet_theta = i->uncorP3().theta();

      cout << "corrmet px=" << pfmet_x << "  met phi=" << pfmet_phi << endl;
    }

    std::cout << "MET:"
	      << "  tcMET="            << tcmet 
	      << "  caloMET="          << calomet 
/* 	      << "  caloMETopt="       << calometopt */
/* 	      << "  caloMEToptnohf="   << calometoptnohf */
/* 	      << "  caloMEToptnohfho=" << calometoptnohfho */
/* 	      << "  caloMEToptho="     << calometoptho */
/* 	      << "  caloMETnohf="      << calometnohf */
/* 	      << "  caloMETnohfho="    << calometnohfho */
/* 	      << "  caloMETho="        << calometho */
	      << "  TypeIMucorrMET="   << cormetmuons
	      << "  pfMET="            << pfmet
	      << std::endl;
    
  }
  
  		      
  void filljets(const edm::Event& iEvent,const edm::EventSetup& iSetup){
    edm::Handle<pat::JetCollection> pfjets,pfjetsmva;

    JetCorrectionUncertainty *jecUnc = NULL;
    if(filljec){   
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    jecUnc = new JetCorrectionUncertainty(JetCorPar);
    }

   iEvent.getByToken(jetsTag_, pfjets);
   iEvent.getByToken(jetsMVATag_, pfjetsmva);
/*    
    if (fillMCTruth == 1){
      cout << "Using Jet energy collection for MC " << endl;

      iEvent.getByToken(jetsTag_, pfjets);
      iEvent.getByToken(jetsMVATag_, pfjetsmva);
      iEvent.getByToken(PuJetMvaMCfullDiscr_,puJetIdMVAMC);
    }
    else{
      iEvent.getByToken(jetsDataTag_, pfjets);
      iEvent.getByToken(jetsMVATag_, pfjetsmva);
      iEvent.getByToken(PuJetMvaDatafullDiscr_,puJetIdMVAData);      
      
    }
*/

    RECO_PFJET_N = pfjets->size();
    cout << "Number of PFJets in the event= " << RECO_PFJET_N << endl;
    
    int index_jets = 0;
    
    for ( pat::JetCollection::const_iterator i=pfjets->begin(); i!=pfjets->end(); i++) {  
      if (index_jets>49) continue;
      
      edm::Ref<pat::JetCollection> pfjetref(pfjets,index_jets);      
      edm::Ref<pat::JetCollection> pfjetrefmva(pfjetsmva,index_jets);
      
      float mva = 0.;
      int pupass = 1;
  
//      mva = i->userFloat("vtxPx");    
      mva = i->userFloat("pileupJetId:fullDiscriminant");
/*
      if (fillMCTruth == 1){
	mva = (*puJetIdMVAMC)[pfjetrefmva];
      }
      else{
	mva = (*puJetIdMVAData)[pfjetrefmva];
      }
*/      
      //      New Selection
      if(i->pt()>20.){
	if(fabs(i->eta())>3.){
	  if(mva<=-0.45) pupass=0;
	}else if(fabs(i->eta())>2.75){
	  if(mva<=-0.55) pupass=0;
	}else if(fabs(i->eta())>2.5){
	  if(mva<=-0.6) pupass=0;
	}else if(mva<=-0.63) pupass=0;
      }else{
	if(fabs(i->eta())>3.){
	  if(mva<=-0.95) pupass=0;
	}else if(fabs(i->eta())>2.75){
	  if(mva<=-0.94) pupass=0;
	}else if(fabs(i->eta())>2.5){
	  if(mva<=-0.96) pupass=0;
	}else if(mva<=-0.95) pupass=0;
      }
        
      // Apply jet ID 
      double nhf = i->neutralHadronEnergyFraction();
      double nef = i->neutralEmEnergyFraction();
      double chf = i->chargedHadronEnergyFraction();
      double cef = i->chargedEmEnergyFraction();
      int nconstituents = i->numberOfDaughters();
      int nch = i->chargedMultiplicity();
      double muf = i->muonEnergyFraction();

      RECO_PFJET_NHF[index_jets] = nhf;
      RECO_PFJET_NEF[index_jets]     = nef;
      RECO_PFJET_CHF[index_jets]     = chf;
      RECO_PFJET_CEF[index_jets]    = cef;
      RECO_PFJET_nconstituents[index_jets] = nconstituents;
      RECO_PFJET_NCH[index_jets] = nch;     
      RECO_PFJET_MUF[index_jets] = muf;      
      cout <<"nconstituents test=" << nconstituents << " another="  << i->chargedMultiplicity()+i->neutralMultiplicity() << endl;

      RECO_PFJET_CHARGE[index_jets] = i->charge();
      RECO_PFJET_ET[index_jets]     = i->et();
      RECO_PFJET_PT[index_jets]     = i->pt();
      RECO_PFJET_ETA[index_jets]    = i->eta();
      RECO_PFJET_PHI[index_jets]    = i->phi();
      RECO_PFJET_PUID[index_jets]     = pupass;
      RECO_PFJET_PUID_MVA[index_jets] = mva;
  
      double unc=0;
      if(filljec){     
      jecUnc->setJetEta(i->eta());
      jecUnc->setJetPt(i->pt());
 
      unc = jecUnc->getUncertainty(true);
      }
      double pt_up = i->pt()*(1+unc);
      double pt_dow = i->pt()*(1-unc);

      RECO_PFJET_PT_UP[index_jets]=pt_up;
      RECO_PFJET_PT_DOW[index_jets]=pt_dow;

      cout 
	<< "PF Jet with ET= " << RECO_PFJET_ET[index_jets]   
	<< " PT="   << RECO_PFJET_PT[index_jets]   
	<< " ETA="  << RECO_PFJET_ETA[index_jets]  
	<< " PHI="  << RECO_PFJET_PHI[index_jets]  
	<< " PUID=" << RECO_PFJET_PUID[index_jets] 
	<< " PUID_MVA=" << RECO_PFJET_PUID_MVA[index_jets]
       //qier test
        << " Uncorrected Pt=" << i->correctedP4("Uncorrected").Pt()
        << " L1FastJet Pt=" << i->correctedP4("L1FastJet").Pt()
        << " L2Relative Pt=" << i->correctedP4("L2Relative").Pt()
        << " L3Absolute Pt=" << i->correctedP4("L3Absolute").Pt()
	<< endl;
      
      index_jets++;
    } // for loop on PFJets jets
    
    edm::Handle<double> rhoHandle;

    iEvent.getByToken(rhojetsTag_,rhoHandle); 
    if (rhoHandle.isValid() ) {
      RHO_mu=*rhoHandle;
      cout << "RHO mu fastjet= " << RHO_mu << endl; 
    }
    else {
      cout << "Not valid RHO mu collection" << endl;
    }
    
    iEvent.getByToken(rhojetsTag_,rhoHandle); 
    if (rhoHandle.isValid() ) {
      RHO_ele=*rhoHandle;
      cout << "RHO ele fastjet= " << RHO_ele << endl; 
    }
    else {
      cout << "Not valid RHO ele collection" << endl;
    } 
    
  }
  


  void fillBTagging(const edm::Event& iEvent){

    // trackCountingHighEffBJetTags
    edm::Handle<pat::JetCollection> bTagHandle;
    iEvent.getByToken(jetsTag_ ,bTagHandle);
    int l=0;
    for (pat::JetCollection::const_iterator btagIter=bTagHandle->begin(); btagIter!=bTagHandle->end();++btagIter) {
      if(l>=99) continue;
      double discrCSV1 = btagIter->bDiscriminator(tCHighEff_bTag_);
	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighEffBJetTags = "<< discrCSV1 
	    << " and jet Pt = "<<btagIter->pt()<<endl;      
	tCHighEff_BTagJet_PT[l]=btagIter->pt();
	tCHighEff_BTagJet_ETA[l]=btagIter->eta();
	tCHighEff_BTagJet_PHI[l]=btagIter->phi();
	tCHighEff_BTagJet_DISCR[l]=discrCSV1;
      

    // trackCountingHighPurBJetTags
     double discrCSV2 = btagIter->bDiscriminator(tCHighPur_bTag_);
	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighPurBJetTags = "<< discrCSV2
	    << " and jet Pt = "<<btagIter->pt()<<endl;      
	tCHighPur_BTagJet_PT[l]=btagIter->pt();
	tCHighPur_BTagJet_ETA[l]=btagIter->eta();
	tCHighPur_BTagJet_PHI[l]=btagIter->phi();
	tCHighPur_BTagJet_DISCR[l]=discrCSV2;

    // combinedSecondaryVertexBJetTags
    double discrCSV3 = btagIter->bDiscriminator(cSV_bTag_);
	cout<<" Jet "<< l
	    <<" has b tag discriminator combinedSecondaryVertexBJetTags = "<< discrCSV3
	    << " and jet Pt = "<<btagIter->pt()<<endl;      
	cSV_BTagJet_PT[l]=btagIter->pt();
	cSV_BTagJet_ETA[l]=btagIter->eta();
	cSV_BTagJet_PHI[l]=btagIter->phi();
	cSV_BTagJet_ET[l]=btagIter->et();
	cSV_BTagJet_DISCR[l]=discrCSV3;
      l++;

      if (fillMCTruth==true){
             if(abs(btagIter->partonFlavour())==5) RECOBOT_MatchingMCTruth[l]= 1;
             if(abs(btagIter->partonFlavour())==4) RECOBOT_MatchingMCTruth[l]= 2;
             if(abs(btagIter->partonFlavour())==0) RECOBOT_MatchingMCTruth[l]= 3;
//             if((abs(btagIter->partonFlavour())>=1&&abs(btagIter->partonFlavour())<=3)||abs(btagIter->partonFlavour())==21) RECOBOT_MatchingMCTruth[l] = 4;
             if(abs(btagIter->partonFlavour())!=5&&abs(btagIter->partonFlavour())!=4&&abs(btagIter->partonFlavour())!=0) RECOBOT_MatchingMCTruth[l] = 5;
     }
    }

  }

  
  
  void fillP3Covariance(const reco::PFCandidate &c, TMatrixDSym &bigCov, int offset) const {
    double dp = PFEnergyResolution().getEnergyResolutionEm(c.energy(), c.eta());
    // In order to produce a 3x3 matrix, we need a jacobian from (p) to (px,py,pz), i.e.
    //            [ Px/P  ]                
    //  C_(3x3) = [ Py/P  ] * sigma^2(P) * [ Px/P Py/P Pz/P  ]
    //            [ Pz/P  ]                
    AlgebraicMatrix31 ptop3;
    ptop3(0,0) = c.px()/c.p();
    ptop3(1,0) = c.py()/c.p();
    ptop3(2,0) = c.pz()/c.p();
    AlgebraicSymMatrix33 mat = ROOT::Math::Similarity(ptop3, AlgebraicSymMatrix11(dp*dp) );
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) {
	bigCov(offset+i,offset+j) = mat(i,j);
      } } 
  }
  
  
  
 private:
 
  std::vector<std::string> module_to_search;
  std::string par_to_search;
  
  // ROOT definition
  TFile *theFile_ ;
  TTree *theTree_ ;
  
  bool useRECOformat;
  std::string inputfileName;

  // Input tags
  std::string decaychannel;
  std::string flaginst;
  std::vector<std::string> flagtags;
  std::string rootFileName;

  // PU
  bool fillPUinfo;
  int num_PU_vertices;
  int PU_BunchCrossing; // bunch crossing for the PU event added
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSrc_;     

  // MC weight
  float MC_weighting,MC_weighting_un[9],PDF_weighting_un;
  edm::EDGetTokenT<GenEventInfoProduct> generator_;
  edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;

  // HLT
  bool fillHLTinfo;
  edm::InputTag HLTInfoFired;
  std::string HLTAnalysisinst;
  std::vector<edm::InputTag> flagHLTnames; 
  std::vector<std::string> HLTFilter_;

  bool fillLHEinfo;
  bool filljec;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent;     

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_,triggerPrescalesL1min_,triggerPrescalesL1max_;

  edm::InputTag triggerMatchObjectEle;
  edm::EDGetTokenT<edm::Association<std::vector<pat::TriggerObjectStandAlone> > > triggerMatchObject;
  std::string triggerFilter,triggerEleFilter,triggerHLTcollection;
  
  
  // SkimEarlyData
  std::string SkimEarlyDataAnalysisinst;
  std::vector<edm::InputTag> flagSkimEarlyDatanames; 
 
  // MC truth
  bool fillMCTruth;
  edm::EDGetTokenT<edm::View<reco::Candidate> > MCcollName;
  // GenParticles;
  edm::EDGetTokenT<vector<reco::GenParticle> > genParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > fourgenleptons_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > digenZ_;

  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
	
  // RECO
  bool useAdditionalRECO;
  bool use2011EA;
  std::vector<edm::InputTag> RECOcollNameBest2e2mu,RECOcollNameBest4e,RECOcollNameBest4mu,
    RECOcollNameBestRestFrame2e2mu,RECOcollNameBestRestFrame4mu,RECOcollNameBestRestFrame4e;
   edm::EDGetTokenT<edm::View<reco::Candidate> > zToMuMu, zToEE,
    zToMuMussmerge, zToEEssmerge, zToCrossLeptons,
    RECOcollNameDiLep_,RECOcollNameMMMM_,RECOcollNameEEEE,RECOcollNameEEMM,
    quadLeptons4Mu, quadLeptons2Mu2E, quadLeptons4E,
    triLeptonsMuMuMu, triLeptonsMuMuE, triLeptonsMuEE, triLeptonsEEE,
    quadLeptons3Mu1E, quadLeptons3E1Mu, quadLeptonsSSOSele,
    quadLeptonsSSOSmu, quadLeptonsSSOSelemu, quadLeptonsSSOSmuele,
    RECOcollNameLLLL;
    std::vector<edm::InputTag> RECOcollNameZ,RECOcollNameZss,
 //   RECOcollNameMMMM,RECOcollNameEEEE,RECOcollNameEEMM,
    RECOcollNameLLLLss,RECOcollNameLLL,RECOcollNameLLLl,RECOcollNameLLLLssos;
  edm::InputTag RECOcollNameDiLep;//RECOcollNameLLLL,RECOcollNameDiLep;
  
  // electron and muon tags
  bool useBestCandidate;
  edm::InputTag BestCandidatesLeptonsTag_;
  edm::InputTag clusterCollectionTag_,gsftrackCollection_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonPFTag_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronEgmTag_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonTag_;


  edm::EDGetTokenT<edm::ValueMap<float> > muonCorrPtErrorMapTag_;
    
  edm::InputTag electronEgmTkMapTag_;
  edm::InputTag electronEgmEcalMapTag_;
  edm::InputTag electronEgmHcalMapTag_;
  edm::EDGetTokenT<edm::View<pat::Electron> > mvaElectronTag_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaTrigV0MapTag_,mvaNonTrigV0MapTag_;
/*
  edm::EDGetTokenT<edm::ValueMap<double> >
    electronPFIsoValueChargedAllTag_,
    electronPFIsoValueChargedTag_,
    electronPFIsoValueNeutralTag_,
    electronPFIsoValueGammaTag_,
    electronPFIsoValuePUTag_;   


  edm::InputTag eleRegressionEnergyErrorTag_,eleRegressionEnergyTag_;

  edm::EDGetTokenT<edm::ValueMap<double> >
    muonPFIsoValueChargedAllTag_,
    muonPFIsoValueChargedTag_,
    muonPFIsoValueNeutralTag_,
    muonPFIsoValueGammaTag_,
    muonPFIsoValuePUTag_;   


  edm::EDGetTokenT<edm::ValueMap<double> >
    photonPFIsoValueChargedAllTag_,
    photonPFIsoValueChargedTag_,
    photonPFIsoValueNeutralTag_,
    photonPFIsoValueGammaTag_,
    photonPFIsoValuePUTag_;   

*/


  // vertexing 3D 
  edm::InputTag muonTag_Vert;
  edm::EDGetTokenT<edm::ValueMap<float> > muonMapTag_Vert,muonMapTag_VertValue,muonMapTag_VertError;

  // vertexing STIP / SLIP
  edm::EDGetTokenT<edm::ValueMap<float> >  electronSTIPMapTag_Vert,muonSTIPMapTag_Vert,
    electronSLIPMapTag_Vert,muonSLIPMapTag_Vert,
    electronSTIPMapTag_VertValue,muonSTIPMapTag_VertValue,
    electronSLIPMapTag_VertValue,muonSLIPMapTag_VertValue,
    electronSTIPMapTag_VertError,muonSTIPMapTag_VertError,
    electronSLIPMapTag_VertError,muonSLIPMapTag_VertError;
  
  // vertexing GD
  edm::InputTag ftsigma_Vert,ftsigmalag_Vert,
                ftsigma_VertMMMM,ftsigmalag_VertMMMM,
                ftsigma_VertEEEE,ftsigmalag_VertEEEE;
  edm::InputTag 
    gdX_Vert,gdY_Vert,gdZ_Vert,
    gdX_VertMMMM,gdY_VertMMMM,gdZ_VertMMMM,
    gdX_VertEEEE,gdY_VertEEEE,gdZ_VertEEEE;
   edm::InputTag 
    gdlagX_Vert,gdlagY_Vert,gdlagZ_Vert,gdlagProb_Vert,gdlagNdof_Vert,
     gdlagX_VertMMMM,gdlagY_VertMMMM,gdlagZ_VertMMMM,gdlagProb_VertMMMM,gdlagNdof_VertMMMM,
     gdlagX_VertEEEE,gdlagY_VertEEEE,gdlagZ_VertEEEE,gdlagProb_VertEEEE,gdlagNdof_VertEEEE;
  //

   // ConstraintVtx
   edm::InputTag 
     StandardFitVertex,StandardFitVertexMMMM,StandardFitVertexEEEE,
     KinematicFitVertex,KinematicFitVertexMMMM,KinematicFitVertexEEEE,
     RefittedMass,RefittedMassMMMM,RefittedMassEEEE;

   edm::InputTag 
     StandardFitVertexEEE,StandardFitVertexMMM,StandardFitVertexMEE,StandardFitVertexMME,StandardFitVertexDiLep;

  //electronID
  std::vector<edm::InputTag> eleIDTag_;
  edm::InputTag EleID_VeryLooseTag_ ;
  edm::InputTag EleID_LooseTag_  ;
  edm::InputTag EleID_MediumTag_  ;
  edm::InputTag EleID_TightTag_ ;

  edm::InputTag EleID_HZZVeryLooseTag_ ;
  edm::InputTag EleID_HZZLooseTag_  ;
  edm::InputTag EleID_HZZMediumTag_  ;
  edm::InputTag EleID_HZZTightTag_ ;
  
  edm::InputTag allelectronsColl, allmuonsColl;
  edm::InputTag isoVarTagElectronsCal,isoVarTag,isoVarTagElectronsTracker,isoVarTagElectronsX;
  edm::InputTag isoVarTagMuonsHCalIso,isoVarTagMuonsECalIso,
    isoVarTagMuonsCalIso,isoVarTagMuonsTracker,isoVarTagMuonsX;
  
  edm::InputTag theECALIsoDepositLabel;    //EM calorimeter Isolation deposit label
  edm::InputTag theHCALIsoDepositLabel;    //Hadron calorimeter Isolation deposit label
  edm::InputTag theHOCALIsoDepositLabel;   //Outer calorimeter Isolation deposit label
  edm::InputTag theTrackerIsoDepositLabel; //Tracker Isolation deposit label 
  
  // Photon, Tracks, Jets, Vertices
//  edm::EDGetTokenT<vector<reco::Track> > tracksTag_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfTag_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonsTag_;
  vector<reco::Vertex> PV;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesTag_;
  edm::EDGetTokenT<double> rhojetsTag_;
  edm::EDGetTokenT<vector<pat::Jet> > jetsTag_,jetsMVATag_;
//  edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaMCfullDiscr_,PuJetMvaMCfullId_;
//  edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaDatafullDiscr_,PuJetMvaDatafullId_;

  // genJet
  edm::EDGetTokenT<vector<reco::GenJet> >  genjetTag_;
  
  // MET
  edm::InputTag trackermetTag_;
  edm::EDGetTokenT<vector<pat::MET> >  pfmetTag_;
  edm::EDGetTokenT<vector<reco::GenMET> > genmetTag_;
  edm::InputTag calometTag_,calometoptTag_,calometoptnohfTag_,calometoptnohfhoTag_;
  edm::InputTag calometopthoTag_,calometnohfTag_,calometnohfhoTag_,calomethoTag_;
  bool useAdditionalMET_;
  edm::InputTag htmetic5Tag_,htmetkt4Tag_,htmetkt6Tag_,htmetsc5Tag_,htmetsc7Tag_;
  edm::InputTag jescormetic5Tag_,jescormetkt4Tag_,jescormetkt6Tag_,jescormetsc5Tag_,jescormetsc7Tag_;
  edm::InputTag cormetMuTag_;

  // Conversion
//  edm::EDGetTokenT<edm::ValueMap<float> > ConvMapDistTag_,ConvMapDcotTag_;

  // Matching
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodElectronMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myElectrons_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodMuonMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myMuons_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodGammaMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myGammas_;

  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodZtoMuMuMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodZtoEEMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToEEMMMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToMMMMMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToEEEEMCMatch_;
  // Beam Spot
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpot_;

  // bTagging
//  edm::EDGetTokenT<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >  
  std::string tCHighEff_bTag_, tCHighPur_bTag_,cSV_bTag_;

//  const std::vector<std::string> bDiscriminators_;

  // counters
  int irun, ievt,ils,nevt;
  float Avginstlumi;

  // HLT
  int RECO_nMuHLTMatch,RECO_nEleHLTMatch;
  float RECOMU_PT_MuHLTMatch[50],RECOMU_ETA_MuHLTMatch[50],RECOMU_PHI_MuHLTMatch[50];
  float RECOELE_PT_EleHLTMatch[50],RECOELE_ETA_EleHLTMatch[50],RECOELE_PHI_EleHLTMatch[50];
  bool RECOMU_sm_MuHLTMatch[50],RECOELE_se_EleHLTMatch[50];
  int RECOMU_dm_MuHLTMatch[50],RECOELE_de_EleHLTMatch[50],RECOJET_JetHLTMatch[50];
  char HLTPathsFired[10000];
  bool dm_trig, sm_trig, de_trig, se_trig;
  bool tri_trig;
  int jet_trig; 

  // MC info
  edm::ESHandle<ParticleDataTable>  pdt_;
  
  // MC truth
  float MC_E[7],MC_PT[7],MC_ETA[7],MC_THETA[7],MC_PHI[7],MC_MASS[7],MC_PDGID[7];

  float MC_LEPT_PT[4],MC_LEPT_ETA[4],MC_LEPT_PHI[4],MC_LEPT_THETA[4],MC_LEPT_PDGID[4];
  float MC_Z_MASS[2][5],MC_Z_PT[2][5],MC_Z_ETA[2][5],MC_Z_PHI[2][5],MC_Z_THETA[2][5],MC_Z_PDGID[2][5];

  // RECO collection
 
  

  int leptonflavor;

  // RECO electrons
  edm::ESHandle<CaloGeometry> theCaloGeom_;  
  float RECOELE_E[50],RECOELE_PT[50],RECOELE_PTError[50],RECOELE_P[50],RECOELE_ETA[50],RECOELE_THETA[50],RECOELE_PHI[50],RECOELE_MASS[50];
  float RECOELE_CHARGE[50], RECOELE_ID[50], RECOELE_PT_uncorr[50];

  bool RECOELE_isEcalDriven[50], RECOELE_isTrackerDriven[50];
  float 
    RECOELE_gsftrack_NPixHits[50],
    RECOELE_gsftrack_NStripHits[50],
    RECOELE_gsftrack_chi2[50], 
    RECOELE_gsftrack_dxyB[50], RECOELE_gsftrack_dxy[50], RECOELE_gsftrack_dxyError[50],
    RECOELE_gsftrack_dzB[50], RECOELE_gsftrack_dz[50],RECOELE_gsftrack_dzError[50];
  float
    RECOELE_scl_E[50],RECOELE_scl_Et[50],RECOELE_scl_Eta[50],RECOELE_scl_Phi[50]; 
  //float RECOELE_eSuperClusterOverP[50], RECOELE_eSeedClusterOverPout[50], RECOELE_deltaEtaSuperClusterTrackAtVtx[50], RECOELE_deltaPhiSuperClusterTrackAtVtx[50];
  float 
    RECOELE_ep[50], RECOELE_eSeedp[50], RECOELE_eSeedpout[50], RECOELE_eElepout[50],
    RECOELE_deltaEtaIn[50],RECOELE_deltaEtaSeed[50],RECOELE_deltaEtaEle[50],RECOELE_deltaPhiIn[50],
    RECOELE_deltaPhiSeed[50],RECOELE_deltaPhiEle[50],RECOELE_ecalEnergy[50];
  int RECOELE_isbarrel[50], RECOELE_isendcap[50], RECOELE_isEBetaGap[50], RECOELE_isEBphiGap[50], RECOELE_isEEdeeGap[50], RECOELE_isEEringGap[50], RECOELE_isGap[50]; 
  float RECOELE_sigmaIetaIeta[50], RECOELE_sigmaEtaEta[50], RECOELE_e15[50], RECOELE_e25max[50], RECOELE_e55[50], RECOELE_he[50], RECOELE_r9[50];
  float RECOELE_mva[50], RECOELE_fbrem[50],RECOELE_fbrem_mean[50],RECOELE_fbrem_mode[50];
  int RECOELE_nbrems[50], RECOELE_Class[50];
  
  float RECOELE_EGMTRACKISO[50],RECOELE_EGMECALISO[50],RECOELE_EGMHCALISO[50],RECOELE_EGMX[50],
    RECOELE_IP[50],RECOELE_SIP[50],RECOELE_IPERROR[50],
    RECOELE_STIP[50],RECOELE_TIP[50],RECOELE_TIPERROR[50],
    RECOELE_SLIP[50],RECOELE_LIP[50],RECOELE_LIPERROR[50];

  double RECOELE_PFchAllPart[50],RECOELE_PFchHad[50],RECOELE_PFneuHad[50],RECOELE_PFphoton[50],
    RECOELE_PFPUchAllPart[50],RECOELE_PFX_dB[50],RECOELE_PFX_rho[50],RECOELE_PF_RingsIsoMVA[50];

  double RECOELE_regEnergy[50],RECOELE_regEnergyError[50];

  float RECOELE_TLE_ParentSC_X[50],RECOELE_TLE_ParentSC_Y[50],RECOELE_TLE_ParentSC_Z[50];
  
  int RECOELE_EEEE_MATCHED[50],RECOELE_EEMM_MATCHED[50],RECOELE_ZEE_MATCHED[50],RECOELE_ZssEE_MATCHED[10],RECOELE_ZEM_MATCHED[50],
    RECOELE_LLL0_MATCHED[50],RECOELE_LLL1_MATCHED[50],RECOELE_LLL2_MATCHED[50],RECOELE_LLL3_MATCHED[50],
    RECOELE_LLLLss0_MATCHED[50],RECOELE_LLLLss1_MATCHED[50],RECOELE_LLLLss2_MATCHED[50],
    RECOELE_LLLl0_MATCHED[50],RECOELE_LLLl1_MATCHED[50],RECOELE_LLLL_MATCHED[50];
  double RECOELE_mvaTrigV0[50],RECOELE_mvaNonTrigV0[50];
  
  double ele_sclRawE[50] ;
  double ele_sclX[50], ele_sclY[50], ele_sclZ[50];
  int ele_seedSubdet1[50];
  double ele_seedDphi1[50], ele_seedDrz1[50];
  int ele_seedSubdet2[50];
  double ele_seedDphi2[50], ele_seedDrz2[50];
  double ele_eidVeryLoose[50], ele_eidLoose[50], ele_eidMedium[50], ele_eidTight[50] ;
  double ele_eidHZZVeryLoose[50], ele_eidHZZLoose[50], ele_eidHZZMedium[50], ele_eidHZZTight[50] ;
  double RECOELE_COV[50][3][3];

   float RECOELE_ecalTrkEnergyErrPostCorr[50],RECOELE_energyScaleValue[50],RECOELE_energySigmaValue[50], RECOELE_energyScaleUp[50], RECOELE_energyScaleDown[50], RECOELE_energyScaleStatUp[50], RECOELE_energyScaleStatDown[50], RECOELE_energyScaleSystUp[50], RECOELE_energyScaleSystDown[50], RECOELE_energyScaleGainUp[50], RECOELE_energyScaleGainDown[50],RECOELE_energyScaleEtUp[50], RECOELE_energyScaleEtDown[50], RECOELE_energySigmaUp[50], RECOELE_energySigmaDown[50], RECOELE_energySigmaPhiUp[50], RECOELE_energySigmaPhiDown[50], RECOELE_energySigmaRhoUp[50], RECOELE_energySigmaRhoDown[50],RECOELE_ecalTrkEnergyPreCorr[50], RECOELE_ecalTrkEnergyErrPreCorr[50];

  // RECO muons
  bool RECOMU_isPFMu[50],RECOMU_isGlobalMu[50],RECOMU_isStandAloneMu[50],RECOMU_isTrackerMu[50],RECOMU_isCaloMu[50],RECOMU_isTrackerHighPtMu[50],RECOMU_isME0Muon[50];
  float RECOMU_E[50],RECOMU_PT[50],RECOMU_P[50],RECOMU_ETA[50],RECOMU_THETA[50],RECOMU_PHI[50],RECOMU_MASS[50],RECOMU_CHARGE[50];
  double RECOMU_COV[50][3][3];

  float 
    RECOMU_TRACKISO[50],RECOMU_TRACKISO_SUMPT[50],RECOMU_ECALISO[50],RECOMU_HCALISO[50], RECOMU_X[50],   
    RECOMU_IP[50],RECOMU_SIP[50],RECOMU_IPERROR[50],
    RECOMU_IP_KF[50],RECOMU_SIP_KF[50],RECOMU_IPERROR_KF[50],
    RECOMU_STIP[50],RECOMU_TIP[50],RECOMU_TIPERROR[50],
    RECOMU_SLIP[50],RECOMU_LIP[50],RECOMU_LIPERROR[50],
    RECOMU_SIP_GD[50], RECOMU_SIP_GDMMMM[50],
    RECOMU_SIP_Std[50], RECOMU_SIP_StdMMMM[50], 
    RECOMU_SIP_Kin[50], RECOMU_SIP_KinMMMM[50];

  double 
    RECOMU_PFchHad[50],RECOMU_PFneuHad[50],RECOMU_PFphoton[50],RECOMU_PFsumPUPt[50],RECOMU_PFX_dB[50],RECOMU_PFPUchAllPart[50],RECOMU_PFchAllPart[50],RECOMU_PFX_rho[50],
    RECOMU_PFchHad42[50],RECOMU_PFneuHad42[50],RECOMU_PFphoton42[50],RECOMU_PFPUchAllPart42[50],RECOMU_PFchAllPart42[50],
    RECOMU_PF_RingsIsoMVA[50],RECOMU_PF_RingsIDMVA[50];

  float
    RECOMU_caloCompatibility[50],RECOMU_segmentCompatibility[50];
  bool RECOMU_glbmuPromptTight[50],RECOMU_isMedium[50];
  float RECOMU_chi2LocalPosition[50], RECOMU_trkKink[50];
 
  int RECOMU_MMMM_MATCHED[50],RECOMU_EEMM_MATCHED[50],
      RECOMU_ZMM_MATCHED[50],RECOMU_ZssMM_MATCHED[50],RECOMU_ZEM_MATCHED[50],
      RECOMU_LLL0_MATCHED[50],RECOMU_LLL1_MATCHED[50],RECOMU_LLL2_MATCHED[50],RECOMU_LLL3_MATCHED[50],
    RECOMU_LLLLss0_MATCHED[50],RECOMU_LLLLss1_MATCHED[50],RECOMU_LLLLss2_MATCHED[50],
    RECOMU_LLLl0_MATCHED[50],RECOMU_LLLl1_MATCHED[50],RECOMU_LLLL_MATCHED[50];
  int RECOMU_numberOfMatches[50],RECOMU_numberOfMatchedStations[50];
  int RECOMU_mubesttrkType[50];
  
  float 
    RECOMU_mutrkPT[50],RECOMU_mutrkPTError[50],
    RECOMU_mutrkDxy[50],RECOMU_mutrkDxyError[50],RECOMU_mutrkDxyB[50],
    RECOMU_mutrkDz[50],RECOMU_mutrkDzError[50],RECOMU_mutrkDzB[50],
    RECOMU_mutrkChi2PerNdof[50],
    RECOMU_mutrktrackerLayersWithMeasurement[50],
    RECOMU_muInnertrkDxy[50],RECOMU_muInnertrkDxyError[50],RECOMU_muInnertrkDxyB[50],
    RECOMU_muInnertrkDz[50],RECOMU_muInnertrkDzError[50],RECOMU_muInnertrkDzB[50],
    RECOMU_mubesttrkDxy[50],RECOMU_mubesttrkDxyError[50],RECOMU_mubesttrkDxyB[50],
    RECOMU_mubesttrkDz[50],RECOMU_mubesttrkDzError[50],RECOMU_mubesttrkDzB[50],
    RECOMU_mubesttrkPTError[50],
    RECOMU_muInnertrkChi2PerNdof[50],
    RECOMU_muInnertrktrackerLayersWithMeasurement[50],RECOMU_muInnertrkPT[50],RECOMU_muInnertrkPTError[50],
    RECOMU_muInnertrkCharge[50],RECOMU_muInnertrkNHits[50],RECOMU_muInnertrkNPixHits[50],RECOMU_muInnertrkNStripHits[50],RECOMU_muInnertrkvalidFraction[50],
    RECOMU_mutrkCharge[50],RECOMU_mutrkNHits[50],RECOMU_mutrkNPixHits[50],RECOMU_mutrkNStripHits[50],RECOMU_mutrkNMuonHits[50];
  bool RECOMU_trkmuArbitration[50],RECOMU_trkmu2DCompatibilityLoose[50],RECOMU_trkmu2DCompatibilityTight[50];
  bool RECOMU_trkmuOneStationLoose[50],RECOMU_trkmuOneStationTight[50];
  bool RECOMU_trkmuLastStationLoose[50],RECOMU_trkmuLastStationTight[50];
  bool RECOMU_trkmuLastStationAngLoose[50],RECOMU_trkmuLastStationAngTight[50];
  bool RECOMU_trkmuOneStationAngLoose[50],RECOMU_trkmuOneStationAngTight[50];
  bool RECOMU_trkmuLastStationOptimizedLowPtLoose[50],RECOMU_trkmuLastStationOptimizedLowPtTight[50];
  
   // Photons
  float RECOPHOT_PT[20],RECOPHOT_ETA[20],RECOPHOT_PHI[20],RECOPHOT_THETA[20],RECOPHOT_TLE_ParentSC_X[20],RECOPHOT_TLE_ParentSC_Y[20],RECOPHOT_TLE_ParentSC_Z[20];
  float RECOPFPHOT_PT[20],RECOPFPHOT_PTError[20],RECOPFPHOT_ETA[20],RECOPFPHOT_PHI[20],RECOPFPHOT_THETA[20];
  double RECOPFPHOT_PFchAllPart[20],RECOPFPHOT_PFchHad[20],RECOPFPHOT_PFneuHad[20],RECOPFPHOT_PFphoton[20],
    RECOPFPHOT_PFPUchAllPart[20],RECOPFPHOT_PFX_rho[20];
  
  
  //Muons Matching
  bool RECOMU_MatchingMCTruth[50];
  float RECOMU_MatchingMCpT[50];
  float RECOMU_MatchingMCEta[50];
  float RECOMU_MatchingMCPhi[50];

  //Bottom Matching
  int RECOBOT_MatchingMCTruth[50];
  float RECOBOT_MatchingMCpT[50];
  float RECOBOT_MatchingMCEta[50];
  float RECOBOT_MatchingMCPhi[50];
  
  //Electrons:
  bool RECOELE_MatchingMCTruth[50];
  float RECOELE_MatchingMCpT[50];
  float RECOELE_MatchingMCEta[50];
  float RECOELE_MatchingMCPhi[50];
  //Gamma:
  bool RECOPHOT_MatchingMCTruth[50];
  float RECOPHOT_MatchingMCpT[50];
  float RECOPHOT_MatchingMCEta[50];
  float RECOPHOT_MatchingMCPhi[50];

 

  // RECO counters
  int RECO_NMU, RECO_NELE, RECO_NTRACK, RECO_NPHOT, RECO_NPFPHOT,RECO_NJET, RECO_NVTX;
  float RECO_TRACK_PT[100], RECO_TRACK_ETA[100], RECO_TRACK_PHI[100],
    RECO_TRACK_CHI2[100],RECO_TRACK_CHI2RED[100],RECO_TRACK_CHI2PROB[100], 
    RECO_TRACK_DXY[100],RECO_TRACK_DXYERR[100], 
    RECO_TRACK_DZ[100],RECO_TRACK_DZERR[100];
  int RECO_TRACK_NHITS[100];
  
  // Primary Vertices
  float RECO_VERTEX_x[15], RECO_VERTEX_y[15], RECO_VERTEX_z[15],RECO_VERTEX_ndof[15],RECO_VERTEX_chi2[15],RECO_VERTEXPROB[15],RECO_VERTEX_TRACK_PT[15][50];
  bool RECO_VERTEX_isValid[15];
  int RECO_VERTEX_ntracks[15];
  
  // RECO JETS
  int RECO_PFJET_N, RECO_PFJET_CHARGE[50],RECO_PFJET_PUID[50];
  int RECO_PFJET_nconstituents[50],RECO_PFJET_NCH[50];
  float RECO_PFJET_ET[50], RECO_PFJET_PT[50], RECO_PFJET_ETA[50], RECO_PFJET_PHI[50],RECO_PFJET_PUID_MVA[50];
  float RECO_PFJET_PT_UP[50], RECO_PFJET_PT_DOW[50];
  float RECO_PFJET_NHF[50],RECO_PFJET_NEF[50],RECO_PFJET_CHF[50],RECO_PFJET_CEF[50],RECO_PFJET_MUF[50];
  double RHO,RHO_ele,RHO_mu;

  // GenJET
  float MC_GENJET_PT[50], MC_GENJET_ETA[50], MC_GENJET_PHI[50];

  // RECO MET
  float genmet;
  float calomet;
    //calometopt,calometoptnohf,calometoptnohfho,calometoptho,calometnohf,calometnohfho,calometho;       
  float pfmet,pfmet_x,pfmet_y,pfmet_phi,pfmet_theta;

    //htmetic5,htmetkt4,htmetkt6,htmetsc5,htmetsc7;        
    float tcmet;
    //jescormetic5,jescormetkt4,jescormetkt6,jescormetsc5,jescormetsc7;    
  float cormetmuons;  
  
  // Beam Spot
  double BeamSpot_X,BeamSpot_Y,BeamSpot_Z;	
  BeamSpot bs;
  
 
  
  float tCHighEff_BTagJet_PT[50],
    tCHighPur_BTagJet_PT[50],
    cSV_BTagJet_PT[50];
  float tCHighEff_BTagJet_ETA[50],
    tCHighPur_BTagJet_ETA[50],
    cSV_BTagJet_ETA[50];
  float tCHighEff_BTagJet_PHI[50],
    tCHighPur_BTagJet_PHI[50],
    cSV_BTagJet_PHI[50];
  float tCHighEff_BTagJet_DISCR[50],
    tCHighPur_BTagJet_DISCR[50],
    cSV_BTagJet_DISCR[50];
  float cSV_BTagJet_ET[50];

  float ConvMapDist[50],ConvMapDcot[50];

  // MVA Ring Isolation
  //MuonMVAEstimator *fMuonIsoMVA;

  // Magnetic Field
  edm::ESHandle<MagneticField> magfield_;

};

#endif



