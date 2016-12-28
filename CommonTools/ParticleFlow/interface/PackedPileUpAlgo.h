#ifndef CommonTools_PFCandProducer_PFPileUpAlgo_
#define CommonTools_PFCandProducer_PFPileUpAlgo_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class PackedPileUpAlgo {
 public:


  typedef std::vector< edm::FwdPtr<pat::PackedCandidate> >  PFCollection;

  PackedPileUpAlgo():checkClosestZVertex_(true), verbose_(false) {;}
    
  PackedPileUpAlgo( bool checkClosestZVertex, bool verbose=false):
    checkClosestZVertex_(checkClosestZVertex), verbose_(verbose) {;}

  ~PackedPileUpAlgo(){;}

  // the last parameter is needed if you want to use the sourceCandidatePtr
  void process(const PFCollection & pfCandidates, 
	       const reco::VertexCollection & vertices)  ;

  inline void setVerbose(bool verbose) { verbose_ = verbose; }

  inline void setCheckClosestZVertex(bool val) { checkClosestZVertex_ = val;}

  const PFCollection & getPFCandidatesFromPU() const {return pfCandidatesFromPU_;}
  
  const PFCollection & getPFCandidatesFromVtx() const {return pfCandidatesFromVtx_;}

  int chargedHadronVertex(const reco::VertexCollection& vertices, 
			const pat::PackedCandidate& pfcand ) const;


 private  :

  /// use the closest z vertex if a track is not in a vertex
  bool   checkClosestZVertex_;
  
  
  /// verbose ?
  bool   verbose_;

  PFCollection pfCandidatesFromVtx_;
  PFCollection pfCandidatesFromPU_;
  
};

#endif
