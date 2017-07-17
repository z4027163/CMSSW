/** \class HZZ4LeptonsHLTAnalysis
 *
 * See header file for documentation
 *
 *
 *  \author Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsHLTAnalysis.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"

// Muons:
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <cassert>

using namespace std;
using namespace edm;
using namespace reco;

//
// constructors and destructor
//
HZZ4LeptonsHLTAnalysis::HZZ4LeptonsHLTAnalysis(const edm::ParameterSet& iConfig)
{
  // get names from module parameters, then derive slot numbers

  inputTag_           = iConfig.getParameter<edm::InputTag> ("TriggerResultsTag");
  muonlabel_          = iConfig.getParameter<edm::InputTag> ("MuonCollectionLabel");
  electronlabel_      = iConfig.getParameter<edm::InputTag> ("ElectronCollectionLabel");
  andOr_              = iConfig.getParameter<bool> ("andOr" );
  n_                  = 0;
  firstevent_         = true;  

	// PG and FRC 06-07-11
	debug	=	iConfig.getUntrackedParameter<bool> ("debug", false);

  HLTPathsByName_= iConfig.getParameter<std::vector<std::string > >("HLTPaths");
  n_=HLTPathsByName_.size();
  HLTPathsByIndex_.resize(n_);

  ntrig.resize(n_);
  boolflag.resize(n_);
  npassed=0;
  
  int loc=0;
  for (unsigned int j = 0; j < n_; ++ j ){
    string collection;
    collection = "flag" + HLTPathsByName_[j];
    loc = collection.find( "_", 0 );
    if (loc > 0) {
      collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
      loc = collection.find( "_", 0 );
      if (loc > 0) collection.erase(loc, 1);
    }
    valias.push_back(collection);
    produces<bool >(valias.at(j)).setBranchAlias( valias.at(j) );
  }

  aliasaccept="flagHLTaccept";
  produces<bool> (aliasaccept).setBranchAlias(aliasaccept);


  // this is a user/analysis filter: it places no product into the event!

}

HZZ4LeptonsHLTAnalysis::~HZZ4LeptonsHLTAnalysis()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HZZ4LeptonsHLTAnalysis::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  
   const string invalid("@@invalid@@");

   // get hold of TriggerResults Object
   Handle<TriggerResults> trh;
   try {iEvent.getByLabel(inputTag_,trh);} catch(...) {;}
   if (trh.isValid()) {
     if(debug) cout << "TriggerResults found, number of HLT paths: " << trh->size() << endl;
     
     // get hold of trigger names - based on TriggerResults object!
     //triggerNames_.init(*trh);
     triggerNames_=iEvent.triggerNames(*trh);
     if (firstevent_) {
       for (unsigned int i=0; i<triggerNames_.size(); i++) {
	 if(debug) cout << "Found the trigger path= " << triggerNames_.triggerName(i) << endl;
       }
     }

     unsigned int n(n_);
     for (unsigned int i=0; i!=n; i++) {
       HLTPathsByIndex_[i]=triggerNames_.triggerIndex(HLTPathsByName_[i]);
     }
     
     // for empty input vectors (n==0), default to all HLT trigger paths!
     if (n==0) {
       n=trh->size();
       HLTPathsByName_.resize(n);
       HLTPathsByIndex_.resize(n);
       for (unsigned int i=0; i!=n; i++) {
	 HLTPathsByName_[i]=triggerNames_.triggerName(i);
	 HLTPathsByIndex_[i]=i;
       }
     }
     
     // report on what is finally used
     if (firstevent_){ 
       if(debug) cout << "HLT trigger paths: " + inputTag_.encode()
	    << "\n - Number requested: " << n
	    << "\n - andOr mode: " << andOr_ << endl; 
       if (n>0) {
	 if(debug) cout << "  HLT trigger paths requested: index, name and valididty:" << endl;
	 for (unsigned int i=0; i!=n; i++) {
	   bool validity ( (HLTPathsByIndex_[i]<trh->size()) && (HLTPathsByName_[i]!=invalid) );
	   
	   if(debug) cout << " " << HLTPathsByIndex_[i]
		<< " " << HLTPathsByName_[i]
		<< " " << validity << endl;
	 }
       }
     }
     
     // count number of requested HLT paths which have fired
     unsigned int fired(0);
     for (unsigned int i=0; i!=n; i++) {
       boolflag[i]=false;
       if (HLTPathsByIndex_[i]<trh->size()) {
	 if (trh->accept(HLTPathsByIndex_[i])) {
	   fired++;
	   if(debug) cout << "Fired HLT path= " << HLTPathsByName_[i] << endl;
	   ntrig[i]++;
	   boolflag[i]=true;
	 }
       }
     }
     

     // Boolean filter result
     const bool accept( ((!andOr_) && (fired==n)) ||
			(( andOr_) && (fired!=0)) );
     if(debug) cout << "Accept = " << accept << endl;

     // Try to study events skipped at trigger level
     if (!accept){
       if(debug) cout << "Perform sanity check analysis for trigger" << endl;

       // Get the muons 
       edm::Handle<edm::View<reco::Muon> >mus;
       iEvent.getByLabel(muonlabel_.label(), mus);
       for ( edm::View<reco::Muon>::const_iterator muons = mus->begin(); muons != mus->end(); ++muons ) {
	 if(debug) cout << "Not triggered event in Run=" << iEvent.id().run() << " and Event=" << iEvent.id().event()
	      << " Muon from not-triggered event with pt= " << muons->pt() 
	      << " and eta= " << muons->eta()
	      << endl;
       }

       // Get the electrons
       edm::Handle<edm::View<reco::GsfElectron> >pTracks;
       iEvent.getByLabel(electronlabel_.label(),pTracks);
       const edm::View<reco::GsfElectron>* eTracks = pTracks.product();       
       for ( edm::View<reco::GsfElectron>::const_iterator electrons = eTracks->begin(); electrons != eTracks->end(); ++electrons ) {
	 if(debug) cout << "Not triggered event in Run=" << iEvent.id().run() << " and Event=" << iEvent.id().event()
	      << " Electron from not-triggered event with pt= " << electrons->pt() 
	      << " and eta= " << electrons->eta()
	      << endl; 
       }

     }
     else {
       npassed++;
     }
     firstevent_=false;

     // store the flags
     for  (unsigned int j = 0; j < n_; ++ j ){
       auto_ptr<bool> flag ( new bool );
       *flag=boolflag[j];
       iEvent.put(flag,valias.at(j));
     }
     auto_ptr<bool> flagaccept ( new bool );
     *flagaccept=accept;
     iEvent.put(flagaccept,aliasaccept);
   }
   
}

void HZZ4LeptonsHLTAnalysis::endJob() {
  
  for (unsigned int i=0; i<ntrig.size(); i++) {
    if(debug) cout << "Triggered paths " << HLTPathsByName_[i] << "= " << ntrig[i] << endl;
  }

  if(debug) cout << "Total passing HLT= " << npassed << endl;

}



