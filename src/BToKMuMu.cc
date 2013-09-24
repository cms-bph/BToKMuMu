// -*- C++ -*-
//
// Package:    BToKMuMu
// Class:      BToKMuMu
// 
/**\class BToKMuMu BToKMuMu.cc BphAna/BToKMuMu/src/BToKMuMu.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Xin Shi <Xin.Shi@cern.ch> 
//         Created:  Tue Sep 24 09:16:20 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


#include <TFile.h>
#include <TTree.h>


using namespace std;

//
// class declaration
//

class BToKMuMu : public edm::EDAnalyzer {
public:
  explicit BToKMuMu(const edm::ParameterSet&);
  ~BToKMuMu();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
				    edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&,
				  edm::EventSetup const&);

  bool buildBuToKMuMu(const edm::Event &); 
  void clearVariables(); 
  bool hasBeamSpot(const edm::Event&);
  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &); 
  bool hasPrimaryVertex(const edm::Event &); 
  void hltReport(const edm::Event&);


  // ----------member data ---------------------------
  // --- begin input from python file --- 
  string OutputFileName_; 
  bool BuildBuToKMuMu_; 

  // labels 
  edm::InputTag TriggerResultsLabel_;
  edm::InputTag BeamSpotLabel_;
  edm::InputTag VertexLabel_;
  edm::InputTag MuonLabel_;

  vector<string> TriggerNames_; 
  vector<string> LastFilterNames_;

  // pre-selection cuts
  double MuonMinPt_; 
  double MuonMaxEta_; 
  double MuonMaxDcaBs_; 

  // --- end input from python file --- 

 
  // Across the event 
  map<string, string> mapTriggerToLastFilter_;
  reco::BeamSpot beamSpot_;  
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;

  // ---- Root Variables ---- 
  TFile* fout_;
  TTree* tree_;
  
  unsigned int run, event, lumiblock, nprivtx; 
  vector<string> *triggernames;
  vector<int> *triggerprescales;

  // B+ and B- 
  int nb; 


  // variables to monitor 
  TDatime t_begin_, t_now_ ;
  int n_processed_, n_selected_; 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BToKMuMu::BToKMuMu(const edm::ParameterSet& iConfig):
  OutputFileName_(iConfig.getParameter<string>("OutputFileName")),
  BuildBuToKMuMu_(iConfig.getUntrackedParameter<bool>("BuildBuToKMuMu")),

  // labels 
  TriggerResultsLabel_(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")),
  MuonLabel_(iConfig.getParameter<edm::InputTag>("MuonLabel")),

  // pre-selection cuts 
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),

  tree_(0), 
  triggernames(0), triggerprescales(0), 
  nb(0) 
{
  //now do what ever initialization is needed
  
  assert(TriggerNames_.size() == LastFilterNames_.size());
  for (size_t i = 0; i < TriggerNames_.size(); ++i)
    mapTriggerToLastFilter_[TriggerNames_[i]] = LastFilterNames_[i];
}


BToKMuMu::~BToKMuMu()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BToKMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  n_processed_ += 1; 
  clearVariables(); 
  run = iEvent.id().run() ;
  event = iEvent.id().event() ;
  lumiblock = iEvent.luminosityBlock(); 
  hltReport(iEvent);

  if ( hasBeamSpot(iEvent) ) {
    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);      
    
    if ( bFieldHandle_.isValid() && hasPrimaryVertex(iEvent) ) { 
      
      if ( BuildBuToKMuMu_ && buildBuToKMuMu(iEvent) ) { 
	tree_->Fill();
	n_selected_ += 1;
      }
    }
  }

  clearVariables(); 
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
BToKMuMu::beginJob()
{
  t_begin_.Set(); 
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;

  fout_ = new TFile(OutputFileName_.c_str(), "RECREATE");
  fout_->cd();

  tree_ = new TTree ("tree", "BToKMuMu");

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/i");
  tree_->Branch("lumiblock", &lumiblock, "lumiblock/i");
  tree_->Branch("nprivtx", &nprivtx, "nprivtx/i");
  tree_->Branch("triggernames", &triggernames);
  tree_->Branch("triggerprescales", &triggerprescales);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
BToKMuMu::endJob() 
{
  fout_->cd();
  tree_->Write();
  fout_->Close();

  t_now_.Set(); 
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();  
  printf(" processed: %i \n selected: %i \n \
 duration: %i sec \n rate: %g evts/sec\n",
	 n_processed_, n_selected_, 
	 t_now_.Convert() - t_begin_.Convert(), 
	 float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );
  
}

// ------------ method called when starting to processes a run  ------------
void 
BToKMuMu::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BToKMuMu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BToKMuMu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BToKMuMu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BToKMuMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void 
BToKMuMu::clearVariables(){
  run = 0; 
  event = 0;
  lumiblock = 0;
  nprivtx = 0; 
  triggernames->clear();
  triggerprescales->clear();
  nb = 0; 
}

void
BToKMuMu::hltReport(const edm::Event& iEvent)
{
  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try {iEvent.getByLabel( TriggerResultsLabel_, hltTriggerResults ); }
  catch ( ... ) { edm::LogInfo("myHLT") 
      << __LINE__ << " : couldn't get handle on HLT Trigger" ; }
  
  HLTConfigProvider hltConfig_;
  if (hltTriggerResults.isValid()) {
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);

    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++){

      // Only consider the triggered case. 
      if ((*hltTriggerResults)[itrig].accept() == 1){ 

	string triggername = triggerNames_.triggerName(itrig);	
	int triggerprescale = hltConfig_.prescaleValue(itrig, triggername);

	// Loop over our interested HLT trigger names to find if this event contains. 
	for (unsigned int it=0; it<TriggerNames_.size(); it++){
	  if (triggername.find(TriggerNames_[it]) != string::npos) {
	    // save the no versioned case 
	    triggernames->push_back(TriggerNames_[it]); 
	    triggerprescales->push_back(triggerprescale); 

	  }}}}}
}


bool
BToKMuMu::hasBeamSpot(const edm::Event& iEvent)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);
  
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ; 
    return false; 
  }
  
  beamSpot_ = *beamSpotHandle; 
  return true; 
}


bool
BToKMuMu::hasPrimaryVertex(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(VertexLabel_, recVtxs);
  nprivtx = recVtxs->size(); 

  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtxs->begin();
       iVertex != recVtxs->end(); iVertex++) { 
    primaryVertex_ = *(iVertex); 
    if (primaryVertex_.isValid()) break; 
  }

  if (!primaryVertex_.isValid()) return false; 
 
  return true; 
}

bool 
BToKMuMu::buildBuToKMuMu(const edm::Event& iEvent)
{
  // init variables 
  edm::Handle< vector<pat::Muon> > patMuonHandle;
  iEvent.getByLabel(MuonLabel_, patMuonHandle);
  if( patMuonHandle->size() < 2 ) return false;
  
  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr;

  // ---------------------------------
  // loop 1: mu-
  // ---------------------------------
  for (vector<pat::Muon>::const_iterator iMuonM = patMuonHandle->begin(); 
       iMuonM != patMuonHandle->end(); iMuonM++){
    
    reco::TrackRef muTrackm = iMuonM->innerTrack(); 
    if ( muTrackm.isNull() ) continue; 
    
    if ( (muTrackm->charge() != -1) ||
	 (muTrackm->pt() < MuonMinPt_) ||
	 (fabs(muTrackm->eta()) > MuonMaxEta_)) continue;
    
    // check mu- DCA to beam spot 
    const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle_));   
    
    if ( ! hasGoodMuonDcaBs(muTrackmTT, DCAmumBS, DCAmumBSErr) ) 
      continue ;
    
    nb++; 
  } // close mu- loop 
  
  if ( nb > 0) {
    edm::LogInfo("myBu") << "Found " << nb << " Bu -> K mu mu."; 
    return true; 
  }

  return false; 
}


bool 
BToKMuMu::hasGoodMuonDcaBs (const reco::TransientTrack muTrackTT, 
			    double &muDcaBs, double &muDcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS = 
    muTrackTT.trajectoryStateClosestToPoint(
     GlobalPoint(beamSpot_.position().x(),
		 beamSpot_.position().y(),beamSpot_.position().z()));
  
  if ( !theDCAXBS.isValid() )  return false; 
  
  muDcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  muDcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( muDcaBs > MuonMaxDcaBs_ )   return false; 
  return true; 
}


//define this as a plug-in
DEFINE_FWK_MODULE(BToKMuMu);
