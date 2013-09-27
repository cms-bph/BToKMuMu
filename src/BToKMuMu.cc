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
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>


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
  bool calClosestApproachTracks(const reco::TransientTrack,
				const reco::TransientTrack, 
				double&, double &, double &);
  void calCosAlpha (double, double, double, double, double,
		    double, double, double, double, double,
		    double, double, double, double,
		    double, double, double, double,
		    double*, double*); 
  void calLS (double, double, double, double, double, double, double, 
	      double, double,  double, double, double, double, double, 
	      double, double, double, double, double*, double*); 

  void clearVariables(); 
  bool hasBeamSpot(const edm::Event&);
  bool hasGoodBuMass(RefCountedKinematicTree, double &); 
  bool hasGoodBuVertex(const reco::TrackRef, const reco::TrackRef, const reco::TrackRef, 
		       double &, double &, RefCountedKinematicTree &); 
  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &); 
  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
			  reco::TransientTrack &, reco::TransientTrack &, 
			  double &, double &, double &, double &, double &,
			  double &, double &, double &);
  bool hasGoodTrack(const edm::Event&, const pat::GenericParticle, double &);
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &); 
  bool hasPrimaryVertex(const edm::Event &); 
  void hltReport(const edm::Event&);
  bool matchMuonTrack (const edm::Event&, const reco::TrackRef);
  void saveDimuVariables(double, double, double, double, double, double, 
			 double, double, double, double, double, double,
			 double, double);
  void saveSoftMuonVariables(pat::Muon, pat::Muon, reco::TrackRef, reco::TrackRef); 

  // ----------member data ---------------------------
  // --- begin input from python file --- 
  string OutputFileName_; 
  bool BuildBuToKMuMu_; 

  // particle properties 
  ParticleMass MuonMass_; 
  float MuonMassErr_; 
  ParticleMass KaonMass_; 
  float KaonMassErr_; 
  double BuMass_; 

  // labels 
  edm::InputTag TriggerResultsLabel_;
  edm::InputTag BeamSpotLabel_;
  edm::InputTag VertexLabel_;
  edm::InputTag MuonLabel_;
  edm::InputTag TrackLabel_;

  vector<string> TriggerNames_; 
  vector<string> LastFilterNames_;

  // pre-selection cuts
  double MuonMinPt_; 
  double MuonMaxEta_; 
  double MuonMaxDcaBs_; 
  double TrkMinPt_; 
  double TrkMaxDcaSigBs_; 
  double TrkMaxR_;
  double TrkMaxZ_; 
  double MuMuMaxDca_; 
  double MuMuMinVtxCl_; 
  double MuMuMinPt_; 
  double MuMuMinInvMass_; 
  double MuMuMaxInvMass_; 
  double MuMuMinLxySigmaBs_; 
  double MuMuMinCosAlphaBs_; 
  double BMinVtxCl_; 
  double BMinMass_; 
  double BMaxMass_; 

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

  // dimuon 
  vector<double> *mumdcabs, *mumdcabserr, *mumpx, *mumpy, *mumpz; 
  vector<double> *mupdcabs, *mupdcabserr, *muppx, *muppy, *muppz; 
  vector<double> *mumutrkr, *mumutrkz , *mumudca; 
  vector<double> *mumuvtxcl, *mumulsbs, *mumulsbserr;  
  vector<double> *mumucosalphabs, *mumucosalphabserr; 
  vector<double> *mumumass, *mumumasserr; 

  // soft muon variables 
  vector<bool>   *mumisgoodmuon, *mupisgoodmuon ; 
  vector<int>    *mumnpixhits, *mupnpixhits, *mumnpixlayers, *mupnpixlayers; 
  vector<int>    *mumntrkhits, *mupntrkhits, *mumntrklayers, *mupntrklayers; 
  vector<double> *mumnormchi2, *mupnormchi2; 
  vector<double> *mumdxyvtx, *mupdxyvtx, *mumdzvtx, *mupdzvtx; 
  vector<string> *mumtriglastfilter, *muptriglastfilter; 
  vector<double> *mumpt, *muppt, *mumeta, *mupeta; 
 
  // kaon track 
  vector<int> *trkchg; // +1 for K+, -1 for K-
  vector<double> *trkpx, *trkpy, *trkpz, *trkpt; 
  vector<double> *trkdcabs, *trkdcabserr; 


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

  // particle properties 
  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
  KaonMass_(iConfig.getUntrackedParameter<double>("KaonMass")),
  KaonMassErr_(iConfig.getUntrackedParameter<double>("KaonMassErr")),
  BuMass_(iConfig.getUntrackedParameter<double>("BuMass")),

  // labels 
  TriggerResultsLabel_(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")),
  MuonLabel_(iConfig.getParameter<edm::InputTag>("MuonLabel")),
  TrackLabel_(iConfig.getParameter<edm::InputTag>("TrackLabel")),

  // pre-selection cuts 
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),

  TrkMinPt_(iConfig.getUntrackedParameter<double>("TrkMinPt")),
  TrkMaxDcaSigBs_(iConfig.getUntrackedParameter<double>("TrkMaxDcaSigBs")),
  TrkMaxR_(iConfig.getUntrackedParameter<double>("TrkMaxR")),
  TrkMaxZ_(iConfig.getUntrackedParameter<double>("TrkMaxZ")),

  MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
  MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
  MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
  MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")),
  MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),
  MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")), 
  MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")), 
  BMinVtxCl_(iConfig.getUntrackedParameter<double>("BMinVtxCl")),
  BMinMass_(iConfig.getUntrackedParameter<double>("BMinMass")),
  BMaxMass_(iConfig.getUntrackedParameter<double>("BMaxMass")),
 
  // root variables  
  tree_(0), 
  triggernames(0), triggerprescales(0), 

  mumdcabs(0), mumdcabserr(0), mumpx(0), mumpy(0), mumpz(0), 
  mupdcabs(0),  mupdcabserr(0), muppx(0),  muppy(0), muppz(0), 
  mumutrkr(0), mumutrkz(0), mumudca(0),  mumuvtxcl(0),  mumulsbs(0), 
  mumulsbserr(0), mumucosalphabs(0),  mumucosalphabserr(0), 
  mumumass(0), mumumasserr(0), 
  mumisgoodmuon(0), mupisgoodmuon(0), 
  mumnpixhits(0), mupnpixhits(0), mumnpixlayers(0), mupnpixlayers(0), 
  mumntrkhits(0), mupntrkhits(0), mumntrklayers(0), mupntrklayers(0), 
  mumnormchi2(0), mupnormchi2(0), mumdxyvtx(0), mupdxyvtx(0),
  mumdzvtx(0), mupdzvtx(0), mumtriglastfilter(0), muptriglastfilter(0), 
  mumpt(0), muppt(0), mumeta(0), mupeta(0), 

  trkchg(0), trkpx(0), trkpy(0), trkpz(0), trkpt(0),  
  trkdcabs(0), trkdcabserr(0), 

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

  tree_->Branch("mumdcabs", &mumdcabs);
  tree_->Branch("mumdcabserr", &mumdcabserr);
  tree_->Branch("mumpx", &mumpx);
  tree_->Branch("mumpy", &mumpy);
  tree_->Branch("mumpz", &mumpz);
  tree_->Branch("mupdcabs", &mupdcabs);
  tree_->Branch("mupdcabserr", &mupdcabserr);
  tree_->Branch("muppx", &muppx);
  tree_->Branch("muppy", &muppy);
  tree_->Branch("muppz", &muppz);
  tree_->Branch("mumutrkr", &mumutrkr);
  tree_->Branch("mumutrkz", &mumutrkz);
  tree_->Branch("mumudca", &mumudca);
  tree_->Branch("mumuvtxcl", &mumuvtxcl);
  tree_->Branch("mumulsbs", &mumulsbs);
  tree_->Branch("mumulsbserr", &mumulsbserr);
  tree_->Branch("mumucosalphabs", &mumucosalphabs);
  tree_->Branch("mumucosalphabserr", &mumucosalphabserr);
  tree_->Branch("mumumass", &mumumass);
  tree_->Branch("mumumasserr", &mumumasserr);
  tree_->Branch("mumisgoodmuon", &mumisgoodmuon);
  tree_->Branch("mupisgoodmuon", &mupisgoodmuon);
  tree_->Branch("mumnpixhits", &mumnpixhits);
  tree_->Branch("mupnpixhits", &mupnpixhits);
  tree_->Branch("mumnpixlayers", &mumnpixlayers);
  tree_->Branch("mupnpixlayers", &mupnpixlayers);
  tree_->Branch("mumntrkhits", &mumntrkhits);
  tree_->Branch("mupntrkhits", &mupntrkhits);
  tree_->Branch("mumntrklayers", &mumntrklayers);
  tree_->Branch("mupntrklayers", &mupntrklayers);
  tree_->Branch("mumnormchi2", &mumnormchi2);
  tree_->Branch("mupnormchi2", &mupnormchi2);
  tree_->Branch("mumdxyvtx", &mumdxyvtx);
  tree_->Branch("mupdxyvtx", &mupdxyvtx);
  tree_->Branch("mumdzvtx", &mumdzvtx);
  tree_->Branch("mupdzvtx", &mupdzvtx);
  tree_->Branch("mumtriglastfilter", &mumtriglastfilter);
  tree_->Branch("muptriglastfilter", &muptriglastfilter);
  tree_->Branch("mumpt", &mumpt);
  tree_->Branch("muppt", &muppt);
  tree_->Branch("mumeta", &mumeta);
  tree_->Branch("mupeta", &mupeta);

  tree_->Branch("trkchg", &trkchg);
  tree_->Branch("trkpx", &trkpx);
  tree_->Branch("trkpy", &trkpy);
  tree_->Branch("trkpz", &trkpz);
  tree_->Branch("trkpt", &trkpt);
  tree_->Branch("trkdcabs", &trkdcabs);
  tree_->Branch("trkdcabserr", &trkdcabserr);

  tree_->Branch("nb", &nb, "nb/I");
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

  mumdcabs->clear();  mumdcabserr->clear();  mumpx->clear();   mumpy->clear();  
  mumpz->clear(); 
  mupdcabs->clear();  mupdcabserr->clear();  muppx->clear();   muppy->clear(); 
  muppz->clear(); 
  mumutrkr->clear(); mumutrkz->clear(); 
  mumudca->clear();  mumuvtxcl->clear();   mumulsbs->clear();  mumulsbserr->clear(); 
  mumucosalphabs->clear();  mumucosalphabserr->clear();
  mumumass->clear(); mumumasserr->clear(); 
  mumisgoodmuon->clear();  mupisgoodmuon->clear(); 
  mumnpixhits->clear();  mupnpixhits->clear();  mumnpixlayers->clear(); 
  mupnpixlayers->clear(); 
  mumntrkhits->clear();  mupntrkhits->clear();  mumntrklayers->clear();
  mupntrklayers->clear(); 
  mumnormchi2->clear(); mupnormchi2->clear(); mumdxyvtx->clear(); mupdxyvtx->clear(); 
  mumdzvtx->clear(); mupdzvtx->clear();  mumtriglastfilter->clear();
  muptriglastfilter->clear(); 
  mumpt->clear(); muppt->clear(); 
  mumeta->clear(); mupeta->clear(); 

  trkchg->clear(); trkpx->clear(); trkpy->clear(); trkpz->clear(); 
  trkpt->clear();
  trkdcabs->clear(); trkdcabserr->clear(); 

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
 
  edm::Handle< vector<pat::GenericParticle> >thePATTrackHandle;
  iEvent.getByLabel(TrackLabel_, thePATTrackHandle);
    
  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr;
  double mumutrk_R, mumutrk_Z, DCAmumu; 
  reco::TransientTrack refitMupTT, refitMumTT; 
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err; 
  double MuMuLSBS, MuMuLSBSErr; 
  double MuMuCosAlphaBS, MuMuCosAlphaBSErr;
  double trk_pt, DCATrkBS, DCATrkBSErr; 
  double b_vtx_chisq, b_vtx_cl, b_mass; 
  RefCountedKinematicTree vertexFitTree; 

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
    
    // ---------------------------------
    // loop 2: mu+ 
    // ---------------------------------
    for (vector<pat::Muon>::const_iterator iMuonP = patMuonHandle->begin(); 
	 iMuonP != patMuonHandle->end(); iMuonP++){

      reco::TrackRef muTrackp = iMuonP->innerTrack(); 
      if ( muTrackp.isNull() || 
	   (muTrackp->charge() != 1) ||
	   (muTrackp->pt() < MuonMinPt_) ||
	   (fabs(muTrackp->eta()) > MuonMaxEta_)) continue;
      
      // check mu+ DCA to beam spot 
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle_)); 
      if ( ! hasGoodMuonDcaBs(muTrackpTT, DCAmupBS, DCAmupBSErr) ) continue; 
      
      if (! calClosestApproachTracks(muTrackpTT, muTrackmTT,
				     mumutrk_R, mumutrk_Z, DCAmumu)) continue ; 
      
      if ( mumutrk_R > TrkMaxR_ ||
	   mumutrk_Z > TrkMaxZ_ || 
	   DCAmumu > MuMuMaxDca_ ) continue; 
      
      // check dimuon vertex 
      if ( ! hasGoodMuMuVertex(muTrackpTT, muTrackmTT, refitMupTT, refitMumTT, 
			       mu_mu_vtx_cl, mu_mu_pt, 
			       mu_mu_mass, mu_mu_mass_err, 
			       MuMuLSBS, MuMuLSBSErr,
			       MuMuCosAlphaBS, MuMuCosAlphaBSErr) ) continue; 
      
      // ---------------------------------
      // loop 3: hadron track 
      // ---------------------------------
      for ( vector<pat::GenericParticle>::const_iterator iTrack
	      = thePATTrackHandle->begin();
	    iTrack != thePATTrackHandle->end(); ++iTrack ) {
	
	if (! hasGoodTrack(iEvent, *iTrack, trk_pt) ) continue; 

	// compute track DCA to beam spot 
	reco::TrackRef theTrack = iTrack->track(); 
	const reco::TransientTrack theTrackTT(theTrack, &(*bFieldHandle_));   
	if (! hasGoodTrackDcaBs(theTrackTT, DCATrkBS, DCATrkBSErr) ) continue ; 

	if (! hasGoodBuVertex(muTrackm, muTrackp, theTrack,
			      b_vtx_chisq, b_vtx_cl, vertexFitTree) ) continue ; 

	if (! hasGoodBuMass(vertexFitTree, b_mass) ) continue;
 
	nb++; 

	// save the tree variables 
	saveDimuVariables(DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr,
			  mumutrk_R, mumutrk_Z, DCAmumu, mu_mu_vtx_cl, 
			  MuMuLSBS, MuMuLSBSErr, 
			  MuMuCosAlphaBS, MuMuCosAlphaBSErr, 
			  mu_mu_mass, mu_mu_mass_err); 

	saveSoftMuonVariables(*iMuonM, *iMuonP, muTrackm, muTrackp); 
	trkpt->push_back(trk_pt); 
	trkdcabs->push_back(DCATrkBS); 
	trkdcabserr->push_back(DCATrkBSErr); 


      } // close track loop
    } // close mu+ loop
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

bool
BToKMuMu::calClosestApproachTracks (const reco::TransientTrack trackpTT, 
				    const reco::TransientTrack trackmTT,
				    double & trk_R, 
				    double & trk_Z, 
				    double & trk_DCA)
{
  ClosestApproachInRPhi ClosestApp; 
  ClosestApp.calculate(trackpTT.initialFreeState(),
		       trackmTT.initialFreeState());
  if (! ClosestApp.status() )  return false ;
  
  GlobalPoint XingPoint = ClosestApp.crossingPoint();
  
  trk_R = sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()); 
  trk_Z = fabs(XingPoint.z()); 

  trk_DCA = ClosestApp.distance();
 
  return true; 
}

bool 
BToKMuMu::hasGoodMuMuVertex (const reco::TransientTrack muTrackpTT, 
			     const reco::TransientTrack muTrackmTT, 
			     reco::TransientTrack &refitMupTT, 
			     reco::TransientTrack &refitMumTT, 
			     double & mu_mu_vtx_cl, double & mu_mu_pt, 
			     double & mu_mu_mass, double & mu_mu_mass_err, 
			     double & MuMuLSBS, double & MuMuLSBSErr, 
			     double & MuMuCosAlphaBS,
			     double & MuMuCosAlphaBSErr)
{
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;   

  vector<RefCountedKinematicParticle> muonParticles;
  double chi = 0.;
  double ndf = 0.;
  muonParticles.push_back(partFactory.particle(muTrackmTT, 
					       MuonMass_,chi,ndf,MuonMassErr_));
  muonParticles.push_back(partFactory.particle(muTrackpTT,
					       MuonMass_,chi,ndf,MuonMassErr_));
  
  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
  
  if ( !mumuVertexFitTree->isValid())  return false;
  
  mumuVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
  RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();
  
  if ( !mumu_KV->vertexIsValid()) return false;
  
  mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
			     int(rint(mumu_KV->degreesOfFreedom())));
  
  if (mu_mu_vtx_cl < MuMuMinVtxCl_)  return false;

  // extract the re-fitted tracks
  mumuVertexFitTree->movePointerToTheTop();
  
  mumuVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
  refitMumTT = refitMum->refittedTransientTrack();
  
  mumuVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
  refitMupTT = refitMup->refittedTransientTrack();
  
  TLorentzVector mymum, mymup, mydimu; 
  
  mymum.SetXYZM(refitMumTT.track().momentum().x(), 
		refitMumTT.track().momentum().y(),
		refitMumTT.track().momentum().z(), MuonMass_); 

  mymup.SetXYZM(refitMupTT.track().momentum().x(),
		refitMupTT.track().momentum().y(),
		refitMupTT.track().momentum().z(), MuonMass_); 
  
  mydimu = mymum + mymup; 
  mu_mu_pt = mydimu.Perp(); 
 
  mu_mu_mass = mumu_KP->currentState().mass(); 
  mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
			matrix()(6,6));

  if ((mu_mu_pt < MuMuMinPt_) || (mu_mu_mass < MuMuMinInvMass_) ||
      (mu_mu_mass > MuMuMaxInvMass_))  return false;

  // compute the distance between mumu vtx and beam spot 
  calLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
	 beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	 mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
	 mumu_KV->error().matrix()(0,1),0.0,0.0,
	 beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	 beamSpot_.covariance()(0,1),0.0,0.0,
	 &MuMuLSBS,&MuMuLSBSErr);
  
  if (MuMuLSBS/MuMuLSBSErr < MuMuMinLxySigmaBs_)  return false;

  calCosAlpha(mumu_KP->currentState().globalMomentum().x(),
	      mumu_KP->currentState().globalMomentum().y(), 
	      0.0,
	      mumu_KV->position().x() - beamSpot_.position().x(),
	      mumu_KV->position().y() - beamSpot_.position().y(),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
	      mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
	      0.0,
	      mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
	      0.0,
	      0.0,
	      mumu_KV->error().cxx() + beamSpot_.covariance()(0,0),
	      mumu_KV->error().cyy() + beamSpot_.covariance()(1,1),
	      0.0,
	      mumu_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
	      0.0,
	      0.0,
	      &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);	  
  
  if (MuMuCosAlphaBS < MuMuMinCosAlphaBs_)  return false;

  return true; 
}

void 
BToKMuMu::calCosAlpha (double Vx, double Vy, double Vz,
		       double Wx, double Wy, double Wz,
		       double VxErr2, double VyErr2, double VzErr2,
		       double VxyCov, double VxzCov, double VyzCov,
		       double WxErr2, double WyErr2, double WzErr2,
		       double WxyCov, double WxzCov, double WyzCov,
		       double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;
  
  if ((Vnorm > 0.) && (Wnorm > 0.)) {
      *cosAlpha = VdotW / (Vnorm * Wnorm);
      *cosAlphaErr = sqrt( (
       (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
       (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
       (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +
       
       (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
       (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
       (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +
			   
			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			    
			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }
}

void 
BToKMuMu::calLS (double Vx, double Vy, double Vz,
		 double Wx, double Wy, double Wz,
		 double VxErr2, double VyErr2, double VzErr2,
		 double VxyCov, double VxzCov, double VyzCov,
		 double WxErr2, double WyErr2, double WzErr2,
		 double WxyCov, double WxzCov, double WyzCov,
		 double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}


bool 
BToKMuMu::hasGoodTrack(const edm::Event& iEvent, 
		       const pat::GenericParticle iTrack, 
		       double & trk_pt)
{
   reco::TrackRef theTrackRef = iTrack.track(); 
   if ( theTrackRef.isNull() ) return false; 

   // veto muon tracks
   if ( matchMuonTrack(iEvent, theTrackRef) ) return false; 
   
   // check the track kinematics
   trk_pt = theTrackRef->pt(); 

   if ( theTrackRef->pt() < TrkMinPt_ ) return false; 

   return true;
}

bool
BToKMuMu::matchMuonTrack (const edm::Event& iEvent, 
			  const reco::TrackRef theTrackRef)
{
  if ( theTrackRef.isNull() ) return false;

  edm::Handle< vector<pat::Muon> > thePATMuonHandle; 
  iEvent.getByLabel(MuonLabel_, thePATMuonHandle);
  
  reco::TrackRef muTrackRef; 
  for (vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin();
       iMuon != thePATMuonHandle->end(); iMuon++){

    muTrackRef = iMuon->innerTrack();
    if ( muTrackRef.isNull() ) continue; 

    if (muTrackRef == theTrackRef) return true; 
  }
  
  return false; 
}


bool 
BToKMuMu::hasGoodTrackDcaBs (const reco::TransientTrack TrackTT, 
			     double &DcaBs, double &DcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS = 
    TrackTT.trajectoryStateClosestToPoint(
     GlobalPoint(beamSpot_.position().x(),
		 beamSpot_.position().y(),beamSpot_.position().z()));
  
  if ( !theDCAXBS.isValid() )  return false; 
  
  DcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  DcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( DcaBs/DcaBsErr > TrkMaxDcaSigBs_ )   return false; 
  return true; 
}

bool 
BToKMuMu::hasGoodBuVertex(const reco::TrackRef mu1Track,
			      const reco::TrackRef mu2Track,
			      const reco::TrackRef theTrack, 
			      double & b_vtx_chisq, double & b_vtx_cl, 
			      RefCountedKinematicTree &vertexFitTree) 
{
  // do vertex fit for Bu
  KinematicParticleFactoryFromTransientTrack pFactory;
  reco::TransientTrack mu1TT(mu1Track, &(*bFieldHandle_) );
  reco::TransientTrack mu2TT(mu2Track, &(*bFieldHandle_) );
  reco::TransientTrack trkTT(theTrack, &(*bFieldHandle_) );

  float chi = 0.;
  float ndf = 0.;
  vector<RefCountedKinematicParticle> vFitMCParticles;
  vFitMCParticles.push_back(pFactory.particle(mu1TT, MuonMass_,
					      chi, ndf, MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(mu2TT, MuonMass_,
					      chi, ndf, MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(trkTT, KaonMass_, chi, 
					      ndf, KaonMassErr_));
  KinematicParticleVertexFitter fitter;   
  vertexFitTree = fitter.fit(vFitMCParticles);
  if (!vertexFitTree->isValid()) return false; 

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
  if ( !bDecayVertexMC->vertexIsValid()) return false; 

  b_vtx_chisq = bDecayVertexMC->chiSquared(); 
  if ( bDecayVertexMC->chiSquared()<0
       || bDecayVertexMC->chiSquared()>1000 ) return false; 

  RefCountedKinematicVertex b_KV = vertexFitTree->currentDecayVertex();
  b_vtx_cl = ChiSquaredProbability((double)(b_KV->chiSquared()),
				   (double)(b_KV->degreesOfFreedom()));

  if ( b_vtx_cl < BMinVtxCl_ ) return false; 
  
  return true; 
}


bool 
BToKMuMu::hasGoodBuMass(RefCountedKinematicTree vertexFitTree, 
			double & b_mass)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  b_mass = b_KP->currentState().mass(); 
  if ( b_mass < BMinMass_ || b_mass > BMaxMass_ ) return false;  
  return true; 
}


void 
BToKMuMu::saveDimuVariables(double DCAmumBS, double DCAmumBSErr, 
			    double DCAmupBS, double DCAmupBSErr,
			    double mumutrk_R, double mumutrk_Z, 
			    double DCAmumu,  double mu_mu_vtx_cl, 
			    double MuMuLSBS, double MuMuLSBSErr, 
			    double MuMuCosAlphaBS, double MuMuCosAlphaBSErr,
			    double mu_mu_mass, double mu_mu_mass_err) 
  
{
  mumdcabs->push_back(DCAmumBS);
  mumdcabserr->push_back(DCAmumBSErr);

  mupdcabs->push_back(DCAmupBS);
  mupdcabserr->push_back(DCAmupBSErr);

  mumutrkr->push_back(mumutrk_R);
  mumutrkz->push_back(mumutrk_Z); 
  mumudca->push_back(DCAmumu);
  mumuvtxcl->push_back(mu_mu_vtx_cl); 
  mumulsbs->push_back(MuMuLSBS);
  mumulsbserr->push_back(MuMuLSBSErr);
  mumucosalphabs->push_back(MuMuCosAlphaBS);
  mumucosalphabserr->push_back(MuMuCosAlphaBSErr); 

  mumumass->push_back(mu_mu_mass); 
  mumumasserr->push_back(mu_mu_mass_err); 
}


void 
BToKMuMu::saveSoftMuonVariables(pat::Muon iMuonM, pat::Muon iMuonP, 
				reco::TrackRef muTrackm, reco::TrackRef muTrackp)
{
  mumisgoodmuon->push_back(muon::isGoodMuon(iMuonM, muon::TMOneStationTight)); 
  mupisgoodmuon->push_back(muon::isGoodMuon(iMuonP, muon::TMOneStationTight)); 
  mumnpixhits->push_back(muTrackm->hitPattern().numberOfValidPixelHits()); 
  mupnpixhits->push_back(muTrackp->hitPattern().numberOfValidPixelHits()); 
  mumnpixlayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement()); 
  mupnpixlayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement()); 

  mumntrkhits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits()); 
  mupntrkhits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits()); 
  mumntrklayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement()); 
  mupntrklayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement()); 

  mumnormchi2->push_back(muTrackm->normalizedChi2()); 
  mupnormchi2->push_back(muTrackp->normalizedChi2()); 

  mumdxyvtx->push_back(muTrackm->dxy(primaryVertex_.position())); 
  mupdxyvtx->push_back(muTrackp->dxy(primaryVertex_.position())); 
	  
  mumdzvtx->push_back(muTrackm->dz(primaryVertex_.position())); 
  mupdzvtx->push_back(muTrackp->dz(primaryVertex_.position())); 

  mumpt->push_back(muTrackm->pt()); 
  muppt->push_back(muTrackp->pt()); 
  
  mumeta->push_back(muTrackm->eta()); 
  mupeta->push_back(muTrackp->eta()); 
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(BToKMuMu);
