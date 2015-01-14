// -----------------------------------------------
// Author: Xin Shi <Xin.Shi@cern.ch>
// Created: <2013-02-26 Tue 12:08>
// -----------------------------------------------
#define SingleBuToKMuMuSelector_cxx

#include <iostream>
#include <sstream>
#include <map>
#include "SingleBuToKMuMuSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TProof.h>
#include <TLorentzVector.h>

// Global Constants
const double K_MASS    = 0.89166; // GeV
const double K_WIDTH   = 0.0508; // GeV
const double MUON_MASS = 0.10565837;
const double KAON_MASS = 0.493677;

// user defined variables
TDatime t_begin_ , t_now_ ;
int n_processed_, n_selected_;
int ww=0;

TTree *tree_;

// Branch variables for new tree
//int     Nb             = -999;

TLorentzVector Mum_4vec   ;  //N.A.
TLorentzVector Mup_4vec   ;  //N.A.
double  Mumumass          = -999;
double  Mumumasserr       = -999;

TLorentzVector Tk_4vec    ;  //N.A.
int     TkChg             = -999;
double  Trkdcasigbs       = -999;  

TLorentzVector B_4vec     ;  //N.A.
double  Bmass             = -999;
int     Bchg              = -999;
double  Bvtxcl            = -999;
double  Blxysig           = -999;
double  Bcosalphabs       = -999;
double  Bcosalphabs2D     = -999;
double  Bctau             = -999;

double  Q2                = -999;
double  CosThetaL         = -999;

// Branches for Generator level information
int     genBChg        = -999;
double  genBPt         = -999;
double  genBEta        = -999;
double  genBPhi        = -999;
double  genBVtxX       = -999;
double  genBVtxY       = -999;
double  genBVtxZ       = -999;
double  genMupPt       = -999;
double  genMupEta      = -999;
double  genMupPhi      = -999;
double  genMumPt       = -999;
double  genMumEta      = -999;
double  genMumPhi      = -999;
int     genTkChg       = -999;
double  genTkPt        = -999;
double  genTkEta       = -999;
double  genTkPhi       = -999;
double  genQ2          = -999;
double  genCosThetaL   = -999;
/*
//TLorentzVector genB_4vec  = -999;  //N.A.

int     genBChg           = -999;
double  genBVtxX          = -999;
double  genBVtxY          = -999;
double  genBVtxZ          = -999;

//TLorentzVector genMum_4vec= -999;  //N.A.
//TLorentzVector genMup_4vec= -999;  //N.A.

//TLorentzVector genTk_4vec = -999;  //N.A.
int     genTkChg          = -999;

double  genQ2             = -999;
double  genCosThetaL      = -999;
*/

void ClearEvent() { 
//	Nb                  = -999;
//	Mum_4vec            ;
//	Mup_4vec            ;
	Mumumass            = -999;
	Mumumasserr         = -999;

//	Tk_4vec             ;
	TkChg               = -999;
	Trkdcasigbs         = -999;   
	
//	B_4vec              ;
	Bmass               = -999;
	Bchg                = -999;
	Bvtxcl              = -999;
	Blxysig             = -999;
	Bcosalphabs         = -999;
	Bcosalphabs2D       = -999; 
	Bctau               = -999;
	
	Q2                  = -999;
	CosThetaL           = -999;

//	GenLevel
	genBChg             = -999;
	genBPt              = -999;
	genBEta             = -999;
	genBPhi             = -999;
	genBVtxX            = -999;
	genBVtxY            = -999;
	genBVtxZ            = -999;
	genMupPt            = -999;
	genMupEta           = -999;
	genMupPhi           = -999;
	genMumPt            = -999;
	genMumEta           = -999;
	genMumPhi           = -999;
	genTkChg            = -999;
	genTkPt             = -999;
	genTkEta            = -999;
	genTkPhi            = -999;
	genQ2               = -999;
	genCosThetaL        = -999;
/*
   genB_4vec           = -999;
	genBChg             = -999;
	genBVtxX            = -999;
	genBVtxY            = -999;
	genBVtxZ            = -999;
	
	genMum_4vec         = -999;
	genMup_4vec         = -999;
	
	genTk_4vec          = -999;
	genTkChg            = -999;
	
	genQ2               = -999;
	genCosThetaL        = -999;
*/
}

void str_replace(std::string& str, const std::string& oldStr, const std::string& newStr) {
	size_t pos = 0;
	while((pos = str.find(oldStr, pos)) != std::string::npos) {
		str.replace(pos, oldStr.length(), newStr);
		pos += newStr.length();
	}
}

string get_option_value(string option, string name) {
	vector<string> args;
	istringstream f(option);
	string s;
	while (getline(f, s, ';')) {
		args.push_back(s);
	}
	
	string value;
	for(vector<string>::iterator it = args.begin(); it != args.end(); ++it) {
		value = *it;
		unsigned found = value.find(name);
		if (found == 0) {
			str_replace(value, name+"=", "");
			break;
		}
	}
	return value;
}


void SingleBuToKMuMuSelector::Begin(TTree * /*tree*/) {
	t_begin_.Set();
	printf("\n ---------- Begin Job ---------- \n");
	t_begin_.Print();
	
	n_processed_   = 0;
	n_selected_    = 0;
}

void SingleBuToKMuMuSelector::SlaveBegin(TTree * /*tree*/) {
	string option = GetOption();
	tree_ = new TTree("tree", "tree");
//	tree_->Branch("Nb      ",        &Nb      ,         "Nb      /I");  //N.A.
	
	tree_->Branch("Mum_4vec"   ,     &Mum_4vec   ,      32000  ,  0 );  //N.A.
	tree_->Branch("Mup_4vec"   ,     &Mup_4vec   ,      32000  ,  0 );  //N.A.
	tree_->Branch("Mumumass"   ,     &Mumumass   ,      "Mumumass/D");
	tree_->Branch("Mumumasserr",     &Mumumasserr,      "Mumumasserr/D");
	
	tree_->Branch("Tk_4vec"    ,     &Tk_4vec    ,      32000  ,  0 );  //N.A.
	tree_->Branch("TkChg"      ,     &TkChg      ,      "TkChg/I");  
	tree_->Branch("Trkdcasigbs",     &Trkdcasigbs,      "Trkdcasigbs/D");  
	
	tree_->Branch("B_4vec"     ,     &B_4vec     ,      32000  ,  0 );  //N.A.
	tree_->Branch("Bmass",           &Bmass,            "Bmass/D");
	tree_->Branch("Bchg",            &Bchg,             "Bchg/I");
	tree_->Branch("Bvtxcl",          &Bvtxcl,           "Bvtxcl/D");
	tree_->Branch("Blxysig",         &Blxysig,          "Blxysig/D");
	tree_->Branch("Bcosalphabs",     &Bcosalphabs,      "Bcosalphabs/D");
	tree_->Branch("Bcosalphabs2D",   &Bcosalphabs2D ,   "Bcosalphabs2D/D"); 
	tree_->Branch("Bctau",           &Bctau,            "Bctau/D");
	
	tree_->Branch("Q2" ,             &Q2 ,              "Q2/D");
	tree_->Branch("CosThetaL" ,      &CosThetaL ,       "CosThetaL/D");
//////////////////////////////////////////////////////////////////////////
	string datatype = get_option_value(option, "datatype");
	std::map<string,int> maptype;
	maptype.insert(std::pair<string,int>("data",1));
	maptype.insert(std::pair<string,int>("mc.lite",2));
	maptype.insert(std::pair<string,int>("mc.AOD",2));
	maptype.insert(std::pair<string,int>("mc.hlt",998));
	maptype.insert(std::pair<string,int>("mc",999));
	switch (maptype[datatype]) {
		case 1: break;
		case 2:
		   tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
			tree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
			tree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
			tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
			tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
			tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
			tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
			tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
			tree_->Branch("genTkChg"     , &genTkChg     , "genTkChg/I");
			tree_->Branch("genTkPt"      , &genTkPt      , "genTkPt/D");
			tree_->Branch("genTkEta"     , &genTkEta     , "genTkEta/D");
			tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
			tree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
		/*
			tree_->Branch("genB_4vec"    , &genB_4vec    ,  32000  ,  0 );  //N.A.
		   tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
			tree_->Branch("genMum_4vec"  , &genMum_4vec  ,  32000  ,  0 );  //N.A.
			tree_->Branch("genMup_4vec"  , &genMup_4vec  ,  32000  ,  0 );  //N.A.
			tree_->Branch("genTk_4vec"   , &genTk_4vec   ,  32000  ,  0 );  //N.A.
		   tree_->Branch("genTkChg"     , &genTkChg     , "genTkChg/I");
			tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
			tree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
		*/	
			break;
		case 998:
		   tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
			tree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
			tree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
			tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
			tree_->Branch("genBVtxX"     , &genBVtxX     , "genBVtxX/D");
			tree_->Branch("genBVtxY"     , &genBVtxY     , "genBVtxY/D");
			tree_->Branch("genBVtxZ"     , &genBVtxZ     , "genBVtxZ/D");
			tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
			tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
			tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
			tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
			tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
			tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");
			tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
			break;
		case 999: 
		   tree_->Branch("genBChg"      , &genBChg      , "genBChg/I");
			tree_->Branch("genBPt"       , &genBPt       , "genBPt/D");
			tree_->Branch("genBEta"      , &genBEta      , "genBEta/D");
			tree_->Branch("genBPhi"      , &genBPhi      , "genBPhi/D");
			tree_->Branch("genBVtxX"     , &genBVtxX     , "genBVtxX/D");
			tree_->Branch("genBVtxY"     , &genBVtxY     , "genBVtxY/D");
			tree_->Branch("genBVtxZ"     , &genBVtxZ     , "genBVtxZ/D");
			tree_->Branch("genMupPt"     , &genMupPt     , "genMupPt/D");
			tree_->Branch("genMupEta"    , &genMupEta    , "genMupEta/D");
			tree_->Branch("genMupPhi"    , &genMupPhi    , "genMupPhi/D");
			tree_->Branch("genMumPt"     , &genMumPt     , "genMumPt/D");
			tree_->Branch("genMumEta"    , &genMumEta    , "genMumEta/D");
			tree_->Branch("genMumPhi"    , &genMumPhi    , "genMumPhi/D");
			tree_->Branch("genTkChg"     , &genTkChg     , "genTkChg/I");
			tree_->Branch("genTkPt"      , &genTkPt      , "genTkPt/D");
			tree_->Branch("genTkEta"     , &genTkEta     , "genTkEta/D");
			tree_->Branch("genTkPhi"     , &genTkPhi     , "genTkPhi/D");
			tree_->Branch("genQ2"        , &genQ2        , "genQ2/D");
			tree_->Branch("genCosThetaL" , &genCosThetaL , "genCosThetaL/D");
         break;
		default:
		printf("No compatible datatype found. Please check use following types...\n\t\t[");
		for (std::map<string,int>::iterator iType = maptype.begin(); iType != maptype.end(); iType++) {
			if (iType->second != 0) printf("%s,",iType->first.c_str());
		}
		printf("]\n");
		break;
	}
//////////////////////////////////////////////////////////////////////  
	fOutput->AddAll(gDirectory->GetList());
}


Bool_t SingleBuToKMuMuSelector::Process(Long64_t entry) {
	ClearEvent(); // 2014-08-13 N.A.
	
	string option = GetOption();
	string datatype = get_option_value(option, "datatype");
	string cut = get_option_value(option, "cut");
	
	GetEntry(entry);
	n_processed_ += 1;
//	Nb = nb;
	
	if (datatype != "data") SaveGen();  ///// 13-08-2014 /////////////////////// for Signal MC ////////////////////////////////////////////////
	
	int i = SelectB(cut);
	if ( i != -1 && ( datatype == "data" || istruebu->at(i) )) {
		n_selected_ += 1;
		SaveEvent(i);
//		if ( datatype != "data" ) SaveGen();    /////  12-04-2014  /////  for J/PsiK & Psi2SK & J/PsiX & K*0MuMu  MC  /////////////////////
//		tree_->Fill();  /////// for data & J/PsiK & Psi2SK & J/PsiX & K*0MuMu  MC
	}
	tree_->Fill();  ////////// for Signal MC
	return kTRUE;
}

void SingleBuToKMuMuSelector::SlaveTerminate(){
}


void SingleBuToKMuMuSelector::Terminate() {
	string option = GetOption();
	TString outfile = get_option_value(option, "outfile");
	
	TFile file(outfile.Data(), "recreate");
	fOutput->Write();
	
	t_now_.Set();
	printf(" \n ---------- End Job ---------- \n" ) ;
	t_now_.Print();
	printf(" processed: %i \n selected: %i \n duration: %i sec \n rate: %g evts/sec\n",
	n_processed_,  n_selected_, 
	t_now_.Convert() - t_begin_.Convert(),
	float (n_processed_) / ( t_now_.Convert() - t_begin_.Convert() ) );

	cout<<"ww"<<ww<<endl;
}


int SingleBuToKMuMuSelector::SelectB(string cut) {
	int best_idx = -1;
	double best_bvtxcl = 0.0;
	
	if (cut == "genonly") {
		best_idx = -1;
	}else if (cut == "nocut") {
		for (int i = 0; i < nb; i++) {
			if (bvtxcl->at(i) > best_bvtxcl) {
				best_bvtxcl = bvtxcl->at(i);
				best_idx = i;
			}
		}
	//	if (nb == 0) best_idx = 0;
	//	for (int i = 0; i < nb; i++) {
	//		best_idx = i;
	//	}
	}else if (cut == "cut0") {
		for (int i = 0; i < nb; i++) {
			if ( ! HasGoodMuons(i) ) continue;
			if ( ! TriggerSelections(i) ) continue;
//			if ( ! KpSelections(i) ) continue; // 2014-05-22 N.A.
			if (bvtxcl->at(i) > best_bvtxcl) {
				best_bvtxcl = bvtxcl->at(i);
				best_idx = i;
			}
		}
	}else if (cut == "cut1") {
		for (int i = 0; i< nb; i++) {
			if ( ! HasGoodMuons(i) ) continue;
			if ( ! TriggerSelections(i) ) continue;
//			if ( ! KpSelections(i) ) continue; // 2014-05-22 N.A.
//			if ( ! OptimizedSelections(i) ) continue;
			if (bvtxcl->at(i) > best_bvtxcl) {
				best_bvtxcl = bvtxcl->at(i);
				best_idx = i;
			}
		}
	}
	return best_idx;
}

bool SingleBuToKMuMuSelector::HasGoodMuons(int i){
  if ( // soft new muon
       mumisgoodmuon->at(i)
       && mupisgoodmuon->at(i)
		 && mumntrklayers->at(i) > 5 // 2014-05-22 N.A.
		 && mupntrklayers->at(i) > 5 // 2014-05-22 N.A.
       && mumntrkhits->at(i) > 5 //10
       && mupntrkhits->at(i) > 5 //10
       && mumnpixlayers->at(i) > 0 //1
       && mupnpixlayers->at(i) > 0 //1
       && mumnormchi2->at(i) < 1.8 //3?
       && mupnormchi2->at(i) < 1.8 //3?
       && fabs(mumdxyvtx->at(i)) < 0.3 //3 // 2014-05-22 N.A.
       && fabs(mupdxyvtx->at(i)) < 0.3 //3 // 2014-05-22 N.A.
       && fabs(mumdzvtx->at(i)) < 20 //30 // 2014-05-23 N.A. 
       && fabs(mupdzvtx->at(i)) < 20 //30 // 2014-05-23 N.A.
        ) return true;
  return false;
}
bool SingleBuToKMuMuSelector::TriggerSelections(int i){
	if ( // trigger
		mumdcabs->at(i) < 0.1  //2
		&& mupdcabs->at(i) < 0.1  //2
		&& mumudca->at(i) < 0.1  //0.5
	//	&& mumuvtxcl->at(i) > 0.09  //0.05    // 27-10-2014 N.A.
	//	&& (mumulsbs->at(i)/mumulsbserr->at(i)) > 3  // 27-10-2014 N.A.
		&& mumucosalphabs->at(i) > 0.999  //0.9
	//	&& fabs( trkdcabs->at(i)/trkdcabserr->at(i) ) > 1.20  // 27-20-2014 N.A.
	//	&& trkpt->at(i) > 1.6  //2.7  ///////////////////////////////////////////////////////////////////////////////////////////////////////
		) return true;
	return false;
}
bool SingleBuToKMuMuSelector::OptimizedSelections(int i) {
/*	if ((mumumass->at(i) > 3.096916-5. *mumumasserr->at(i)) && (mumumass->at(i) < 3.096916+3.5*mumumasserr->at(i))
		) return false;
	else if (( mumumass->at(i) > 3.686109-3.5*mumumasserr->at(i) ) && (mumumass->at(i) < 3.686109+3.5*mumumasserr->at(i)) 
		) return false;
*/		
	if ( // Optimized Selections after cut0
		//   bmass->at(i) >= 5.23
		 trkpt->at(i) > 1.6  //2.7
/*		&& bvtxcl->at(i) > 0.09  //0.12
		&& ( blsbs->at(i) / blsbserr->at(i) ) > 8.2
		&& bcosalphabs2D->at(i) > 0.9997
		&& bmass->at(i) > 5.0 && bmass->at(i) < 5.56
		&& fabs( (trkdcabs->at(i))/(trkdcabserr->at(i)) )> 1.4  //1.2
	//	&& bctau->at(i) < 2.  // N.A. 12-04 
		&& ( ( mumumass->at(i) < 3.096916-5. *mumumasserr->at(i) ) || ( ( mumumass->at(i) > 3.096916+3.5*mumumasserr->at(i) ) && ( mumumass->at(i) < 3.686109-3.5*mumumasserr->at(i) ) )  || (mumumass->at(i) > 3.686109+3.5*mumumasserr->at(i)) ) // J/PSi && PSi(2S) cut
	//	&& fabs(bmass->at(i) - 5.279) < 0.060  /////////////////////////////////////////////////////////////////////////////
	//	&& fabs( bmass->at(i) - mumumass->at(i) - 2.182 ) > 0.18  //0.12 // CDF cut for JPSi
	//	&& fabs( bmass->at(i) - mumumass->at(i) - 1.593 ) > 0.09  //0.08 // CDF cut for PSi(2S)
*/		) return true;
/*	else if (
		   bmass->at(i) < 5.23
		&& trkpt->at(i) > 1.6  //2.7
		&& bvtxcl->at(i) > 0.09  //0.12
		&& ( blsbs->at(i) / blsbserr->at(i) ) > 8.2
		&& bcosalphabs2D->at(i) > 0.9997
		&& bmass->at(i) > 5.0 && bmass->at(i) < 5.56
		&& fabs( (trkdcabs->at(i))/(trkdcabserr->at(i)) )> 1.4  //1.2
	//	&& bctau->at(i) < 2.  // N.A. 12-04 
		&& ( ( mumumass->at(i) < 3.096916-5. *mumumasserr->at(i) ) || ( ( mumumass->at(i) > 3.096916+3.5*mumumasserr->at(i) ) && ( mumumass->at(i) < 3.686109-3.5*mumumasserr->at(i) ) )  || (mumumass->at(i) > 3.686109+3.5*mumumasserr->at(i)) ) // J/PSi && PSi(2S) cut
		
		&& fabs( bmass->at(i) - mumumass->at(i) - 2.182 ) > 0.14  //0.12 // CDF cut for JPSi
		&& fabs( bmass->at(i) - mumumass->at(i) - 1.593 ) > 0.08  //0.08 // CDF cut for PSi(2S)
		) return true;
*/	return false;
}
bool SingleBuToKMuMuSelector::KpSelections(int i){
	if ( // Kaon
		trkpt->at(i) > 2.3
		&& trkchg->at(i) == +1
      ) return true;
	return false;
}
void SingleBuToKMuMuSelector::SaveEvent(int i)
{
	TLorentzVector buff1, buff2, buff3;
	Mup_4vec.SetXYZM(muppx->at(i),muppy->at(i),muppz->at(i),MUON_MASS);
	Mum_4vec.SetXYZM(mumpx->at(i),mumpy->at(i),mumpz->at(i),MUON_MASS);
	Tk_4vec.SetXYZM(trkpx->at(i),trkpy->at(i),trkpz->at(i),KAON_MASS);
	
	buff2 = Mup_4vec + Mum_4vec;
	buff1 = buff2 + Tk_4vec;
	B_4vec = buff1;
	Bmass = buff1.M();
	Bchg = bchg->at(i);
	TkChg = trkchg->at(i);
	Bvtxcl = bvtxcl->at(i);
	Blxysig = (blsbs->at(i)/blsbserr->at(i));
	Bcosalphabs = bcosalphabs->at(i);
	Bcosalphabs2D = bcosalphabs2D->at(i); 
	Bctau = bctau->at(i);
	if (Tk_4vec.Pt() != trkpt->at(i)) { cout<<"Tk_4vec.Pt()="<<Tk_4vec.Pt()<<" TrkPt="<<trkpt->at(i)<<" trkpx->at(i)="<<trkpx->at(i)<<" Tk_4vec.Px()="<<Tk_4vec.Px()<<" trkpy->at(i)="<<trkpy->at(i)<<" Tk_4vec.Py()="<<Tk_4vec.Py()<<" trkpz->at(i)="<<trkpz->at(i)<<" Tk_4vec.Pz()="<<Tk_4vec.Pz()<<endl; ww++; }
//	if (B_4vec->Pt() < 8.0 ) 
//	if (bmass->at(i) == buff1.M()) { cout<<buff1.M()<<"    Bmass   "<<bmass->at(i)<<endl;}
//	if (trkpt->at(i) != Tk_4vec.Pt()) cout<<trkpt->at(i)<<" TrkP   "<<Tk_4vec.Pt()<<endl;
//	if (mumumass->at(i) != buff2.M()) cout<<mumumass->at(i)<<" mumumass   "<<buff2.M()<<endl;
	
	Mumumass = buff2.M();
	Mumumasserr = mumumasserr->at(i); 
	Trkdcasigbs = fabs( trkdcabs->at(i)/trkdcabserr->at(i) ); 
	Q2 = pow(buff2.M(),2);
//	2014-05-23 N.A.
	buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
	if ( Bchg > 0) {
		buff3 = Mum_4vec;//Take mu- to avoid extra minus sign.
	} else {
		buff3 = Mup_4vec;
	}
	buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
	CosThetaL = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();
}
void SingleBuToKMuMuSelector::SaveGen()
{
	TLorentzVector genB_4vec,  genMup_4vec, genMum_4vec, genTk_4vec, buff1, buff2, buff3;
	genMup_4vec.SetXYZM(genmuppx,genmuppy,genmuppz,MUON_MASS);
	genMum_4vec.SetXYZM(genmumpx,genmumpy,genmumpz,MUON_MASS);
	genTk_4vec.SetXYZM(gentrkpx,gentrkpy,gentrkpz,KAON_MASS);

	buff2 = genMup_4vec + genMum_4vec;
	buff1 = buff2 + genTk_4vec;
	genB_4vec = buff1;

	genBChg      = genbchg;
	genBPt       = genB_4vec.Pt();
	genBEta      = genB_4vec.Eta();
	genBPhi      = genB_4vec.Phi();
	genBVtxX     = 0;//Should be at PV?
	genBVtxY     = 0;
	genBVtxZ     = 0;
	genMupPt     = genMup_4vec.Pt();
	genMupEta    = genMup_4vec.Eta();
	genMupPhi    = genMup_4vec.Phi();
	genMumPt     = genMum_4vec.Pt();
	genMumEta    = genMum_4vec.Eta();
	genMumPhi    = genMum_4vec.Phi();
	genTkChg     = gentrkchg;
	genTkPt      = genTk_4vec.Pt();
	genTkEta     = genTk_4vec.Eta();
	genTkPhi     = genTk_4vec.Phi();
	genQ2        = (genMup_4vec+genMum_4vec).Mag2();

	buff1.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
	if (genBChg > 0){
		buff3 = genMum_4vec;//Take mu- to avoid extra minus sign.
	}else{
		buff3 = genMup_4vec;
	}
	buff3.Boost(-buff2.X()/buff2.T(),-buff2.Y()/buff2.T(),-buff2.Z()/buff2.T());
	genCosThetaL = buff1.Vect().Dot(buff3.Vect())/buff1.Vect().Mag()/buff3.Vect().Mag();

}

#ifndef __CINT__
#include <algorithm>

char* get_option(char ** begin, char ** end, const std::string & option) {
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end) return *itr;
	return 0;
} 

bool option_exists(char** begin, char** end, const std::string& option) {
	return std::find(begin, end, option) != end;
}

void print_usage() {
	cerr << "Usage: SingleBuToKMuMuSelector datatype cut infile outfile [-n] [-s] [-j] [-h]\n"
	     << " datatype: data, mc, mc.AOD, mc.lite, mc.hlt\n"  // 2014-05-23 N.A.
		  << " cut : genonly, nocut, cut0, cut1.\n"  // 2014-05-23 N.A.
		  << "Options: \n"
		  << " -h \t\tPrint this info\n"
		  << " -n \t\tNumber of entries\n"
		  << " -s \t\tStarting run number.\n" // 2014-05-23 N.A.
		  << " -j \t\tNumber of workers\n"
		  << endl;
}

int main(int argc, char** argv) {
	if ( (argc < 3) or option_exists(argv, argv+argc, "-h") ) {
		print_usage() ;
		return -1;
	}
	TString datatype = argv[1];
	TString cut = argv[2];
	TString infile = argv[3];
	TString outfile = argv[4];
	
	Printf("   datatype: '%s'", datatype.Data());
	Printf("        cut: '%s'", cut.Data());
	Printf(" input file: '%s'", infile.Data());
	Printf("output file: '%s'", outfile.Data());
   
	TString option;
	option.Form("datatype=%s;cut=%s;outfile=%s", datatype.Data(), cut.Data(), outfile.Data());
	
	TChain *ch = new TChain("tree");
	ch->Add(infile.Data());
	
	char * j = get_option(argv, argv + argc, "-j");
	if (j) {
		TProof::Open(Form("workers=%s", j));
		ch->SetProof();
	}
	
	Long64_t nentries = 1000000000;
	char * n = get_option(argv, argv+argc, "-n");
	if (n) nentries = atoi(n);
	
	int     iStart = 0;
	char *s = get_option(argv, argv+argc, "-s");
	if (s) {
		iStart = atoi(s);
		if (iStart > ch->GetEntries()) {
			printf("ERROR: Number of entries is %lld.\n",ch->GetEntries());
			return -1;
		}
	}
	
	// It's not allowed to run with fat trees!
	if (datatype.Data() == "mc" && (!(s) || !(n))) { 
		printf("WARNING: You must specify #entries(-n) and start run(-s) for datatype '%s'.\n",datatype.Data());
		return -1;
	}
	
	ch->Process("SingleBuToKMuMuSelector.cc+", option, nentries, iStart);
	gSystem->Exit(0);
	
	return 0 ;
}

#endif


