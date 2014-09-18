//////////////////////////////////////////////////////////
// [2013-07-18 Thu 10:52] Update based on BuToKMuMu_v2
//
// This class was automatically generated on
// Tue Feb 26 11:42:08 2013 by ROOT version 5.27/06b
// from TTree tree/BToKMuMu
// found on file: root://eoscms//eos/cms/store/user/xshi/dat/ntp/mc/BuToKJPsi_7TeV_5E5_v1_run2011v0_2/BToKstarMuMu_merged_1.root
//////////////////////////////////////////////////////////

#ifndef SingleBuToKMuMuSelector_h
#define SingleBuToKMuMuSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

using namespace std;

class SingleBuToKMuMuSelector : public TSelector {
public :

  TTree *fChain; //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t run;
   UInt_t event;
   UInt_t lumiblock;
   UInt_t nprivtx;
   vector<string> *triggernames;
   vector<int> *triggerprescales;
   vector<double> *mumdcabs;
   vector<double> *mumdcabserr;
   vector<double> *mumpx;
   vector<double> *mumpy;
   vector<double> *mumpz;
   vector<double> *mupdcabs;
   vector<double> *mupdcabserr;
   vector<double> *muppx;
   vector<double> *muppy;
   vector<double> *muppz;
   vector<double> *mumutrkr;
   vector<double> *mumutrkz;
   vector<double> *mumudca;
   vector<double> *mumuvtxcl;
   vector<double> *mumulsbs;
   vector<double> *mumulsbserr;
   vector<double> *mumucosalphabs;
   vector<double> *mumucosalphabserr;
   vector<double> *mumumass;
   vector<double> *mumumasserr;
   vector<bool> *mumisgoodmuon;
   vector<bool> *mupisgoodmuon;
   vector<int> *mumnpixhits;
   vector<int> *mupnpixhits;
   vector<int> *mumnpixlayers;
   vector<int> *mupnpixlayers;
   vector<int> *mumntrkhits;
   vector<int> *mupntrkhits;
   vector<int> *mumntrklayers;
   vector<int> *mupntrklayers;
   vector<double> *mumnormchi2;
   vector<double> *mupnormchi2;
   vector<double> *mumdxyvtx;
   vector<double> *mupdxyvtx;
   vector<double> *mumdzvtx;
   vector<double> *mupdzvtx;
   vector<string> *mumtriglastfilter;
   vector<string> *muptriglastfilter;
   vector<double> *mumpt;
   vector<double> *muppt;
   vector<double> *mumeta;
   vector<double> *mupeta;
   vector<int> *trkchg;
   vector<double> *trkpx;
   vector<double> *trkpy;
   vector<double> *trkpz;
   vector<double> *trkpt;
   vector<double> *trkdcabs;
   vector<double> *trkdcabserr;
   Int_t nb;
   vector<int> *bchg;
   vector<double> *bpx;
   vector<double> *bpxerr;
   vector<double> *bpy;
   vector<double> *bpyerr;
   vector<double> *bpz;
   vector<double> *bpzerr;
   vector<double> *bmass;
   vector<double> *bmasserr;
   vector<double> *bvtxcl;
   vector<double> *bvtxx;
   vector<double> *bvtxxerr;
   vector<double> *bvtxy;
   vector<double> *bvtxyerr;
   vector<double> *bvtxz;
   vector<double> *bvtxzerr;
   vector<double> *bcosalphabs;
   vector<double> *bcosalphabserr;
   vector<double> *bcosalphabs2D;    // 2014-05-22 N.A.
	vector<double> *bcosalphabs2Derr; // 2014-05-22 N.A.
	vector<double> *blsbs;
   vector<double> *blsbserr;
   vector<double> *bctau;
   vector<double> *bctauerr; 
   
	Int_t genbchg;
   Double_t genbpx;
   Double_t genbpy;
   Double_t genbpz;
   Int_t gentrkchg;
   Double_t gentrkpx;
   Double_t gentrkpy;
   Double_t gentrkpz;
   Double_t genmumpx;
   Double_t genmumpy;
   Double_t genmumpz;
   Double_t genmuppx;
   Double_t genmuppy;
   Double_t genmuppz;

   string *decname;
   vector<bool> *istruemum;
   vector<bool> *istruemup;
   vector<bool> *istruetrk;
   vector<bool> *istruebu;

   // List of branches
   TBranch *b_run; //!
   TBranch *b_event; //!
   TBranch *b_lumiblock; //!
   TBranch *b_nprivtx; //!
   TBranch *b_triggernames; //!
   TBranch *b_triggerprescales; //!
   TBranch *b_mumdcabs; //!
   TBranch *b_mumdcabserr; //!
   TBranch *b_mumpx; //!
   TBranch *b_mumpy; //!
   TBranch *b_mumpz; //!
   TBranch *b_mupdcabs; //!
   TBranch *b_mupdcabserr; //!
   TBranch *b_muppx; //!
   TBranch *b_muppy; //!
   TBranch *b_muppz; //!
   TBranch *b_mumutrkr; //!
   TBranch *b_mumutrkz; //!
   TBranch *b_mumudca; //!
   TBranch *b_mumuvtxcl; //!
   TBranch *b_mumulsbs; //!
   TBranch *b_mumulsbserr; //!
   TBranch *b_mumucosalphabs; //!
   TBranch *b_mumucosalphabserr; //!
   TBranch *b_mumumass; //!
   TBranch *b_mumumasserr; //!
   TBranch *b_mumisgoodmuon; //!
   TBranch *b_mupisgoodmuon; //!
   TBranch *b_mumnpixhits; //!
   TBranch *b_mupnpixhits; //!
   TBranch *b_mumnpixlayers; //!
   TBranch *b_mupnpixlayers; //!
   TBranch *b_mumntrkhits; //!
   TBranch *b_mupntrkhits; //!
   TBranch *b_mumntrklayers; //!
   TBranch *b_mupntrklayers; //!
   TBranch *b_mumnormchi2; //!
   TBranch *b_mupnormchi2; //!
   TBranch *b_mumdxyvtx; //!
   TBranch *b_mupdxyvtx; //!
   TBranch *b_mumdzvtx; //!
   TBranch *b_mupdzvtx; //!
   TBranch *b_mumtriglastfilter; //!
   TBranch *b_muptriglastfilter; //!
   TBranch *b_mumpt; //!
   TBranch *b_muppt; //!
   TBranch *b_mumeta; //!
   TBranch *b_mupeta; //!
   TBranch *b_trkchg; //!
   TBranch *b_trkpx; //!
   TBranch *b_trkpy; //!
   TBranch *b_trkpz; //!
   TBranch *b_trkpt; //!
   TBranch *b_trkdcabs; //!
   TBranch *b_trkdcabserr; //!
   TBranch *b_nb; //!
   TBranch *b_bchg; //!
   TBranch *b_bpx; //!
   TBranch *b_bpxerr; //!
   TBranch *b_bpy; //!
   TBranch *b_bpyerr; //!
   TBranch *b_bpz; //!
   TBranch *b_bpzerr; //!
   TBranch *b_bmass; //!
   TBranch *b_bmasserr; //!
   TBranch *b_bvtxcl; //!
   TBranch *b_bvtxx; //!
   TBranch *b_bvtxxerr; //!
   TBranch *b_bvtxy; //!
   TBranch *b_bvtxyerr; //!
   TBranch *b_bvtxz; //!
   TBranch *b_bvtxzerr; //!
   TBranch *b_bcosalphabs; //!
   TBranch *b_bcosalphabserr; //!
   TBranch *b_bcosalphabs2D; //!    // 2014-05-22 N.A.
	TBranch *b_bcosalphabs2Derr; //! // 2014-05-22 N.A.
	TBranch *b_blsbs; //!
   TBranch *b_blsbserr; //!
   TBranch *b_bctau; //!
   TBranch *b_bctauerr; //!

   TBranch *b_genbchg; //!
   TBranch *b_genbpx; //!
   TBranch *b_genbpy; //!
   TBranch *b_genbpz; //!
   TBranch *b_gentrkchg; //!
   TBranch *b_gentrkpx; //!
   TBranch *b_gentrkpy; //!
   TBranch *b_gentrkpz; //!
   TBranch *b_genmumpx; //!
   TBranch *b_genmumpy; //!
   TBranch *b_genmumpz; //!
   TBranch *b_genmuppx; //!
   TBranch *b_genmuppy; //!
   TBranch *b_genmuppz; //!

   TBranch *b_decname; //!
   TBranch *b_istruemum; //!
   TBranch *b_istruemup; //!
   TBranch *b_istruetrk; //!
   TBranch *b_istruebu; //!


   SingleBuToKMuMuSelector(TTree * /*tree*/ =0) { }
   virtual ~SingleBuToKMuMuSelector() { }
   virtual Int_t Version() const { return 2; }
   virtual void Begin(TTree *tree);
   virtual void SlaveBegin(TTree *tree);
   virtual void Init(TTree *tree);
   virtual Bool_t Notify();
   virtual Bool_t Process(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void SetOption(const char *option) { fOption = option; }
   virtual void SetObject(TObject *obj) { fObject = obj; }
   virtual void SetInputList(TList *input) { fInput = input; }
   virtual TList *GetOutputList() const { return fOutput; }
   virtual void SlaveTerminate();
   virtual void Terminate();

   // ClassDef(SingleBuToKMuMuSelector,0);
   int SelectB(string);
   bool HasGoodMuons(int);
   void SaveEvent(int);
	bool TriggerSelections(int);
	bool KpSelections(int);
	bool OptimizedSelections(int);
	void SaveGen();
};

#endif

#ifdef SingleBuToKMuMuSelector_cxx
void SingleBuToKMuMuSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggernames = 0;
   triggerprescales = 0;
   mumdcabs = 0;
   mumdcabserr = 0;
   mumpx = 0;
   mumpy = 0;
   mumpz = 0;
   mupdcabs = 0;
   mupdcabserr = 0;
   muppx = 0;
   muppy = 0;
   muppz = 0;
   mumutrkr = 0;
   mumutrkz = 0;
   mumudca = 0;
   mumuvtxcl = 0;
   mumulsbs = 0;
   mumulsbserr = 0;
   mumucosalphabs = 0;
   mumucosalphabserr = 0;
   mumumass = 0;
   mumumasserr = 0;
   mumisgoodmuon = 0;
   mupisgoodmuon = 0;
   mumnpixhits = 0;
   mupnpixhits = 0;
   mumnpixlayers = 0;
   mupnpixlayers = 0;
   mumntrkhits = 0;
   mupntrkhits = 0;
   mumntrklayers = 0;
   mupntrklayers = 0;
   mumnormchi2 = 0;
   mupnormchi2 = 0;
   mumdxyvtx = 0;
   mupdxyvtx = 0;
   mumdzvtx = 0;
   mupdzvtx = 0;
   mumtriglastfilter = 0;
   muptriglastfilter = 0;
   mumpt = 0;
   muppt = 0;
   mumeta = 0;
   mupeta = 0;
   trkchg = 0;
   trkpx = 0;
   trkpy = 0;
   trkpz = 0;
   trkpt = 0;
   trkdcabs = 0;
   trkdcabserr = 0;
   bchg = 0;
   bpx = 0;
   bpxerr = 0;
   bpy = 0;
   bpyerr = 0;
   bpz = 0;
   bpzerr = 0;
   bmass = 0;
   bmasserr = 0;
   bvtxcl = 0;
   bvtxx = 0;
   bvtxxerr = 0;
   bvtxy = 0;
   bvtxyerr = 0;
   bvtxz = 0;
   bvtxzerr = 0;
   bcosalphabs = 0;
   bcosalphabserr = 0;
   bcosalphabs2D = 0;     // 2014-05-22 N.A.
   bcosalphabs2Derr = 0;  // 2014-05-22 N.A.
   blsbs = 0;
   blsbserr = 0;
   bctau = 0;
   bctauerr = 0;
   decname = 0;
   istruemum = 0;
   istruemup = 0;
   istruetrk = 0;
   istruebu = 0;
  // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("nprivtx", &nprivtx, &b_nprivtx);
   fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
   fChain->SetBranchAddress("triggerprescales", &triggerprescales, &b_triggerprescales);
   fChain->SetBranchAddress("mumdcabs", &mumdcabs, &b_mumdcabs);
   fChain->SetBranchAddress("mumdcabserr", &mumdcabserr, &b_mumdcabserr);
   fChain->SetBranchAddress("mumpx", &mumpx, &b_mumpx);
   fChain->SetBranchAddress("mumpy", &mumpy, &b_mumpy);
   fChain->SetBranchAddress("mumpz", &mumpz, &b_mumpz);
   fChain->SetBranchAddress("mupdcabs", &mupdcabs, &b_mupdcabs);
   fChain->SetBranchAddress("mupdcabserr", &mupdcabserr, &b_mupdcabserr);
   fChain->SetBranchAddress("muppx", &muppx, &b_muppx);
   fChain->SetBranchAddress("muppy", &muppy, &b_muppy);
   fChain->SetBranchAddress("muppz", &muppz, &b_muppz);
   fChain->SetBranchAddress("mumutrkr", &mumutrkr, &b_mumutrkr);
   fChain->SetBranchAddress("mumutrkz", &mumutrkz, &b_mumutrkz);
   fChain->SetBranchAddress("mumudca", &mumudca, &b_mumudca);
   fChain->SetBranchAddress("mumuvtxcl", &mumuvtxcl, &b_mumuvtxcl);
   fChain->SetBranchAddress("mumulsbs", &mumulsbs, &b_mumulsbs);
   fChain->SetBranchAddress("mumulsbserr", &mumulsbserr, &b_mumulsbserr);
   fChain->SetBranchAddress("mumucosalphabs", &mumucosalphabs, &b_mumucosalphabs);
   fChain->SetBranchAddress("mumucosalphabserr", &mumucosalphabserr, &b_mumucosalphabserr);
   fChain->SetBranchAddress("mumumass", &mumumass, &b_mumumass);
   fChain->SetBranchAddress("mumumasserr", &mumumasserr, &b_mumumasserr);
   fChain->SetBranchAddress("mumisgoodmuon", &mumisgoodmuon, &b_mumisgoodmuon);
   fChain->SetBranchAddress("mupisgoodmuon", &mupisgoodmuon, &b_mupisgoodmuon);
   fChain->SetBranchAddress("mumnpixhits", &mumnpixhits, &b_mumnpixhits);
   fChain->SetBranchAddress("mupnpixhits", &mupnpixhits, &b_mupnpixhits);
   fChain->SetBranchAddress("mumnpixlayers", &mumnpixlayers, &b_mumnpixlayers);
   fChain->SetBranchAddress("mupnpixlayers", &mupnpixlayers, &b_mupnpixlayers);
   fChain->SetBranchAddress("mumntrkhits", &mumntrkhits, &b_mumntrkhits);
   fChain->SetBranchAddress("mupntrkhits", &mupntrkhits, &b_mupntrkhits);
   fChain->SetBranchAddress("mumntrklayers", &mumntrklayers, &b_mumntrklayers);
   fChain->SetBranchAddress("mupntrklayers", &mupntrklayers, &b_mupntrklayers);
   fChain->SetBranchAddress("mumnormchi2", &mumnormchi2, &b_mumnormchi2);
   fChain->SetBranchAddress("mupnormchi2", &mupnormchi2, &b_mupnormchi2);
   fChain->SetBranchAddress("mumdxyvtx", &mumdxyvtx, &b_mumdxyvtx);
   fChain->SetBranchAddress("mupdxyvtx", &mupdxyvtx, &b_mupdxyvtx);
   fChain->SetBranchAddress("mumdzvtx", &mumdzvtx, &b_mumdzvtx);
   fChain->SetBranchAddress("mupdzvtx", &mupdzvtx, &b_mupdzvtx);
   fChain->SetBranchAddress("mumtriglastfilter", &mumtriglastfilter, &b_mumtriglastfilter);
   fChain->SetBranchAddress("muptriglastfilter", &muptriglastfilter, &b_muptriglastfilter);
   fChain->SetBranchAddress("mumpt", &mumpt, &b_mumpt);
   fChain->SetBranchAddress("muppt", &muppt, &b_muppt);
   fChain->SetBranchAddress("mumeta", &mumeta, &b_mumeta);
   fChain->SetBranchAddress("mupeta", &mupeta, &b_mupeta);
   fChain->SetBranchAddress("trkchg", &trkchg, &b_trkchg);
   fChain->SetBranchAddress("trkpx", &trkpx, &b_trkpx);
   fChain->SetBranchAddress("trkpy", &trkpy, &b_trkpy);
   fChain->SetBranchAddress("trkpz", &trkpz, &b_trkpz);
   fChain->SetBranchAddress("trkpt", &trkpt, &b_trkpt);
   fChain->SetBranchAddress("trkdcabs", &trkdcabs, &b_trkdcabs);
   fChain->SetBranchAddress("trkdcabserr", &trkdcabserr, &b_trkdcabserr);
   fChain->SetBranchAddress("nb", &nb, &b_nb);
   fChain->SetBranchAddress("bchg", &bchg, &b_bchg);
   fChain->SetBranchAddress("bpx", &bpx, &b_bpx);
   fChain->SetBranchAddress("bpxerr", &bpxerr, &b_bpxerr);
   fChain->SetBranchAddress("bpy", &bpy, &b_bpy);
   fChain->SetBranchAddress("bpyerr", &bpyerr, &b_bpyerr);
   fChain->SetBranchAddress("bpz", &bpz, &b_bpz);
   fChain->SetBranchAddress("bpzerr", &bpzerr, &b_bpzerr);
   fChain->SetBranchAddress("bmass", &bmass, &b_bmass);
   fChain->SetBranchAddress("bmasserr", &bmasserr, &b_bmasserr);
   fChain->SetBranchAddress("bvtxcl", &bvtxcl, &b_bvtxcl);
   fChain->SetBranchAddress("bvtxx", &bvtxx, &b_bvtxx);
   fChain->SetBranchAddress("bvtxxerr", &bvtxxerr, &b_bvtxxerr);
   fChain->SetBranchAddress("bvtxy", &bvtxy, &b_bvtxy);
   fChain->SetBranchAddress("bvtxyerr", &bvtxyerr, &b_bvtxyerr);
   fChain->SetBranchAddress("bvtxz", &bvtxz, &b_bvtxz);
   fChain->SetBranchAddress("bvtxzerr", &bvtxzerr, &b_bvtxzerr);
   fChain->SetBranchAddress("bcosalphabs", &bcosalphabs, &b_bcosalphabs);
   fChain->SetBranchAddress("bcosalphabserr", &bcosalphabserr, &b_bcosalphabserr);
   fChain->SetBranchAddress("bcosalphabs2D", &bcosalphabs2D, &b_bcosalphabs2D);           // 2014-05-22 N.A.  
   fChain->SetBranchAddress("bcosalphabs2Derr", &bcosalphabs2Derr, &b_bcosalphabs2Derr);  // 2014-05-22 N.A.
   fChain->SetBranchAddress("blsbs", &blsbs, &b_blsbs);
   fChain->SetBranchAddress("blsbserr", &blsbserr, &b_blsbserr);
   fChain->SetBranchAddress("bctau", &bctau, &b_bctau);
   fChain->SetBranchAddress("bctauerr", &bctauerr, &b_bctauerr);

   fChain->SetBranchAddress("genbchg", &genbchg, &b_genbchg);
   fChain->SetBranchAddress("genbpx", &genbpx, &b_genbpx);
   fChain->SetBranchAddress("genbpy", &genbpy, &b_genbpy);
   fChain->SetBranchAddress("genbpz", &genbpz, &b_genbpz);
   fChain->SetBranchAddress("gentrkchg", &gentrkchg, &b_gentrkchg);
   fChain->SetBranchAddress("gentrkpx", &gentrkpx, &b_gentrkpx);
   fChain->SetBranchAddress("gentrkpy", &gentrkpy, &b_gentrkpy);
   fChain->SetBranchAddress("gentrkpz", &gentrkpz, &b_gentrkpz);
   fChain->SetBranchAddress("genmumpx", &genmumpx, &b_genmumpx);
   fChain->SetBranchAddress("genmumpy", &genmumpy, &b_genmumpy);
   fChain->SetBranchAddress("genmumpz", &genmumpz, &b_genmumpz);
   fChain->SetBranchAddress("genmuppx", &genmuppx, &b_genmuppx);
   fChain->SetBranchAddress("genmuppy", &genmuppy, &b_genmuppy);
   fChain->SetBranchAddress("genmuppz", &genmuppz, &b_genmuppz);

	fChain->SetBranchAddress("decname", &decname, &b_decname);
   fChain->SetBranchAddress("istruemum", &istruemum, &b_istruemum);
   fChain->SetBranchAddress("istruemup", &istruemup, &b_istruemup);
   fChain->SetBranchAddress("istruetrk", &istruetrk, &b_istruetrk);
   fChain->SetBranchAddress("istruebu", &istruebu, &b_istruebu);

}

Bool_t SingleBuToKMuMuSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef SingleBuToKMuMuSelector_cxx
