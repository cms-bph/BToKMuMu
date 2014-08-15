// -----------------------------------------------
// Author: Xin Shi <Xin.Shi@cern.ch>
// Created: [2013-04-04 Thu 08:42]
// -----------------------------------------------
#include <TSystem.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TROOT.h>

#include "tools.cc"

using namespace std;

void bmass(TString datatype, TString label, TString cut, TString outfile){
  // plot the BuToKJPsi part
  TH1F* h_bmass = new TH1F("h_bmass", "B^{+/-} mass; M(B^{+/-}) [GeV/c^{2}]",
100, 5, 5.6);
  TChain* ch = add_chain(datatype, label, cut);
  if (ch == NULL) gSystem->Exit(0);

  double Bmass = 0;
  ch->SetBranchAddress("Bmass", &Bmass);

  Int_t nentries = (Int_t)ch->GetEntries();
  if (nentries == 0) {
    cerr << "No entries found!" << endl;
    gSystem->Exit(0);
  }
  
  for (Int_t i=0;i<nentries;i++) {
    ch->GetEntry(i);
    h_bmass->Fill(Bmass);
  }

  TCanvas* c = new TCanvas("c","c", 400, 400);
  set_root_style();
  c->UseCurrentStyle() ;

  h_bmass->SetMinimum(0);
  h_bmass->Draw();
  c->Print(outfile.Data());
  h_bmass->Delete();
  delete c;
}


int main(int argc, char** argv) {
  TString func = argv[1];
  TString datatype = argv[2];
  TString label = argv[3];
  TString cut = argv[4];
  TString outfile = argv[5];
  
  if (func == "bmass")
    bmass(datatype, label, cut, outfile);
  else
    cerr << "No function available for: " << func.Data() << endl;

  gSystem->Exit(0);
  return 0 ;
}



