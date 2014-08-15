// -----------------------------------------------
// Author: Xin Shi <Xin.Shi@cern.ch>
// Created: [2013-08-15 Thu 14:54]
// -----------------------------------------------
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <math.h>

#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TPad.h> 
#include <TCanvas.h>
#include <TChain.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>

#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
#include <RooChebychev.h>
#include <RooGenericPdf.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAddition.h>

#include "tools.cc"

using namespace std;
using namespace RooFit ;

void bmass(TString datatype, TString label, TString cut, TString outfile)
{
//	bool test = true;
	bool test = false;
	
//	Importing a TTree into a RooDataSet with cuts
//	--------------------------------------------------------------------------
	TChain* ch = add_chain(datatype, label, cut);
	if (ch == NULL) gSystem->Exit(0);
	
	RooRealVar     Bmass("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.27926-0.28, 5.27926+0.28);
	RooDataSet     data("data", "data", RooArgSet(Bmass), Import(*ch));
	data.Print();
	
//	Create model and dataset
//	-----------------------------------------------
//	Gaussian signal
/*
//	Single Gauss: // 15-08
	RooRealVar     mean("mean","mean of gaussians", 5.27925, 5.27908, 5.27942);
	RooRealVar     sigma("sigma","width of gaussians", 0.065, 0, 0.1);
	RooGaussian    sig("sig","Signal component", x, mean, sigma);
*/
//	Double Gauss: // 15-08
	RooRealVar     mean("mean","mean of gaussians", 5.27, 5.23, 5.32);
	RooRealVar     sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05);
	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.065, 0.05, 0.1);
	RooRealVar     sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
	RooGaussian    sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1);
	RooGaussian    sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2);
	RooAddPdf      sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));
	
//	Build Chebychev polynomial p.d.f.
	RooRealVar     a0("a0", "constant", 0.5, -1, 1.);
	RooRealVar     a1("a1", "linear", 0.6, -1, 1);
	RooRealVar     a2("a2", "quadratic", 0.1, -1, 1);
	RooChebychev   bkg("bkg", "Background", Bmass, RooArgSet(a0, a1, a2));
//	Construct signal+background PDF
	RooRealVar     nsig("nsig", "number of signal events", 4648, 0, 1E8);
	RooRealVar     nbkg("nbkg", "number of background events", 21472, 0, 1E8);
	RooAddPdf      model("model", "g+c", RooArgList(bkg, sig), RooArgList(nbkg, nsig));
	
//	Print structure of composite p.d.f.
	model.Print("t");
	
//	Fit model to data, save fitresult
//	-----------------------------------------------------------------
	RooFitResult  *fitres;
	if ( !test) {
		fitres = model.fitTo(data, Extended(true), Save(true));
		fitres->Print("v");
	}
	
//	Plot model
//	---------------------------------------------------------
	TString title = "B^{+/-} mass";
	int nbins = 50;
	RooPlot *xframe = Bmass.frame(Title(title), Bins(nbins));
	data.plotOn(xframe);
	model.plotOn(xframe);

//	Overlay the background component of model with a dashed line
	model.plotOn(xframe,Components("bkg"), LineStyle(kDashed));
	
//	Draw the frame on the canvas
	TCanvas* c = new TCanvas("c", "c", 800, 600);
	set_root_style();
	c->UseCurrentStyle();
	
	gPad->SetLeftMargin(0.15);
	xframe->GetYaxis()->SetTitleOffset(1.7);
	xframe->Draw();

//	15-08 N.A.
	TPaveText* paveText = new TPaveText(0.17, 0.70, 0.41, 0.88, "NDC");
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("nsig = %.0f #pm %.0f " , nsig.getVal()  , nsig.getError()));
	paveText->AddText(Form("nbkg = %.0f #pm %.0f " , nbkg.getVal()  , nbkg.getError()));
	paveText->AddText(Form("mean = %.3f #pm %.3f " , mean.getVal()  , mean.getError()));
	paveText->AddText(Form("sigma1 = %.3f #pm %.3f ", sigma1.getVal(), sigma1.getError()));
	paveText->AddText(Form("sigma2 = %.3f #pm %.3f ", sigma2.getVal(), sigma2.getError()));
	paveText->AddText(Form("frac = %.3f #pm %.3f ", sigM_frac.getVal(), sigM_frac.getError()));
	paveText->Draw();

	TString pdffile = "./fitresults/" + outfile + ".pdf";
	TString pngfile = "./fitresults/" + outfile + ".png";
	c->Print(pdffile);
	c->Print(pngfile);
	
//	Persist fit result in root file
//	-------------------------------------------------------------
	TString resfile = "./fitresults/" + outfile + ".root";
	TFile resf(resfile, "RECREATE");
	gPad->Write("plot");
	if (!test) fitres->Write("fitres");
	resf.Close() ;
	
//	In a clean ROOT session retrieve the persisted fit result as follows:
//	RooFitResult* r = gDirectory->Get("fitres");
	
	delete paveText;
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

