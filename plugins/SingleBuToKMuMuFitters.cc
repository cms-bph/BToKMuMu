//vim: sw=4 ts=4 fdm=marker et:

// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
//#include <regex>
#include "Math/MinimizerOptions.h"
#include "TROOT.h"

#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TPad.h> 
#include <TLegend.h> 
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
#include <RooMinuit.h>

#include "tools.cc" 

using namespace std; 
using namespace RooFit;

// Tags configration
bool is7TeVCheck = false; // Using 2011 efficiency map.
TChain *ch=new TChain("tree");

// Lumis.
const float Lumi_Sig = 3295.08494, Lumi_Jpsi = 18.53167, Lumi_Psi = 209.13758, Lumi_Data = 20.47;   // fb-1
const float Lumi_Scale = 2.41281;

//Constants, Fit results for efficiency, etc.. //{{{
char genQ2range[11][200] = {"genQ2 <=  2.00 && genQ2 >  1.00",
                           "genQ2 <=  4.30 && genQ2 >  2.00",
                           "genQ2 <=  8.68 && genQ2 >  4.30",
                           "genQ2 <= 10.09 && genQ2 >  8.68",
									"genQ2 <= 12.86 && genQ2 > 10.09",
                           "genQ2 <= 14.18 && genQ2 > 12.86",
									"genQ2 <= 16.00 && genQ2 > 14.18",
                           "genQ2 <= 18.00 && genQ2 > 16.00",
                           "genQ2 <= 22.00 && genQ2 > 18.00",
                           "genQ2 <=  6.00 && genQ2 >  1.00",
									"(genQ2 > 1.  && genQ2 < 8.68) || (genQ2 >= 14.18  && genQ2 < 22.) || (genQ2 >= 10.09  && genQ2 < 12.86)"};
								//	"genQ2 <= 22.00 && genQ2 >  1.00"};
char Q2range[11][200] = {"Q2 <=  2.00 && Q2 >  1.00",
                        "Q2 <=  4.30 && Q2 >  2.00",
                        "Q2 <=  8.68 && Q2 >  4.30",
                        "Q2 <= 10.09 && Q2 > 8.68",
								"Q2 <= 12.86 && Q2 > 10.09",
                        "Q2 <= 14.18 && Q2 > 12.86",
								"Q2 <= 16.00 && Q2 > 14.18",
                        "Q2 <= 18.00 && Q2 > 16.00",
                        "Q2 <= 22.00 && Q2 > 18.00",
                        "Q2 <=  6.00 && Q2 >  1.00",
								"(Q2 > 1.  && Q2 < 8.68) || (Q2 >= 14.18 && Q2 < 22.) || (Q2 >= 10.09  && Q2 < 12.86)"};
							//	"Q2 <= 22.00 && Q2 >  1.00"};
double Q2rangedn[11] = {1.00 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 1.00 ,  1.00};
double Q2rangeup[11] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 22.00 , 6.00 , 22.00};
char mumuMassWindow[5][300] = { 
	" Mumumass > 0 ",
	" (Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)", // No JPsi && Psi(2S)
//	" (Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.5*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)",
//	&& (fabs(Bmass - Mumumass - 2.182) > 0.14) && (fabs(Bmass - Mumumass - 1.593) > 0.09)",
	" (Mumumass < 3.096916+3.*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr) \
	  || (Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr)",  // Only JPsi && Psi(2S)
	" Mumumass < 3.096916+3*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr",  // Only JPsi
	" Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr"}; //Only Psi(2S)
double genAfb[9]={0.00, 0.07, -0.02, -0.03, -0.01, -0.09, 0.02, 0.02, -0.01};   // LHCb results.....
double genFh [9]={0.00, 0.14,  0.04,  0.11,  0.08,  0.14, 0.14, 0.05,  0.02};
std::string f_accXrecoEff_ord0[8] = {
    "76.637435*((exp(-0.5*((CosThetaL-(-0.001125))/0.421561)**2)*(-0.027922*CosThetaL**2-0.000022*CosThetaL+0.026002))*(0.010083-0.000846*CosThetaK+0.014753*CosThetaK**2-0.001424*CosThetaK**3-0.013006*CosThetaK**4-0.001910*CosThetaK**5-0.003559*CosThetaK**6))",
    "80.036712*((0.020509+0.000112*CosThetaL-0.046503*CosThetaL**2+0.001045*CosThetaL**3+0.037179*CosThetaL**4-0.001041*CosThetaL**5-0.011163*CosThetaL**6)*(0.009374-0.000138*CosThetaK+0.023015*CosThetaK**2-0.004101*CosThetaK**3-0.036478*CosThetaK**4+0.000861*CosThetaK**5+0.016003*CosThetaK**6))",
    "110.528399*((0.012585+0.000951*CosThetaL-0.015603*CosThetaL**2-0.001576*CosThetaL**3+0.005612*CosThetaL**4+0.001148*CosThetaL**5-0.001163*CosThetaL**6)*(0.007571-0.000567*CosThetaK+0.011154*CosThetaK**2+0.000546*CosThetaK**3-0.018860*CosThetaK**4-0.002906*CosThetaK**5+0.009528*CosThetaK**6))",
    "461.476571*((0.002469+0.000225*CosThetaL-0.000676*CosThetaL**2-0.000874*CosThetaL**3-0.000779*CosThetaL**4+0.001016*CosThetaL**5-0.000323*CosThetaL**6)*(0.001929-0.000503*CosThetaK+0.001211*CosThetaK**2+0.000547*CosThetaK**3-0.001014*CosThetaK**4-0.000494*CosThetaK**5+0.000002*CosThetaK**6))",
    "75.454746*((0.014286+0.001027*CosThetaL-0.003358*CosThetaL**2-0.002559*CosThetaL**3+0.001474*CosThetaL**4+0.002134*CosThetaL**5-0.003858*CosThetaL**6)*(0.012601-0.000826*CosThetaK+0.003595*CosThetaK**2-0.001026*CosThetaK**3-0.003101*CosThetaK**4-0.000333*CosThetaK**5+0.000221*CosThetaK**6))",
    "551.279776*((0.001797+0.000110*CosThetaL+0.000791*CosThetaL**2+0.000274*CosThetaL**3-0.003945*CosThetaL**4-0.000380*CosThetaL**5+0.003388*CosThetaL**6)*(0.001864+0.000467*CosThetaK-0.001804*CosThetaK**2-0.002188*CosThetaK**3+0.005723*CosThetaK**4+0.001518*CosThetaK**5-0.004205*CosThetaK**6))",
    "65.806692*((0.014983+0.000288*CosThetaL+0.000503*CosThetaL**2+0.001627*CosThetaL**3-0.003861*CosThetaL**4-0.002529*CosThetaL**5+0.004376*CosThetaL**6)*(0.014860-0.002541*CosThetaK+0.001603*CosThetaK**2+0.004161*CosThetaK**3-0.001351*CosThetaK**4-0.003963*CosThetaK**5+0.000169*CosThetaK**6))",
    "48.998722*((0.020047-0.001508*CosThetaL+0.001947*CosThetaL**2+0.005356*CosThetaL**3-0.004953*CosThetaL**4-0.004908*CosThetaL**5+0.005197*CosThetaL**6)*(0.020431-0.001608*CosThetaK-0.003238*CosThetaK**2-0.001358*CosThetaK**3+0.012846*CosThetaK**4+0.001144*CosThetaK**5-0.011104*CosThetaK**6))"
    };
//}}}
//

double readParam(int iBin, const char parName[], int iColumn)
{//{{{
    std::vector<double> output;
    char lineBuff[512];
    char *valBuff;
    memset(lineBuff,' ',512*sizeof(char));
    FILE *fp = fopen(TString::Format("./fitParameters/fitParameters%d.txt",iBin),"r");
    while(fgets(lineBuff,512,fp) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            printf("INFO: readParam, matched %15s!\n",valBuff);
            valBuff = strtok(NULL," ");
            while(valBuff != NULL){
                //output.push_back(stof(valBuff));//stof if c++11 function, use other function
                output.push_back(std::atof(valBuff));
                valBuff = strtok(NULL," ");
            }
            break;
        }
        memset(lineBuff,' ',512*sizeof(char));
    }
    fclose(fp);
    
    if (iColumn < output.size() ){
		 cout<<"INFO: readParam, get "<<parName<<"["<<iColumn<<"] = "<<output.at(iColumn)<<endl;
      //  printf("INFO: readParam, get %s[%d]= %f \n",parName,iColumn, output.at(iColumn));
        return output.at(iColumn);
    }else{
        printf("ERROR: readParam, empty column! Return 0.\n");
        return 0.;
    }
}//}}}
void writeParam(int iBin, const char parName[], double *val, int nVal=2, bool overwrite=true)
{//{{{
    struct stat fiBuff;
    FILE *fi = 0;
    if (stat(TString::Format("./fitParameters/fitParameters%d.txt",iBin),&fiBuff) == 0){
        rename(TString::Format("./fitParameters/fitParameters%d.txt",iBin),TString::Format("./fitParameters/fitParameters%d.txt.temp",iBin));
        fi = fopen(TString::Format("./fitParameters/fitParameters%d.txt.temp",iBin),"r");
    }else{
        fi = fopen(TString::Format("./fitParameters/fitParameters%d.txt.temp",iBin),"w");
    }
    
    bool parExist = false;
    char lineBuff[512];
    char *valBuff = 0;
    memset(lineBuff,' ',512*sizeof(char));
    FILE *fp = fopen(TString::Format("./fitParameters/fitParameters%d.txt",iBin),"w");
    while(fgets(lineBuff,512,fi) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            fprintf(fp,"%15s",parName);
            int iVal = 0;
            while(iVal < nVal){
                fprintf(fp," %20.15f",val[iVal]);
                iVal++;
            }
            fprintf(fp,"\n");
            parExist = true;
        }else{
            fprintf(fp,"%15s",lineBuff);
            valBuff = strtok(NULL," ");
            while( valBuff != NULL ){
                fprintf(fp," %15s",valBuff);
                valBuff = strtok(NULL," ");
            }
        }
        memset(lineBuff,' ',512*sizeof(char));
    }
    if (parExist == false){
        fprintf(fp,"%15s",parName);
        int iVal = 0;
        while(iVal < nVal){
            fprintf(fp," %20.15f",val[iVal]);
            iVal++;
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    fclose(fi);
    remove(TString::Format("./fitParameters/fitParameters%d.txt.temp",iBin));
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
/*
TF2  *f2_fcn = NULL;
double model_2D(double *x, double *par)
{//{{{
    double xx = x[0];
    double yy = x[1];
    for (int i = 0; i < f2_fcn->GetNpar(); i++) f2_fcn->SetParameter(i,par[i]);
    return f2_fcn->Eval(xx,yy);
}//}}}

TH1F *h1_fcn = NULL;
// nParameters, ???, return fcn value, parameter array, strategy
void fcn_binnedChi2_2D(int &npar, double *gin, double &f, double *par, int ifhag)
{//{{{
    f=0;
    for (int i = 1; i <= h1_fcn->GetNbinsX(); i++) {
        for (int j = 1; j <= h1_fcn->GetNbinsY(); j++) {
            int gBin = h1_fcn->GetBin(i);
            double x[2] = {h1_fcn->GetXaxis()->GetBinCenter(i),h1_fcn->GetYaxis()->GetBinCenter(j)};
            double measure  = h1_fcn->GetBinContent(gBin);
            double error    = h1_fcn->GetBinError(gBin);
            
            //// Naively test using center value
            //double func     = model_2D(x, par);//Take center value
            //double delta    = (measure-func)/error;
            //if (measure != 0) 
            //    f+=delta*delta;
            
            //// Real run using integral
            for (int i = 0; i < f2_fcn->GetNpar(); i++){//nPar MUST be the same value as f2_fcn
                f2_fcn->SetParameter(i,par[i]);
            }
            double xi = h1_fcn->GetXaxis()->GetBinLowEdge(i);
            double xf = h1_fcn->GetXaxis()->GetBinUpEdge(i);
            double yi = h1_fcn->GetYaxis()->GetBinLowEdge(j);
            double yf = h1_fcn->GetYaxis()->GetBinUpEdge(j);
            f += pow( (f2_fcn->Integral(xi,xf,yi,yf)/(xf-xi)/(yf-yi)-measure)/error,2);

            //f2_square = 0;
            //delete f2_square;
        }
    }
    //printf("FCN in calls = %f\n",f);
    //printf("npar=%d ",npar);
}//}}}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////
//_________________________________________________________________________________

void bmass(int iBin, const char outfile[] = "bmass")
{//{{{
	bool test = false; 
	
	RooRealVar     Bmass("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.27925-0.28, 5.27925+0.28);
	RooRealVar     Q2("Q2","q^{2}",1.0,22.);
	RooRealVar     Mumumass("Mumumass","M^{#mu#mu}",1.,10.);
	RooRealVar     Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
//	RooDataSet     *data = new RooDataSet("data", "data", RooArgSet(Bmass), Import(*ch));
	RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr),
//	TString::Format("(%s) && (%s) ",Q2range[iBin],mumuMassWindow[1]),0);	// data_v2
	TString::Format("(%s) && (%s) ",Q2range[iBin],mumuMassWindow[0]),0);	// data_v1

	data->Print();
	
//	Create model and dataset
//	-------------------------------------------------------------------------
//	Gaussian signal 
	RooRealVar     mean("mean","mean of gaussians", 5.27925, 5.23, 5.32);
	RooRealVar     sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05);
	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.065, 0.05, 0.35);
	RooRealVar     sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
	RooGaussian    sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1);
	RooGaussian    sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2);
	RooAddPdf      sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));
	
//	Build Chebychev polynomial p.d.f.  
	RooRealVar     a0("a0", "constant", 0.5, -1, 1);
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
//	------------------------------------------------------------------------
	RooFitResult* fitres; 
	if (! test) {
		fitres = model.fitTo(*data, Extended(kTRUE), Save(kTRUE));
		fitres->Print("v"); 
	}
	
//	Plot model 
//	---------------------------------------------------------
	TString title = "B^{+/-} mass";
	int nbins = 20; 
	RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
	data->plotOn(frame);
	model.plotOn(frame, LineColor(1));
	
//	Overlay the background component of model with a dashed line
	model.plotOn(frame,Components("bkg"), LineStyle(kDashed), LineColor(2));
//	Overlay the signal component of model with a blue line
	model.plotOn(frame,Components("sig"), LineStyle(1), LineColor(4));
	
//	Draw the frame on the canvas
	TCanvas *c = new TCanvas("c", "c", 800, 600); 
	set_root_style(); 
	c->UseCurrentStyle();
	
	gPad->SetLeftMargin(0.15);
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->Draw();
	
//	double chi2Val=0;
//	chi2Val = model.GetChisquare();
//	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
	TPaveText* paveText = new TPaveText( 0.17, 0.70, 0.41, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
	paveText->AddText(Form("nbkg   = %.0f #pm %.0f ", nbkg.getVal(), nbkg.getError())); 
	paveText->AddText(Form("mean   = %.5f #pm %.5f ", mean.getVal(), mean.getError())); 
	paveText->AddText(Form("sigma1 = %.5f #pm %.5f ", sigma1.getVal(), sigma1.getError())); 
	paveText->AddText(Form("sigma2 = %.5f #pm %.5f ", sigma2.getVal(), sigma2.getError())); 
	paveText->AddText(Form("frac   = %.5f #pm %.5f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->Draw(); 
	
//	c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_bin%d.png",outfile,iBin));

  // Persist fit result in root file 
  // -------------------------------------------------------------
  //TFile resf(TString::Format("./plots/%s.root",outfile), "RECREATE") ;
  //gPad->Write("plot"); 
  //if (! test) fitres->Write("fitres") ;
  //resf.Close() ;

  // In a clean ROOT session retrieve the persisted fit result as follows:
  // RooFitResult* r = gDirectory->Get("fitres") ;
   
	delete paveText; 
	delete c;

}//}}}

std::vector<double> angular_gen_bin(int iBin, const char outfile[] = "angular_gen")
{//{{{
	RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{L}", -1., 1.);
	RooRealVar genQ2("genQ2","q^{2}",1.0,22.);
	RooRealVar fh("fh", "F_{H}", genFh[iBin], 0., 1.);
	RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -1., 1.);
	
	RooRealVar nsig("nsig","nsig",1E6,1E2,1E9);
//	RooRealVar nbkg("nbkg","nbkg",10,0.1,1E4);
	
	RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL", RooArgSet(genCosThetaL,fh,afb));
	RooExtendPdf f("f","",f_sig,nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaL,genQ2),genQ2range[iBin],0);
	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(1));
	
//	Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = -0.;
	
	RooPlot* framecosl = genCosThetaL.frame(); 
	data->plotOn(framecosl,Binning(100)); 
	f.plotOn(framecosl); 
	
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->Draw();
	framecosl->GetYaxis()->SetLabelFont(22);
	framecosl->GetYaxis()->SetLabelSize(0.04);
	framecosl->GetYaxis()->SetTitleSize(0.04);
	framecosl->GetYaxis()->SetTitleOffset(1.2);
	framecosl->GetYaxis()->SetTitleFont(22);
	framecosl->GetXaxis()->SetLabelFont(22);
	framecosl->GetXaxis()->SetLabelSize(0.04);
	framecosl->GetXaxis()->SetTitleSize(0.04);
	framecosl->GetXaxis()->SetTitleOffset(1.15);
	framecosl->GetXaxis()->SetTitleFont(22);
	
	fixNDC = -0.5;
//	if (iBin > 4) fixNDC = 0.;
	if (iBin == 10) { 
		t1->DrawLatex(.30,.75+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.30,.75+fixNDC,TString::Format("%s",genQ2range[iBin]));
	t1->DrawLatex(.15,.69+fixNDC,TString::Format("F_{H}  =%8.5f#pm%8.5f",fh.getVal(),fh.getError()));
	t1->DrawLatex(.51,.69+fixNDC,TString::Format("A_{FB} =%8.5f#pm%8.5f",afb.getVal(),afb.getError()));
	c->Update();
//	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
	
//	clear
	delete t1;
	delete c;
	delete data;
	
	// write output
	double val[3]={0,0,0};
	val[0] = fh.getVal();val[1] = fh.getError();
	writeParam(iBin, "genfh", val);
	val[0] = afb.getVal();val[1] = afb.getError();
	writeParam(iBin, "genafb",val);
	
	printf("genAfb[%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"genafb",0), fabs(readParam(iBin,"genafb",1)));
	printf("genFh [%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"genfh",0),  fabs(readParam(iBin,"genfh",1)));

//	write output
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;

}//}}}

void angular_gen( const char outfile[] = "angular_gen")
{//{{{

	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yerrfh[9], yafb[9], yerrafb[9]; 
			
//	Check input data
	for(int i = 0, ibin = 0; i < 9, ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		yafb[i]     = readParam(ibin,"genafb",0);
		yerrafb[i]  = fabs(readParam(ibin,"genafb",1));
		yfh[i]      = readParam(ibin,"genfh",0);
		yerrfh[i]   = fabs(readParam(ibin,"genfh",1));
		printf("genAfb[%d] =%12.6f +- %12.6f     ",i,yafb[i],yerrafb[i]);
		printf("genFh[%d]  =%12.6f +- %12.6f\n",i,yfh[i],yerrfh[i]);
	}
//	plotting
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->GetYaxis()->SetLabelFont(22);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleFont(22);
	frame->GetXaxis()->SetLabelFont(22);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleOffset(1.15);
	frame->GetXaxis()->SetTitleFont(22);
//	Fh
	frame->SetTitle("");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.02,0.08,"Y");
	frame->Draw();
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yerrfh,yerrfh);
	g_fh->SetMarkerColor(4);
	g_fh->SetMarkerStyle(20);
	g_fh->SetFillColor(2);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
//	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Print(TString::Format("./plots/%s_fh.png",outfile));
	c->Clear();
//	Afb	
	frame->SetTitle("");
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetAxisRange(-0.02,0.02,"Y");
	frame->Draw();
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yerrafb,yerrafb);
	g_afb->SetMarkerColor(4);
	g_afb->SetMarkerStyle(20);
	g_afb->SetFillColor(2);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
//	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
	c->Print(TString::Format("./plots/%s_afb.png",outfile));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
/*
//  04-09-2014
//_________________________________________________________________________________
//Fit parameters of acceptance and efficiency using RooFit.

std::vector<double> acceptance(int iBin) // acceptance, just for check...
{//{{{
	TH1::SetDefaultSumw2();
	double accUpperBound = 0.09;
	double gQ2 = 0;
	double gCosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;

	
	ch->SetBranchStatus("*",0);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("genCosThetaL"  , 1);
	ch->SetBranchStatus("genMu*"        , 1);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Fill histograms
	int nbinsL = 20;
	TH1F h1_ngen("h1_ngen","h1_ngen",nbinsL,-1.,1);  // #events in gen
	TH1F h1_nacc("h1_nacc","h1_nacc",nbinsL,-1.,1);  // #events in acceptance
//	Read data
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
		h1_ngen.Fill(gCosThetaL);
		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ) h1_nacc.Fill(gCosThetaL);
	}  
	
//	Calculate acceptance
	TH1F h1_acc("h1_acc","",nbinsL,-1,1);
	h1_acc.SetAxisRange(0.,1.,"y");
	for (int i = 1; i <= nbinsL; i++) {
	//	Generate toy model
	//	double ii = -1. + i*2./nbinsL - 1./nbinsL;
	//	double jj = -1. + j*2./nbinsK - 1./nbinsK;
	//	TF2 f2_gen("f2_gen","(0.06+0.02*(3*y**2-1)/2)-(0.03+0.03*(3*y**2-1)/2)*x**2+(0.005+0.003*(3*y**2-1)/2)*x**4",-1.,1.,-1.,1.);
	//	h1_ngen.SetBinContent(i,400);
	//	h1_nacc.SetBinContent(i,400*f2_gen.Eval(iij));
		
	//	Fill acceptance
		if (h1_ngen.GetBinContent(i) == 0) {
			printf("WARNING: Acceptance(%d)=%f/%f\n",i,h1_nacc.GetBinContent(i),h1_ngen.GetBinContent(i));
			h1_acc.SetBinContent(i,0.);
			h1_acc.SetBinError(i,1.);
		}else{
			h1_acc.SetBinContent(i,h1_nacc.GetBinContent(i)/h1_ngen.GetBinContent(i));
			if (h1_nacc.GetBinContent(i) != 0){
				h1_acc.SetBinError(i,sqrt(h1_acc.GetBinContent(i)*(1.-h1_acc.GetBinContent(i))/h1_ngen.GetBinContent(i)));
			}else{
				h1_acc.SetBinError(i,sqrt(0.05/h1_ngen.GetBinContent(i)));
			}
		//	printf("INFO: Angular bin(%d,%d)= %f +- %f ( %f / %f).\n",i,h1_acc.GetBinContent(i),h1_acc.GetBinError(i),h1_nacc.GetBinContent(i),h1_ngen.GetBinContent(i));
		}	
	}
	printf("INFO: h1_acc built.\n");
	
//	TF1 f1_model("f1_model","([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6",-1.,1.,-1.,1.);
//	TF1 f1_model("f1_model","[0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2",-1.,1.);
	TF1 f1_model("f1_model","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4",-1.,1.);
	f1_model.SetParameter(0,0.01);
	f1_model.SetParameter(1,0.01);
	f1_model.SetParameter(2,0.01);
	f1_model.SetParameter(3,0.01);
	f1_model.SetParameter(4,0.01);
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
	
//	Draw efficiency
	h1_acc.SetStats(0);
	h1_acc.SetMinimum(0.);
//	h1_acc.SetMaximum(accUpperBound);
	h1_acc.SetTitleOffset(1.3,"XY");
	h1_acc.SetMarkerStyle(20);
	h1_acc.SetMarkerSize(1.2);
	h1_acc.SetXTitle("genCosThetaL");
	h1_acc.SetYTitle("Acceptance");
	h1_acc.Draw("P E1");
	latex->DrawLatexNDC(0.35,0.95,TString::Format("Acceptance in Bin%d",iBin));
	h1_acc.Fit("f1_model"); //// 09-09
//	Draw FitResult
	f1_model.SetTitle("");
//	f1_model.SetMaximum(accUpperBound);
	f1_model.SetLineWidth(1);
//////////////////////////////////////////////////////////////////////////
// 09-09
	int nPar = 5;
	double arrPar[5], arrParErr[5], chi2Val;
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model.GetParameter(iPar);
		arrParErr[iPar] = f1_model.GetParError(iPar);
		chi2Val         = f1_model.GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
/////////////////////////////////////////////////////////////////////////////////////	
//	latex->DrawLatexNDC(0.01,0.90,TString::Format("DoF = %d",nbinsK*nbinsL-gMinuit->GetNumFreePars()));
	f1_model.Draw("SAME ");
	canvas.Print(TString::Format("./plots/acceptance_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/acceptance_bin%d.png",iBin));
	
//	Draw compare
	TH1F h1_compFit("h1_compFit","",nbinsL,-1,1);
	for (int i = 1; i <= nbinsL; i++) {//thetaL
		if (h1_acc.GetBinContent(i) != 0){
			h1_compFit.SetBinContent(i,f1_model.Eval(h1_acc.GetXaxis()->GetBinCenter(i))/h1_acc.GetBinContent(i));
		}else{
			h1_compFit.SetBinContent(i,0.);
		}		
	}
	h1_compFit.SetMinimum(0.);
	h1_compFit.SetStats(0);
	h1_compFit.SetTitleOffset(1.3,"XY");
	h1_compFit.SetXTitle("genCosThetaL");
	h1_compFit.SetMarkerStyle(20);
	h1_compFit.SetMarkerSize(1.2);
	h1_compFit.Draw("PE1");
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.30,0.95,TString::Format("acceptance_{measured} / acceptance_{fit} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/acceptance_compFit_bin%d.pdf",iBin));

//	Draw significance
	TH1F h_pull("Deviation/Error","",15,-3.,3.);
	h_pull.SetXTitle("Significance of deviation");
	h_pull.SetYTitle("Angular bins");
	for (int i = 1; i <= nbinsL; i++) {//thetaL
		double _xlo = h1_acc.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h1_acc.GetXaxis()->GetBinUpEdge(i);
		if (h1_nacc.GetBinContent(i) != 0){
			h_pull.Fill((f1_model.Integral(_xlo,_xhi)/(_xhi-_xlo)-h1_acc.GetBinContent(i))/h1_acc.GetBinError(i));
		}
	}
	h_pull.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/acceptance_sigma_bin%d.pdf",iBin));

//	Clear
	delete latex;
	
//	prepare output
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
	return output;
}//}}}

////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> recoEff(int iBin) // reconstruction efficiency for check!
{//{{{
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	double effUpperBound = 0.02;
	double BMass         = 0;
	double Mumumass      = 0;
	double Mumumasserr   = 0;
	double gQ2           = 0;
	double gCosThetaL    = 0;
	double gmuppt        = 0;
	double gmupeta       = 0;
	double gmumpt        = 0;
	double gmumeta       = 0;
	
	ch->SetBranchStatus("*",0);
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("genCosThetaL"  , 1);
	ch->SetBranchStatus("genMu*"        , 1);
	ch->SetBranchAddress("Bmass"        , &BMass);
	ch->SetBranchAddress("Mumumass"     , &Mumumass);
	ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Fill histograms
	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
	TH1F h1_nacc("h1_nacc" ,"h1_nacc" ,6,thetaLBins); 
	TH1F h1_nreco("h1_nreco","h1_nreco",6,thetaLBins);
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ) h1_nacc.Fill(gCosThetaL);
		if (BMass != 0 && ((Mumumass > 3.096916+3.*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) 
			&& (Mumumass > 3.686109+3.*Mumumasserr || Mumumass < 3.686109- 3.*Mumumasserr)) ) {
			h1_nreco.Fill(gCosThetaL);
		}
	}
	
//	Calculate efficiency
	TH1F h1_rec("h1_rec","",6,thetaLBins);
	h1_rec.SetMinimum(0.);
	h1_rec.SetTitleOffset(1.3,"XY");
	h1_rec.SetXTitle("genCosThetaL");
	h1_rec.SetYTitle("RecoEfficiency");
	for (int i = 1; i <= 6; i++) {
	//	Build from MC samples
		if (h1_nacc.GetBinContent(i) == 0 || h1_nreco.GetBinContent(i) == 0) {
			printf("WARNING: Efficiency(%d)=0, set error to be 1.\n",i);
			h1_rec.SetBinContent(i,0.);
			h1_rec.SetBinError(i,1.);
		}else{
			h1_rec.SetBinContent(i,h1_nreco.GetBinContent(i)/h1_nacc.GetBinContent(i));
			h1_rec.SetBinError(i,sqrt(h1_rec.GetBinContent(i)*(1-h1_rec.GetBinContent(i))/h1_nacc.GetBinContent(i)));
		}
	//	Build a toy sample
	//	h1_rec.SetBinContent(i,fabs(sin((thetaLBins[i]+thetaLBins[i-1])*3.1415926/4)/(20.-2*j))+0.1);
	//	h1_rec.SetBinError(i,sqrt(h1_rec.GetBinContent(i)*(1-h1_rec.GetBinContent(i))/20));
	}
//	Use Legendre polynomial for better convergance
//	1,x,(3x^2-1)/2,(5x^3-3x)/2
	TF1 f1_model("f1_model","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4",-1.,1.);
//	TF2 f2_model("f2_model","([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6",-1.,1.,-1.,1.);

//	Prepare draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
	
//	Draw efficiency
	h1_rec.SetStats(0);
//	h1_rec.SetMaximum(effUpperBound);
	h1_rec.Draw("PE1");
//	h1_rec.Fit("f1_model");
	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon_{RECO} in Bin%d",iBin));
	
//	Draw FitResult
	f1_model.SetParameter(0,0.01);
	f1_model.SetParameter(1,0.01);
	f1_model.SetParameter(2,0.01);
	f1_model.SetParameter(3,0.01);
	f1_model.SetParameter(4,0.01);
	f1_model.SetTitle("");
	
	h1_rec.Fit("f1_model"); //// 09-09
	
//	f1_model.SetMaximum(effUpperBound);
	f1_model.SetLineWidth(1);
	f1_model.Draw("SAME");
	canvas.Print(TString::Format("./plots/recoEff_bin%d.pdf",iBin));
	
//	Draw compare
	int nPar = 5;
	double arrPar[5], arrParErr[5], chi2Val;
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model.GetParameter(iPar);
		arrParErr[iPar] = f1_model.GetParError(iPar);
		chi2Val         = f1_model.GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
	TH1F h1_compFit("h1_compFit","",6,-1.,1.);
	h1_compFit.SetTitleOffset(2,"XY");
	h1_compFit.SetXTitle("genCosThetaL");
	for (int i = 1; i <= 6; i++) {//thetaL
		if (h1_rec.GetBinContent(i) != 0){
			h1_compFit.SetBinContent(i,f1_model.Eval(h1_rec.GetXaxis()->GetBinCenter(i))/h1_rec.GetBinContent(i));
		}else{
			h1_compFit.SetBinContent(i,0.);
		}
	}
	h1_compFit.SetMinimum(0.);
	h1_compFit.SetStats(0);
	h1_compFit.Draw("PE1");
	h1_compFit.SetMarkerStyle(20);
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.3,0.95,TString::Format("#varepsilon_{RECO,fit} / #varepsilon_{RECO,measured} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEff_compFit_bin%d.pdf",iBin));
	
//	Draw significance of deviation
	TH1F h_pull("Deviation/Error","",15,-3.,3.);
	h_pull.SetXTitle("Significance of deviation");
	h_pull.SetYTitle("Angular bins");
	for (int i = 1; i <= 6; i++) {//thetaL
		double _xlo = h1_rec.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h1_rec.GetXaxis()->GetBinUpEdge(i);
		if (h1_rec.GetBinContent(i) != 0){
			h_pull.Fill((f1_model.Integral(_xlo,_xhi)/(_xhi-_xlo)-h1_rec.GetBinContent(i))/h1_rec.GetBinError(i));
		}
	}
	h_pull.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEff_sigma_bin%d.pdf",iBin));
	
//	Clear
	delete latex;
	
//	prepare output
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
	return output;
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////////
*/
void createAccptanceHist() // create acceptance histogram from UNFILTERED GEN.
{//{{{
	double accUpperBound = 0.05;
	double gQ2 = 0;
	double gCosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	
//	int gen = 0;   //////////////////////////////////////////////////////////////////////////////////////////////////////////////// 12-23
//	int reco = 0;
	
	TChain *treein=new TChain("tree");
	treein->Add("../RootFiles/MC_GENOnly/BToKMuMu_GENOnly_8TeV_genonly_v3-1.root");
	if (treein == NULL) gSystem->Exit(0);
	treein->SetBranchStatus("*",0);
	treein->SetBranchStatus("genQ2"         , 1);
	treein->SetBranchStatus("genCosThetaL"  , 1);
	treein->SetBranchStatus("genMu*"        , 1);
	treein->SetBranchAddress("genQ2"        , &gQ2);
	treein->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	treein->SetBranchAddress("genMupPt"     , &gmuppt);
	treein->SetBranchAddress("genMupEta"    , &gmupeta);
	treein->SetBranchAddress("genMumPt"     , &gmumpt);
	treein->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Create histograms
	TFile *fout = new TFile("./RootFiles/acceptance_8TeV.root","RECREATE");
//	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
//	TH1F *h1_ngen[11];
//	TH1F *h1_nacc[11];
//	TH1F *h1_acc[11];
	TH1F *h1_ngen_fine[11];
	TH1F *h1_nacc_fine[11];
	TH1F *h1_acc_fine[11];
	for(int iBin = 0; iBin < 11; iBin++){
	//	if (iBin == 3 || iBin == 5) continue;    ///////////////////////////////////////////////////////////// 12-23
//		h1_ngen[iBin] = new TH1F(TString::Format("h1_ngen_bin%d",iBin),"h1_ngen",6,thetaLBins);
//		h1_nacc[iBin] = new TH1F(TString::Format("h1_nacc_bin%d",iBin) ,"h1_nacc" ,6,thetaLBins); 
//		h1_acc [iBin] = new TH1F(TString::Format("h1_acc_bin%d",iBin),"",6,thetaLBins);
		h1_ngen_fine[iBin] = new TH1F(TString::Format("h1_ngen_fine_bin%d",iBin),"h1_ngen",20,-1,1);
		h1_nacc_fine[iBin] = new TH1F(TString::Format("h1_nacc_fine_bin%d",iBin) ,"h1_nacc" ,20,-1,1); 
		h1_acc_fine[iBin]  = new TH1F(TString::Format("h1_acc_fine_bin%d",iBin),"",20,-1,1);
//		h1_ngen[iBin]->SetTitleOffset(1.3,"XY");
//		h1_ngen[iBin]->SetXTitle("genCosThetaL");
//		h1_ngen[iBin]->SetYTitle("Generated events");
//		h1_nacc[iBin]->SetTitleOffset(1.3,"XY");
//		h1_nacc[iBin]->SetXTitle("genCosThetaL");
//		h1_nacc[iBin]->SetYTitle("Events in acceptance");
//		h1_acc [iBin]->SetStats(0);
//		h1_acc [iBin]->SetMinimum(0.);
//		h1_acc [iBin]->SetMaximum(accUpperBound);
//		h1_acc [iBin]->SetTitleOffset(1.3,"XY");
//		h1_acc [iBin]->SetXTitle("genCosThetaL");
//		h1_acc [iBin]->SetYTitle("Acceptance");
		h1_ngen_fine[iBin]->SetTitleOffset(1.3,"XY");
		h1_ngen_fine[iBin]->SetXTitle("genCosThetaL");
		h1_ngen_fine[iBin]->SetYTitle("Generated events");
		h1_nacc_fine[iBin]->SetTitleOffset(1.3,"XY");
		h1_nacc_fine[iBin]->SetXTitle("genCosThetaL");
		h1_nacc_fine[iBin]->SetYTitle("Events in acceptance");
		h1_acc_fine [iBin]->SetStats(0);
		h1_acc_fine [iBin]->SetMinimum(0.);
		h1_acc_fine [iBin]->SetMaximum(accUpperBound);
		h1_acc_fine [iBin]->SetTitleOffset(1.3,"XY");
		h1_acc_fine [iBin]->SetXTitle("genCosThetaL");
		h1_acc_fine [iBin]->SetYTitle("Acceptance");
	}

//	Fill histograms
	// Read data
	for (int entry = 0; entry < treein->GetEntries(); entry++) {
		treein->GetEntry(entry);
		for(int iBin = 0; iBin < 11; iBin++){
		//	if (iBin == 3 || iBin == 5) continue;    //////////////////////////////////// 12-23
			if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
//			h1_ngen[iBin]->Fill(gCosThetaL);
			h1_ngen_fine[iBin]->Fill(gCosThetaL);
//			if (iBin != 10 && iBin != 9) gen++;   ///////////////////////////////////////////////////////////////////////////////////////////// 12-23
			if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
//				h1_nacc[iBin]->Fill(gCosThetaL);
				h1_nacc_fine[iBin]->Fill(gCosThetaL);
//				if (iBin != 10 && iBin != 9) reco++; //////////////////////////////////////////////////////////////////////////////////////////////////////// 12-23
			}
		}
	}
	for(int iBin = 0; iBin < 11; iBin++){
	//	if (iBin == 3 || iBin == 5) continue;  /////////////////////////////////////////////////////////////// 12-23
	//	Calculate acceptance
	//	h1_acc[iBin]->SetAxisRange(0.,1.,"Y");
//		for (int i = 1; i <= 6; i++) {
//		//	Fill acceptance
//			if (h1_ngen[iBin]->GetBinContent(i) == 0) {
//				printf("WARNING: Acceptance(%d)=%f/%f\n",i,h1_nacc[iBin]->GetBinContent(i),h1_ngen[iBin]->GetBinContent(i));
//				h1_acc[iBin]->SetBinContent(i,0.);
//				h1_acc[iBin]->SetBinError(i,1.);
//			}else{
//				h1_acc[iBin]->SetBinContent(i,h1_nacc[iBin]->GetBinContent(i)/h1_ngen[iBin]->GetBinContent(i));
//				if (h1_nacc[iBin]->GetBinContent(i) != 0){
//					h1_acc[iBin]->SetBinError(i,sqrt(h1_acc[iBin]->GetBinContent(i)*(1.-h1_acc[iBin]->GetBinContent(i))/h1_ngen[iBin]->GetBinContent(i)));
//				}else{
//					h1_acc[iBin]->SetBinError(i,0.);
//				}
//			}
//		}
//		printf("INFO: h1_acc_bin%d built.\n",iBin);
//		
//	//	h1_acc_fine[iBin]->SetAxisRange(0.,1.,"Y");
		for (int i = 1; i <= 20; i++) {//L
		//	Fill acceptance
			if (h1_nacc_fine[iBin]->GetBinContent(i) == 0 || h1_ngen_fine[iBin]->GetBinContent(i) == 0) {
				h1_acc_fine[iBin]->SetBinContent(i,0.);
				h1_acc_fine[iBin]->SetBinError(i,1.);
			}else{
				h1_acc_fine[iBin]->SetBinContent(i,h1_nacc_fine[iBin]->GetBinContent(i)/h1_ngen_fine[iBin]->GetBinContent(i));
				if (h1_nacc_fine[iBin]->GetBinContent(i) != 0){
					h1_acc_fine[iBin]->SetBinError(i,sqrt(h1_acc_fine[iBin]->GetBinContent(i)*(1.-h1_acc_fine[iBin]->GetBinContent(i))/h1_ngen_fine[iBin]->GetBinContent(i)));
				}else{
					h1_acc_fine[iBin]->SetBinError(i,0.);
				}
			}
	//		if ( i == 20 && iBin == 0) { cout<<h1_acc_fine[iBin]->GetBinContent(i)<<endl<<h1_nacc_fine[iBin]->GetBinContent(i)<<endl<<h1_ngen_fine[iBin]->GetBinContent(i)<<endl; }
		}
		printf("INFO: h1_acc_fine_bin%d built.\n",iBin);
//	Draw FitResult
	double up , dn ;
	up = 1.; dn = -1.;
//	if (iBin == 0) {up = 0.85; dn = -0.85;}
//	else if (iBin ==1 || iBin == 4 || iBin == 6 || iBin == 9) { up = 0.95; dn = -0.95;}
//	else {up = 1.; dn = -1.;}
//	TString f1_model_format_0 = "[0] + [1]*x + [2]*x**2 + [3]*x**3 + [4]*x**4";
	TString f1_model_format_0 = "[0] + [1]*exp(-0.5*((x-[2])/[3])^2) ";
	const int nPar = 4;
	TF1 *f1_model = new TF1 ("f1_model", f1_model_format_0, dn, up);
	f1_model->FixParameter(0,0.);
//	f1_model->SetParameter(1,0.01);
//	f1_model->SetParameter(2,0.);
	f1_model->SetParameter(3,10);
//	f1_model->SetParameter(4,0.01);  // f1_model_format_0
		
	TCanvas canvas("canvas");
	h1_acc_fine[iBin]->SetMinimum(0.);
	h1_acc_fine[iBin]->SetTitleOffset(1.3,"XY");
	h1_acc_fine[iBin]->SetXTitle("genCosThetaL");
	h1_acc_fine[iBin]->SetYTitle("Acceptance");
	h1_acc_fine[iBin]->SetStats(0);
	h1_acc_fine[iBin]->GetYaxis()->SetLabelFont(22);
	h1_acc_fine[iBin]->GetYaxis()->SetLabelSize(0.04);
	h1_acc_fine[iBin]->GetYaxis()->SetTitleSize(0.04);
	h1_acc_fine[iBin]->GetYaxis()->SetTitleOffset(1.2);
	h1_acc_fine[iBin]->GetYaxis()->SetTitleFont(22);
	h1_acc_fine[iBin]->GetXaxis()->SetLabelFont(22);
	h1_acc_fine[iBin]->GetXaxis()->SetLabelSize(0.04);
	h1_acc_fine[iBin]->GetXaxis()->SetTitleSize(0.04);
	h1_acc_fine[iBin]->GetXaxis()->SetTitleOffset(1.15);
	h1_acc_fine[iBin]->GetXaxis()->SetTitleFont(22);
	//	h1_acc_fine[iBin]->SetMaximum(effUpperBound2[iBin]);
	//	if (iBin == 0) h1_acc_fine[iBin]->SetMaximum(0.01);
	h1_acc_fine[iBin]->Draw("PE1");
	//	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
	//	h1_acc_fine[iBin]->Draw();
	h1_acc_fine[iBin]->Fit(f1_model,"R"); //// 09-09
	
	f1_model->SetTitle("");
//	f1_model->SetMaximum(effUpperBound[iBin]); //03-11
	f1_model->SetLineWidth(1);
//	f1_model->SetRange(-0.89,0.79);
//	if (iBin == 1) f1_model->SetRange(-0.89,0.89);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
//	Draw compare
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
	}
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
//	return output;
	writeParam(iBin,"acc",   arrPar,   4);  // f1_model_format_0
	writeParam(iBin,"accErr",arrParErr,4);  // f1_model_format_0
//	return output.c_str();	
		canvas.Update();
//		canvas.Print(TString::Format("./plots/accXrecoEff_accL_fine_bin%d.pdf",iBin));
		canvas.Print(TString::Format("./plots/accXrecoEff_accL_fine_bin%d.png",iBin));
	}
//	cout<<" gen = "<<gen<<endl;   ///////////////////////////////////////////////////////////////////////////////// 12-23
//	cout<<"reco = "<<reco<<endl;
	fout->Write();
	fout->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

void accXrecoEff(int iBin)
//std::vector<double> accXrecoEff(int iBin)
{//{{{
//	TH1::SetDefaultSumw2();
//	ROOT::Math::Minimizer::SetDefaultMaxFunctionCalls(10000);
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	double effUpperBound[11]  = {0.5e-3, 0.5e-3, 0.5e-3, 0.5e-3, 0.5e-3, 0.5e-3, 0.8e-3, 1.5e-3, 2.0e-3, 0.5e-3, 0.5e-3};
	double effUpperBound2[11] = {  0.18,   0.08,   0.04,   0.05,  0.035,   0.05,  0.035,   0.05,   0.06,   0.06,   0.04};
//	TLorentzVector B_4vec;
//	double Bctau = 0;
	double BMass = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double Q2 = 0;
	double gCosThetaL = 0;
	double CosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	
//	ch->SetBranchStatus("*",0);
//	ch->SetBranchStatus("B_4vec"        , 1);
//	ch->SetBranchStatus("Bctau"         , 1);
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("genCosThetaL"  , 1);
	ch->SetBranchStatus("CosThetaL"     , 1);
	ch->SetBranchStatus("genMu*"        , 1);
//	ch->SetBranchStatus("B_4vec"        , &B_4vec);
//	ch->SetBranchStatus("Bctau"         , &Bctau);
	ch->SetBranchAddress("Bmass"        , &BMass);
	ch->SetBranchAddress("Mumumass"     , &Mumumass);
	ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
	ch->SetBranchAddress("Q2"           , &Q2);
	ch->SetBranchAddress("CosThetaL"    , &CosThetaL);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Load acceptance
	TFile f_acc("./RootFiles/acceptance_8TeV.root");
//	TH1F *h1_acc       = (TH1F*)f_acc.Get(TString::Format("h1_acc_bin%d",      iBin));
	TH1F *h1_acc_fine  = (TH1F*)f_acc.Get(TString::Format("h1_acc_fine_bin%d", iBin));
//	TH1F *h1_ngen      = (TH1F*)f_acc.Get(TString::Format("h1_ngen_bin%d",     iBin));
//	TH1F *h1_ngen_fine = (TH1F*)f_acc.Get(TString::Format("h1_ngen_fine_bin%d",iBin));
	
//	Fill histograms
//	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
//	TH1F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins); 
//	TH1F h2_nreco("h2_nreco","h2_nreco",6,thetaLBins);
	int nLBins = 20;// The same value as h_acc
	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,-1,1); 
	TH1F h2_nreco_fine("h2_nreco_fine","h2_nreco_fine",nLBins,-1,1);
	h2_nacc_fine.SetStats(0);
	h2_nacc_fine.SetMinimum(0.);
	h2_nacc_fine.SetXTitle("CosThetaL");
	h2_nacc_fine.SetYTitle("Acceptated Events/0.2");
	h2_nreco_fine.SetStats(0);
	h2_nreco_fine.SetMinimum(0.);
	h2_nreco_fine.SetXTitle("CosThetaL");
	h2_nreco_fine.SetYTitle("Reconstructed Events/0.2");
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 <= Q2rangedn[iBin]) continue;
		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
//			h2_nacc.Fill(gCosThetaL);
			h2_nacc_fine.Fill(gCosThetaL);
		}
	//	if (BMass != -999 && ((Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)) ){
		if (BMass != -999 ){    ////////////////////////   12-10 N.A.
//			h2_nreco.Fill(CosThetaL);
			h2_nreco_fine.Fill(CosThetaL);
		}
	}
	
//	Calculate efficiency
//	TH1F h2_eff("h2_eff","",6,thetaLBins);
//	for (int i = 1; i <= 6; i++) {//L
//	//	Build from MC samples
//		if (h2_nacc.GetBinContent(i) == 0 || h2_nreco.GetBinContent(i) == 0) {
//			printf("WARNING: Efficiency(%d)0, set error to be 1.\n",i);
//			h2_eff.SetBinContent(i,0.);
//			h2_eff.SetBinError(i,1.);
//		}else{
//			h2_eff.SetBinContent(i,h2_nreco.GetBinContent(i)/h2_nacc.GetBinContent(i) * h1_acc->GetBinContent(i));
//			h2_eff.SetBinError(i,h2_eff.GetBinContent(i)*sqrt(-1./h2_nacc.GetBinContent(i)+1./h2_nreco.GetBinContent(i)+pow(h1_acc->GetBinError(i)/h1_acc->GetBinContent(i),2)));
//			printf("INFO: Efficiency(%d) : %f +- %f.\n",i,h2_eff.GetBinContent(i),h2_eff.GetBinError(i));
//		}
//	}
//	
	TH1F h2_reco_fine("h2_reco_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
	//	if (!(Selection[1])) continue;   //////////////////////////////  12-08
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_reco_fine.SetBinContent(i,0.);
			h2_reco_fine.SetBinError(i,1.);
		//	h2_eff_fine.SetBinError(i,0.001e-3); ////////////////////////////////////////////////////////////////// test
		}else{
			h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
			printf("INFO: recoEfficiency_fine(%d)=%f +- %f.\n",i,h2_reco_fine.GetBinContent(i),h2_reco_fine.GetBinError(i));
		}
	//	if ( i == 20 && iBin == 0) { cout<<h2_reco_fine.GetBinContent(i)<<endl<<h2_nreco_fine.GetBinContent(i)<<endl<<h2_nacc_fine.GetBinContent(i)<<endl; }
	}
	
	TH1F h2_eff_fine("h2_eff_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
	//	if (!(Selection[1])) continue;   //////////////////////////////  12-08
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_eff_fine.SetBinContent(i,0.);
			h2_eff_fine.SetBinError(i,1.);
		//	h2_eff_fine.SetBinError(i,0.001e-3); ////////////////////////////////////////////////////////////////// test
		}else{
			h2_eff_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i) * h1_acc_fine->GetBinContent(i));
			h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*sqrt(-1./h2_nacc_fine.GetBinContent(i)+1./h2_nreco_fine.GetBinContent(i)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
		//	h2_eff_fine.SetBinContent(i,h2_reco_fine.GetBinContent(i)*h1_acc_fine->GetBinContent(i));
		//	h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i));
			printf("INFO: Efficiency_fine(%d)=%f +- %f.\n",i,h2_eff_fine.GetBinContent(i),h2_eff_fine.GetBinError(i));
		}
	}
	TString f1_model_format_1 ;
	TString f1_model_format_2 ;
	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
//	if (iBin == 0 || iBin ==2 || iBin == 4 || iBin == 8 ) { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_1 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8. + [6]*(231*x**6-315*x**4+105*x**2-5)/16."; 
		f1_model_format_2 = "( [0]*exp(-0.5*((x-[1])/[2])**2) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
	//	f1_model_format_2 = "( [0]*exp(-0.5*((x-[1])/[2])**2) ) * ( [3]+[4]*x+[5]*(3*x**2-1)/2.+[6]*(5*x**3-3*x)/2.+[7]*(35*x**4-30*x**2+3)/8.+[8]*(63*x**5-70*x**3+15*x)/8. +[9]*(231*x**6-315*x**4+105*x**2-5)/16. )";
	//	f1_model_format_2 = "( ([0]*exp(-0.5*((x-[1])/[2])**2) ) + ([3]*exp(-0.5*((x-[4])/[5])**2)) + ( [6]*exp(-0.5*((x-[7])/[8])**2)) ) + [9]";
	}
	if (iBin ==20 ) { 
		f1_model_format_1 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8. + [6]*(231*x**6-315*x**4+105*x**2-5)/16."; 
	//	f1_model_format_1 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8. + [6]*x**6"; 
	//	f1_model_format_2 = "( [0] + [1]*exp(-0.5*((x-[2])/[3])^2) ) * ( [4]+[5]*x+[6]*(3*x**2-1)/2.+[7]*(5*x**3-3*x)/2.+[8]*(35*x**4-30*x**2+3)/8.+[9]*(63*x**5-70*x**3+15*x)/8. +[10]*(231*x**6-315*x**4+105*x**2-5)/16.)";
		f1_model_format_2 = "( ( [0]*exp(-0.5*((x-[1])/[2])**2) ) + ([3]*exp(-0.5*((x-[4])/[5])**2)) + ( [6]*exp(-0.5*((x-[7])/[8])**2)) ) + [9]";
	}else { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
		f1_model_format_2 = "( [0]*exp(-0.5*((x-[1])/[2])**2) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
	}
//	TString f1_model_format_5 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8.";
//	TString f1_model_format_6 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8.+[6]*(231*x**6-315*x**4+105*x**2-5)/16.";
//	TString f1_model_format_7 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8.+[6]*(231*x**6-315*x**4+105*x**2-5)/16.+[7]*(427*x**7-693*x**5+315*x**3-35*x)/16.";
//	TString f1_model_format_8 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8.+[6]*(231*x**6-315*x**4+105*x**2-5)/16.+[7]*(427*x**7-693*x**5+315*x**3-35*x)/16.+[8]*(6435*x**8-12012*x**6+6930*x**4-1260*x**2+35)/128.+[9]*(12155*x**9-25740*x**7+18010*x**5-4620*x**3+315*x)/128.";
//	Draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
//	double chi2Val=0;
//	fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
//	printf("Chi2(Bin center)=%f \n",chi2Val);

//	Draw #events in acceptance
	h2_nacc_fine.Draw();
	h2_nacc_fine.SetLabelFont(22,"XY");
	h2_nacc_fine.SetLabelSize(0.04,"XY");
	h2_nacc_fine.SetTitleSize(0.04,"XY");
	h2_nacc_fine.SetTitleFont(22,"XY");
	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff_naccL_fine_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_naccL_fine_bin%d.png",iBin));
	
//	Draw #events pass all cuts
	h2_nreco_fine.Draw();
	h2_nreco_fine.SetLabelFont(22,"XY");
	h2_nreco_fine.SetLabelSize(0.04,"XY");
	h2_nreco_fine.SetTitleSize(0.04,"XY");
	h2_nreco_fine.SetTitleFont(22,"XY");
	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff_nrecoL_fine_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_nrecoL_fine_bin%d.png",iBin));

//	Draw FitResult for recoEfficiency
	double up_r , dn_r ;
	up_r = 1.; dn_r = -1.;
	const int nPar_r = 7;
	TF1 *f1_model_r = new TF1 ("f1_model_r", f1_model_format_1, dn_r, up_r);
	f1_model_r->SetParameter(0,0.);
	f1_model_r->SetParameter(1,0.01);
	f1_model_r->SetParameter(2,0.01);
	f1_model_r->SetParameter(3,0.01);
	f1_model_r->SetParameter(4,0.01);
	f1_model_r->SetParameter(5,0.01);
	f1_model_r->SetParameter(6,0.01);  // f1_model_format_1
		
	h2_reco_fine.Fit(f1_model_r,"R"); 
	
	f1_model_r->SetTitle("");
	f1_model_r->SetMaximum(effUpperBound[iBin]); 
	f1_model_r->SetLineWidth(1);
	f1_model_r->SetLineColor(2);
	f1_model_r->Draw(" SAME ");
	
	h2_reco_fine.SetMinimum(0.);
	h2_reco_fine.SetTitleOffset(1.3,"XY");
	h2_reco_fine.SetXTitle("CosThetaL");
	h2_reco_fine.SetYTitle("recoEfficiency");
	h2_reco_fine.SetStats(0);
	h2_reco_fine.SetMaximum(effUpperBound2[iBin]);
	h2_reco_fine.Draw("PE1");
	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
	h2_reco_fine.Draw();
	h2_reco_fine.SetLabelFont(22,"XY");
	h2_reco_fine.SetLabelSize(0.04,"XY");
	h2_reco_fine.SetTitleSize(0.04,"XY");
	h2_reco_fine.SetTitleFont(22,"XY");
	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff_recoL_fine_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_recoL_fine_bin%d.png",iBin));

//	Save fitting results
	double chi2Val_r=0;
	double arrPar_r[nPar_r], arrParErr_r[nPar_r];
	for (int iPar = 0; iPar < nPar_r; iPar++) {
		arrPar_r[iPar]    = f1_model_r->GetParameter(iPar);
		arrParErr_r[iPar] = f1_model_r->GetParError(iPar);
		chi2Val_r         = f1_model_r->GetChisquare();
	}
	std::vector<double> output_r;
	for (int iPar = 0; iPar < nPar_r; iPar++){
		output_r.push_back(arrPar_r[iPar]);
		output_r.push_back(arrParErr_r[iPar]);
		printf("%18.15f,",arrPar_r[iPar]);
		if (iPar+1 >= nPar_r) printf("\n");
	}
	for (int i = 0; i < output_r.size(); i=i+2) {
		printf("%18.15f,",output_r[i+1]);
		if (i+2 >= output_r.size()) printf("\n");
	}
	writeParam(iBin,"reco",   arrPar_r,   nPar_r);  // f1_model_format_2
	writeParam(iBin,"recoErr",arrParErr_r,nPar_r);  // f1_model_format_2
/*
	TH1F h_theoL("h_theoL" ,"h_theoL" ,nLBins,-1,1); 
	h_theoL.SetStats(0);
	h_theoL.SetMinimum(0.);
	h_theoL.SetXTitle("CosThetaL");
	h_theoL.SetYTitle("#Events / 0.2");
	for (int lBin = 1; lBin <= nLBins; lBin++) {
		h_theoL.SetBinContent(lBin,h_ngenL->GetBinContent(lBin)*h_effL.GetBinContent(lBin));
		if (h_effL.GetBinContent(lBin) != 0){
			h_theoL.SetBinError(lBin,h_theoL.GetBinContent(lBin)*h_effL.GetBinError(lBin)/h_effL.GetBinContent(lBin));
		}else{
			h_theoL.SetBinError(lBin,sqrt(h_theoL.GetBinContent(lBin)));
		}
	}
	h_theoL.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_theoL_bin%d.pdf",iBin));
*/	
//	Draw efficiency
//	h2_eff.SetMinimum(0.);
//	h2_eff.SetTitleOffset(1.3,"XY");
//	h2_eff.SetXTitle("genCosThetaL");
//	h2_eff.SetYTitle("Efficiency");
//	h2_eff.SetStats(0);
//	h2_eff.SetMaximum(effUpperBound[iBin]);  //03-11
//	h2_eff.Draw("PE1");
//	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
///////////////////////////////////////////////////////////////////////////////////////////////

//	Draw FitResult for Total Efficiency
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
	double up , dn ;
	up = 1.; dn = -1.;
	const int nPar = 10;
//	TF1 f1_model("f1_model", f1_model_format_7, dn, up);
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, dn, up);
	f1_model->SetParameter(0, readParam(iBin,"acc", 1));
	f1_model->SetParameter(1, readParam(iBin,"acc", 2));
	f1_model->SetParameter(2, readParam(iBin,"acc", 3));
	f1_model->SetParameter(3, readParam(iBin,"reco", 0));
	f1_model->SetParameter(4, readParam(iBin,"reco", 1));
	f1_model->SetParameter(5, readParam(iBin,"reco", 2));
	f1_model->SetParameter(6, readParam(iBin,"reco", 3));
	f1_model->SetParameter(7, readParam(iBin,"reco", 4));
	f1_model->SetParameter(8, readParam(iBin,"reco", 5));
	f1_model->SetParameter(9, readParam(iBin,"reco", 6));
	if ( iBin != 2 && iBin != 4 && iBin != 7 ) {
		f1_model->FixParameter(0, readParam(iBin,"acc", 1));
		f1_model->FixParameter(1, readParam(iBin,"acc", 2));
		f1_model->FixParameter(2, readParam(iBin,"acc", 3));
		f1_model->FixParameter(3, readParam(iBin,"reco", 0));
		f1_model->FixParameter(4, readParam(iBin,"reco", 1));
		f1_model->FixParameter(5, readParam(iBin,"reco", 2));
		f1_model->FixParameter(6, readParam(iBin,"reco", 3));
		f1_model->FixParameter(7, readParam(iBin,"reco", 4));
		f1_model->FixParameter(8, readParam(iBin,"reco", 5));
		f1_model->FixParameter(9, readParam(iBin,"reco", 6));
	}

	h2_eff_fine.SetStats(0);
	h2_eff_fine.SetMinimum(0.);
	h2_eff_fine.SetTitleOffset(1.3,"XY");
	h2_eff_fine.SetXTitle("CosThetaL");
	h2_eff_fine.SetYTitle("Efficiency");
	h2_eff_fine.SetMaximum(effUpperBound[iBin]); //03-11
//	h2_eff_fine.Draw("TEXT");
	h2_eff_fine.Draw("PE1");
	h2_eff_fine.SetLabelFont(22,"XY");
	h2_eff_fine.SetLabelSize(0.04,"XY");
	h2_eff_fine.SetTitleSize(0.04,"XY");
	h2_eff_fine.SetTitleFont(22,"XY");
	
//	f1_model->SetDefaultMaxFunctionCalls(10000);
//	ROOT::Math::MinimizerOptions::SetDefaultErrorDef(10000);
//	h2_eff_fine.Fit(f1_model,"WL S R"); //// 09-09
	h2_eff_fine.Fit(f1_model,"R"); //// 09-09
//	TFitResultPtr r = h2_eff_fine.Fit(f1_model,"WL S R");
//	r->Print();
	
	f1_model->SetTitle("");
	f1_model->SetMaximum(effUpperBound[iBin]); //03-11
	f1_model->SetLineWidth(1);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d.png",iBin));
/////////////////////////////////////////////////////////////////////////////////////////////////
	
//	Save Fitting results
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
//	Draw compare
	TH1F h2_compFit("h2_compFit","",nLBins,-1.,1.);
	h2_compFit.SetTitleOffset(1.3,"XY");
	h2_compFit.SetXTitle("genCosThetaL");
	TH1F h2_pullFit("h2_pullFit","",nLBins,-1.,1.);
	h2_pullFit.SetTitleOffset(1.3,"XY");
	h2_pullFit.SetXTitle("genCosThetaL");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		if (h2_eff_fine.GetBinContent(i) != 0){
			h2_compFit.SetBinContent(i,f1_model->Eval(h2_eff_fine.GetXaxis()->GetBinCenter(i))/h2_eff_fine.GetBinContent(i));
			double _xlo = h2_eff_fine.GetXaxis()->GetBinLowEdge(i);
			double _xhi = h2_eff_fine.GetXaxis()->GetBinUpEdge(i);
			h2_pullFit.SetBinContent(i,(f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff_fine.GetBinContent(i))/h2_eff_fine.GetBinError(i));
		}else{
			h2_compFit.SetBinContent(i,0.);
			h2_pullFit.SetBinContent(i,0.);
		}
	}
	h2_compFit.SetMinimum(0.);
	h2_compFit.SetStats(0);
	h2_compFit.SetMarkerStyle(20);
	h2_compFit.SetMarkerSize(1.0);
	h2_compFit.Draw("PE1");
	h2_compFit.SetTitleOffset(1.3,"XY");
	h2_compFit.SetLabelFont(22,"XY");
	h2_compFit.SetLabelSize(0.04,"XY");
	h2_compFit.SetTitleSize(0.04,"XY");
	h2_compFit.SetTitleFont(22,"XY");
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.3,0.95,TString::Format("#varepsilon_{fit} / #varepsilon_{measured} in Bin%d",iBin));
	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff_compFit_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_compFit_bin%d.png",iBin));
	
	h2_pullFit.SetStats(0);
	h2_pullFit.Draw(" HIST TEXT");
	h2_pullFit.SetTitleOffset(1.3,"XY");
	h2_pullFit.SetLabelFont(22,"XY");
	h2_pullFit.SetLabelSize(0.04,"XY");
	h2_pullFit.SetTitleSize(0.04,"XY");
	h2_pullFit.SetTitleFont(22,"XY");
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.3,0.95,TString::Format("(#varepsilon_{fit} - #varepsilon_{measured})/Error in Bin%d",iBin));
	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff_pullFit_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_pullFit_bin%d.png",iBin));
	
//	Draw significance of deviation
	TH1F h2_pull("Deviation/Error","",15,-3.,3.);
	h2_pull.SetXTitle("Significance of deviation");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		double _xlo = h2_eff_fine.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h2_eff_fine.GetXaxis()->GetBinUpEdge(i);
		if (h2_eff_fine.GetBinContent(i) != 0){
			h2_pull.Fill((f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff_fine.GetBinContent(i))/h2_eff_fine.GetBinError(i));
		}
	}
	h2_pull.Draw("HIST E1");
	h2_pull.SetTitleOffset(1.3,"XY");
	h2_pull.SetLabelFont(22,"XY");
	h2_pull.SetLabelSize(0.04,"XY");
	h2_pull.SetTitleSize(0.04,"XY");
	h2_pull.SetTitleFont(22,"XY");
	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff_sigma_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_sigma_bin%d.png",iBin));
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	delete latex;
	
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
//	return output;
	writeParam(iBin,"accXrecoEff",   arrPar,   nPar);  // f1_model_format_2
	writeParam(iBin,"accXrecoEffErr",arrParErr,nPar);  // f1_model_format_2
//	return output.c_str();	
}//}}}

std::vector<double> angular_reco_bin(int iBin, const char outfile[] = "angular_reco")
{//{{{
	double up , dn ;
//	if (iBin == 0) {up = 0.80; dn = -0.80;} // 19-09
//	else if (iBin == 1) { up = 0.85; dn = -0.85;}
//	else if (iBin == 9) { up = 0.95; dn = -0.95;}
//	else { up = 1.; dn = -1.;}
	up = 1., dn = -1.;
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", dn, up);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar fh("fh", "F_{H}", genFh[iBin], 0., 1.);
	RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -1., 1.);
	
	RooRealVar nsig("nsig","nsig",1E6,1E1,1E9);
//	RooRealVar nbkg("nbkg","nbkg",10,0.1,1E4);
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	Acceptance
	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));   // f1_model_format_1
	accP0.setError(readParam(iBin,"accErr", 0));
	accP1.setError(readParam(iBin,"accErr", 1));
	accP2.setError(readParam(iBin,"accErr", 2));
	accP3.setError(readParam(iBin,"accErr", 3));  // f1_model_format_1
//	reco Efficiency
	RooRealVar recoP0("recoP0","recoP0",readParam(iBin,"reco", 0));
	RooRealVar recoP1("recoP1","recoP1",readParam(iBin,"reco", 1));
	RooRealVar recoP2("recoP2","recoP2",readParam(iBin,"reco", 2));
	RooRealVar recoP3("recoP3","recoP3",readParam(iBin,"reco", 3));   // f1_model_format_1
	RooRealVar recoP4("recoP4","recoP4",readParam(iBin,"reco", 4));   // f1_model_format_2  // f1_model_format_4
	RooRealVar recoP5("recoP5","recoP5",readParam(iBin,"reco", 5));   // f1_model_format_2  // f1_model_format_5  
	RooRealVar recoP6("recoP6","recoP6",readParam(iBin,"reco", 6));   // f1_model_format_2  // f1_model_format_6
	recoP0.setError(readParam(iBin,"recoErr", 0));
	recoP1.setError(readParam(iBin,"recoErr", 1));
	recoP2.setError(readParam(iBin,"recoErr", 2));
	recoP3.setError(readParam(iBin,"recoErr", 3));  // f1_model_format_1
	recoP4.setError(readParam(iBin,"recoErr", 4));  // f1_model_format_2  // f1_model_format_4
	recoP5.setError(readParam(iBin,"recoErr", 5));  // f1_model_format_2  // f1_model_format_5 
	recoP6.setError(readParam(iBin,"recoErr", 6));  // f1_model_format_2  // f1_model_format_6
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////

/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
//	Total Efficiency
	RooRealVar effP0("effP0","effP0",readParam(iBin,"accXrecoEff", 0));
	RooRealVar effP1("effP1","effP1",readParam(iBin,"accXrecoEff", 1));
	RooRealVar effP2("effP2","effP2",readParam(iBin,"accXrecoEff", 2));
	RooRealVar effP3("effP3","effP3",readParam(iBin,"accXrecoEff", 3));   // f1_model_format_1
	RooRealVar effP4("effP4","effP4",readParam(iBin,"accXrecoEff", 4));   // f1_model_format_2  // f1_model_format_4
	RooRealVar effP5("effP5","effP5",readParam(iBin,"accXrecoEff", 5));   // f1_model_format_2  // f1_model_format_5  
	RooRealVar effP6("effP6","effP6",readParam(iBin,"accXrecoEff", 6));   // f1_model_format_2  // f1_model_format_6
	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff", 7));   // f1_model_format_3  // f1_model_format_7
	RooRealVar effP8("effP8","effP8",readParam(iBin,"accXrecoEff", 8));   // f1_model_format_8
	RooRealVar effP9("effP9","effP9",readParam(iBin,"accXrecoEff", 9));   // f1_model_format_8
//	RooRealVar effP10("effP10","effP10",readParam(iBin,"accXrecoEff", 10));   // f1_model_format_2
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  // f1_model_format_1
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  // f1_model_format_2  // f1_model_format_4
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5));  // f1_model_format_2  // f1_model_format_5 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6));  // f1_model_format_2  // f1_model_format_6
	effP7.setError(readParam(iBin,"accXrecoEffErr", 7));  // f1_model_format_3  // f1_model_format_7
	effP8.setError(readParam(iBin,"accXrecoEffErr", 8));  // f1_model_format_8
	effP9.setError(readParam(iBin,"accXrecoEffErr", 9));  // f1_model_format_8
//	effP10.setError(readParam(iBin,"accXrecoEffErr", 10));  // f1_model_format_2

//	RooArgSet f_effA_argset(CosThetaL);
//	f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));
//	f_effA_argset.add(RooArgSet(effP8, effP9));      // f1_model_format8_
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	RooArgSet f_accA_argset(CosThetaL);
//	f_accA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
//	RooArgSet f_recoA_argset(CosThetaL);
//	f_recoA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
	
//	RooGenericPdf f_acc("f_acc", "accP0 + accP1 *exp(-0.5*((CosThetaL-accP2)/accP3)^2) ", f_accA_argset);     // f1_model_format_0
//	RooGenericPdf f_reco("f_reco", "recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6", f_recoA_argset);     // f1_model_format_0
//	RooGenericPdf f_sigA("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb));
//	//RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb),TString::Format("fabs(%s) <= (%s)/2.",afb,fh);
	
//	RooProdPdf    f_eff("f_eff","", f_acc, f_reco);
//	//RooProdPdf    f_effXsig("f_effXsig","", f_acc, f_reco, f_sig);
//	RooProdPdf    f_sig("f_effXsig","", f_eff, f_sigA); 
//	RooExtendPdf  f("f","", f_sig, nsig);
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////

/////////////	For bin 0, 1, 9 ////////////////////////////////////////////////////////////////////////
//	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
//	if (iBin == 0 || iBin ==1 || iBin == 9 || iBin ==2 || iBin == 4 || iBin == 7 ) { 
	if (iBin == 2 || iBin == 4 || iBin == 7 ) { 
		RooArgSet f_effA_argset(CosThetaL, fh, afb);
		f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));
		f_effA_argset.add(RooArgSet(effP8, effP9));      // f1_model_format8_
	//	RooGenericPdf f_sig("f_sig", "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL ) * (( effP0*exp(-0.5*((CosThetaL-effP1)/effP2)**2) ) + ( effP3*exp(-0.5*((CosThetaL-effP4)/effP5)**2)) + ( effP6*exp(-0.5*((CosThetaL-effP7)/effP8)**2) ) + effP9)", f_effA_argset); 
	//	RooGenericPdf f_eff("f_eff", "effP0+effP1*CosThetaL+effP2*(3*CosThetaL**2-1)/2+effP3*(5*CosThetaL**3-3*CosThetaL)/2+effP4*(35*CosThetaL**4-30*CosThetaL**2+3)/8.+effP5*(63*CosThetaL**5-70*CosThetaL**3+15*CosThetaL)/8.+effP6*(231*CosThetaL**6-315*CosThetaL**4+105*CosThetaL**2-5)/16.+effP7*(429*CosThetaL**7-693*CosThetaL**5+315*CosThetaL**3-35*CosThetaL)/16.", f_effA_argset); // f1_model_format_7
	//	RooGenericPdf f_sig("f_sig", "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL ) *( effP0 + effP1*exp(-0.5*((CosThetaL-effP2)/effP3)^2) ) * (effP4+effP5*CosThetaL+effP6*(3*CosThetaL**2-1)/2+effP7*(5*CosThetaL**3-3*CosThetaL)/2+effP8*(35*CosThetaL**4-30*CosThetaL**2+3)/8.+effP9*(63*CosThetaL**5-70*CosThetaL**3+15*CosThetaL)/8.+effP10*(231*CosThetaL**6-315*CosThetaL**4+105*CosThetaL**2-5)/16. )", f_effA_argset); 
	//	RooGenericPdf f_sig("f_sig", "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL ) *( effP0 *exp(-0.5*((CosThetaL-effP1)/effP2)**2) ) * (effP3+effP4*CosThetaL+effP5*(3*CosThetaL**2-1)/2+effP6*(5*CosThetaL**3-3*CosThetaL)/2+effP7*(35*CosThetaL**4-30*CosThetaL**2+3)/8.+effP8*(63*CosThetaL**5-70*CosThetaL**3+15*CosThetaL)/8.+effP9*(231*CosThetaL**6-315*CosThetaL**4+105*CosThetaL**2-5)/16. )", f_effA_argset); 
		RooGenericPdf f_sig("f_sig", "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL ) * (effP0*exp(-0.5*((CosThetaL-effP1)/effP2)**2) ) * ( effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**4+effP8*CosThetaL**5+effP9*CosThetaL**6 ) ", f_effA_argset);
		RooExtendPdf  f("f","", f_sig, nsig);
		RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);    // 12-08
	//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
	//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Hesse(0),Minos(1),Strategy(2));
	//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(2));
		RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(1));
	//	Draw the frame on the canvas
		TCanvas* c = new TCanvas("c");
		TLatex *t1 = new TLatex();
		t1->SetNDC();
		double fixNDC = -0.;
		double UpperBound[11] ={350,800,1400,0.,1000,0.,800,1000,1200,1600,6000};
		
		RooPlot* framecosl = CosThetaL.frame(); 
		data->plotOn(framecosl,Binning(100)); 
	//	data->plotOn(framecosl,Binning(85)); 
		f.plotOn(framecosl); 
		framecosl->SetTitle("");
		framecosl->SetMinimum(0);
		framecosl->SetMaximum(UpperBound[iBin]);  //11-03
		framecosl->Draw();
		framecosl->GetYaxis()->SetLabelFont(22);
		framecosl->GetYaxis()->SetLabelSize(0.04);
		framecosl->GetYaxis()->SetTitleSize(0.04);
		framecosl->GetYaxis()->SetTitleOffset(1.2);
		framecosl->GetYaxis()->SetTitleFont(22);
		framecosl->GetXaxis()->SetLabelFont(22);
		framecosl->GetXaxis()->SetLabelSize(0.04);
		framecosl->GetXaxis()->SetTitleSize(0.04);
		framecosl->GetXaxis()->SetTitleOffset(1.15);
		framecosl->GetXaxis()->SetTitleFont(22);
		
	//	fixNDC = -0.5;
	//	if (iBin > 4) fixNDC = 0.;
		if (iBin == 10) { 
			t1->DrawLatex(.30,.85+fixNDC,TString::Format(" Total Signal Region"));
		} else t1->DrawLatex(.30,.85+fixNDC,TString::Format("%s",Q2range[iBin]));
		t1->DrawLatex(.10,.79+fixNDC,TString::Format("F_{H}  =%9.5f #pm%9.5f",fh.getVal(),fh.getError()));
		t1->DrawLatex(.50,.79+fixNDC,TString::Format("A_{FB} =%9.5f #pm%9.5f",afb.getVal(),afb.getError()));
		c->Update();
	//	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
		c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
	//	clear
		delete t1;
		delete c;
		delete data;
	//	write output
		double val[3]={0,0,0};
		val[0] = fh.getVal();val[1] = fh.getError();
		writeParam(iBin, "recofh", val);
		val[0] = afb.getVal();val[1] = afb.getError();
		writeParam(iBin, "recoafb",val);
		std::vector<double> output;
		output.push_back(fh.getVal());
		output.push_back(fh.getError());
		output.push_back(afb.getVal());
		output.push_back(afb.getError());
		return output;
	}
//////////////////  For other bins  /////////////////////////////////////////////////////////////////////////////////		
//	if (iBin != 0 && iBin !=1 && iBin != 9 ) { 
//	if (iBin != 0 && iBin !=1 && iBin != 9 && iBin != 4 && iBin != 2 && iBin != 7) { 
//	if (iBin != 0 && iBin != 4 && iBin != 2 && iBin != 7) { 
	else { 
		RooArgSet f_accA_argset(CosThetaL);
		f_accA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		RooArgSet f_recoA_argset(CosThetaL);
		f_recoA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		RooGenericPdf f_acc("f_acc", "accP0 + accP1 *exp(-0.5*((CosThetaL-accP2)/accP3)^2) ", f_accA_argset);     // f1_model_format_0
		RooGenericPdf f_reco("f_reco", "recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6", f_recoA_argset);     // f1_model_format_0
		RooGenericPdf f_sigA("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb));
		//	RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb),TString::Format("fabs(%s) <= (%s)/2.",afb,fh);
		RooProdPdf    f_eff("f_eff","", f_acc, f_reco);
		RooProdPdf    f_sig("f_effXsig","", f_eff, f_sigA); 
		RooExtendPdf  f("f","", f_sig, nsig);
	//	RooGenericPdf f_sig("f_sig", "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL ) * (effP0*exp(-0.5*((CosThetaL-effP1)/effP2)**2) ) * ( effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**4+effP8*CosThetaL**5+effP9*CosThetaL**6 ) ", f_effA_argset);
	//	RooExtendPdf  f("f","", f_sig, nsig);
		RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);    // 12-08
	//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),TString::Format("(%s) && (%s)",Q2range[iBin],Selection[1]),0);    ////// 12-08
	//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
		RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(2));
	//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(1));
		if ( iBin == 0 || iBin == 9 ) {
			RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(1));
		}
	//	Draw the frame on the canvas
		TCanvas* c = new TCanvas("c");
		TLatex *t1 = new TLatex();
		t1->SetNDC();
		double fixNDC = -0.;
		double UpperBound[11] ={350,800,1400,0.,1000,0.,800,1000,1200,1600,6000};
	
		RooPlot* framecosl = CosThetaL.frame(); 
		data->plotOn(framecosl,Binning(100)); 
	//	data->plotOn(framecosl,Binning(85)); 
		f.plotOn(framecosl); 
		framecosl->SetTitle("");
		framecosl->SetMinimum(0);
		framecosl->SetMaximum(UpperBound[iBin]);  //11-03
		framecosl->Draw();
		framecosl->GetYaxis()->SetLabelFont(22);
		framecosl->GetYaxis()->SetLabelSize(0.04);
		framecosl->GetYaxis()->SetTitleSize(0.04);
		framecosl->GetYaxis()->SetTitleOffset(1.2);
		framecosl->GetYaxis()->SetTitleFont(22);
		framecosl->GetXaxis()->SetLabelFont(22);
		framecosl->GetXaxis()->SetLabelSize(0.04);
		framecosl->GetXaxis()->SetTitleSize(0.04);
		framecosl->GetXaxis()->SetTitleOffset(1.15);
		framecosl->GetXaxis()->SetTitleFont(22);
		
	//	fixNDC = -0.5;
	//	if (iBin > 4) fixNDC = 0.;
		if (iBin == 10) { 
			t1->DrawLatex(.30,.85+fixNDC,TString::Format(" Total Signal Region"));
		} else t1->DrawLatex(.30,.85+fixNDC,TString::Format("%s",Q2range[iBin]));
		t1->DrawLatex(.10,.79+fixNDC,TString::Format("F_{H}  =%9.5f #pm%9.5f",fh.getVal(),fh.getError()));
		t1->DrawLatex(.50,.79+fixNDC,TString::Format("A_{FB} =%9.5f #pm%9.5f",afb.getVal(),afb.getError()));
		c->Update();
	//	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
		c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
	//	clear
		delete t1;
		delete c;
		delete data;
	//	write output
		double val[3]={0,0,0};
		val[0] = fh.getVal();val[1] = fh.getError();
		writeParam(iBin, "recofh", val);
		val[0] = afb.getVal();val[1] = afb.getError();
		writeParam(iBin, "recoafb",val);
		
		std::vector<double> output;
		output.push_back(fh.getVal());
		output.push_back(fh.getError());
		output.push_back(afb.getVal());
		output.push_back(afb.getError());
		return output;
	}
	
	printf("recoAfb[%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"recoafb",0), fabs(readParam(iBin,"recoafb",1)));
	printf("recoFh [%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"recofh",0),  fabs(readParam(iBin,"recofh",1)));
}//}}}

void angular_reco(const char outfile[] = "angular_reco")
{//{{{
//	bool refit = false; // Turn to true if you want to fit again.
	
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double ygenfh[9], ygenuerrfh[9], ygenderrfh[9], ygenafb[9], ygenerrafb[9]; 
	double yfh[9], yuerrfh[9], yderrfh[9], yafb[9], yerrafb[9]; 

//	Check input 
	for(int i = 0, ibin = 0; i < 9, ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
	//	reco
		yafb[i]      = readParam(ibin,"recoafb",0);
		yerrafb[i]  = fabs(readParam(ibin,"recoafb",1));
		yfh[i]       = readParam(ibin,"recofh",0);
		yuerrfh[i]   = fabs(readParam(ibin,"recofh",1));
		if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		else { yderrfh[i] = yuerrfh[i]; }
		printf("recoAfb[%d]=%6.4f +- %6.4f\n",i,yafb[i],yerrafb[i]);
		printf("recoFh [%d]=%6.4f +- %6.4f\n",i,yfh[i],yuerrfh[i]);
	//	gen
		ygenafb[i]      = readParam(ibin,"genafb",0);
		ygenerrafb[i]  = fabs(readParam(ibin,"genafb",1));
		ygenfh[i]       = readParam(ibin,"genfh",0);
		ygenuerrfh[i]   = fabs(readParam(ibin,"genfh",1));
		if (ygenuerrfh[i] > fabs(ygenfh[i])) { ygenderrfh[i] = fabs(ygenfh[i]);}
		else { ygenderrfh[i] = ygenuerrfh[i]; }
		printf("genAfb[%d]=%6.4f +- %6.4f\n",i,ygenafb[i],ygenerrafb[i]);
		printf("genFh [%d]=%6.4f +- %6.4f\n",i,ygenfh[i],ygenuerrfh[i]);
	}
//	plotting
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->GetYaxis()->SetLabelFont(22);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleFont(22);
	frame->GetXaxis()->SetLabelFont(22);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleOffset(1.15);
	frame->GetXaxis()->SetTitleFont(22);
	frame->SetTitle("");
	frame->Draw();
//	reco
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.02,0.08,"Y");  //11-03
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	g_fh->SetMarkerColor(4);
	g_fh->SetMarkerStyle(20);	
	g_fh->SetFillColor(2);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
//	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Print(TString::Format("./plots/%s_fh_reco.png",outfile));
// reco and gen
	TGraphAsymmErrors *gen_fh  = new TGraphAsymmErrors(7,x,ygenfh,xerr,xerr,ygenderrfh,ygenuerrfh);
	gen_fh->SetMarkerColor(3);
	gen_fh->SetMarkerStyle(24);
	
	gen_fh->SetFillColor(1);
	gen_fh->SetFillStyle(3005);
	gen_fh->Draw("2");
	gen_fh->Draw("P");
	c->Print(TString::Format("./plots/%s_fh.png",outfile));
	c->Clear();
	
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetAxisRange(-0.02,0.02,"Y"); //11-03
	frame->Draw();
//	reco
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yerrafb,yerrafb);
	g_afb->SetMarkerColor(4);
	g_afb->SetMarkerStyle(20);
	g_afb->SetFillColor(2);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
//	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
	c->Print(TString::Format("./plots/%s_afb_reco.png",outfile));
	TGraphAsymmErrors *gen_afb = new TGraphAsymmErrors(7,x,ygenafb,xerr,xerr,ygenerrafb,ygenerrafb);
	gen_afb->SetMarkerColor(3);
	gen_afb->SetMarkerStyle(24);
	gen_afb->SetFillColor(1);
	gen_afb->SetFillStyle(3005);
	gen_afb->Draw("2");
	gen_afb->Draw("P");
	c->Print(TString::Format("./plots/%s_afb.png",outfile));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

//  04-09-2014
///////////////////////////////////////////////////////////////////////////////////////////
/*   //// 17-08-2014 N.A.
void getToyFromUnfilterGen(int iBin)
{//{{{
    bool fhatModel = true;

    // Generate random toy using unfiltered input and efficiency function
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double Q2 = 0;
    double CosThetaK = 0;
    double CosThetaL = 0;
    double Mumumass = 0;
    double Mumumasserr = 0;

    //TChain *treein=new TChain("tree");
    //treein->Add("./data/2012/sel_BuToKstarMuMu_NoGenFilter_8TeV_part*_mc.lite.root");//Contact po-hsun.chen@cern.ch to get these files.
    TChain *treein=ch;
    if (treein == NULL) {
        printf("Unfiltered MC sample is missing. Please contact pchen@cern.ch to get it.");
        return;
    }
    treein->SetBranchAddress("genQ2"        , &gQ2);
    treein->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    treein->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    treein->SetBranchAddress("Q2"           , &Q2);
    treein->SetBranchAddress("CosThetaK"    , &CosThetaK);
    treein->SetBranchAddress("CosThetaL"    , &CosThetaL);
    treein->SetBranchAddress("Mumumass"     , &Mumumass);
    treein->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
    
    // Get efficiency map, y=cosThetaK, x=cosThetaL
    std::string f2_model_format;
    if (fhatModel){
        f2_model_format = TString::Format("0.01+0.*x+0.*y");
    }else if (is7TeVCheck){
        f2_model_format = TString::Format("([0]+[1]*CosThetaK+[2]*CosThetaK**2+[3]*CosThetaK**3)+([4]+[5]*CosThetaK+[6]*CosThetaK**2+[7]*CosThetaK**3)*CosThetaL**2+([8]+[9]*CosThetaK+[10]*CosThetaK**2+[11]*CosThetaK**3)*CosThetaL**3+([12]+[13]*CosThetaK+[14]*CosThetaK**2+[15]*CosThetaK**3)*CosThetaL**4+([16]+[17]*CosThetaK+[18]*CosThetaK**2+[19]*CosThetaK**3)*CosThetaL**6");
    }else{
        f2_model_format = TString::Format("%s+([0]+[1]*CosThetaK+[2]*(3*CosThetaK**2-1)/2+[3]*(5*CosThetaK**3-3*CosThetaK)/2)+([4]+[5]*CosThetaK+[6]*(3*CosThetaK**2-1)/2+[7]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2+([8]+[9]*CosThetaK+[10]*(3*CosThetaK**2-1)/2+[11]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3+([12]+[13]*CosThetaK+[14]*(3*CosThetaK**2-1)/2+[15]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4+([16]+[17]*CosThetaK+[18]*(3*CosThetaK**2-1)/2+[19]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6",f_accXrecoEff_ord0[iBin].data());
    }
    f2_model_format = regex_replace(f2_model_format,regex("CosThetaK"),"y");
    f2_model_format = regex_replace(f2_model_format,regex("CosThetaL"),"x");
    TF2 f2_model("f2_model",f2_model_format.c_str(),-1.,1.,-1.,1.);
    if (fhatModel){
        //Set parameters here.
    }else if (is7TeVCheck){
        f2_model.SetParameters(arrRecPar2011[iBin]);
        f2_model.SetParErrors(arrRecParErr2011[iBin]);
    }else{
        for (int i = 0; i < 20; i++) {
            f2_model.SetParameter(i,readParam(iBin,"accXrecoEff",i));
            f2_model.SetParError(i,readParam(iBin,"accXrecoEffErr",i));
        }
    }

    // 
    TFile fout(TString::Format("./rndToy_Bin%d.root",iBin), "RECREATE") ;
    TTree *treeout = treein->CloneTree(0);
    int _count = 0;//number of accepted events
    int _entry = 0;
    TRandom3 *rndGenerator = new TRandom3();
    do {
        treein->GetEntry(_entry); _entry++;
        if (gQ2 > Q2rangedn[iBin] && gQ2 < Q2rangeup[iBin]) {
            if (rndGenerator->Rndm() < f2_model.Eval(gCosThetaL,gCosThetaK)) {
                Q2 = gQ2;
                CosThetaL   = gCosThetaL;
                CosThetaK   = gCosThetaK;
                Mumumass    = sqrt(gQ2);
                Mumumasserr = 0.001;
                treeout->Fill();
                _count++;
            }
        }
    } while ( _entry < treein->GetEntries() );
    fout.Write();
    fout.Close();
}//}}}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 12-10-2014 N.A.

//void angular2D_1a_Sm(int iBin, const char outfile[] = "angular2D_1a_Sm", bool keepParam = false)
void angular2D_1a_Sm(int iBin, const char outfile[] = "angular2D_1a_Sm", bool keepParam = true)
{//{{{
	// Fit to signal simulation by YsSm+YcCm to determine Sm
	RooRealVar Bmass("Bmass","M_{K^{+/-}#Mu#Mu}",5.0, 5.56);
//	RooRealVar Bmass("Bmass","M_{K^{+/-}#Mu#Mu}",5.27925-0.28, 5.27925+0.28);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	
	// Create parameters and PDFs
	   // Signal double gaussian
	RooRealVar sigGauss_mean("sigGauss_mean","M_{K#Mu#Mu}",5.27925,5.23,5.32);
	RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",.028,0.01,0.05);
//	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",.065,0.05,0.1);
	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",.065,0.01,0.15);
	RooRealVar sigM_frac("sigM_frac","sigM_frac",0.5,0.,1.);
	
	// Create signal distribution
	   // mass distro of signal
	RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
	RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
	RooAddPdf f_sigM("f_sigM","f_sigM", f_sigMGauss1, f_sigMGauss2, sigM_frac);
   
/*	// Create combinatorial background distributio
 	RooRealVar bkgCombM_c("bkgCombM_c","c1",0,-30,50);
	RooRealVar offset("offset","offset",-5.);
	RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
	RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
*/		
	RooRealVar nsig("nsig","nsig",0,1E8);
//	RooRealVar nbkg("nbkg","nbkg",0,1E4);
//	RooAddPdf f("f", "f",RooArgList(f_sigM,f_bkgCombM),RooArgList(nsig,nbkg));
	RooAddPdf f("f", "f",RooArgList(f_sigM),RooArgList(nsig));
	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),Q2range[iBin],0);
	RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));
//	RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Extended(kTRUE));  // 12-10-2014

	// Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	RooPlot* frame = Bmass.frame(); 
	data->plotOn(frame,Binning(20)); 
	f.plotOn(frame,LineColor(1)); 
	f.plotOn(frame,Components(f_sigM),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
	f.plotOn(frame,Components(f_sigM),LineStyle(2),LineColor(4),LineWidth(2));
//	f.plotOn(frame,Components(f_bkgCombM),LineColor(2),LineStyle(2),LineWidth(2));

	frame->SetTitle("");
	frame->SetMinimum(0);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->Draw();
	frame->GetYaxis()->SetLabelFont(22);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleFont(22);
	frame->GetXaxis()->SetLabelFont(22);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleOffset(1.15);
	frame->GetXaxis()->SetTitleFont(22);
	TPaveText* paveText = new TPaveText( 0.13, 0.70, 0.41, 0.88, "NDC" ); 
//	TPaveText* paveText = new TPaveText( 0.17, 0.70, 0.41, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("   nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
//	paveText->AddText(Form(" nbkg   = %.0f #pm %.0f ", nbkg.getVal(), nbkg.getError())); 
	paveText->AddText(Form(" sigmean   = %.3f #pm %.3f ", sigGauss_mean.getVal(), sigGauss_mean.getError())); 
	paveText->AddText(Form("sig_sigma1 = %.3f #pm %.3f ", sigGauss1_sigma.getVal(), sigGauss1_sigma.getError())); 
	paveText->AddText(Form("sig_sigma2 = %.3f #pm %.3f ", sigGauss2_sigma.getVal(), sigGauss2_sigma.getError())); 
	paveText->AddText(Form("    frac   = %.3f #pm %.3f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->AddText(Form("#Chi^{2}   = %.2f  ", frame->chiSquare())); 
	paveText->Draw(); 
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = 0.05;
	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
//	t1->DrawLatex(.35,.10,TString::Format("nbkg=%5.3f#pm%5.3f",nbkg.getVal(),nbkg.getError()));
//	c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_bin%d.png",outfile,iBin));
	
	// clear
	delete t1;
	delete c;
	delete data;
	// Prepare datacard
	if (keepParam){
		double val[3]={0,0,0};
		writeParam(iBin, "iBin", new double((double)iBin), 1);
		if (is7TeVCheck){
			writeParam(iBin, "mode", new double(2011), 1);
		}else{
			writeParam(iBin, "mode", new double(2012), 1);
		}
		val[0]=sigGauss1_sigma.getVal();val[1]=sigGauss1_sigma.getError();
		writeParam(iBin, "sigGauss1_sigma", val);
		val[0]=sigGauss2_sigma.getVal();val[1]=sigGauss2_sigma.getError();
		writeParam(iBin, "sigGauss2_sigma", val);
		val[0]=sigM_frac.getVal();val[1]=sigM_frac.getError();
		writeParam(iBin, "sigM_frac", val);
	}
}//}}}

//void angular2D_1b_YpPm(int iBin, const char outfile[] = "angular2D_1b_YpPm", bool keepParam = false)
void angular2D_1b_YpPm(int iBin, const char outfile[] = "angular2D_1b_YpPm", bool keepParam = true)
{//{{{
//	if (iBin ==0 || iBin%2 == 1 || iBin/2 > 3){    ///////////   Peaking bkg. Only for bin 2,4,6
	if (iBin != 10 && (iBin ==0 || iBin%2 == 1 || iBin/2 > 3)){    ///////////   Peaking bkg. Only for bin 2,4,6,9,10
//	if (iBin != 9 && iBin != 10 && (iBin ==0 || iBin%2 == 1 || iBin/2 > 3)){    ///////////   Peaking bkg. Only for bin 2,4,6,9,10
		if (keepParam){
			double val[3]={1,0,0};
			writeParam(iBin, "bkgGauss1_mean1", val);
			writeParam(iBin, "bkgGauss1_mean2", val);
			writeParam(iBin, "bkgGauss1_sigma1", val);
			writeParam(iBin, "bkgGauss1_sigma2", val);
			writeParam(iBin, "bkgM_frac1", val);
			writeParam(iBin, "bkgGauss2_mean1", val);
			writeParam(iBin, "bkgGauss2_mean2", val);
			writeParam(iBin, "bkgGauss2_sigma1", val);
			writeParam(iBin, "bkgGauss2_sigma2", val);
			writeParam(iBin, "bkgM_frac2", val);
			writeParam(iBin, "bkgM_frac12", val);
			val[0]=0;
			writeParam(iBin, "nbkgPeak", val);
		}
		return;
	}
	// Fit to control channel simulations by YpPm to determine Yp,Pm.
	RooRealVar Bmass("Bmass","M_{K^{+/-}#Mu#Mu}",5.,5.56);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	
	// Create peak background distribution
	RooRealVar bkgGauss1_mean1("bkgGauss1_mean1","M_{K#Mu#Mu}",5.05,5.,5.2);
	RooRealVar bkgGauss1_mean2("bkgGauss1_mean2","M_{K#Mu#Mu}",5.27,5.,5.35);
	RooRealVar bkgGauss2_mean1("bkgGauss2_mean1","M_{K#Mu#Mu}",5.4,5.33,5.53);
	RooRealVar bkgGauss2_mean2("bkgGauss2_mean2","M_{K#Mu#Mu}",5.27,5.17,5.37);
	RooRealVar bkgGauss1_sigma1("bkgGauss1_sigma1","#sigma_{11}",.03,.01,.08);
	RooRealVar bkgGauss1_sigma2("bkgGauss1_sigma2","#sigma_{12}",.12,.08,.50);
	RooRealVar bkgGauss2_sigma1("bkgGauss2_sigma1","#sigma_{21}",.03,.01,.05);
	RooRealVar bkgGauss2_sigma2("bkgGauss2_sigma2","#sigma_{22}",.12,.05,.50);
	RooRealVar bkgM_frac1("bkgM_frac1","bkgM_frac1",1.,0.,1.);
	RooRealVar bkgM_frac2("bkgM_frac2","bkgM_frac2",1.,0.,1.);
	RooRealVar bkgM_frac12("bkgM_frac12","bkgM_frac12",0.,0.,1.);
	RooGaussian f_bkgPeakMGauss11("f_bkgPeakMGauss11","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean1, bkgGauss1_sigma1);
	RooGaussian f_bkgPeakMGauss12("f_bkgPeakMGauss12","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean2, bkgGauss1_sigma2);
	RooGaussian f_bkgPeakMGauss21("f_bkgPeakMGauss21","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean1, bkgGauss2_sigma1);
	RooGaussian f_bkgPeakMGauss22("f_bkgPeakMGauss22","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean2, bkgGauss2_sigma2);
	RooAddPdf f_bkgPeakM1("f_bkgPeakM1","f_bkgPeakM1", RooArgList(f_bkgPeakMGauss11, f_bkgPeakMGauss12), bkgM_frac1);
	RooAddPdf f_bkgPeakM2("f_bkgPeakM2","f_bkgPeakM2", RooArgList(f_bkgPeakMGauss21, f_bkgPeakMGauss22), bkgM_frac2);
	RooAddPdf f_bkgPeakM12("f_bkgPeakM12","f_bkgPeakM12", RooArgList(f_bkgPeakM1,f_bkgPeakM2), bkgM_frac12);
	
	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",1E2,1,1E7);
	RooExtendPdf *f = 0;
	switch (iBin) {
		case 2:
		//1 double guassian ,4+4 deg. ploy
		      f = new RooExtendPdf("f","f",f_bkgPeakM1,nbkgPeak);
				break;
		case 4:
		//2 double guassian ,4+4 deg. ploy
		      //f = new RooExtendPdf("f","f",f_bkgPeakMGauss22,nbkgPeak);
				f = new RooExtendPdf("f","f",f_bkgPeakM12,nbkgPeak);
				break;
		case 6:
		//1 guassian ,2+2 deg. ploy
		      f = new RooExtendPdf("f","f",f_bkgPeakM2,nbkgPeak);
				break;
	/*	case 9:
		//1 double guassian ,4+4 deg. ploy
		      f = new RooExtendPdf("f","f",f_bkgPeakM1,nbkgPeak);
				break;
	*/	case 10:
		//2 double guassian ,4+4 deg. ploy
				f = new RooExtendPdf("f","f",f_bkgPeakM12,nbkgPeak);
				break;
		default:
		      break;
	}	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),Q2range[iBin],0);
	RooFitResult *f_fitresult = f->fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Extended());
	
	// Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	RooPlot* frame = Bmass.frame(); 
	data->plotOn(frame,Binning(20)); 
	f->plotOn(frame); 
	frame->SetTitle("");
	frame->SetMinimum(0);
	frame->Draw();
	frame->GetYaxis()->SetLabelFont(22);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleFont(22);
	frame->GetXaxis()->SetLabelFont(22);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleOffset(1.15);
	frame->GetXaxis()->SetTitleFont(22);
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = 0.05;
	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
//	c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_bin%d.png",outfile,iBin));
	
	// clear
	delete t1;
	delete c;
	delete data;
	
	if (keepParam){
		double val[3]={0,0,0};
		val[0] = bkgGauss1_mean1.getVal();val[1] = bkgGauss1_mean1.getError();
		writeParam(iBin, "bkgGauss1_mean1", val);
		val[0] = bkgGauss1_mean2.getVal();val[1] = bkgGauss1_mean2.getError();
		writeParam(iBin, "bkgGauss1_mean2", val);
		val[0] = bkgGauss1_sigma1.getVal();val[1] = bkgGauss1_sigma1.getError();
		writeParam(iBin, "bkgGauss1_sigma1", val);
		val[0] = bkgGauss1_sigma2.getVal();val[1] = bkgGauss1_sigma2.getError();
		writeParam(iBin, "bkgGauss1_sigma2", val);
		val[0] = bkgM_frac1.getVal();val[1] = bkgM_frac1.getError();
		writeParam(iBin, "bkgM_frac1", val);
		val[0] = bkgGauss2_mean1.getVal();val[1] = bkgGauss2_mean1.getError();
		writeParam(iBin, "bkgGauss2_mean1", val);
		val[0] = bkgGauss2_mean2.getVal();val[1] = bkgGauss2_mean2.getError();
		writeParam(iBin, "bkgGauss2_mean2", val);
		val[0] = bkgGauss2_sigma1.getVal();val[1] = bkgGauss2_sigma1.getError();
		writeParam(iBin, "bkgGauss2_sigma1", val);
		val[0] = bkgGauss2_sigma2.getVal();val[1] = bkgGauss2_sigma2.getError();
		writeParam(iBin, "bkgGauss2_sigma2", val);
		val[0] = bkgM_frac2.getVal();val[1] = bkgM_frac2.getError();
		writeParam(iBin, "bkgM_frac2", val);
		val[0] = bkgM_frac12.getVal();val[1] = bkgM_frac12.getError();
		writeParam(iBin, "bkgM_frac12", val);
		if (is7TeVCheck){
			switch (iBin) {
				case 2:
				   val[0]=470;val[1]=11;
					break;
				case 4:
				   val[0]=155;val[1]=6.9;
					break;
				case 6:
				   val[0]=6.8;val[1]=1.4;
					break;
				default:
				   val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
			}
		}else{
			switch (iBin) {
				case 2:
				   val[0]=1410;val[1]=33;
					break;
				case 4:
				   val[0]=465;val[1]=20.7;
					break;
				case 6:
				   val[0]=20.4;val[1]=4.2;
					break;
				default:
				   val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
			}
			val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
		}
		writeParam(iBin, "nbkgPeak", val);
	}
}//}}}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 13-10-2014

//void angular2D_2a_PkPl(int iBin, const char outfile[] = "angular2D_2a_PkPl", bool keepParam = false)
void angular2D_2a_PkPl(int iBin, const char outfile[] = "angular2D_2a_PkPl", bool keepParam = true)
{//{{{
	// Gaussian constraint on yields and mass is needed.
//	if (iBin ==0 || iBin%2 == 1 || iBin/2 > 3) {    /////   Peaking bkg. CosThetaL...
	if (iBin != 10 && (iBin ==0 || iBin%2 == 1 || iBin/2 > 3)){    ///////////   Peaking bkg. Only for bin 2,4,6,9,10
//	if (iBin != 9 && iBin != 10 && (iBin ==0 || iBin%2 == 1 || iBin/2 > 3)){    ///////////   Peaking bkg. Only for bin 2,4,6,9,10
		// Pm is fhat(and the yield is 0) for bins other than 2,4,6
		if (keepParam){
			double val[3]={0,0,0};
			writeParam(iBin, "bkgPeakL_c1", val);
			writeParam(iBin, "bkgPeakL_c2", val);
			writeParam(iBin, "bkgPeakL_c3", val);
			writeParam(iBin, "bkgPeakL_c4", val);
		}
		return;
	}
	
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
	RooArgSet f_bkgPeakL_argset;
	RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",0,-5,5);
	RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",0,-5,5);
	RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",0,-5,5);
	RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",0,-5,5);
	switch (iBin) {
		case 2:
		   //1 double guassian ,4+4 deg. ploy
			f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
			break;
		case 4:
		   //2 double guassian ,4+4 deg. ploy
			f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
			break;
		case 6:
		   //1 guassian ,2+2 deg. ploy
			//f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2));
			f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
			break;
	/*	case 9:
		   //1 double guassian ,4+4 deg. ploy
			f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
			break;
	*/	case 10:
		   //1 double guassian ,4+4 deg. ploy
			f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
			break;
		default:
		   break;
	}
	RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",50,0.,1E4);
	RooExtendPdf f_bkgPeakL_ext("f_bkgPeakL_ext","f_bkgPeakL_ext",f_bkgPeakL,nbkgPeak);
	
	// Gaussian Constraint
	RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,"nbkgPeak", 0)),RooConst(readParam(iBin, "nbkgPeak", 1)));
	// Get data
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);
	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*data,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));
	
	// Draw CosThetaL
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = 0.05;
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,Binning(20)); 
	f_bkgPeakL.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->Draw();
	framecosl->GetYaxis()->SetLabelFont(22);
	framecosl->GetYaxis()->SetLabelSize(0.04);
	framecosl->GetYaxis()->SetTitleSize(0.04);
	framecosl->GetYaxis()->SetTitleOffset(1.2);
	framecosl->GetYaxis()->SetTitleFont(22);
	framecosl->GetXaxis()->SetLabelFont(22);
	framecosl->GetXaxis()->SetLabelSize(0.04);
	framecosl->GetXaxis()->SetTitleSize(0.04);
	framecosl->GetXaxis()->SetTitleOffset(1.15);
	framecosl->GetXaxis()->SetTitleFont(22);
	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
	c->Update();
//	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
	
	// clear
	delete t1;
	delete c;
	delete data;
	
	if (keepParam){
		double val[3] = {0,0,0};
		val[0] = bkgPeakL_c1.getVal();val[1] = bkgPeakL_c1.getError();
		writeParam(iBin, "bkgPeakL_c1", val);
		val[0] = bkgPeakL_c2.getVal();val[1] = bkgPeakL_c2.getError();
		writeParam(iBin, "bkgPeakL_c2", val);
		val[0] = bkgPeakL_c3.getVal();val[1] = bkgPeakL_c3.getError();
		writeParam(iBin, "bkgPeakL_c3", val);
		val[0] = bkgPeakL_c4.getVal();val[1] = bkgPeakL_c4.getError();
		writeParam(iBin, "bkgPeakL_c4", val);
	}
}//}}}

void angular2D_prior(int iBin, const char outfile[] = "angular2D_prior", bool keepParam = true)
{//{{{
	// Fit to signal simulation by YsSm+YcCm to determine Sm
	RooRealVar Bmass("Bmass","M_{K^{+/-}#Mu#Mu}",5.,5.56);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
	
	// Create combinatorial background distribution
	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",0.,-3.,3.);
	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",0.,-3.,3.);
	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",0.,-3.,3.);
	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",0.,-3.,3.);
	RooRealVar bkgCombL_c4("bkgCombL_c4","c4",0.,-3.,3.);
	RooRealVar bkgCombL_c5("bkgCombL_c5","c5",0.,-3.,3.);
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
     /*   case 7:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            bkgCombL_c5.setConstant(kTRUE);
            break;
        case 0:
        case 1:
        case 4:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2));
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            bkgCombL_c5.setConstant(kTRUE);
            break;
        case 2:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3));
            bkgCombL_c4.setConstant(kTRUE);
            bkgCombL_c5.setConstant(kTRUE);
            break;
       
        case 0:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            bkgCombL_c5.setConstant(kTRUE);
            break;
      */  
	   default:
		f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4,bkgCombL_c5));
	//	f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4,bkgCombL_c5));
	//	f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
		f_bkgCombL_argset.add(RooArgSet(bkgCombL_c0));
		break;
	}
	RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
	
	// Get data and apply unbinned fit
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.38 || Bmass < 5.18)",Q2range[iBin]),0);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.474 || Bmass < 5.184)",Q2range[iBin]),0);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.384 || Bmass < 5.174)",Q2range[iBin]),0);
	RooFitResult *f_fitresult = f_bkgCombL.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));
	
	// Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = 0.05;
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,Binning(20)); 
	f_bkgCombL.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->Draw();
	framecosl->GetYaxis()->SetLabelFont(22);
	framecosl->GetYaxis()->SetLabelSize(0.04);
	framecosl->GetYaxis()->SetTitleSize(0.04);
	framecosl->GetYaxis()->SetTitleOffset(1.2);
	framecosl->GetYaxis()->SetTitleFont(22);
	framecosl->GetXaxis()->SetLabelFont(22);
	framecosl->GetXaxis()->SetLabelSize(0.04);
	framecosl->GetXaxis()->SetTitleSize(0.04);
	framecosl->GetXaxis()->SetTitleOffset(1.15);
	framecosl->GetXaxis()->SetTitleFont(22);
	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
	c->Update();
//	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
	
	// clear
	delete t1;
	delete c;
	delete data;
	// Prepare datacard
	if (keepParam){
		double val[3] = {0,0,0};
		val[0] = bkgCombL_c0.getVal();val[1] = bkgCombL_c0.getError();
		writeParam(iBin, "bkgCombL_c0", val);
		val[0] = bkgCombL_c1.getVal();val[1] = bkgCombL_c1.getError();
		writeParam(iBin, "bkgCombL_c1", val);
		val[0] = bkgCombL_c2.getVal();val[1] = bkgCombL_c2.getError();
		writeParam(iBin, "bkgCombL_c2", val);
		val[0] = bkgCombL_c3.getVal();val[1] = bkgCombL_c3.getError();
		writeParam(iBin, "bkgCombL_c3", val);
		val[0] = bkgCombL_c4.getVal();val[1] = bkgCombL_c4.getError();
		writeParam(iBin, "bkgCombL_c4", val);
		val[0] = bkgCombL_c5.getVal();val[1] = bkgCombL_c5.getError();
		writeParam(iBin, "bkgCombL_c5", val);
	}
}//}}}

//void angular2D_bin(int iBin, const char outfile[] = "angular2D")
std::vector<double> angular2D_bin(int iBin, const char outfile[] = "angular2D")
//std::vector<double> angular_reco_bin(int iBin, const char outfile[] = "angular_reco")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    // Need some modification for accXrecoEff.
    
    cout<<"iBin = "<<iBin<<endl<<endl; //////////////////////////////// 12-26
	double up, dn;
	if (iBin == 0) {up = 0.80; dn = -0.80;} // 19-09
//	else if (iBin == 1) {up = 0.95; dn = -0.85;} // 19-09
	else if (iBin == 9) { up = 0.95; dn = -0.95;}
	else if (iBin == 1) { up = 0.85; dn = -0.85;}
//	else if (iBin == 9) { up = 0.89; dn = -0.89;}
	else { up = 1.; dn = -1.;}
	 // Read data
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", dn, up);
    //RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar Bmass("Bmass","M_{K^{+/-}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",1.0,22.);

    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{+/-}#Mu#Mu}",5.28,5.26,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0));
    sigGauss1_sigma.setError(readParam(iBin,"sigGauss1_sigma",1));
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0));
    sigGauss2_sigma.setError(readParam(iBin,"sigGauss2_sigma",1));
    RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0));
    sigM_frac.setError(readParam(iBin,"sigM_frac",1));
        // Angular parameters
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fh("fh", "F_{H}", 0.2, 0., 1.);
    if (iBin != 3 && iBin != 5){
        // 2011 cross check

        // read parameter from datacard
    }
//	Efficiency
	RooRealVar effP0("effP0","effP0",readParam(iBin,"accXrecoEff", 0));
	RooRealVar effP1("effP1","effP1",readParam(iBin,"accXrecoEff", 1));
	RooRealVar effP2("effP2","effP2",readParam(iBin,"accXrecoEff", 2));
	RooRealVar effP3("effP3","effP3",readParam(iBin,"accXrecoEff", 3));   // f1_model_format_1
	RooRealVar effP4("effP4","effP4",readParam(iBin,"accXrecoEff", 4));   // f1_model_format_2  // f1_model_format_4
	RooRealVar effP5("effP5","effP5",readParam(iBin,"accXrecoEff", 5));   // f1_model_format_2  // f1_model_format_5  
	RooRealVar effP6("effP6","effP6",readParam(iBin,"accXrecoEff", 6));   // f1_model_format_2  // f1_model_format_6
	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff", 7));   // f1_model_format_3  // f1_model_format_7
	RooRealVar effP8("effP8","effP8",readParam(iBin,"accXrecoEff", 8));   // f1_model_format_3  // f1_model_format_8
	RooRealVar effP9("effP9","effP9",readParam(iBin,"accXrecoEff", 9));   // f1_model_format_3  // f1_model_format_9
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  // f1_model_format_1
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  // f1_model_format_2  // f1_model_format_4
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5));  // f1_model_format_2  // f1_model_format_5 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6));  // f1_model_format_2  // f1_model_format_6
	effP7.setError(readParam(iBin,"accXrecoEffErr", 7));  // f1_model_format_3  // f1_model_format_7
	effP8.setError(readParam(iBin,"accXrecoEffErr", 8));  // f1_model_format_3  // f1_model_format_8
	effP9.setError(readParam(iBin,"accXrecoEffErr", 9));  // f1_model_format_3  // f1_model_format_9
        // Efficiency and acceptance
    RooArgSet f_sigA_argset(CosThetaL);
    f_sigA_argset.add(RooArgSet(fh,afb));
	f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));       // f1_model_format_3 // f1_model_format_7
	f_sigA_argset.add(RooArgSet( effP8, effP9));       // f1_model_format_3 // f1_model_format_7
    TString f_sigA_format;
    TString f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
//    TString f_ang_format = "0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL";
//	 TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*fh*CosThetaK**2*(1-CosThetaL**2)+1/2*(1-fh)*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*afb*(1-CosThetaK**2)*CosThetaL))";
    TString f_rec_format; // f_rec_L0, f_rec_L2, f_rec_L3, f_rec_L4, f_rec_L6;
    if (is7TeVCheck){
    //    f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*CosThetaK**2+recK3L0*CosThetaK**3)";
    //    f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*CosThetaK**2+recK3L2*CosThetaK**3)*CosThetaL**2";
    //    f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*CosThetaK**2+recK3L3*CosThetaK**3)*CosThetaL**3";
    //    f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*CosThetaK**2+recK3L4*CosThetaK**3)*CosThetaL**4";
    //    f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*CosThetaK**2+recK3L6*CosThetaK**3)*CosThetaL**6";
    }else if (iBin == 0 || iBin == 1 || iBin == 9){
		 f_rec_format = "(effP0*exp(-0.5*((CosThetaL-effP1)/effP2)**2) ) * ( effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**4+effP8*CosThetaL**5+effP9*CosThetaL**6 )";
		 //f_rec_format = "effP0+effP1*CosThetaL+effP2*(3*CosThetaL**2-1)/2+effP3*(5*CosThetaL**3-3*CosThetaL)/2+effP4*(35*CosThetaL**4-30*CosThetaL**2+3)/8.+effP5*(63*CosThetaL**5-70*CosThetaL**3+15*CosThetaL)/8.+effP6*(231*CosThetaL**6-315*CosThetaL**4+105*CosThetaL**2-5)/16.+effP7*(429*CosThetaL**7-693*CosThetaL**5+315*CosThetaL**3-35*CosThetaL)/16.";
    }else if (iBin == 20) {
		 f_rec_format = "(( effP0*exp(-0.5*((CosThetaL-effP1)/effP2)**2) ) + ( effP3*exp(-0.5*((CosThetaL-effP4)/effP5)**2)) + ( effP6*exp(-0.5*((CosThetaL-effP7)/effP8)**2) ) + effP9)";
	 }else if (iBin !=0 && iBin != 1 && iBin != 9) {
		 f_rec_format = "(effP0*exp(-0.5*((CosThetaL-effP1)/effP2)**2) ) * ( effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**4+effP8*CosThetaL**5+effP9*CosThetaL**6 )";
	 }

	 ///////////////////    02-11-2014   N.A.
   // if (iBin == 0) {
   //     f_sigA_argset.add(RooArgSet(f_rec_format));
   //     f_sigA_format = TString::Format("(%s+%s+%s+%s)*%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data(),f_ang_format.Data());
   // }else if (iBin == 1) {
   //     f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
   //     f_sigA_format = TString::Format("(%s+%s+%s)*%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_ang_format.Data());
   // }else if (iBin > 1 && iBin < 6) {
   //     f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
   //     f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
   //     f_sigA_format = TString::Format("(%s+%s+%s+%s)*%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data(),f_ang_format.Data());
   //}else{
        f_sigA_argset.add(RooArgSet(f_rec_format));
        f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
   // }
        // angular map of signal
    RooGenericPdf f_sigA("f_sigA", f_sigA_format,f_sigA_argset);

    // Create signal distribution
        // mass distro of signal
    RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
    RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
    RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
    RooProdPdf f_sig("f_sig","f_sig",f_sigM,f_sigA);
    printf("INFO: f_sig prepared.\n");


    // Create combinatorial background distribution
    RooRealVar bkgCombM_c("bkgCombM_c","c",-0.1,-10.,1.);
    RooRealVar offset("offset","offset",-5.);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
    RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
//    RooRealVar     a0("a0", "constant", -0.5, -1, 1);
//	 RooRealVar     a1("a1", "linear", 0., -1, 1);
//	 RooRealVar     a2("a2", "linear", -0.1, -1, 1);
//	 RooChebychev   f_bkgCombM("f_bkgCombM", "Background", Bmass, RooArgSet(a0, a1, a2));

	 RooRealVar bkgCombL_c0("bkgCombL_c0","c0",readParam(iBin,"bkgCombL_c0",0),-3.,3.);
	 RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0),-3.,3.);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0),-3.,3.);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0),-3.,3.);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0),-3.,3.);
    RooRealVar bkgCombL_c5("bkgCombL_c5","c5",readParam(iBin,"bkgCombL_c5",0),-3.,3.);
    RooArgSet f_bkgCombL_argset;
    switch (iBin) {
    /*    case 0:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            bkgCombL_c5.setVal(0.);
            bkgCombL_c5.setConstant(kTRUE);
            break;
    */    default:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4,bkgCombL_c5));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c0));
            bkgCombL_c0.setConstant(kTRUE);
            bkgCombL_c1.setConstant(kTRUE);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            bkgCombL_c5.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
    RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombL,f_bkgCombM);
    printf("INFO: f_bkgComb prepared.\n");


    // Create peak background distribution
    RooRealVar bkgGauss1_mean1("bkgGauss1_mean1","M_{K+#Mu#Mu}",readParam(iBin,"bkgGauss1_mean1",0));
    bkgGauss1_mean1.setError(readParam(iBin,"bkgGauss1_mean1",1));
    RooRealVar bkgGauss1_mean2("bkgGauss1_mean2","M_{K+#Mu#Mu}",readParam(iBin,"bkgGauss1_mean2",0));
    bkgGauss1_mean2.setError(readParam(iBin,"bkgGauss1_mean2",1));
    RooRealVar bkgGauss1_sigma1("bkgGauss1_sigma1","#sigma_{11}",readParam(iBin,"bkgGauss1_sigma1",0));
    bkgGauss1_sigma1.setError(readParam(iBin,"bkgGauss1_sigma1",1));
    RooRealVar bkgGauss1_sigma2("bkgGauss1_sigma2","#sigma_{12}",readParam(iBin,"bkgGauss1_sigma2",0));
    bkgGauss1_sigma2.setError(readParam(iBin,"bkgGauss1_sigma2",1));
    RooRealVar bkgM_frac1("bkgM_frac1","bkgM_frac1",readParam(iBin,"bkgM_frac1",0));
    bkgM_frac1.setError(readParam(iBin,"bkgM_frac1",1));

    RooRealVar bkgGauss2_mean1("bkgGauss2_mean1","M_{K+#Mu#Mu}",readParam(iBin,"bkgGauss2_mean1",0));
    bkgGauss2_mean1.setError(readParam(iBin,"bkgGauss2_mean1",1));
    RooRealVar bkgGauss2_mean2("bkgGauss2_mean2","M_{K+#Mu#Mu}",readParam(iBin,"bkgGauss2_mean2",0));
    bkgGauss2_mean2.setError(readParam(iBin,"bkgGauss2_mean2",1));
    RooRealVar bkgGauss2_sigma1("bkgGauss2_sigma1","#sigma_{21}",readParam(iBin,"bkgGauss2_sigma1",0));
    bkgGauss2_sigma1.setError(readParam(iBin,"bkgGauss2_sigma1",1));
    RooRealVar bkgGauss2_sigma2("bkgGauss2_sigma2","#sigma_{22}",readParam(iBin,"bkgGauss2_sigma2",0));
    bkgGauss2_sigma2.setError(readParam(iBin,"bkgGauss2_sigma2",1));
    RooRealVar bkgM_frac2("bkgM_frac2","bkgM_frac2",readParam(iBin,"bkgM_frac2",0));
    bkgM_frac2.setError(readParam(iBin,"bkgM_frac2",1));

    RooRealVar bkgM_frac12("bkgM_frac12","bkgM_frac12",readParam(iBin,"bkgM_frac12",0));
    bkgM_frac12.setError(readParam(iBin,"bkgM_frac12",1));
    
	 RooGaussian f_bkgPeakMGauss11("f_bkgPeakMGauss11","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean1, bkgGauss1_sigma1);
    RooGaussian f_bkgPeakMGauss12("f_bkgPeakMGauss12","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean2, bkgGauss1_sigma2);
    RooGaussian f_bkgPeakMGauss21("f_bkgPeakMGauss21","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean1, bkgGauss2_sigma1);
    RooGaussian f_bkgPeakMGauss22("f_bkgPeakMGauss22","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean2, bkgGauss2_sigma2);
    RooAddPdf f_bkgPeakM1("f_bkgPeakM1","f_bkgPeakM1", RooArgList(f_bkgPeakMGauss11, f_bkgPeakMGauss12), bkgM_frac1);
    RooAddPdf f_bkgPeakM2("f_bkgPeakM2","f_bkgPeakM2", RooArgList(f_bkgPeakMGauss21, f_bkgPeakMGauss22), bkgM_frac2);
    RooAddPdf f_bkgPeakM12("f_bkgPeakM12","f_bkgPeakM12", RooArgList(f_bkgPeakM1, f_bkgPeakM2), bkgM_frac12);
    
	 RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",readParam(iBin,"bkgPeakL_c1",0));
    bkgPeakL_c1.setError(readParam(iBin,"bkgPeakL_c1",1));
    RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",readParam(iBin,"bkgPeakL_c2",0));
    bkgPeakL_c2.setError(readParam(iBin,"bkgPeakL_c2",1));
    RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",readParam(iBin,"bkgPeakL_c3",0));
    bkgPeakL_c3.setError(readParam(iBin,"bkgPeakL_c3",1));
    RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",readParam(iBin,"bkgPeakL_c4",0));
    bkgPeakL_c4.setError(readParam(iBin,"bkgPeakL_c4",1));
    RooArgSet f_bkgPeakL_argset(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4);
    switch (iBin) {// Should be fixed constants already.
        case 2:
            //1 double guassian ,4+4 deg. ploy
            bkgM_frac12.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            //bkgM_frac12.setConstant(kTRUE);
            //bkgM_frac2.setConstant(kTRUE);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            //bkgPeakL_c3.setConstant(kTRUE);
            //bkgPeakL_c4.setConstant(kTRUE);
            bkgM_frac12.setVal(1.);
            //bkgM_frac1.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            //bkgM_frac2.setConstant(kTRUE);
            break;
     /*   case 9:
            //1 double guassian ,4+4 deg. ploy
            bkgM_frac12.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            bkgM_frac12.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
        bkgPeakL_c1.setMin(NULL,0.);
        bkgPeakL_c1.setVal(0.);
        bkgPeakL_c2.setMin(NULL,0.);
        bkgPeakL_c2.setVal(0.);
        bkgPeakL_c3.setMin(NULL,0.);
        bkgPeakL_c3.setVal(0.);
        bkgPeakL_c4.setMin(NULL,0.);
        bkgPeakL_c4.setVal(0.);
				break;
     */   case 10:
            //2 double guassian ,4+4 deg. ploy
            //bkgM_frac12.setConstant(kTRUE);
            //bkgM_frac2.setConstant(kTRUE);
            break;
        default:
            bkgPeakL_c1.setConstant(kTRUE);
            bkgPeakL_c2.setConstant(kTRUE);
            bkgPeakL_c3.setConstant(kTRUE);
            bkgPeakL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
    RooProdPdf f_bkgPeak("f_bkgPeak", "f_bkgPeak",f_bkgPeakL,f_bkgPeakM12);
    printf("INFO: f_bkgPeak prepared.\n");


    // Observed spectrum = model*fullEfficiency
    RooRealVar nsig("nsig","nsig",10,0,5E3);
    RooRealVar nbkgComb("nbkgComb","nbkgComb",20,0,1E4);
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",readParam(iBin,"nbkgPeak",0) / Lumi_Scale);   //////////////////////////// 12-30
    //nbkgPeak.setError(readParam(iBin,"nbkgPeak",1) / Lumi_Scale);
	 nbkgPeak.setConstant(kTRUE);
    //RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",50,0.,1E4);
    //if (iBin !=10 && ( iBin == 0 || iBin %2 == 1 || iBin/2 > 3) ){
    //if (iBin !=9 && iBin !=10 && ( iBin == 0 || iBin %2 == 1 || iBin/2 > 3) ){
    if (iBin !=10 && ( iBin == 0 || iBin %2 == 1 || iBin/2 > 3) ){
        nbkgPeak.setMin(NULL,0.);
        nbkgPeak.setVal(0.);
        nbkgPeak.setConstant(kTRUE);
    }
    //RooAddPdf kernel("kernel","kernel",RooArgList(f_sig,f_bkgComb,f_bkgPeak),RooArgList(nsig,nbkgComb,nbkgPeak));
    RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_bkgPeak,f_sig),RooArgList(nbkgComb,nbkgPeak,nsig));// no penalty term

    // Extra penalty term to confine As, Fs, Fh, Afb.
    //RooRealVar t_penalty("t_penalty","t",0.01);
    //RooGenericPdf f_penaltyAfb("f_penaltyAfb","(1-TMath::Erf((afb-0.75*(1-fh))/(1.5*t_penalty*(1-fh))))*(1-TMath::Erf((-afb-0.75*(1-fh))/(1.5*t_penalty*(1-fh))))",RooArgSet(afb,fh,t_penalty));
    //RooGenericPdf f_penaltyAs("f_penaltyAfb","(1-TMath::Erf((afb-2*(1-fh)/3)/(1.5*t_penalty*(1-fh))))*(1-TMath::Erf((-afb-0.75*(1-fh))/(1.5*t_penalty*(1-fh))))",RooArgSet(afb,fh,t_penalty));
    //RooProdPdf f_penalty("f_penalty","f_penalty",f_penaltyAfb,f_penaltyAs);
    //RooProdPdf f("f","f",f_model,f_penalty);
    printf("INFO: f_penalty NOT prepared.\n");

    // Gaussian constraints
    RooGaussian gaus_sigGauss1_sigma("gaus_sigGauss1_sigma","gaus_sigGauss1_sigma",sigGauss1_sigma,RooConst(readParam(iBin,"sigGauss1_sigma",0)),RooConst(readParam(iBin,"sigGauss1_sigma",1)));
    RooGaussian gaus_sigGauss2_sigma("gaus_sigGauss2_sigma","gaus_sigGauss2_sigma",sigGauss2_sigma,RooConst(readParam(iBin,"sigGauss2_sigma",0)),RooConst(readParam(iBin,"sigGauss2_sigma",1)));
    RooGaussian gaus_sigM_frac("gaus_sigM_frac","gaus_sigM_frac",sigM_frac,RooConst(readParam(iBin,"sigM_frac",0)),RooConst(readParam(iBin,"sigM_frac",1)));

   // bkgPeak    //////  Lumi. Scale
	 //RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,"nbkgPeak",0) ),RooConst(readParam(iBin,"nbkgPeak",1) ));  
	 RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,"nbkgPeak",0)/Lumi_Scale ),RooConst(readParam(iBin,"nbkgPeak",1)/Lumi_Scale ));  
    
	 RooGaussian gaus_bkgGauss1_mean1("gaus_bkgGauss1_mean1","gaus_bkgGauss1_mean1",bkgGauss1_mean1,RooConst(readParam(iBin,"bkgGauss1_mean1",0)),RooConst(readParam(iBin,"bkgGauss1_mean1",1)));
    RooGaussian gaus_bkgGauss1_mean2("gaus_bkgGauss1_mean2","gaus_bkgGauss1_mean2",bkgGauss1_mean2,RooConst(readParam(iBin,"bkgGauss1_mean2",0)),RooConst(readParam(iBin,"bkgGauss1_mean2",1)));
    RooGaussian gaus_bkgGauss1_sigma1("gaus_bkgGauss1_sigma1","gaus_bkgGauss1_sigma1",bkgGauss1_sigma1,RooConst(readParam(iBin,"bkgGauss1_sigma1",0)),RooConst(readParam(iBin,"bkgGauss1_sigma1",1)));
    RooGaussian gaus_bkgGauss1_sigma2("gaus_bkgGauss1_sigma2","gaus_bkgGauss1_sigma2",bkgGauss1_sigma2,RooConst(readParam(iBin,"bkgGauss1_sigma2",0)),RooConst(readParam(iBin,"bkgGauss1_sigma2",1)));
    RooGaussian gaus_bkgM_frac1("gaus_bkgM_frac1","gaus_bkgM_frac1",bkgM_frac1,RooConst(readParam(iBin,"bkgM_frac1",0)),RooConst(readParam(iBin,"bkgM_frac1",1)));
    RooGaussian gaus_bkgGauss2_mean1("gaus_bkgGauss2_mean1","gaus_bkgGauss2_mean1",bkgGauss2_mean1,RooConst(readParam(iBin,"bkgGauss2_mean1",0)),RooConst(readParam(iBin,"bkgGauss2_mean1",1)));
    RooGaussian gaus_bkgGauss2_mean2("gaus_bkgGauss2_mean2","gaus_bkgGauss2_mean2",bkgGauss2_mean2,RooConst(readParam(iBin,"bkgGauss2_mean2",0)),RooConst(readParam(iBin,"bkgGauss2_mean2",1)));
    RooGaussian gaus_bkgGauss2_sigma1("gaus_bkgGauss2_sigma1","gaus_bkgGauss2_sigma1",bkgGauss2_sigma1,RooConst(readParam(iBin,"bkgGauss2_sigma1",0)),RooConst(readParam(iBin,"bkgGauss2_sigma1",1)));
    RooGaussian gaus_bkgGauss2_sigma2("gaus_bkgGauss2_sigma2","gaus_bkgGauss2_sigma2",bkgGauss2_sigma2,RooConst(readParam(iBin,"bkgGauss2_sigma2",0)),RooConst(readParam(iBin,"bkgGauss2_sigma2",1)));
    RooGaussian gaus_bkgM_frac2("gaus_bkgM_frac2","gaus_bkgM_frac2",bkgM_frac2,RooConst(readParam(iBin,"bkgM_frac2",0)),RooConst(readParam(iBin,"bkgM_frac2",1)));
    RooGaussian gaus_bkgM_frac12("gaus_bkgM_frac12","gaus_bkgM_frac12",bkgM_frac12,RooConst(readParam(iBin,"bkgM_frac12",0)),RooConst(readParam(iBin,"bkgM_frac12",1)));
    
    RooGaussian gaus_bkgPeakL_c1("gaus_bkgPeakL_c1","gaus_bkgPeakL_c1",bkgPeakL_c1,RooConst(readParam(iBin,"bkgPeakL_c1",0)),RooConst(readParam(iBin,"bkgPeakL_c1",1)));
    RooGaussian gaus_bkgPeakL_c2("gaus_bkgPeakL_c2","gaus_bkgPeakL_c2",bkgPeakL_c2,RooConst(readParam(iBin,"bkgPeakL_c2",0)),RooConst(readParam(iBin,"bkgPeakL_c2",1)));
    RooGaussian gaus_bkgPeakL_c3("gaus_bkgPeakL_c3","gaus_bkgPeakL_c3",bkgPeakL_c3,RooConst(readParam(iBin,"bkgPeakL_c3",0)),RooConst(readParam(iBin,"bkgPeakL_c3",1)));
    RooGaussian gaus_bkgPeakL_c4("gaus_bkgPeakL_c4","gaus_bkgPeakL_c4",bkgPeakL_c4,RooConst(readParam(iBin,"bkgPeakL_c4",0)),RooConst(readParam(iBin,"bkgPeakL_c4",1)));
    
    RooArgSet gausConstraints(gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac);
    switch (iBin) {
        case 2:
            //1 double guassian ,4+4 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            gausConstraints.add(gaus_nbkgPeak);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_mean2,gaus_bkgGauss2_sigma1,gaus_bkgGauss2_sigma2,gaus_bkgM_frac2));
            gausConstraints.add(gaus_bkgM_frac12);
            gausConstraints.add(gaus_nbkgPeak);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            //gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2));
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_mean2,gaus_bkgGauss2_sigma1,gaus_bkgGauss2_sigma2,gaus_bkgM_frac2));
            //gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_sigma1));
            gausConstraints.add(gaus_nbkgPeak);
            break;
     /*   case 9:
            //1 double guassian ,4+4 deg. ploy
            //gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            //gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            //gausConstraints.add(gaus_nbkgPeak);
            break;
     */   case 10:
            //2 double guassian ,4+4 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_mean2,gaus_bkgGauss2_sigma1,gaus_bkgGauss2_sigma2,gaus_bkgM_frac2));
            gausConstraints.add(gaus_bkgM_frac12);
            gausConstraints.add(gaus_nbkgPeak);
            break;
    }
    
    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaL,Q2),Q2range[iBin],0);
    f.fitTo(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framemass = Bmass.frame();
	 data->plotOn(framemass,RooFit::Name("data"), Binning(20)); 
    f.plotOn(framemass,RooFit::Name("pdf"), LineColor(1)); 
    //f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
    f.plotOn(framemass, Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
    f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(4),LineWidth(2));
    f.plotOn(framemass,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    //f.plotOn(framemass,RooFit::Name("bkgPeak"), Components(f_bkgPeak),FillStyle(3004),FillColor(3),VLines(), DrawOption("F"));
    f.plotOn(framemass, Components(f_bkgPeak),FillStyle(3004),FillColor(3),VLines(), DrawOption("F"));
    f.plotOn(framemass,RooFit::Name("bkgPeak"), Components(f_bkgPeak),LineColor(3),LineWidth(2),LineStyle(2));

    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();
	framemass->GetYaxis()->SetLabelFont(22);
	framemass->GetYaxis()->SetLabelSize(0.04);
	framemass->GetYaxis()->SetTitleSize(0.04);
	framemass->GetYaxis()->SetTitleOffset(1.2);
	framemass->GetYaxis()->SetTitleFont(22);
	framemass->GetXaxis()->SetLabelFont(22);
	framemass->GetXaxis()->SetLabelSize(0.04);
	framemass->GetXaxis()->SetTitleSize(0.04);
	framemass->GetXaxis()->SetTitleOffset(1.15);
	framemass->GetXaxis()->SetTitleFont(22);
    
	TLegend *leg =new TLegend(0.74,0.70,0.89,0.89,NULL,"brNDC");
	leg->AddEntry("data"," Data "," PE ");
	leg->AddEntry("pdf"," Total P.d.f. "," L ");
	leg->AddEntry("sig"," Signal "," L ");
	leg->AddEntry("bkgComb"," Comb. bkg. "," L");
	leg->AddEntry("bkgPeak"," Peak. bkg. ","L ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
	leg->Draw();

	TPaveText* paveText = new TPaveText( 0.13, 0.70, 0.41, 0.88, "NDC" ); 
//	TPaveText* paveText = new TPaveText( 0.17, 0.70, 0.41, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("F_{H}  =%9.5f #pm%9.5f ", fh.getVal(),fh.getError())); 
	paveText->AddText(Form("A_{FB} =%9.5f #pm%9.5f ", afb.getVal(),afb.getError())); 
	paveText->AddText(Form("   nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
	paveText->AddText(Form(" sigmean   = %.3f #pm %.3f ", sigGauss_mean.getVal(), sigGauss_mean.getError())); 
	paveText->AddText(Form("#Chi^{2}_{Bmass}   = %.2f  ", framemass->chiSquare())); 
	paveText->Draw(); 
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0.05;
	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Update();
//    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
    c->Print(TString::Format("./plots/%s_bin%d.png",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,RooFit::Name("data"), Binning(20)); 
    f.plotOn(framecosl,RooFit::Name("pdf"), LineColor(1)); 
    //f.plotOn(framecosl,Components(f_sig),LineColor(4),LineWidth(2));
    //f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
    f.plotOn(framecosl, Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
    f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(4),LineWidth(2));
    f.plotOn(framecosl,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    //f.plotOn(framecosl,RooFit::Name("bkgPeak"), Components(f_bkgPeak),FillStyle(3004),FillColor(3),VLines(), DrawOption("F"));
    f.plotOn(framecosl, Components(f_bkgPeak),FillStyle(3004),FillColor(3),VLines(), DrawOption("F"));
    f.plotOn(framecosl,RooFit::Name("bkgPeak"), Components(f_bkgPeak),LineColor(3),LineWidth(2),LineStyle(2));

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();
	framecosl->GetYaxis()->SetLabelFont(22);
	framecosl->GetYaxis()->SetLabelSize(0.04);
	framecosl->GetYaxis()->SetTitleSize(0.04);
	framecosl->GetYaxis()->SetTitleOffset(1.2);
	framecosl->GetYaxis()->SetTitleFont(22);
	framecosl->GetXaxis()->SetLabelFont(22);
	framecosl->GetXaxis()->SetLabelSize(0.04);
	framecosl->GetXaxis()->SetTitleSize(0.04);
	framecosl->GetXaxis()->SetTitleOffset(1.15);
	framecosl->GetXaxis()->SetTitleFont(22);

/*	TLegend *leg =new TLegend(0.74,0.70,0.89,0.89,NULL,"brNDC");
	leg->AddEntry("data"," Data "," P ");
	leg->AddEntry("pdf"," Fit Pdf "," L ");
	leg->AddEntry("sig"," Signal "," FL ");
	leg->AddEntry("bkgComb"," Comb. bkg. "," L");
	leg->AddEntry("bkgPeak"," Peak. bkg. ","FL ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
*/	leg->Draw();


	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
//	t1->DrawLatex(.10,.81+fixNDC,TString::Format("F_{H}  =%9.5f #pm%9.5f",fh.getVal(),fh.getError()));
//	t1->DrawLatex(.10,.76+fixNDC,TString::Format("A_{FB} =%9.5f #pm%9.5f",afb.getVal(),afb.getError()));
	paveText->AddText(Form("#Chi^{2}_{cos#theta}   = %.2f  ", framecosl->chiSquare())); 
	paveText->Draw(); 
    c->Update();
//    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
    c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
    delete c;
    delete t1;
    delete data;
    
	 // write output
    double val[3]={0,0,0};
    val[0] = fh.getVal();val[1] = fh.getError();
    writeParam(iBin, "fh", val);
    val[0] = afb.getVal();val[1] = afb.getError();
    writeParam(iBin, "afb",val);
    
	 std::vector<double> output;
    output.push_back(fh.getVal());
    output.push_back(fh.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    return output;
}//}}}

void angular2D( const char outfile[] = "angular2D")
{//{{{
//	doFit = false;
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yuerrfh[9], yderrfh[9], yafb[9], yuerrafb[9], yderrafb[9];
	for (int i = 0; i < 9; i++) {
		yfh[i]  = 0; yuerrfh[i]  = 0; yderrfh[i]  = 0;
		yafb[i] = 0; yuerrafb[i] = 0; yderrafb[i] = 0;
	}
// Checkout input data
	for(int i = 0, ibin = 0; i < 9, ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
		yafb[i]      = readParam(ibin,"afb",0);
		yuerrafb[i]  = fabs(readParam(ibin,"afb",1));
		yderrafb[i]  = yuerrafb[i];
		yfh[i]       = readParam(ibin,"fh",0);
		yuerrfh[i]   = fabs(readParam(ibin,"fh",1));
		if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		else { yderrfh[i] = yuerrfh[i]; }
		printf("Afb[%d]=%6.4f +- %6.4f\n",i,yafb[i],yuerrafb[i]);
		printf("Fh [%d]=%6.4f +- %6.4f\n",i,yfh[i],yuerrfh[i]);
	}
//	Draw
	TCanvas *c = new TCanvas("c");
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->GetYaxis()->SetLabelFont(22);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleFont(22);
	frame->GetXaxis()->SetLabelFont(22);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleOffset(1.15);
	frame->GetXaxis()->SetTitleFont(22);
	frame->SetTitle("");
	frame->Draw();
	
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,1.,"Y");
	TGraphAsymmErrors *d_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	d_fh->SetMarkerColor(1);
	d_fh->SetMarkerStyle(20);
//	d_fh->SetLineWidth(0.005);	
	d_fh->SetFillColor(2);
	d_fh->SetFillStyle(3001);
	d_fh->Draw("2");
	d_fh->Draw("P");
//	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Print(TString::Format("./plots/%s_fh.png",outfile));
	c->Clear();
	
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("A_{FB}");
	frame->SetAxisRange(-1.,1.,"Y");
	TGraphAsymmErrors *d_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	d_afb->SetMarkerColor(1);
	d_afb->SetMarkerStyle(20);
//	d_afb->SetLineWidth(0.005);	
	d_afb->SetFillColor(2);
	d_afb->SetFillStyle(3001);
	d_afb->Draw("2");
	d_afb->Draw("P");
//	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
	c->Print(TString::Format("./plots/%s_afb.png",outfile));
}//}}}

//////////////////////////////////////////////////////////////////////////
//_________________________________________________________________________________
//_________________________________________________________________________________
int main(int argc, char** argv) {
//	Tags
	is7TeVCheck = false;   
//	Help message
	if (argc <= 2) {
		printf("Usage       : ./fit Function infile binID\n");
		printf("Functions   :\n");
		printf("    0. bmass               Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg.\n");
		printf("    1. angular_gen         Derive F_{H} and A_{FB} from cosThetaL distribution at GEN level.\n");
		printf("       acceptance          Get acceptance map from unfiltered signal GEN, |Mu pT| > 2.8 GeV, |Mu eta| < 2.3.\n");
		printf("       recoEff             Get reconstruction efficiency map from signal simulation.\n");
		printf("    2. accXrecoEff         Get efficiency map from signal simulation.\n");
		printf("    3. angular_reco        Derive F_{H} and A_{FB} from cosThetaL distribution at RECO level.\n");
		printf("    4. angular2D_1a_Sm     Leading step1 to angular2D, determine signal shape from simulation.\n");
		printf("    5. angular2D_1b_YpPm   Leading step2 to angular2D, determine mass spectrum of peaking bkg from simulation.\n");
		printf("    6. angular2D_2a_PkPl   Leading step3 to angular2D, determine angular dist. of peaking bkg from simulation.\n");
		printf("    7. angular2D_prior     Leading step4 to angular2D, fit to data sideband to get initial values of combinatorial bkg.\n");
		printf("    8. angular2D           Derive F_{H} and A_{FB} by fitting to mass and angular distribution.\n");
		printf("Remark      :\n");
		printf("    1. Outputs will be stored in ./plots, please keep the directory.\n");
		printf("    2. Fitted parameters will be stored in ./fitParameters/*.txt, please keep the directory.\n");
	//	printf("    3. Wildcard is allowed for infile. But you must quote infile like \"inputData_Run*.root\"!\n");
		return 0;
	}
//	main
	if (argc != 4){
		printf("./fit func infile binID\n");
		for (int i = 0; i < 11; i++) {
		//	if (i == 3 || i == 5) continue;
			printf("    Bin %d : %s\n",i,Q2range[i]);
		}
		return 0;
	}
	TString func    = argv[1];
	TString infile  = argv[2];
	int iBin        = atoi(argv[3]);
	
	if (func == "bmass") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="bmass";
		bmass(iBin, outfile); 
	}else if (func == "angular_gen"){
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			std::vector<double> vbin;
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			vbin = angular_gen_bin(iBin);
		}else if (iBin == 999) {
			const char outfile[]="angular_gen";
			angular_gen(outfile);
		}else { 
			cout<<"Refit gen level,iBin counts from 0 to 10, or 999!"<<endl;
			return 0; 
		}
/*	}else if (func == "acceptance") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		for (int iBin = 0; iBin < 11; iBin++) {
			if (iBin == 3 || iBin == 5) continue;
			acceptance(iBin);
		}
	}else if (func == "recoEff") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		for (int iBin = 0; iBin < 11; iBin++) {
			if (iBin == 3 || iBin == 5) continue;
			recoEff(iBin);
		}
*/	}else if (func == "accXrecoEff") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			accXrecoEff(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	}else if (func == "angular_reco"){
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			std::vector<double> vbin;
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			vbin = angular_reco_bin(iBin);
		}else if (iBin == 999) {
			const char outfile[]="angular_reco";
			angular_reco( outfile);
		}else { 
			cout<<"Refit reco level,iBin counts from 0 to 10; or 999 to plot the results!"<<endl;
			return 0; 
		}
	}else if (func == "angular2D_1a_Sm" || func == "angular2D_1b_YpPm" || func == "angular2D_2a_PkPl" || func == "angular2D_prior"){
		void (*fx)(int, const char*, bool);
		if ( func == "angular2D_1a_Sm" ){
			fx = angular2D_1a_Sm;
			}else if (func == "angular2D_1b_YpPm"){
				fx = angular2D_1b_YpPm;
			}else if (func == "angular2D_2a_PkPl"){
				fx = angular2D_2a_PkPl;
			}else{
				fx = angular2D_prior;
			}
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			fx(iBin,func,true);// By default overwrite exist parameters.
	}else if (func == "angular2D"){
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			std::vector<double> vbin;
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			vbin = angular2D_bin(iBin);
		}else if (iBin == 999) {
			const char outfile[]="angular2D";
			angular2D(outfile);
		}else { 
			cout<<"Refit data, iBin counts from 0 to 10; or 999 to plot the results!"<<endl;
			return 0; 
		}
	} 
// 04-09-2014
/////////////////////////////////////////////////////////////////////////////////////////////
	 else if (func == "test"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
       // const char outfile[]="test";
        for (int iBin = 0; iBin < 8; iBin++) {
            //getToyFromUnfilterGen(iBin);
            //createRecoEffHist(iBin);
            //accXrecoEff(iBin);
            //angular2D_bin(iBin);
        }
        //createAccptanceHist();
        //angular2D_bin(0);
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    printf("%lld entries processed.\n",ch->GetEntries());
    gSystem->Exit(0);

    return 0 ;
}
