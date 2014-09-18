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

#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
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
using namespace RooFit;

// Tags configration
bool is7TeVCheck = false; // Using 2011 efficiency map.
TChain *ch=new TChain("tree");
//Constants, Fit results for efficiency, etc.. //{{{
char genQ2range[11][32] = {"genQ2 <=  2.00 && genQ2 >  0.05",
                           "genQ2 <=  4.30 && genQ2 >  2.00",
                           "genQ2 <=  8.68 && genQ2 >  4.30",
                           "genQ2 <= 10.09 && genQ2 >  8.68",
									"genQ2 <= 12.86 && genQ2 > 10.09",
                           "genQ2 <= 14.18 && genQ2 > 12.86",
									"genQ2 <= 16.00 && genQ2 > 14.18",
                           "genQ2 <= 18.00 && genQ2 > 16.00",
                           "genQ2 <= 22.00 && genQ2 > 18.00",
                           "genQ2 <=  6.00 && genQ2 >  1.00",
									"genQ2 <= 22.00 && genQ2 >  1.00"};
char Q2range[11][32] = {"Q2 <=  2.00 && Q2 >  0.05",
                        "Q2 <=  4.30 && Q2 >  2.00",
                        "Q2 <=  8.68 && Q2 >  4.30",
                        "Q2 <= 10.09 && Q2 > 8.68",
								"Q2 <= 12.86 && Q2 > 10.09",
                        "Q2 <= 14.18 && Q2 > 12.86",
								"Q2 <= 16.00 && Q2 > 14.18",
                        "Q2 <= 18.00 && Q2 > 16.00",
                        "Q2 <= 22.00 && Q2 > 18.00",
                        "Q2 <=  6.00 && Q2 >  1.00",
								"Q2 <= 22.00 && Q2 >  1.00"};
double Q2rangedn[11] = {0.05 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 1.00 ,  1.00};
double Q2rangeup[11] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 22.00 , 6.00 , 22.00};
char mumuMassWindow[5][300] = { 
	" Mumumass > 0 ",
	" (Mumumass > 3.096916+3.*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) && (Mumumass > 3.686109+3.*Mumumasserr || Mumumass < 3.686109-3.*Mumumasserr)", // No JPsi && Psi(2S)
//	" (Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.5*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)",
//	&& (fabs(Bmass - Mumumass - 2.182) > 0.14) && (fabs(Bmass - Mumumass - 1.593) > 0.09)",
	" (Mumumass < 3.096916+3.*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr) \
	  || (Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr)",  // Only JPsi && Psi(2S)
	" Mumumass < 3.096916+3*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr",  // Only JPsi
	" Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr"}; //Only Psi(2S)
double genAfb[9]={0.00, 0.07, -0.02, -0.03, -0.01, -0.09, 0.02, 0.02, -0.01};
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
	RooRealVar     Q2("Q2","q^{2}",0.5,20.);
	RooRealVar     Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
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
	model.plotOn(frame);
	
//	Overlay the background component of model with a dashed line
	model.plotOn(frame,Components("bkg"), LineStyle(kDashed), LineColor(2));
	
//	Draw the frame on the canvas
	TCanvas *c = new TCanvas("c", "c", 800, 600); 
	set_root_style(); 
	c->UseCurrentStyle();
	
	gPad->SetLeftMargin(0.15);
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->Draw();
	
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
	
	c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
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

//_________________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
std::vector<double> fh_gen_bin(int iBin, const char outfile[] = "fh_gen")
{//{{{
  // From fomula (2) in LHCb 2012 PRL108, 181806(2012)
  // integrated over theta_l and phi: 
  // 
  // 1/Gamma * d Gamma/d cos(theta_L) =
  // 3/4(1-F_H)(1-cos^2theta_L) + 1/2 * F_H + A_FB * cos(theta_L)
  // |A_FB| <= F_H / 2
  //

  bool test = false;

  RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{L}", -1, 1);
  RooRealVar genQ2("genQ2","q^{2}",0.05,22.);
  RooRealVar Q2("Q2","q^{2}",0.05,22.);
  RooRealVar fh("fh", "F_{H}", 0.02, 0, 0.5);
  RooRealVar afb("afb", "A_{FB}", -0.01, -0.2, 0.5);

  RooGenericPdf f("f", "0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL", RooArgSet(genCosThetaL,fh,afb));
  RooDataSet* data;
  
  if (test){
      fh.setVal(0.5);
      data = f.generate(RooArgSet(genCosThetaL,Q2), 10000);
  }else{
      data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaL,genQ2),genQ2range[iBin],0);
  }
  
  //f.fitTo(*data,Extended(kTRUE)); 
  f.fitTo(*data); 

  RooPlot* framecosl = genCosThetaL.frame(); 
  data->plotOn(framecosl); 
  f.plotOn(framecosl); 

  // Draw the frame on the canvas
  TCanvas *c = new TCanvas("c"); 
  framecosl->SetTitle("");
  framecosl->Draw();

  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->DrawLatex(.30,.75,TString::Format("%s",genQ2range[iBin]));
  t1->DrawLatex(.20,.69,TString::Format("F_{H}  =%8.4f#pm%8.4f",fh.getVal(),fh.getError()));
  t1->DrawLatex(.50,.69,TString::Format("A_{FB} =%8.4f#pm%8.4f",afb.getVal(),afb.getError()));

  c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
  c->Print(TString::Format("./plots/%s_bin%d.png",outfile,iBin));

  delete c;
  delete t1;
  delete data;

  std::vector<double> outvect;
  outvect.push_back(fh.getVal());
  outvect.push_back(fh.getError());
  outvect.push_back(afb.getVal());
  outvect.push_back(afb.getError());
  return outvect;
}//}}}

void fh_gen(const char outfile[] = "fh_gen")
{//{{{

    TCanvas *c = new TCanvas();
    TH2F *frame = new TH2F("frame","",22,0.,22,9,-0.1,0.8);
    frame->SetStats(kFALSE);
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{H}");
    frame->Draw();

    double x[11]   ={1.025, 3.15, 6.49, 9.385, 11.475, 13.52, 15.09, 17.0, 20.0, 3.5, 11.5};
    double xerr[11]={0.975, 1.15, 2.09, 0.705,  1.385,  0.66,  0.91,  1.0,  2.0, 2.5, 10.5};
    double yfh[11],yerrfh[11];

    std::vector<double> vbin;
    for(int ibin = 0; ibin < 11; ibin++){
        vbin = fh_gen_bin(ibin);
        yfh[ibin]       =vbin.at(0);
        yerrfh[ibin]    =vbin.at(1);
    }
    
    // Check input data
    for(int ibin = 0; ibin < 11; ibin++){
        printf("yfh [%d]=%10.6f +- %10.6f\n",ibin,yfh[ibin],yerrfh[ibin]);
    }

    TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(9,x,yfh,xerr,xerr,yerrfh,yerrfh);
    g_fh->Draw("P*");
    c->Print(TString::Format("./plots/%s.pdf",outfile));
    c->Print(TString::Format("./plots/%s.png",outfile));


    delete g_fh;
    delete frame;
    delete c;
    
	 TCanvas *c1 = new TCanvas();
    TH2F *frame1 = new TH2F("frame1","",22,0.,22,11,-0.3,0.8);
    frame1->SetStats(kFALSE);
    frame1->SetXTitle("q^{2} [(GeV)^{2}]");
    frame1->SetYTitle("A_{FB}");
    frame1->Draw();

//    double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
//    double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
    double yafb[11],yerrafb[11];

    std::vector<double> vbin1;
    for(int ibin = 0; ibin < 11; ibin++){
        vbin1 = fh_gen_bin(ibin);
        yafb[ibin]       =vbin1.at(2);
        yerrafb[ibin]    =vbin1.at(3);
    }
    
    // Check input data
    for(int ibin = 0; ibin < 11; ibin++){
        printf("yafb [%d]=%10.6f +- %10.6f    yfh [%d]=%10.6f +- %10.6f\n",
		         ibin,yafb[ibin],yerrafb[ibin],ibin,yfh[ibin],yerrfh[ibin]);
    }

    TGraphAsymmErrors *g_afb  = new TGraphAsymmErrors(9,x,yafb,xerr,xerr,yerrafb,yerrafb);
    g_afb->Draw("P*");
    c1->Print(TString::Format("./plots/afb_gen.pdf"));
    c1->Print(TString::Format("./plots/afb_gen.png"));

    delete g_afb;
    delete frame1;
    delete c1;
}//}}}
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//_________________________________________________________________________________

std::vector<double> angular_gen_bin(int iBin, const char outfile[] = "angular_gen")
{//{{{
	RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{L}", -1., 1.);
	RooRealVar genQ2("genQ2","q^{2}",0.05,22.);
	RooRealVar fh("fh", "F_{H}", genFh[iBin], 0.0, 0.9);
	RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -0.3, 0.5);
	
	RooRealVar nsig("nsig","nsig",1E6,1E2,1E9);
//	RooRealVar nbkg("nbkg","nbkg",10,0.1,1E4);
	
	RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL", RooArgSet(genCosThetaL,fh,afb));
	RooExtendPdf f("f","",f_sig,nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaL,genQ2),genQ2range[iBin],0);
	
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
	
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
	
	fixNDC = -0.5;
//	if (iBin > 4) fixNDC = 0.;
	t1->DrawLatex(.30,.75+fixNDC,TString::Format("%s",genQ2range[iBin]));
	t1->DrawLatex(.15,.69+fixNDC,TString::Format("F_{H}  =%8.5f#pm%8.5f",fh.getVal(),fh.getError()));
	t1->DrawLatex(.51,.69+fixNDC,TString::Format("A_{FB} =%8.5f#pm%8.5f",afb.getVal(),afb.getError()));
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
	
//	clear
	delete t1;
	delete c;
	delete data;
	
//	write output
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;

}//}}}

void angular_gen(const char outfile[] = "angular_gen")
{//{{{
//	bool refit = false; // Turn to true if you want to fit again.
	bool refit = true; // Turn to true if you want to fit again.
	
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
//	TH2F *frame = new TH2F("frame","",18,1,19,10,-0.5,1);
	frame->SetStats(kFALSE);
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.02,0.16,"Y");
	frame->Draw();
	
	double x[9]   ={1.25, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.75, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yerrfh[9], yafb[9], yerrafb[9]; 
/*	double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
 	double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
	double yfh[8]       ={0.705,0.791,0.649,0.524,0.454,0.399,0.369,0.341};
	double yerrfh[8]    ={0.000565,0.000410,0.000273,0.000428,0.000296,0.000413,0.000367,0.000359};
	double yafb[8]      ={-0.160,-0.066,0.182,0.317,0.374,0.412,0.421,0.376};
	double yerrafb[8]   ={0.000432,0.000284,0.000224,0.000390,0.000286,0.000419,0.000393,0.420};
*/
   if (refit){ // Turn to true if you want to fit again.
		std::vector<double> vbin;
		for(int ibin = 0, i = 0; i < 11; i++){
			if (i == 3 || i == 5) continue;
			vbin = angular_gen_bin(i);
			yfh[ibin]       =vbin.at(0);
			yerrfh[ibin]    =vbin.at(1);
			yafb[ibin]      =vbin.at(2);
			yerrafb[ibin]   =vbin.at(3);
			ibin++;
		}
	}
//	Check input data
	for(int ibin = 0, i = 0; i < 11; i++){
		if (i == 3 || i == 5) continue;
		printf("  yafb[%d]=%12.6f +- %12.6f     ",ibin,yafb[ibin],yerrafb[ibin]);
		printf("yfh[%d]=%12.6f +- %12.6f\n",ibin,yfh[ibin],yerrfh[ibin]);
		printf("genAfb[%d]=%12.6f      >>>>>>     ",ibin,genAfb[ibin]);
		printf("genFh[%d]=%12.6f\n",ibin,genFh[ibin]);
		ibin++;
	}

//	plotting
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yerrfh,yerrfh);
	g_fh->SetMarkerColor(4);
	g_fh->SetMarkerStyle(20);
	
	g_fh->SetFillColor(2);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Print(TString::Format("./plots/%s_fh.png",outfile));
	c->Clear();
	
	frame->SetTitle("");
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetAxisRange(-0.04,0.04,"Y");
	frame->Draw();
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yerrafb,yerrafb);
	g_afb->SetMarkerColor(4);
	g_afb->SetMarkerStyle(20);
	
	g_afb->SetFillColor(2);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
	c->Print(TString::Format("./plots/%s_afb.png",outfile));
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

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
		if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ) h1_nacc.Fill(gCosThetaL);
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
std::vector<double> recoEff(int iBin) // reconstruction efficiency
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
		if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ) h1_nacc.Fill(gCosThetaL);
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

void createAccptanceHist() // create acceptance histogram from UNFILTERED GEN.
{//{{{
	double accUpperBound = 0.09;
	double gQ2 = 0;
	double gCosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	
	TChain *treein=new TChain("tree");
	treein->Add("../RootFiles/MC_GENOnly/MC_genonly/BToKMuMu_GENOnly_8TeV_genonly_v3.root");
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
	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
	TH1F *h1_ngen[11];
	TH1F *h1_nacc[11];
	TH1F *h1_acc[11];
	TH1F *h1_ngen_fine[11];
	TH1F *h1_nacc_fine[11];
	TH1F *h1_acc_fine[11];
	for(int iBin = 0; iBin < 11; iBin++){
		if (iBin == 3 || iBin == 5) continue;
		h1_ngen[iBin] = new TH1F(TString::Format("h1_ngen_bin%d",iBin),"h1_ngen",6,thetaLBins);
		h1_nacc[iBin] = new TH1F(TString::Format("h1_nacc_bin%d",iBin) ,"h1_nacc" ,6,thetaLBins); 
		h1_acc [iBin] = new TH1F(TString::Format("h1_acc_bin%d",iBin),"",6,thetaLBins);
		h1_ngen_fine[iBin] = new TH1F(TString::Format("h1_ngen_fine_bin%d",iBin),"h1_ngen",20,-1,1);
		h1_nacc_fine[iBin] = new TH1F(TString::Format("h1_nacc_fine_bin%d",iBin) ,"h1_nacc" ,20,-1,1); 
		h1_acc_fine[iBin]  = new TH1F(TString::Format("h1_acc_fine_bin%d",iBin),"",20,-1,1);
		h1_ngen[iBin]->SetTitleOffset(1.3,"XY");
		h1_ngen[iBin]->SetXTitle("genCosThetaL");
		h1_ngen[iBin]->SetYTitle("Generated events");
		h1_nacc[iBin]->SetTitleOffset(1.3,"XY");
		h1_nacc[iBin]->SetXTitle("genCosThetaL");
		h1_nacc[iBin]->SetYTitle("Events in acceptance");
		h1_acc [iBin]->SetStats(0);
		h1_acc [iBin]->SetMinimum(0.);
		h1_acc [iBin]->SetMaximum(accUpperBound);
		h1_acc [iBin]->SetTitleOffset(1.3,"XY");
		h1_acc [iBin]->SetXTitle("genCosThetaL");
		h1_acc [iBin]->SetYTitle("Acceptance");
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
			if (iBin == 3 || iBin == 5) continue;
			if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
			h1_ngen[iBin]->Fill(gCosThetaL);
			h1_ngen_fine[iBin]->Fill(gCosThetaL);
			if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ){
				h1_nacc[iBin]->Fill(gCosThetaL);
				h1_nacc_fine[iBin]->Fill(gCosThetaL);
			}
		}
	}
	for(int iBin = 0; iBin < 11; iBin++){
		if (iBin == 3 || iBin == 5) continue;
	//	Calculate acceptance
	//	h1_acc[iBin]->SetAxisRange(0.,1.,"Y");
		for (int i = 1; i <= 6; i++) {
		//	Fill acceptance
			if (h1_ngen[iBin]->GetBinContent(i) == 0) {
				printf("WARNING: Acceptance(%d)=%f/%f\n",i,h1_nacc[iBin]->GetBinContent(i),h1_ngen[iBin]->GetBinContent(i));
				h1_acc[iBin]->SetBinContent(i,0.);
				h1_acc[iBin]->SetBinError(i,1.);
			}else{
				h1_acc[iBin]->SetBinContent(i,h1_nacc[iBin]->GetBinContent(i)/h1_ngen[iBin]->GetBinContent(i));
				if (h1_nacc[iBin]->GetBinContent(i) != 0){
					h1_acc[iBin]->SetBinError(i,sqrt(h1_acc[iBin]->GetBinContent(i)*(1.-h1_acc[iBin]->GetBinContent(i))/h1_ngen[iBin]->GetBinContent(i)));
				}else{
					h1_acc[iBin]->SetBinError(i,0.);
				}
			}
		}
		printf("INFO: h1_acc_bin%d built.\n",iBin);
		
	//	h1_acc_fine[iBin]->SetAxisRange(0.,1.,"Y");
		for (int i = 1; i <= 20; i++) {//L
		//	Fill acceptance
			if (h1_ngen_fine[iBin]->GetBinContent(i) == 0) {
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
		}
		printf("INFO: h1_acc_fine_bin%d built.\n",iBin);
	}
	fout->Write();
	fout->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

void createRecoEffHist(int iBin)
{//{{{
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
//	double effUpperBound = 0.03;
	double BMass = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double gCosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	
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
	const int nLBins = 20;
//	float thetaLBins[nLBins]={-1,-0.7,-0.3,0.,0.3,0.7,1};//nLBins=6
//	TH1F h1_nacc("h1_nacc" ,"h1_nacc" ,6,thetaLBins,5,thetaKBins); 
//	TH1F h1_nreco("h1_nreco","h1_nreco",6,thetaLBins,5,thetaKBins);
	TH1F h1_nacc("h1_nacc","h1_nacc",nLBins,-1,1);
	TH1F h1_nreco("h1_nreco","h1_nreco",nLBins,-1,1);
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (iBin == 3 || iBin == 5) continue;
		if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
		if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ) h1_nacc.Fill(gCosThetaL);
		if (BMass != 0 && ((Mumumass > 3.096916+3.*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) 
			 && (Mumumass > 3.686109+3.*Mumumasserr || Mumumass < 3.686109-     3.*Mumumasserr)) ){
			h1_nreco.Fill(gCosThetaL);
		}
	}
	
//	Calculate efficiency
//	TH1F h1_rec("h1_rec","",6,thetaLBins,5,thetaKBins);
	TH1F h1_rec("h1_rec","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
	//	Build from MC samples
		if (h1_nacc.GetBinContent(i) == 0 || h1_nreco.GetBinContent(i) == 0) {
			printf("WARNING: RecoEfficiency(%d)=0, set error to be 1.\n",i);
			h1_rec.SetBinContent(i,0.);
			h1_rec.SetBinError(i,1.);
		}else{
			h1_rec.SetBinContent(i,h1_nreco.GetBinContent(i)/h1_nacc.GetBinContent(i));
			h1_rec.SetBinError(i,sqrt(h1_rec.GetBinContent(i)*(1-h1_rec.GetBinContent(i))/h1_nreco.GetBinContent(i)));
			printf("INFO: RecoEfficiency(%d)=%f +- %f.\n",i,h1_rec.GetBinContent(i),h1_rec.GetBinError(i));
		}
	}
	h1_rec.SetTitleOffset(1.3,"XY");
	h1_rec.SetXTitle("CosThetaL");
	h1_rec.SetStats(0);
	h1_rec.SetMinimum(0.);
	
	TH1F h_recL("h_recL","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		double nacc = 0;
		double nreco = 0;
		nacc+= h1_nacc.GetBinContent(i);
		nreco+= h1_nreco.GetBinContent(i);
		if (nacc !=0 ){
			h_recL.SetBinContent(i,nreco/nacc);
			h_recL.SetBinError(i,sqrt(h_recL.GetBinContent(i)*(1-h_recL.GetBinContent(i))/nacc));
		}else{
			h_recL.SetBinContent(i,0);
			h_recL.SetBinError(i,1);
		}
	}
	h_recL.SetStats(0);
	h_recL.SetMinimum(0.);
	h_recL.SetXTitle("CosThetaL");
	
//	Print
	TCanvas canvas("canvas");
	h1_rec.Draw("PE1 TEXT");
	canvas.Print(TString::Format("./plots/recoEff_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/recoEff_bin%d.png",iBin));
	h_recL.Draw("PE1 TEXT");
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEff_cosl_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/recoEff_cosl_bin%d.png",iBin));
}//}}}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////   2D
/*   
std::vector<double> accXrecoEff(int iBin) // acceptance*reconstruction efficiency
{//{{{
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
//	double effUpperBound = 0.03;
	double BMass = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double gCosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	
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
//	Load acceptance
	TFile f_acc("./RootFiles/acceptance_8TeV.root");
	TH1F *h1_acc = (TH1F*)f_acc.Get(TString::Format("h1_acc_bin%d",iBin));
//	TH1F *h1_acc_fine = (TH1F*)f_acc_fine.Get(TString::Format("h1_acc_fine_bin%d",iBin));
//	Fill histograms
	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
	TH1F h1_nacc("h1_nacc" ,"h1_nacc" ,6,thetaLBins);
	TH1F h1_nreco("h1_nreco","h1_nreco",6,thetaLBins);
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
		if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ) h1_nacc.Fill(gCosThetaL);
		if (BMass != 0 && ((Mumumass > 3.096916+3.*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) 
			 && (Mumumass > 3.686109+3.*Mumumasserr || Mumumass < 3.686109-3.*Mumumasserr)) ) {
			h1_nreco.Fill(gCosThetaL);
		}
	}
//	Calculate efficiency
	TH1F h1_rec("h1_rec","",6,thetaLBins);
	for (int i = 1; i <= 6; i++) {
	//	Build from MC samples
		if (h1_nacc.GetBinContent(i) == 0 || h1_nreco.GetBinContent(i) == 0) {
			printf("WARNING: Efficiency(%d)=0, set error to be 1.\n",i);
			h1_rec.SetBinContent(i,0.);
			h1_rec.SetBinError(i,1.);
		}else{
			h1_rec.SetBinContent(i,h1_nreco.GetBinContent(i)/h1_nacc.GetBinContent(i)*h1_acc->GetBinContent(i));
			h1_rec.SetBinError(i,h1_rec.GetBinContent(i)*sqrt(-1./h1_nacc.GetBinContent(i)+1./h1_nreco.GetBinContent(i)+pow(h1_acc->GetBinError(i)/h1_acc->GetBinContent(i),2)));
			printf("INFO: Efficiency(%d)=%f +- %f.\n",i,h1_rec.GetBinContent(i),h1_rec.GetBinError(i));
		}
	}
	
//	Use Legendre polynomial for better convergance
//	1,x,(3x^2-1)/2,(5x^3-3x)/2,(35x^4-30x^2+3)/8
//	TString f2_model_format = "([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6";
	TString f1_model_format = "[0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2";
//	if (is7TeVCheck) f2_model_format = "([0]+[1]*y+[2]*y**2+[3]*y**3)+([4]+[5]*y+[6]*y**2+[7]*y**3)*x**2+([8]+[9]*y+[10]*y**2+[11]*y**3)*x**3+([12]+[13]*y+[14]*y**2+[15]*y**3)*x**4+([16]+[17]*y+[18]*y**2+[19]*y**3)*x**6";

//	Prepare draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
//	Draw efficiency
	TF1 f1_model("f1_model", f1_model_format, -1., 1.);
	h1_rec.SetMinimum(0.);
	h1_rec.SetTitleOffset(1.3,"XY");
	h1_rec.SetXTitle("genCosThetaL");
	h1_rec.SetStats(0);
//	h1_rec.SetMaximum(effUpperBound);
	h1_rec.Draw("PE1");
	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));

//	Draw FitResult
	f1_model.SetParameter(0,0.01);
	f1_model.SetParameter(1,0.01);
	f1_model.SetParameter(2,0.01);
	f1_model.SetParameter(3,0.01);
//	f1_model.SetParameter(4,0.01);
	
	h1_rec.Fit("f1_model"); //// 09-09
	
//	f1_model.SetTitle("");
//	f1_model.SetMaximum(effUpperBound);
//	f1_model.SetLineWidth(1);
//	f1_model.Draw(" SAME ");
	canvas.Print(TString::Format("./plots/accXrecoEff_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_bin%d.png",iBin));
	
//	Draw compare
	double chi2Val=0;
	int nPar = 4;
	double arrPar[4], arrParErr[4];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model.GetParameter(iPar);
		arrParErr[iPar] = f1_model.GetParError(iPar);
		chi2Val         = f1_model.GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
//	fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
//	printf("Chi2(Bin center)=%f \n",chi2Val);
	
	TH1F h1_compFit("h1_compFit","",6,thetaLBins);
	h1_compFit.SetTitleOffset(1.3,"XY");
	h1_compFit.SetXTitle("genCosThetaL");
	h1_compFit.SetYTitle(" ");
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
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.3,0.95,TString::Format("#varepsilon_{fit} / #varepsilon_{measured} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_compFit_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_compFit_bin%d.png",iBin));
	
//	Draw significance of deviation
	TH1F h_pull("Deviation/Error","",15,-3.,3.);
	h_pull.SetXTitle("Significance of deviation");
	for (int i = 1; i <= 6; i++) {//thetaL
		double _xlo = h1_rec.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h1_rec.GetXaxis()->GetBinUpEdge(i);
		if (h1_rec.GetBinContent(i) != 0){
			h_pull.Fill((f1_model.Integral(_xlo,_xhi)/(_xhi-_xlo)-h1_rec.GetBinContent(i))/h1_rec.GetBinError(i));
		}
	}
	h_pull.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_sigma_bin%d.pdf",iBin));
	
//	Draw projection to cosThetaK
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
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void accXrecoEff2(int iBin)
//std::vector<double> accXrecoEff2(int iBin)
{//{{{
	TH1::SetDefaultSumw2();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	double effUpperBound[11] = {0.14e-3, 0.25e-3, 0.30e-3, 0., 0.45e-3, 0., 0.65e-3, 1.0e-3, 1.30e-3, 0.25e-3, 0.40e-3};
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
	
	ch->SetBranchStatus("*",0);
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("genCosTheta*"  , 1);
	ch->SetBranchStatus("CosTheta*"     , 1);
	ch->SetBranchStatus("genMu*"        , 1);
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
	TH1F *h1_acc       = (TH1F*)f_acc.Get(TString::Format("h1_acc_bin%d",      iBin));
	TH1F *h1_acc_fine  = (TH1F*)f_acc.Get(TString::Format("h1_acc_fine_bin%d", iBin));
//	TH1F *h1_ngen      = (TH1F*)f_acc.Get(TString::Format("h1_ngen_bin%d",     iBin));
//	TH1F *h1_ngen_fine = (TH1F*)f_acc.Get(TString::Format("h1_ngen_fine_bin%d",iBin));
	
//	Fill histograms
	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
	TH1F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins); 
	TH1F h2_nreco("h2_nreco","h2_nreco",6,thetaLBins);
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
		if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
		if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ){
			h2_nacc.Fill(gCosThetaL);
			h2_nacc_fine.Fill(gCosThetaL);
		}
		if (BMass != 0 && ((Mumumass > 3.096916+3.*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) && (Mumumass > 3.686109+3.*Mumumasserr || Mumumass < 3.686109-3.*Mumumasserr)) ){
			h2_nreco.Fill(CosThetaL);
			h2_nreco_fine.Fill(CosThetaL);
		}
	}
	
//	Calculate efficiency
	TH1F h2_eff("h2_eff","",6,thetaLBins);
	for (int i = 1; i <= 6; i++) {//L
	//	Build from MC samples
		if (h2_nacc.GetBinContent(i) == 0 || h2_nreco.GetBinContent(i) == 0) {
			printf("WARNING: Efficiency(%d)0, set error to be 1.\n",i);
			h2_eff.SetBinContent(i,0.);
			h2_eff.SetBinError(i,1.);
		}else{
			h2_eff.SetBinContent(i,h2_nreco.GetBinContent(i)/h2_nacc.GetBinContent(i) * h1_acc->GetBinContent(i));
			h2_eff.SetBinError(i,h2_eff.GetBinContent(i)*sqrt(-1./h2_nacc.GetBinContent(i)+1./h2_nreco.GetBinContent(i)+pow(h1_acc->GetBinError(i)/h1_acc->GetBinContent(i),2)));
			printf("INFO: Efficiency(%d) : %f +- %f.\n",i,h2_eff.GetBinContent(i),h2_eff.GetBinError(i));
		}
	}
/*	
	TH1F h2_reco_fine("h2_reco_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
		h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_nreco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
	}
*/	
	TH1F h2_eff_fine("h2_eff_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_eff_fine.SetBinContent(i,0.);
			h2_eff_fine.SetBinError(i,1.);
		}else{
			h2_eff_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i) * h1_acc_fine->GetBinContent(i));
			h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*sqrt(-1./h2_nacc_fine.GetBinContent(i)+1./h2_nreco_fine.GetBinContent(i)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
		//	h2_eff_fine.SetBinContent(i,h2_reco_fine.GetBinContent(i)*h1_acc_fine->GetBinContent(i));
		//	h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i));
			printf("INFO: Efficiency_fine(%d)=%f +- %f.\n",i,h2_eff_fine.GetBinContent(i),h2_eff_fine.GetBinError(i));
		}
	}
/*	
//	Quick check 0th order (decoupled)
	TF1 *f_effL_ord0 = 0;
	if (iBin < 1){
		f_effL_ord0 = new TF1("f_effL_ord0","exp(-0.5*((x-[0])/[1])**2)*([2]*x**2+[3]*x+[4])",-1,1);//x
		f_effL_ord0->SetParameter(1,1.);
	}else{
		f_effL_ord0 = new TF1("f_effL_ord0","pol6",-1,1);//x
	}
	h2_eff_fine.Fit("f_effL_ord0","S");
*/

//	Use Legendre polynomial for better convergance
//	1,x,(3x^2-1)/2,(5x^3-3x)/2,(35x^4-30x^2+3)/8
//	TString f2_model_format_ord0 = TString::Format("(exp(-0.5*((x-(%f))/%f)**2)*(%f*x**2%+f*x%+f))*(%f%+f*y%+f*y**2%+f*y**3%+f*y**4%+f*y**5%+f*y**6)",f_effL_ord0->GetParameter(0),f_effL_ord0->GetParameter(1),f_effL_ord0->GetParameter(2),f_effL_ord0->GetParameter(3),f_effL_ord0->GetParameter(4),f_effK_ord0->GetParameter(0),f_effK_ord0->GetParameter(1),f_effK_ord0->GetParameter(2),f_effK_ord0->GetParameter(3),f_effK_ord0->GetParameter(4),f_effK_ord0->GetParameter(5),f_effK_ord0->GetParameter(6));
//	TString f1_model_format_ord0 = TString::Format("(exp(-0.5*((x-(%f))/%f)**2)*(%f*x**2%+f*x%+f))*(%f%+f*y%+f*y**2%+f*y**3%+f*y**4%+f*y**5%+f*y**6)",f_effL_ord0->GetParameter(0),f_effL_ord0->GetParameter(1),f_effL_ord0->GetParameter(2),f_effL_ord0->GetParameter(3),f_effL_ord0->GetParameter(4),f_effK_ord0->GetParameter(0),f_effK_ord0->GetParameter(1),f_effK_ord0->GetParameter(2),f_effK_ord0->GetParameter(3),f_effK_ord0->GetParameter(4),f_effK_ord0->GetParameter(5),f_effK_ord0->GetParameter(6));
//	if (iBin > 0) f2_model_format_ord0 = TString::Format("(%f%+f*x%+f*x**2%+f*x**3%+f*x**4%+f*x**5%+f*x**6)*(%f%+f*y%+f*y**2%+f*y**3%+f*y**4%+f*y**5%+f*y**6)",f_effL_ord0->GetParameter(0),f_effL_ord0->GetParameter(1),f_effL_ord0->GetParameter(2),f_effL_ord0->GetParameter(3),f_effL_ord0->GetParameter(4),f_effL_ord0->GetParameter(5),f_effL_ord0->GetParameter(6),f_effK_ord0->GetParameter(0),f_effK_ord0->GetParameter(1),f_effK_ord0->GetParameter(2),f_effK_ord0->GetParameter(3),f_effK_ord0->GetParameter(4),f_effK_ord0->GetParameter(5),f_effK_ord0->GetParameter(6));
//	printf("%s\n",f2_model_format_ord0.Data());
//	TString f2_model_format_ord1 = "([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6";
//	
//
//	TString f1_model_format_0 = "[0]+[1]*x+[2]*x**2";
//	TString f1_model_format_1 = "[0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2";
//	TString f1_model_format_2 = "[0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2+[4]*x**4+[5]*x**5+[6]*x**6";
	TString f1_model_format_3 = "[0]+[1]*x+[2]*(3*x**2-1)/2+[3]*(5*x**3-3*x)/2+[4]*x**4+[5]*x**5+[6]*x**6+[7]*x**7";

//	Draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
//	double chi2Val=0;
//	fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
//	printf("Chi2(Bin center)=%f \n",chi2Val);

//	Draw 1-D
	h2_nacc_fine.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff2_naccL_fine_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff2_naccL_fine_bin%d.png",iBin));
	
	h2_nreco_fine.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff2_nrecoL_fine_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff2_nrecoL_fine_bin%d.png",iBin));
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
	canvas.Print(TString::Format("./plots/accXrecoEff2_theoL_bin%d.pdf",iBin));
*/	
//	Draw efficiency
	h2_eff.SetMinimum(0.);
	h2_eff.SetTitleOffset(1.3,"XY");
	h2_eff.SetXTitle("genCosThetaL");
	h2_eff.SetYTitle("Efficiency");
	h2_eff.SetStats(0);
	h2_eff.SetMaximum(effUpperBound[iBin]);
	h2_eff.Draw("PE1");
	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
///////////////////////////////////////////////////////////////////////////////////////////////

//	Draw FitResult
	double up , dn ;
	if (iBin == 0) {up = 0.7; dn = -0.7;}
	else if (iBin == 1) { up = 0.8; dn = -0.8;}
	else {up = 1.; dn = -1.;}
	const int nPar = 8;
//	TF1 f1_model("f1_model", f1_model_format_0, dn, up);
//	TF1 f1_model("f1_model", f1_model_format_1, dn, up);
//	TF1 f1_model("f1_model", f1_model_format_2, dn, up);
	TF1 f1_model("f1_model", f1_model_format_3, dn, up);
	f1_model.SetParameter(0,0.01);
	f1_model.SetParameter(1,0.01);
	f1_model.SetParameter(2,0.01);
	f1_model.SetParameter(3,0.01);  // f1_model_format_1
	f1_model.SetParameter(4,0.01);  // f1_model_format_2
	f1_model.SetParameter(5,0.01);  // f1_model_format_2
	f1_model.SetParameter(6,0.01);  // f1_model_format_2
	f1_model.SetParameter(7,0.01);  // f1_model_format_3

	h2_eff.Fit("f1_model"); //// 09-09
	
	f1_model.SetTitle("");
	f1_model.SetMaximum(effUpperBound[iBin]);
	f1_model.SetLineColor(2);
	f1_model.SetLineWidth(1);
	f1_model.Draw(" SAME ");
	
/////////////////////////////////////////////////////////////	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff2_Eff_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff2_Eff_bin%d.png",iBin));
	
	h2_eff_fine.SetStats(0);
	h2_eff_fine.SetMinimum(0.);
	h2_eff_fine.SetTitleOffset(1.3,"XY");
	h2_eff_fine.SetXTitle("CosThetaL");
	h2_eff_fine.SetYTitle("Efficiency");
	h2_eff_fine.SetMaximum(effUpperBound[iBin]);
//	h2_eff_fine.Draw("TEXT");
	h2_eff_fine.Draw("PE1");
//	canvas.Update();
//	canvas.Print(TString::Format("./plots/accXrecoEff2_Eff_fine_bin%d.pdf",iBin));
/////////////////////////////////////////////////////////////////////////////////////////
	
//	Draw FitResult
//	TF1 f1_model("f1_model", f1_model_format, -1., 1.);
	f1_model.SetParameter(0,0.01);
	f1_model.SetParameter(1,0.01);
	f1_model.SetParameter(2,0.01);
	f1_model.SetParameter(3,0.01); // f1_model_format_1
	f1_model.SetParameter(4,0.01);  // f1_model_format_2
	f1_model.SetParameter(5,0.01);  // f1_model_format_2
	f1_model.SetParameter(6,0.01);  // f1_model_format_2
	f1_model.SetParameter(7,0.01);  // f1_model_format_3

	h2_eff_fine.Fit("f1_model"); //// 09-09
	
	f1_model.SetTitle("");
	f1_model.SetMaximum(effUpperBound[iBin]);
	f1_model.SetLineWidth(1);
	f1_model.SetLineColor(2);
	f1_model.Draw(" SAME ");
	
////////////////////////////////////////////////////////////////////
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d.png",iBin));
/////////////////////////////////////////////////////////////////////////////////////////////////
	
//	Draw compare
	double chi2Val=0;
//	const int nPar = 6;
	double arrPar[nPar], arrParErr[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model.GetParameter(iPar);
		arrParErr[iPar] = f1_model.GetParError(iPar);
		chi2Val         = f1_model.GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
//	Draw compare
	TH1F h2_compFit("h2_compFit","",6,thetaLBins);
	h2_compFit.SetTitleOffset(1.3,"XY");
	h2_compFit.SetXTitle("genCosThetaL");
	TH1F h2_pullFit("h2_pullFit","",6,thetaLBins);
	h2_pullFit.SetTitleOffset(1.3,"XY");
	h2_pullFit.SetXTitle("genCosThetaL");
	for (int i = 1; i <= 6; i++) {//thetaL
		if (h2_eff.GetBinContent(i) != 0){
			h2_compFit.SetBinContent(i,f1_model.Eval(h2_eff.GetXaxis()->GetBinCenter(i))/h2_eff.GetBinContent(i));
			double _xlo = h2_eff.GetXaxis()->GetBinLowEdge(i);
			double _xhi = h2_eff.GetXaxis()->GetBinUpEdge(i);
			h2_pullFit.SetBinContent(i,(f1_model.Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff.GetBinContent(i))/h2_eff.GetBinError(i));
		}else{
			h2_compFit.SetBinContent(i,0.);
			h2_pullFit.SetBinContent(i,0.);
		}
	}
	h2_compFit.SetMinimum(0.);
	h2_compFit.SetStats(0);
	h2_compFit.Draw("PE1");
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.3,0.95,TString::Format("#varepsilon_{fit} / #varepsilon_{measured} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff2_compFit_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff2_compFit_bin%d.png",iBin));
	
	h2_pullFit.SetStats(0);
	h2_pullFit.Draw(" TEXT");
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.3,0.95,TString::Format("(#varepsilon_{fit} - #varepsilon_{measured})/Error in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff2_pullFit_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff2_pullFit_bin%d.png",iBin));
	
//	Draw significance of deviation
	TH1F h2_pull("Deviation/Error","",15,-3.,3.);
	h2_pull.SetXTitle("Significance of deviation");
	for (int i = 1; i <= 6; i++) {//thetaL
		double _xlo = h2_eff.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h2_eff.GetXaxis()->GetBinUpEdge(i);
		if (h2_eff.GetBinContent(i) != 0){
			h2_pull.Fill((f1_model.Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff.GetBinContent(i))/h2_eff.GetBinError(i));
		}
	}
	h2_pull.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff2_sigma_bin%d.pdf",iBin));
	canvas.Print(TString::Format("./plots/accXrecoEff2_sigma_bin%d.png",iBin));
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Draw projection to cosThetaK
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
	
/*	
//	prepare output
	string output;
	output = TString::Format("%f*(%s)",arrPar[20],f1_model_format_ord0.Data());
//	output = TString::Format("%f*(%s)+(%f%+f*y%+f*(3*y**2-1)/2%+f*(5*y**3-3*y)/2)+(%f%+f*y%+f*(3*y**2-1)/2%+f*(5*y**3-3*y)/2)*x**2+(%f%+f*y%+f*(3*y**2-1)/2%+f*(5*y**3-3*y)/2)*x**3+(%f%+f*y%+f*(3*y**2-1)/2%+f*(5*y**3-3*y)/2)*x**4+(%f%+f*y%+f*(3*y**2-1)/2%+f*(5*y**3-3*y)/2)*x**6",arrPar[20],f2_model_format_ord0.Data(),arrPar[0],arrPar[1],arrPar[2],arrPar[3],arrPar[4],arrPar[5],arrPar[6],arrPar[7],arrPar[8],arrPar[9],arrPar[10],arrPar[11],arrPar[12],arrPar[13],arrPar[14],arrPar[15],arrPar[16],arrPar[17],arrPar[18],arrPar[19]);
//	printf("\"%s\",\n",output.c_str());
*/
	writeParam(iBin,"accXrecoEff2",arrPar,8);        // f1_model_format_3
	writeParam(iBin,"accXrecoEff2Err",arrParErr,8);  // f1_model_format_3
//	return output.c_str();	
}//}}}

std::vector<double> angular_reco_bin(int iBin, const char outfile[] = "angular_reco")
{//{{{
	double up , dn ;
	if (iBin == 0) {up = 0.89; dn = -0.89;}
	else if (iBin == 1) { up = 0.89; dn = -0.89;}
	else if (iBin == 2) { up = 0.89; dn = -0.89;}
	else { up = 1.; dn = -1.;}
//	double up = 1., dn = -1.;
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", dn, up);
	RooRealVar Q2("Q2","q^{2}",0.05,22.);
	RooRealVar fh("fh", "F_{H}", genFh[iBin], 0.0, 0.9);
	RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -0.3, 0.5);
	
	RooRealVar nsig("nsig","nsig",1E6,1E1,1E9);
//	RooRealVar nbkg("nbkg","nbkg",10,0.1,1E4);
	
//	Efficiency
	RooRealVar effP0("effP0","effP0",readParam(iBin,"accXrecoEff2", 0));
	RooRealVar effP1("effP1","effP1",readParam(iBin,"accXrecoEff2", 1));
	RooRealVar effP2("effP2","effP2",readParam(iBin,"accXrecoEff2", 2));
	RooRealVar effP3("effP3","effP3",readParam(iBin,"accXrecoEff2", 3));   // f1_model_format_1
	RooRealVar effP4("effP4","effP4",readParam(iBin,"accXrecoEff2", 4));   // f1_model_format_2
	RooRealVar effP5("effP5","effP5",readParam(iBin,"accXrecoEff2", 5));   // f1_model_format_2
	RooRealVar effP6("effP6","effP6",readParam(iBin,"accXrecoEff2", 6));   // f1_model_format_2
	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff2", 7));   // f1_model_format_3
	effP0.setError(readParam(iBin,"accXrecoEff2Err", 0));
	effP1.setError(readParam(iBin,"accXrecoEff2Err", 1));
	effP2.setError(readParam(iBin,"accXrecoEff2Err", 2));
	effP3.setError(readParam(iBin,"accXrecoEff2Err", 3));  // f1_model_format_1
	effP4.setError(readParam(iBin,"accXrecoEff2Err", 4));  // f1_model_format_2
	effP5.setError(readParam(iBin,"accXrecoEff2Err", 5));  // f1_model_format_2
	effP6.setError(readParam(iBin,"accXrecoEff2Err", 6));  // f1_model_format_2
	effP6.setError(readParam(iBin,"accXrecoEff2Err", 7));  // f1_model_format_3

/*  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	if (iBin == 0) {
		RooRealVar effP0("effP0","effP0", 9.93199e-05, 9.93199e-05, 9.93199e-05);
		RooRealVar effP1("effP1","effP1",-3.17861e-06,-3.17861e-06,-3.17861e-06);
		RooRealVar effP2("effP2","effP2",-1.28917e-04,-1.28917e-04,-1.28917e-04);
		effP0.setVal( 9.93199e-05);
		effP1.setVal(-3.17861e-06);
		effP2.setVal(-1.28917e-04);
		effP0.setError(2.10395e-06);
		effP1.setError(1.89847e-06);
		effP2.setError(4.15964e-06);
	}else if (iBin == 1) {
		RooRealVar effP0("effP0","effP0", 2.05111e-04, 2.05111e-04, 2.05111e-04);
		RooRealVar effP1("effP1","effP1",-1.20106e-06,-1.20106e-06,-1.20106e-06);
		RooRealVar effP2("effP2","effP2",-2.18593e-04,-2.18593e-04,-2.18593e-04);
		effP0.setError(2.98695e-06);
		effP1.setError(3.66300e-06);
		effP2.setError(6.71177e-06);
	}else if (iBin == 2) {
		RooRealVar effP0("effP0","effP0", 2.04444e-04, 2.04444e-04, 2.04444e-04);
		RooRealVar effP1("effP1","effP1", 1.57368e-05, 1.57368e-05, 1.57368e-05);
		RooRealVar effP2("effP2","effP2",-7.06549e-05,-7.06549e-05,-7.06549e-05);
		effP0.setError(2.36211e-06);
		effP1.setError(4.11057e-06);
		effP2.setError(8.25056e-06);
	}else if (iBin == 4) {
		RooRealVar effP0("effP0","effP0", 2.58065e-04, 2.58065e-04, 2.58065e-04);
		RooRealVar effP1("effP1","effP1", 2.83624e-05, 2.83624e-05, 2.83624e-05);
		RooRealVar effP2("effP2","effP2", 5.67763e-05, 5.67763e-05, 5.67763e-05);
		effP0.setError(3.71257e-06);
		effP1.setError(6.77705e-06);
		effP2.setError(1.41492e-05);
	}else if (iBin == 6) {
		RooRealVar effP0("effP0","effP0", 4.79717e-04, 4.79717e-04, 4.79717e-04);
		RooRealVar effP1("effP1","effP1", 4.82819e-05, 4.82819e-05, 4.82819e-05);
		RooRealVar effP2("effP2","effP2",-6.42579e-05,-6.42579e-05,-6.42579e-05);
		effP0.setError(7.33057e-06);
		effP1.setError(1.20176e-05);
		effP2.setError(2.46917e-05);
	}else if (iBin == 7) {
		RooRealVar effP0("effP0","effP0", 7.59667e-04, 7.59667e-04, 7.59667e-04);
		RooRealVar effP1("effP1","effP1", 9.87008e-05, 9.87008e-05, 9.87008e-05);
		RooRealVar effP2("effP2","effP2",-1.92149e-04,-1.92149e-04,-1.92149e-04);
		effP0.setError(1.06633e-05);
		effP1.setError(1.64954e-05);
		effP2.setError(3.47913e-05);
	}else if (iBin == 8) {
		RooRealVar effP0("effP0","effP0", 1.05255e-03, 1.05255e-03, 1.05255e-03);
		RooRealVar effP1("effP1","effP1", 1.07757e-04, 1.07757e-04, 1.07757e-04);
		RooRealVar effP2("effP2","effP2",-4.36530e-04,-4.36530e-04,-4.36530e-04);
		effP0.setError(1.29244e-05);
		effP1.setError(1.85828e-05);
		effP2.setError(3.88124e-05);
	}else if (iBin == 9) {
		RooRealVar effP0("effP0","effP0", 2.01912e-04, 2.01912e-04, 2.01912e-04);
		RooRealVar effP1("effP1","effP1", 3.57268e-06, 3.57268e-06, 3.57268e-06);
		RooRealVar effP2("effP2","effP2",-1.95027e-04,-1.95027e-04,-1.95027e-04);
		effP0.setError(2.06010e-06);
		effP1.setError(2.88150e-06);
		effP2.setError(5.43976e-06);
	}else if (iBin == 10) {
		RooRealVar effP0("effP0","effP0", 3.18287e-04, 3.18287e-04, 3.18287e-04);
		RooRealVar effP1("effP1","effP1", 3.18666e-05, 3.18666e-05, 3.18666e-05);
		RooRealVar effP2("effP2","effP2",-1.08468e-04,-1.08468e-04,-1.08468e-04);
		effP0.setError(1.60132e-06);
		effP1.setError(2.62222e-06);
		effP2.setError(5.34281e-06);
	}
*/	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	RooArgSet f_effA_argset(CosThetaL);
//	f_sigA_argset.add(RooArgSet(fh,afb));
//	f_effA_argset.add(RooArgSet(effP0, effP1, effP2));              // f1_model_format_0
//	f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3));       // f1_model_format_1
//	f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));       // f1_model_format_2
	f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));       // f1_model_format_3
//	TString f_sigA_format;
//	TString f_ang_format = "effP0 + effP1 * CosThetaL + effP2 * CosThetaL**2";
//	RooGenericPdf f_sigA("f_sigA", f_sigA_format, f_sigA_argset);
//
//
//	RooGenericPdf f_eff("f_eff", "effP0 + effP1 * CosThetaL + effP2 * CosThetaL**2", f_effA_argset);     // f1_model_format_0
//	RooGenericPdf f_eff("f_eff", "effP0+effP1*CosThetaL+effP2*(3*CosThetaL**2-1)/2+effP3*(5*CosThetaL**3-3*CosThetaL)/2", f_effA_argset); // f1_model_format_1
//	RooGenericPdf f_eff("f_eff", "effP0+effP1*CosThetaL+effP2*(3*CosThetaL**2-1)/2+effP3*(5*CosThetaL**3-3*CosThetaL)/2+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6", f_effA_argset); // f1_model_format_2
	RooGenericPdf f_eff("f_eff", "effP0+effP1*CosThetaL+effP2*(3*CosThetaL**2-1)/2+effP3*(5*CosThetaL**3-3*CosThetaL)/2+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6+effP7*CosThetaL**7", f_effA_argset); // f1_model_format_3
	RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb));
//	RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb),TString::Format("fabs(%s) <= (%s)/2.",afb,fh);
	
	RooProdPdf    f_effXsig("f_effXsig","", f_eff, f_sig);
	RooExtendPdf  f("f","", f_effXsig, nsig);
//	RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
//	RooExtendPdf f_ext("f_ext","f_ext",f_sigA,nsig);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Mumumass, Mumumasserr, CosThetaL),TString::Format("(%s) && (%s)",Q2range[iBin],mumuMassWindow[1]),0);
//	RooFitResult *f_fitresult = f_ext.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
	
//	Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = -0.;
	double UpperBound[11] ={110,220,400,0.,340,0.,340,430,640,430,2300};
	
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,Binning(100)); 
	f.plotOn(framecosl); 
	
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetMaximum(UpperBound[iBin]);
	framecosl->Draw();
	
//	fixNDC = -0.5;
//	if (iBin > 4) fixNDC = 0.;
	t1->DrawLatex(.30,.85+fixNDC,TString::Format("%s",Q2range[iBin]));
	t1->DrawLatex(.10,.79+fixNDC,TString::Format("F_{H}  =%9.5f #pm%9.5f",fh.getVal(),fh.getError()));
	t1->DrawLatex(.50,.79+fixNDC,TString::Format("A_{FB} =%9.5f #pm%9.5f",afb.getVal(),afb.getError()));
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_cosl_bin%d.png",outfile,iBin));
	
//	clear
	delete t1;
	delete c;
	delete data;
	
//	write output
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;

}//}}}

void angular_reco(const char outfile[] = "angular_reco")
{//{{{
//	bool refit = false; // Turn to true if you want to fit again.
	bool refit = true; // Turn to true if you want to fit again.
	
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
//	TH2F *frame = new TH2F("frame","",18,1,19,10,-0.5,1);
	frame->SetStats(kFALSE);
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.02,0.16,"Y");
	frame->Draw();
	
	double x[9]   ={1.25, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.75, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yuerrfh[9], yderrfh[9], yafb[9], yerrafb[9]; 
/*	double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
 	double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
	double yfh[8]       ={0.705,0.791,0.649,0.524,0.454,0.399,0.369,0.341};
	double yerrfh[8]    ={0.000565,0.000410,0.000273,0.000428,0.000296,0.000413,0.000367,0.000359};
	double yafb[8]      ={-0.160,-0.066,0.182,0.317,0.374,0.412,0.421,0.376};
	double yerrafb[8]   ={0.000432,0.000284,0.000224,0.000390,0.000286,0.000419,0.000393,0.420};
*/
   if (refit){ // Turn to true if you want to fit again.
		std::vector<double> vbin;
		for(int ibin = 0, i = 0; i < 11; i++){
			if (i == 3 || i == 5) continue;
			vbin = angular_reco_bin(i);
			yfh[ibin]       =vbin.at(0);
			yuerrfh[ibin]    =vbin.at(1);
			if (yuerrfh[ibin] > yfh[ibin]) { yderrfh[ibin] = yfh[ibin];}
			else yderrfh[ibin] = yuerrfh[ibin];
			yafb[ibin]      =vbin.at(2);
			yerrafb[ibin]   =vbin.at(3);
			ibin++;
		}
	}
//	Check input data
	for(int ibin = 0, i = 0; i < 11; i++){
		if (i == 3 || i == 5) continue;
		printf("  yafb[%d]=%12.6f +- %12.6f     ",ibin,yafb[ibin],yerrafb[ibin]);
		printf("yfh[%d]=%12.6f +- %12.6f\n",ibin,yfh[ibin],yuerrfh[ibin]);
		printf("genAfb[%d]=%6.4f      >>>>>>    ",ibin,genAfb[ibin]);
		printf("genFh[%d]=%6.4f\n\n",ibin,genFh[ibin]);
		ibin++;
	}

//	plotting
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	g_fh->SetMarkerColor(4);
	g_fh->SetMarkerStyle(20);
	
	g_fh->SetFillColor(2);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Print(TString::Format("./plots/%s_fh.png",outfile));
	c->Clear();
	
	frame->SetTitle("");
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetAxisRange(-0.04,0.06,"Y");
	frame->Draw();
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yerrafb,yerrafb);
	g_afb->SetMarkerColor(4);
	g_afb->SetMarkerStyle(20);
	
	g_afb->SetFillColor(2);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
	c->Print(TString::Format("./plots/%s_afb.png",outfile));
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
            f2_model.SetParameter(i,readParam(iBin,"accXrecoEff2",i));
            f2_model.SetParError(i,readParam(iBin,"accXrecoEff2Err",i));
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




/*   // 04-09-2014
//_________________________________________________________________________________
void angular2D_bin(int iBin, const char outfile[] = "angular2D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fitis adopted by Mauro, just follow!!

    // Read data
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar fh("fh", "F_{H}", genFh[iBin], 0., 1.);
    RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -1., 1.);

    // Efficiency
    RooRealVar recK0L0("recK0L0","recK0L0",readParam(iBin,"accXrecoEff2", 0));
    RooRealVar recK1L0("recK1L0","recK1L0",readParam(iBin,"accXrecoEff2", 1));
    RooRealVar recK2L0("recK2L0","recK2L0",readParam(iBin,"accXrecoEff2", 2));
    RooRealVar recK3L0("recK3L0","recK3L0",readParam(iBin,"accXrecoEff2", 3));
    RooRealVar recK0L2("recK0L2","recK0L2",readParam(iBin,"accXrecoEff2", 4));
    RooRealVar recK1L2("recK1L2","recK1L2",readParam(iBin,"accXrecoEff2", 5));
    RooRealVar recK2L2("recK2L2","recK2L2",readParam(iBin,"accXrecoEff2", 6));
    RooRealVar recK3L2("recK3L2","recK3L2",readParam(iBin,"accXrecoEff2", 7));
    RooRealVar recK0L3("recK0L3","recK0L3",readParam(iBin,"accXrecoEff2", 8));
    RooRealVar recK1L3("recK1L3","recK1L3",readParam(iBin,"accXrecoEff2", 9));
    RooRealVar recK2L3("recK2L3","recK2L3",readParam(iBin,"accXrecoEff2",10));
    RooRealVar recK3L3("recK3L3","recK3L3",readParam(iBin,"accXrecoEff2",11));
    RooRealVar recK0L4("recK0L4","recK0L4",readParam(iBin,"accXrecoEff2",12));
    RooRealVar recK1L4("recK1L4","recK1L4",readParam(iBin,"accXrecoEff2",13));
    RooRealVar recK2L4("recK2L4","recK2L4",readParam(iBin,"accXrecoEff2",14));
    RooRealVar recK3L4("recK3L4","recK3L4",readParam(iBin,"accXrecoEff2",15));
    RooRealVar recK0L6("recK0L6","recK0L6",readParam(iBin,"accXrecoEff2",16));
    RooRealVar recK1L6("recK1L6","recK1L6",readParam(iBin,"accXrecoEff2",17));
    RooRealVar recK2L6("recK2L6","recK2L6",readParam(iBin,"accXrecoEff2",18));
    RooRealVar recK3L6("recK3L6","recK3L6",readParam(iBin,"accXrecoEff2",19));
    recK0L0.setError(readParam(iBin,"accXrecoEff2Err", 0));
    recK1L0.setError(readParam(iBin,"accXrecoEff2Err", 1));
    recK2L0.setError(readParam(iBin,"accXrecoEff2Err", 2));
    recK3L0.setError(readParam(iBin,"accXrecoEff2Err", 3));
    recK0L2.setError(readParam(iBin,"accXrecoEff2Err", 4));
    recK1L2.setError(readParam(iBin,"accXrecoEff2Err", 5));
    recK2L2.setError(readParam(iBin,"accXrecoEff2Err", 6));
    recK3L2.setError(readParam(iBin,"accXrecoEff2Err", 7));
    recK0L3.setError(readParam(iBin,"accXrecoEff2Err", 8));
    recK1L3.setError(readParam(iBin,"accXrecoEff2Err", 9));
    recK2L3.setError(readParam(iBin,"accXrecoEff2Err",10));
    recK3L3.setError(readParam(iBin,"accXrecoEff2Err",11));
    recK0L4.setError(readParam(iBin,"accXrecoEff2Err",12));
    recK1L4.setError(readParam(iBin,"accXrecoEff2Err",13));
    recK2L4.setError(readParam(iBin,"accXrecoEff2Err",14));
    recK3L4.setError(readParam(iBin,"accXrecoEff2Err",15));
    recK0L6.setError(readParam(iBin,"accXrecoEff2Err",16));
    recK1L6.setError(readParam(iBin,"accXrecoEff2Err",17));
    recK2L6.setError(readParam(iBin,"accXrecoEff2Err",18));
    recK3L6.setError(readParam(iBin,"accXrecoEff2Err",19));
    RooArgSet f_sigA_argset(CosThetaL,CosThetaK);
    f_sigA_argset.add(RooArgSet(fh,afb,fs,as));
    f_sigA_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
    f_sigA_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
    TString f_sigA_format;
    TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*fh*CosThetaK**2*(1-CosThetaL**2)+1/2*(1-fh)*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*afb*(1-CosThetaK**2)*CosThetaL))";
    TString f_rec_ord0 = f_accXrecoEff_ord0[iBin];
    TString f_rec_format, f_rec_L0, f_rec_L2, f_rec_L3, f_rec_L4, f_rec_L6;
    if (is7TeVCheck){
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*CosThetaK**2+recK3L0*CosThetaK**3)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*CosThetaK**2+recK3L2*CosThetaK**3)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*CosThetaK**2+recK3L3*CosThetaK**3)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*CosThetaK**2+recK3L4*CosThetaK**3)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*CosThetaK**2+recK3L6*CosThetaK**3)*CosThetaL**6";
    }else{
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
    }

    if (iBin == 0) {
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
        f_sigA_format = TString::Format("(%s+%s+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data(),f_ang_format.Data());
    }else if (iBin == 1) {
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else if (iBin > 1 && iBin < 6) {
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else{
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*%s",f_rec_ord0.Data(),f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_ang_format.Data());
    }
        // angular map of signal
    RooGenericPdf f_sigA("f_sigA", f_sigA_format, f_sigA_argset);
    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooExtendPdf f_ext("f_ext","f_ext",f_sigA,nsig);
    
    // Get data and apply unbinned fit
    //fs.setConstant(kTRUE); as.setConstant(kTRUE);
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Mumumass, Mumumasserr, CosThetaK,CosThetaL),TString::Format("(%s) && (%s)",Q2range[iBin],mumuMassWindow[1]),0);
    RooFitResult *f_fitresult = f_ext.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f_ext.plotOn(framecosk); 
    if (true) {
        double buffFh = fh.getVal();
        double buffAfb = afb.getVal();
        fh.setVal(genFh[iBin]); afb.setVal(genAfb[iBin]);
        f_ext.plotOn(framecosk, LineColor(2),LineWidth(2),LineStyle(2)); 
        fh.setVal(buffFh); afb.setVal(buffAfb);
    }
    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();
    
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    if (iBin > 3) fixNDC = -0.5;
    t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{H}=%5.3f#pm%8.6f",fh.getVal(),fh.getError()));
    t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB}=%5.3f#pm%8.6f",afb.getVal(),afb.getError()));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f_ext.plotOn(framecosl); 
    if (true) {
        double buffFh = fh.getVal();
        double buffAfb = afb.getVal();
        fh.setVal(genFh[iBin]); afb.setVal(genAfb[iBin]);
        f_ext.plotOn(framecosl, LineColor(2),LineWidth(2),LineStyle(2)); 
        fh.setVal(buffFh); afb.setVal(buffAfb);
    }
    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    fixNDC = 0.;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("CosThetaL,CosThetaK", 6, 5);
    h1->SetXTitle("CosThetaL");
    h1->Draw("LEGO2");
    c->Update();
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    //write output
    double output[4] = {0,0};
    output[0] = fh.getVal();
    output[1] = fh.getError();
    writeParam(iBin,"fh",output);
    output[0] = afb.getVal();
    output[1] = afb.getError();
    writeParam(iBin,"afb",output);
    return;
}//}}}




void angular3D_1a_Sm(int iBin, const char outfile[] = "angular3D_1a_Sm", bool keepParam = false)
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    
    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.28,5.25,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",.03,.01,.05);
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",.08,.05,.35);
    RooRealVar sigM_frac("sigM_frac","sigM_frac",.5,0.,1.);
    
    // Create signal distribution
        // mass distro of signal
    RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
    RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
    RooAddPdf f_sigM("f_sigM","f_sigM", f_sigMGauss1, f_sigMGauss2, sigM_frac);
    
    // Create combinatorial background distribution
    RooRealVar bkgCombM_c("bkgCombM_c","c1",0,-30,50);
    RooRealVar offset("offset","offset",-5.);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
    RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
    
    RooRealVar nsig("nsig","nsig",0,1E8);
    RooRealVar nbkg("nbkg","nbkg",0,1E8);
    RooAddPdf f("f", "f",RooArgList(f_sigM,f_bkgCombM),RooArgList(nsig,nbkg));

    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),Q2range[iBin],0);
    RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* frame = Bmass.frame(); 
    data->plotOn(frame,Binning(20)); 
    f.plotOn(frame); 
    f.plotOn(frame,Components(f_sigM),LineColor(2),LineWidth(2));
    f.plotOn(frame,Components(f_bkgCombM),LineColor(3),LineStyle(2),LineWidth(2));

    frame->SetTitle("");
    frame->SetMinimum(0);
    frame->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    //t1->DrawLatex(.35,.10,TString::Format("nbkg=%5.3f#pm%5.3f",nbkg.getVal(),nbkg.getError()));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

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
void angular3D_1b_YpPm(int iBin, const char outfile[] = "angular3D_1b_YpPm", bool keepParam = false)
{//{{{
    if (iBin ==0 || iBin%2 == 1){
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
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    
    // Create peak background distribution
    RooRealVar bkgGauss1_mean1("bkgGauss1_mean1","M_{K*#Mu#Mu}",5.05,5.,5.12);
    RooRealVar bkgGauss1_mean2("bkgGauss1_mean2","M_{K*#Mu#Mu}",5.0,4.8,5.20);
    RooRealVar bkgGauss2_mean1("bkgGauss2_mean1","M_{K*#Mu#Mu}",5.40,5.35,5.45);
    RooRealVar bkgGauss2_mean2("bkgGauss2_mean2","M_{K*#Mu#Mu}",5.40,5.35,5.45);
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
    
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",1E1,1,1E7);
    RooExtendPdf *f = 0;
    switch (iBin) {
        case 2:
            //1 double guassian ,4+4 deg. ploy
            f = new RooExtendPdf("f","f",f_bkgPeakM1,nbkgPeak);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            f = new RooExtendPdf("f","f",f_bkgPeakM12,nbkgPeak);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            f = new RooExtendPdf("f","f",f_bkgPeakMGauss21,nbkgPeak);
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

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

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
void angular3D_2a_PkPl(int iBin, const char outfile[] = "angular3D_2a_PkPl", bool keepParam = false)
{//{{{
    // Gaussian constraint on yields and mass is needed.
    if (iBin ==0 || iBin%2 == 1){
        // Pm is fhat(and the yield is 0) for bins other than 2,4,6
        if (keepParam){
            double val[3]={0,0,0};
            writeParam(iBin, "bkgPeakL_c1", val);
            writeParam(iBin, "bkgPeakL_c2", val);
            writeParam(iBin, "bkgPeakL_c3", val);
            writeParam(iBin, "bkgPeakL_c4", val);
            writeParam(iBin, "bkgPeakK_c1", val);
            writeParam(iBin, "bkgPeakK_c2", val);
            writeParam(iBin, "bkgPeakK_c3", val);
            writeParam(iBin, "bkgPeakK_c4", val);
        }
        return;
    }

    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    
    RooArgSet f_bkgPeakL_argset;
    RooArgSet f_bkgPeakK_argset;
    RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",0,-5,5);
    RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",0,-5,5);
    RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",0,-5,5);
    RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",0,-5,5);
    RooRealVar bkgPeakK_c1("bkgPeakK_c1","c1",0,-5,5);
    RooRealVar bkgPeakK_c2("bkgPeakK_c2","c2",0,-5,5);
    RooRealVar bkgPeakK_c3("bkgPeakK_c3","c3",0,-5,5);
    RooRealVar bkgPeakK_c4("bkgPeakK_c4","c4",0,-5,5);
    switch (iBin) {
        case 2:
            //1 double guassian ,4+4 deg. ploy
            f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
            f_bkgPeakK_argset.add(RooArgSet(bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3,bkgPeakK_c4));
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4));
            f_bkgPeakK_argset.add(RooArgSet(bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3,bkgPeakK_c4));
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c1,bkgPeakL_c2));
            f_bkgPeakK_argset.add(RooArgSet(bkgPeakK_c1,bkgPeakK_c2));
            break;
        default:
            break;
    }
    RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
    RooPolynomial f_bkgPeakK("f_bkgPeakK","f_bkgPeakK",CosThetaK,f_bkgPeakK_argset);
    RooProdPdf f_bkgPeakA("f_bkgPeakA", "f_bckPeakA",f_bkgPeakK,f_bkgPeakL);
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",50,0.,1E4);
    RooExtendPdf f_bkgPeakA_ext("f_bkgPeakA_ext","f_bkgPeakA_ext",f_bkgPeakA,nbkgPeak);

    // Gaussian Constraint
    RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,"nbkgPeak", 0)),RooConst(readParam(iBin, "nbkgPeak", 1)));
    
    // Get data
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK,CosThetaL,Q2),Q2range[iBin],0);
    RooFitResult *f_fitresult = f_bkgPeakA_ext.fitTo(*data,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));

    // Draw CosThetaK
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f_bkgPeakK.plotOn(framecosk); 

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));
    
    // Draw CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f_bkgPeakL.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("CosThetaL,CosThetaK", 6, 5);
    h1->Draw("LEGO2");
    c->Update();
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

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
        val[0] = bkgPeakK_c1.getVal();val[1] = bkgPeakK_c1.getError();
        writeParam(iBin, "bkgPeakK_c1", val);
        val[0] = bkgPeakK_c2.getVal();val[1] = bkgPeakK_c2.getError();
        writeParam(iBin, "bkgPeakK_c2", val);
        val[0] = bkgPeakK_c3.getVal();val[1] = bkgPeakK_c3.getError();
        writeParam(iBin, "bkgPeakK_c3", val);
        val[0] = bkgPeakK_c4.getVal();val[1] = bkgPeakK_c4.getError();
        writeParam(iBin, "bkgPeakK_c4", val);
    }
}//}}}
void angular3D_prior(int iBin, const char outfile[] = "angular3D_prior", bool keepParam = false)
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    
    // Create combinatorial background distribution
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",0.,-2.5,2.5);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",0.,-2.5,2.5);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",0.,-2.5,2.5);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",0.,-2.5,2.5);
    RooArgSet f_bkgCombL_argset;
    switch (iBin) {
        case 7:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
        case 4:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2));
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 2:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3));
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 3:
        case 5:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            break;
        default:
            bkgCombL_c1.setConstant(kTRUE);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",0.,-2.5,2.5);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",0.,-2.5,2.5);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",0.,-5,5);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",0.,-5,5);
    RooArgSet f_bkgCombK_argset;
    switch (iBin) {
        case 2:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2));
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 3:
        case 4:
        case 5:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4));
            break;
        default:
            bkgCombK_c1.setConstant(kTRUE);
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombK("f_bkgCombK","f_bkgCombK",CosThetaK,f_bkgCombK_argset);
    RooProdPdf f_bkgCombA("f_bkgCombA", "f_bckCombA",f_bkgCombK,f_bkgCombL);
    
    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaK,CosThetaL),TString::Format("%s && (Bmass > 5.38 || Bmass < 5.18)",Q2range[iBin]),0);
    RooFitResult *f_fitresult = f_bkgCombA.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f_bkgCombA.plotOn(framecosk); 

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));
    
    // 
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f_bkgCombA.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    // Prepare datacard
    if (keepParam){
        double val[3] = {0,0,0};
        val[0] = bkgCombL_c1.getVal();val[1] = bkgCombL_c1.getError();
        writeParam(iBin, "bkgCombL_c1", val);
        val[0] = bkgCombL_c2.getVal();val[1] = bkgCombL_c2.getError();
        writeParam(iBin, "bkgCombL_c2", val);
        val[0] = bkgCombL_c3.getVal();val[1] = bkgCombL_c3.getError();
        writeParam(iBin, "bkgCombL_c3", val);
        val[0] = bkgCombL_c4.getVal();val[1] = bkgCombL_c4.getError();
        writeParam(iBin, "bkgCombL_c4", val);
        val[0] = bkgCombK_c1.getVal();val[1] = bkgCombK_c1.getError();
        writeParam(iBin, "bkgCombK_c1", val);
        val[0] = bkgCombK_c2.getVal();val[1] = bkgCombK_c2.getError();
        writeParam(iBin, "bkgCombK_c2", val);
        val[0] = bkgCombK_c3.getVal();val[1] = bkgCombK_c3.getError();
        writeParam(iBin, "bkgCombK_c3", val);
        val[0] = bkgCombK_c4.getVal();val[1] = bkgCombK_c4.getError();
        writeParam(iBin, "bkgCombK_c4", val);
    }
}//}}}

std::vector<double> angular3D_bin(int iBin, const char outfile[] = "angular3D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    // Need some modification for accXrecoEff2.
    
    // Read data
    RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.,5.56);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);

    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.28,5.26,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0));
    sigGauss1_sigma.setError(readParam(iBin,"sigGauss1_sigma",1));
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0));
    sigGauss2_sigma.setError(readParam(iBin,"sigGauss2_sigma",1));
    RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0));
    sigM_frac.setError(readParam(iBin,"sigM_frac",1));
        // Angular parameters
    RooRealVar afb("afb", "A_{FB}", 0., -1., 1.);
    RooRealVar fh("fh", "F_{H}", 0.8, 0., 1.);
    RooRealVar fs("fs","F_{S}",0.01,0.,1.);//Derive from B0ToKstarJpsi, Bin3
    RooRealVar as("as","A_{S}",-0.1,-1.,1.);//Derive from B0ToKstarJpsi, Bin3
    if (iBin != 3 && iBin != 5){
        // 2011 cross check
        fs.setVal(0.0129254);
        fs.setAsymError(-0.00898344,0.0101371);
        as.setVal(-0.0975919);
        as.setAsymError(-0.00490805,0.0049092);

        // read parameter from datacard
        //fs.setVal(readParam(3,"fs",0));
        //fs.setAsymError(readParam(3,"fs",1),readParam(3,"fs",2));
        //as.setVal(readParam(3,"as",0));
        //fs.setAsymError(readParam(3,"as",1),readParam(3,"as",2));
    }
        // Efficiency and acceptance
    RooRealVar recK0L0("recK0L0","recK0L0",readParam(iBin,"accXrecoEff2", 0));
    RooRealVar recK1L0("recK1L0","recK1L0",readParam(iBin,"accXrecoEff2", 1));
    RooRealVar recK2L0("recK2L0","recK2L0",readParam(iBin,"accXrecoEff2", 2));
    RooRealVar recK3L0("recK3L0","recK3L0",readParam(iBin,"accXrecoEff2", 3));
    RooRealVar recK0L2("recK0L2","recK0L2",readParam(iBin,"accXrecoEff2", 4));
    RooRealVar recK1L2("recK1L2","recK1L2",readParam(iBin,"accXrecoEff2", 5));
    RooRealVar recK2L2("recK2L2","recK2L2",readParam(iBin,"accXrecoEff2", 6));
    RooRealVar recK3L2("recK3L2","recK3L2",readParam(iBin,"accXrecoEff2", 7));
    RooRealVar recK0L3("recK0L3","recK0L3",readParam(iBin,"accXrecoEff2", 8));
    RooRealVar recK1L3("recK1L3","recK1L3",readParam(iBin,"accXrecoEff2", 9));
    RooRealVar recK2L3("recK2L3","recK2L3",readParam(iBin,"accXrecoEff2",10));
    RooRealVar recK3L3("recK3L3","recK3L3",readParam(iBin,"accXrecoEff2",11));
    RooRealVar recK0L4("recK0L4","recK0L4",readParam(iBin,"accXrecoEff2",12));
    RooRealVar recK1L4("recK1L4","recK1L4",readParam(iBin,"accXrecoEff2",13));
    RooRealVar recK2L4("recK2L4","recK2L4",readParam(iBin,"accXrecoEff2",14));
    RooRealVar recK3L4("recK3L4","recK3L4",readParam(iBin,"accXrecoEff2",15));
    RooRealVar recK0L6("recK0L6","recK0L6",readParam(iBin,"accXrecoEff2",16));
    RooRealVar recK1L6("recK1L6","recK1L6",readParam(iBin,"accXrecoEff2",17));
    RooRealVar recK2L6("recK2L6","recK2L6",readParam(iBin,"accXrecoEff2",18));
    RooRealVar recK3L6("recK3L6","recK3L6",readParam(iBin,"accXrecoEff2",19));
    recK0L0.setError(readParam(iBin,"accXrecoEff2Err", 0));
    recK1L0.setError(readParam(iBin,"accXrecoEff2Err", 1));
    recK2L0.setError(readParam(iBin,"accXrecoEff2Err", 2));
    recK3L0.setError(readParam(iBin,"accXrecoEff2Err", 3));
    recK0L2.setError(readParam(iBin,"accXrecoEff2Err", 4));
    recK1L2.setError(readParam(iBin,"accXrecoEff2Err", 5));
    recK2L2.setError(readParam(iBin,"accXrecoEff2Err", 6));
    recK3L2.setError(readParam(iBin,"accXrecoEff2Err", 7));
    recK0L3.setError(readParam(iBin,"accXrecoEff2Err", 8));
    recK1L3.setError(readParam(iBin,"accXrecoEff2Err", 9));
    recK2L3.setError(readParam(iBin,"accXrecoEff2Err",10));
    recK3L3.setError(readParam(iBin,"accXrecoEff2Err",11));
    recK0L4.setError(readParam(iBin,"accXrecoEff2Err",12));
    recK1L4.setError(readParam(iBin,"accXrecoEff2Err",13));
    recK2L4.setError(readParam(iBin,"accXrecoEff2Err",14));
    recK3L4.setError(readParam(iBin,"accXrecoEff2Err",15));
    recK0L6.setError(readParam(iBin,"accXrecoEff2Err",16));
    recK1L6.setError(readParam(iBin,"accXrecoEff2Err",17));
    recK2L6.setError(readParam(iBin,"accXrecoEff2Err",18));
    recK3L6.setError(readParam(iBin,"accXrecoEff2Err",19));
    RooArgSet f_sigA_argset(CosThetaL,CosThetaK);
    f_sigA_argset.add(RooArgSet(fh,afb,fs,as));
    f_sigA_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
    f_sigA_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
    TString f_sigA_format;
    TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*fh*CosThetaK**2*(1-CosThetaL**2)+1/2*(1-fh)*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*afb*(1-CosThetaK**2)*CosThetaL))";
    TString f_rec_format, f_rec_L0, f_rec_L2, f_rec_L3, f_rec_L4, f_rec_L6;
    if (is7TeVCheck){
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*CosThetaK**2+recK3L0*CosThetaK**3)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*CosThetaK**2+recK3L2*CosThetaK**3)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*CosThetaK**2+recK3L3*CosThetaK**3)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*CosThetaK**2+recK3L4*CosThetaK**3)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*CosThetaK**2+recK3L6*CosThetaK**3)*CosThetaL**6";
    }else{
        f_rec_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_rec_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_rec_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_rec_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_rec_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";
    }

    if (iBin == 0) {
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_rec_L6.Data(),f_ang_format.Data());
    }else if (iBin == 1) {
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s)*%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else if (iBin > 1 && iBin < 6) {
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
        f_sigA_format = TString::Format("(%s+%s+%s+%s)*%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_rec_L4.Data(),f_ang_format.Data());
    }else{
        f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
        f_sigA_format = TString::Format("(%s+%s+%s)*%s",f_rec_L0.Data(),f_rec_L2.Data(),f_rec_L3.Data(),f_ang_format.Data());
    }
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
    RooRealVar bkgCombM_c("bkgCombM_c","c1",0.,-20,1);
    RooRealVar offset("offset","offset",-5.);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
    RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0),-2.5,2.5);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0),-2.5,2.5);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0),-2.5,2.5);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0),-2.5,2.5);
    RooArgSet f_bkgCombL_argset;
    switch (iBin) {
        case 7:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            bkgCombL_c2.setVal(0.);
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
        case 4:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2));
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 2:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3));
            bkgCombL_c4.setVal(0.);
            bkgCombL_c4.setConstant(kTRUE);
            break;
        case 3:
        case 5:
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            break;
        default:
            bkgCombL_c1.setVal(0.);
            bkgCombL_c2.setVal(0.);
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.);
            bkgCombL_c1.setConstant(kTRUE);
            bkgCombL_c2.setConstant(kTRUE);
            bkgCombL_c3.setConstant(kTRUE);
            bkgCombL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",readParam(iBin,"bkgCombK_c1",0),-2.5,2.5);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",readParam(iBin,"bkgCombK_c2",0),-2.5,2.5);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",readParam(iBin,"bkgCombK_c3",0),-2.5,2.5);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",readParam(iBin,"bkgCombK_c4",0),-2.5,2.5);
    RooArgSet f_bkgCombK_argset;
    switch (iBin) {
        case 2:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            bkgCombK_c2.setVal(0.);
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 0:
        case 1:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2));
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
        case 3:
        case 4:
        case 5:
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3,bkgCombK_c4));
            break;
        default:
            bkgCombK_c1.setVal(0.);
            bkgCombK_c2.setVal(0.);
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
            bkgCombK_c1.setConstant(kTRUE);
            bkgCombK_c2.setConstant(kTRUE);
            bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgCombK("f_bkgCombK","f_bkgCombK",CosThetaK,f_bkgCombK_argset);
    RooProdPdf f_bkgCombA("f_bkgCombA", "f_bckCombA",f_bkgCombK,f_bkgCombL);
    RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombA,f_bkgCombM);
    printf("INFO: f_bkgComb prepared.\n");
    
    // Create peak background distribution
    RooRealVar bkgGauss1_mean1("bkgGauss1_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss1_mean1",0));
    bkgGauss1_mean1.setError(readParam(iBin,"bkgGauss1_mean1",1));
    RooRealVar bkgGauss1_mean2("bkgGauss1_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss1_mean2",0));
    bkgGauss1_mean2.setError(readParam(iBin,"bkgGauss1_mean2",1));
    RooRealVar bkgGauss1_sigma1("bkgGauss1_sigma1","#sigma_{11}",readParam(iBin,"bkgGauss1_sigma1",0));
    bkgGauss1_sigma1.setError(readParam(iBin,"bkgGauss1_sigma1",1));
    RooRealVar bkgGauss1_sigma2("bkgGauss1_sigma2","#sigma_{12}",readParam(iBin,"bkgGauss1_sigma2",0));
    bkgGauss1_sigma2.setError(readParam(iBin,"bkgGauss1_sigma2",1));
    RooRealVar bkgM_frac1("bkgM_frac1","bkgM_frac1",readParam(iBin,"bkgM_frac1",0));
    bkgM_frac1.setError(readParam(iBin,"bkgM_frac1",1));
    RooRealVar bkgGauss2_mean1("bkgGauss2_mean1","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss2_mean1",0));
    bkgGauss2_mean1.setError(readParam(iBin,"bkgGauss2_mean1",1));
    RooRealVar bkgGauss2_mean2("bkgGauss2_mean2","M_{K*#Mu#Mu}",readParam(iBin,"bkgGauss2_mean2",0));
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
    RooRealVar bkgPeakK_c1("bkgPeakK_c1","c1",readParam(iBin,"bkgPeakK_c1",0));
    bkgPeakL_c1.setError(readParam(iBin,"bkgPeakK_c1",1));
    RooRealVar bkgPeakK_c2("bkgPeakK_c2","c2",readParam(iBin,"bkgPeakK_c2",0));
    bkgPeakL_c2.setError(readParam(iBin,"bkgPeakK_c2",1));
    RooRealVar bkgPeakK_c3("bkgPeakK_c3","c3",readParam(iBin,"bkgPeakK_c3",0));
    bkgPeakL_c3.setError(readParam(iBin,"bkgPeakK_c3",1));
    RooRealVar bkgPeakK_c4("bkgPeakK_c4","c4",readParam(iBin,"bkgPeakK_c4",0));
    bkgPeakL_c4.setError(readParam(iBin,"bkgPeakK_c4",1));
    RooArgSet f_bkgPeakL_argset(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakL_c4);
    RooArgSet f_bkgPeakK_argset(bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3,bkgPeakK_c4);
    switch (iBin) {// Should be fixed constants already.
        case 2:
            //1 double guassian ,4+4 deg. ploy
            bkgM_frac12.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            bkgM_frac12.setConstant(kTRUE);
            bkgM_frac2.setConstant(kTRUE);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            bkgPeakK_c3.setConstant(kTRUE);
            bkgPeakL_c3.setConstant(kTRUE);
            bkgPeakK_c4.setConstant(kTRUE);
            bkgPeakL_c4.setConstant(kTRUE);
            bkgM_frac12.setVal(0.);
            bkgM_frac1.setVal(1.);
            bkgM_frac12.setConstant(kTRUE);
            bkgM_frac2.setConstant(kTRUE);
            break;
        case 3:
        case 5:
            bkgPeakK_c1.setConstant(kFALSE);
            bkgPeakL_c1.setConstant(kFALSE);
            bkgPeakK_c2.setConstant(kFALSE);
            bkgPeakL_c2.setConstant(kFALSE);
            bkgPeakK_c3.setConstant(kFALSE);
            bkgPeakL_c3.setConstant(kFALSE);
            bkgPeakK_c4.setConstant(kFALSE);
            bkgPeakL_c4.setConstant(kFALSE);
        default:
            bkgPeakK_c1.setConstant(kTRUE);
            bkgPeakL_c1.setConstant(kTRUE);
            bkgPeakK_c2.setConstant(kTRUE);
            bkgPeakL_c2.setConstant(kTRUE);
            bkgPeakK_c3.setConstant(kTRUE);
            bkgPeakL_c3.setConstant(kTRUE);
            bkgPeakK_c4.setConstant(kTRUE);
            bkgPeakL_c4.setConstant(kTRUE);
            break;
    }
    RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
    RooPolynomial f_bkgPeakK("f_bkgPeakK","f_bkgPeakK",CosThetaK,f_bkgPeakK_argset);
    RooProdPdf f_bkgPeakA("f_bkgPeakA", "f_bckPeakA",f_bkgPeakK,f_bkgPeakL);
    RooProdPdf f_bkgPeak("f_bkgPeak", "f_bkgPeak",f_bkgPeakA,f_bkgPeakM12);
    printf("INFO: f_bkgPeak prepared.\n");

    // Observed spectrum = model*fullEfficiency
    RooRealVar nsig("nsig","nsig",10,0,5E3);
    RooRealVar nbkgComb("nbkgComb","nbkgComb",20,0,1E4);
    //RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",readParam(iBin,"nbkgPeak",0));
    //nbkgPeak.setError(readParam(iBin,"nbkgPeak",1));
    RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",50,0.,1E4);
    if (iBin == 0 || iBin %2 == 1){
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

    RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,"nbkgPeak",0)),RooConst(readParam(iBin,"nbkgPeak",1)));
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
    RooGaussian gaus_bkgPeakK_c1("gaus_bkgPeakK_c1","gaus_bkgPeakK_c1",bkgPeakK_c1,RooConst(readParam(iBin,"bkgPeakK_c1",0)),RooConst(readParam(iBin,"bkgPeakK_c1",1)));
    RooGaussian gaus_bkgPeakK_c2("gaus_bkgPeakK_c2","gaus_bkgPeakK_c2",bkgPeakK_c2,RooConst(readParam(iBin,"bkgPeakK_c2",0)),RooConst(readParam(iBin,"bkgPeakK_c2",1)));
    RooGaussian gaus_bkgPeakK_c3("gaus_bkgPeakK_c3","gaus_bkgPeakK_c3",bkgPeakK_c3,RooConst(readParam(iBin,"bkgPeakK_c3",0)),RooConst(readParam(iBin,"bkgPeakK_c3",1)));
    RooGaussian gaus_bkgPeakK_c4("gaus_bkgPeakK_c4","gaus_bkgPeakK_c4",bkgPeakK_c4,RooConst(readParam(iBin,"bkgPeakK_c4",0)),RooConst(readParam(iBin,"bkgPeakK_c4",1)));
    //RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(readParam(iBin,"fs",0)),RooConst(readParam(iBin,"fs",1)),RooConst(readParam(iBin,"fs",2)));
    //RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(readParam(iBin,"as",0)),RooConst(readParam(iBin,"as",1)),RooConst(readParam(iBin,"as",2)));
    RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(0.0129254),RooConst(0.00898344),RooConst(0.0101371));// 2011 result
    RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(-0.0975919),RooConst(0.00490805),RooConst(0.0049092));
    
    RooArgSet gausConstraints(gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac);
    switch (iBin) {
        case 2:
            //1 double guassian ,4+4 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgPeakK_c1,gaus_bkgPeakK_c2,gaus_bkgPeakK_c3,gaus_bkgPeakK_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            gausConstraints.add(gaus_nbkgPeak);
            break;
        case 4:
            //2 double guassian ,4+4 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2,gaus_bkgPeakL_c3,gaus_bkgPeakL_c4));
            gausConstraints.add(RooArgSet(gaus_bkgPeakK_c1,gaus_bkgPeakK_c2,gaus_bkgPeakK_c3,gaus_bkgPeakK_c4));
            gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1,gaus_bkgGauss1_mean2,gaus_bkgGauss1_sigma1,gaus_bkgGauss1_sigma2,gaus_bkgM_frac1));
            gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_mean2,gaus_bkgGauss2_sigma1,gaus_bkgGauss2_sigma2,gaus_bkgM_frac2));
            gausConstraints.add(gaus_bkgM_frac12);
            gausConstraints.add(gaus_nbkgPeak);
            break;
        case 6:
            //1 guassian ,2+2 deg. ploy
            gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1,gaus_bkgPeakL_c2));
            gausConstraints.add(RooArgSet(gaus_bkgPeakK_c1,gaus_bkgPeakK_c2));
            gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1,gaus_bkgGauss2_sigma1));
            gausConstraints.add(gaus_nbkgPeak);
            break;
    }
    if (iBin == 3 || iBin == 5) gausConstraints.add(RooArgSet(gaus_fs,gaus_as));
    
    // Get data and apply unbinned fit
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaK,CosThetaL,Q2),Q2range[iBin],0);
    f.fitTo(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framemass = Bmass.frame();
    data->plotOn(framemass,Binning(20));
    f.plotOn(framemass,LineColor(1));
    f.plotOn(framemass,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framemass,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framemass,Components(f_bkgPeak),LineColor(6),LineWidth(2),LineStyle(2));

    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaK
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f.plotOn(framecosk,LineColor(1)); 
    f.plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosk,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framecosk,Components(f_bkgPeak),LineColor(6),LineWidth(2),LineStyle(2));

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosk_bin%d.pdf",outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f.plotOn(framecosl,LineColor(1)); 
    f.plotOn(framecosl,Components(f_sig),LineColor(4),LineWidth(2));
    f.plotOn(framecosl,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f.plotOn(framecosl,Components(f_bkgPeak),LineColor(6),LineWidth(2),LineStyle(2));

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
    c->Update();
    c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    // write output
    double val[3]={0,0,0};
    val[0] = fh.getVal();val[1] = fh.getError();
    writeParam(iBin, "fh", val);
    val[0] = afb.getVal();val[1] = afb.getError();
    writeParam(iBin, "afb",val);
    val[0] = fs.getVal();val[1] = fs.getErrorHi();val[2] = fs.getErrorLo();
    writeParam(iBin, "fs", val, 3);
    val[0] = as.getVal();val[1] = as.getErrorHi();val[2] = fs.getErrorLo();
    writeParam(iBin, "as", val, 3);

    std::vector<double> output;
    output.push_back(fh.getVal());
    output.push_back(fh.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    output.push_back(fs.getVal());
    output.push_back(fs.getError());
    output.push_back(as.getVal());
    output.push_back(as.getError());
    return output;
}//}}}

void angular(const char outfile[] = "angular", bool doFit = true)
{//{{{

    double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
    double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
    double yafb[8],yerrafb[8],yfh[8],yerrfh[8];

    if (doFit){
        angular3D_bin(3);
        angular3D_bin(5);

        angular3D_bin(0);
        angular3D_bin(1);
        angular3D_bin(2);
        angular3D_bin(4);
        angular3D_bin(6);
        angular3D_bin(7);
    }

    // Checkout input data
    for(int ibin = 0; ibin < 8; ibin++){
        yfh[ibin]       = readParam(ibin,"fh",0);
        yerrfh[ibin]    = readParam(ibin,"fh",1);
        yafb[ibin]      = readParam(ibin,"afb",0);
        yerrafb[ibin]   = readParam(ibin,"afb",1);
        printf("yafb[%d]=%6.4f +- %6.4f\n",ibin,yafb[ibin],yerrafb[ibin]);
        printf("yfh [%d]=%6.4f +- %6.4f\n",ibin,yfh[ibin],yerrfh[ibin]);
    }
    
    // Draw
    TCanvas *c = new TCanvas("c");
    TH1F *frame = new TH1F("frame","",18,1,19,10,-1,1);
    frame->SetStats(kFALSE);

    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{H}");
    frame->SetAxisRange(0,1,"Y");
    frame->Draw();
    TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(8,x,yfh,xerr,xerr,yerrfh,yerrfh);
    g_fh->SetFillColor(2);
    g_fh->SetFillStyle(3001);
    g_fh->Draw("P2");
    c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
    c->Clear();

    frame->SetTitle("");
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("A_{FB}");
    frame->SetAxisRange(-1,1,"Y");
    frame->Draw();
    TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(8,x,yafb,xerr,xerr,yerrafb,yerrafb);
    g_afb->SetFillColor(2);
    g_afb->SetFillStyle(3001);
    g_afb->Draw("P2");
    c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}
*/  // 04-09-2014
//////////////////////////////////////////////////////////////////////////
//_________________________________________________________________________________
//_________________________________________________________________________________
int main(int argc, char** argv) {
//	Tags
	is7TeVCheck = false;   
//	Help message
	if (argc <= 2) {
		printf("Usage       : ./fit Function infile\n");
		printf("Functions   :\n");
		printf("    bmass               Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg.\n");
	//	printf("    fh_gen              Derive F_{H} and A_{FB} from cosThetaL distribution at GEN level.\n");
		printf("    angular_gen         Derive F_{H} and A_{FB} from cosThetaL distribution at GEN level.\n");
		printf("    acceptance          Get acceptance map from unfiltered signal GEN, |Mu pT| > 2.8 GeV, |Mu eta| < 2.3.\n");
		printf("    recoEff             Get reconstruction efficiency map from signal simulation.\n");
		printf("    createRecoEffHist \n");
		printf("    accXrecoEff         Get efficiency map from signal simulation.\n");
		printf("    angular_reco        Derive F_{H} and A_{FB} from cosThetaL distribution at RECO level.\n");
	//	printf("    angular2D           Same as angular_gen, but fit to data with efficiency correction, bkg component is NOT considered.\n");
	//	printf("    angular3D_1a_Sm     Leading step1 to angular3D, determine signal shape from simulation.\n");
	//	printf("    angular3D_1b_YpPm   Leading step2 to angular3D, determine mass spectrum of peaking bkg from simulation.\n");
	//	printf("    angular3D_2a_PkPl   Leading step3 to angular3D, determine angular dist. of peaking bkg from simulation.\n");
	//	printf("    angular3D_prior     Leading step4 to angular3D, fit to data sideband to get initial values of combinatorial bkg.\n");
	//	printf("    angular3D           Derive F_{H} and A_{FB} by fitting to mass and angular distribution.\n");
		printf("Remark      :\n");
		printf("    1. Outputs will be stored in ./plots, please keep the directory.\n");
		printf("    2. Wildcard is allowed for infile. But you must quote infile like \"inputData_Run*.root\"!\n");
		return 0;
	}
//	main
	TString func    = argv[1];
	TString infile  = argv[2];
	
	if (func == "bmass") {
		if (argc != 4){
			printf("./fit bmass infile binID\n");
			for (int i = 0; i < 11; i++) {
				if (i == 3 || i == 5) continue;
				printf("    Bin %d : %s\n",i,Q2range[i]);
			}
			return 0;
		}
		int iBin = atoi(argv[3]);
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="bmass";
		bmass(iBin, outfile); 
/*	}else if (func == "fh_gen"){
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="fh_gen";
		fh_gen(outfile);
*/	}else if (func == "angular_gen"){
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="angular_gen";
		angular_gen(outfile);
	}else if (func == "acceptance") {
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
	}else if (func == "createRecoEffHist") {
		ch->Add(infile.Data());
		createAccptanceHist();
		if (ch == NULL) gSystem->Exit(0);
		for (int iBin = 0; iBin < 11; iBin++) {
			if (iBin == 3 || iBin == 5) continue;
			createRecoEffHist(iBin);
		}
	}else if (func == "accXrecoEff") {
		ch->Add(infile.Data());
		createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
	//	ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		for (int iBin = 0; iBin < 11; iBin++) {
			if (iBin == 3 || iBin == 5) continue;
		//	if (iBin != 0) continue;  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
		//	accXrecoEff(iBin);//old style efficiency fitter
			accXrecoEff2(iBin);
		}
	}else if (func == "angular_reco"){
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="angular_reco";
		angular_reco(outfile);
	}
	 
	 
/*	 
	 else if (func == "angular2D"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="angular2D";
        for (int iBin = 0; iBin < 8; iBin++) {
            if (iBin == 3 || iBin == 5) continue;
            angular2D_bin(iBin);
        }
    }else if (func == "angular3D_1a_Sm" || func == "angular3D_1b_YpPm" || func == "angular3D_2a_PkPl" || func == "angular3D_prior"){
        void (*fx)(int, const char*, bool);
        if ( func == "angular3D_1a_Sm" ){
            fx = angular3D_1a_Sm;
        }else if (func == "angular3D_1b_YpPm"){
            fx = angular3D_1b_YpPm;
        }else if (func == "angular3D_2a_PkPl"){
            fx = angular3D_2a_PkPl;
        }else{
            fx = angular3D_prior;
        }
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        for (int iBin = 0; iBin < 8; iBin++) {
            fx(iBin,func,true);// By default overwrite exist parameters.
        }
    }else if (func == "angular3D"){
        if (argc != 4){
            printf("./fit angular3D infile doFit\n");
            printf("    If dofit is non-zero, perform fit and write results first.\n");
            printf("    Else, directly read fitParameters?.txt and make plots.\n");
            return 0;
        }
        bool doFit = false;
        if ( atoi(argv[3]) != 0 ) doFit = true;

        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        const char outfile[]="angular3D";
        angular(outfile, doFit);
    }
*/	 
// 04-09-2014
/////////////////////////////////////////////////////////////////////////////////////////////
	 else if (func == "test"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
       // const char outfile[]="test";
        for (int iBin = 0; iBin < 8; iBin++) {
            //getToyFromUnfilterGen(iBin);
            //createRecoEffHist(iBin);
            //accXrecoEff2(iBin);
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
