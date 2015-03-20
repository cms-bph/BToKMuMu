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
#include "tdrstyle.C"

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
//double genAfb[11]={-0.015, 0.003, -0.005, 0.02, 0.004,  0.02, 0.004,  0.003, -0.003, -0.003,  0.003};   // initi.....
//double genFh [11]={ 0.068, 0.016,  0.040, 0.05, 0.011,  0.05, 0.010,  0.006,  0.020,  0.013,  0.032};
//double genAfb[11]={0.01, 0.02,  0.02, 0.02,  0.02, 0.02,  0.02,  0.01,  0.02, 0.02,  0.02};   // diff68.....
//double genFh [11]={0.10, 0.05,  0.05, 0.05,  0.05, 0.05,  0.14,  0.01,  0.05, 0.05,  0.05};
//double genAfb[11]={0.02, 0.02,  0.02, 0.02,  0.02, 0.02,  0.02,  0.02, 0.02, 0.02,  0.02};   // diffall.....
//double genFh [11]={0.05, 0.05,  0.05, 0.05,  0.05, 0.05,  0.05,  0.05, 0.05, 0.05,  0.05};
double genAfb[11]={0.002, 0.002,  0.002, 0.002,  0.002, 0.002,  0.002,  0.002, 0.002, 0.002,  0.002};   // newgen.....
double genFh [11]={ 0.05,  0.05,   0.05,  0.05,   0.05,  0.05,   0.05,   0.05,  0.05,  0.05,   0.05};
// Gen Level results.....
//double iniAfb[11]={-0.00222, -0.00238, -0.00088, 0.02, -0.00020, 0.02, 0.00034, 0.00220, 0.00058, -0.00181, -0.00049};  
//double iniFh [11]={ 0.06656,  0.03288,  0.01399, 0.05,  0.00894, 0.05, 0.00728, 0.00733, 0.00750,  0.03436,  0.01763};
double iniAfb[11]={-0.00100, -0.00200, -0.00100, 0.02, -0.00100, 0.02, -0.00100, 2.60209e-03, -0.00100, -0.00100, -0.00100};  
double iniFh [11]={ 0.09000,  0.02000,  0.02000, 0.05,  0.02000, 0.05,  0.02000, 6.29467e-03,  0.02000,  0.02000,  0.02000};
//double iniAfb[11]={-0.00000, -0.00200, -0.01000, 0.02, -0.01000, 0.02, -0.01000,  0.02200, -0.01000, -0.01000, -0.01000};  
//double iniFh [11]={ 0.06656,  0.02000,  0.02000, 0.05,  0.02000, 0.05,  0.02000,  0.05000,  0.02000,  0.02000,  0.02000};
//double genAfb[9]={0.02, 0.07, -0.02, -0.03, -0.01, -0.09, 0.02, 0.02, -0.01};   // LHCb results.....
//double genFh [9]={0.05, 0.14,  0.04,  0.11,  0.08,  0.14, 0.14, 0.05,  0.02};
//double dataAfb[11]={0.002, 0.002,  0.002, 0.002,  0.002, 0.002,  0.002,  0.002, 0.002, 0.002,  0.002};   // newgen.....
//double dataFh [11]={ 0.05,  0.05,   0.05,  0.05,   0.05,  0.05,   0.05,   0.05,  0.05,  0.05,   0.05};
double dataAfb[11]={ 3.62745e-02, 6.99728e-02, 5.53337e-02, 0.02, -0.00100, 0.02, -2.82606e-02, -0.00100, 3.96597e-02, 0.00100, 6.47723e-02};  
double dataFh [11]={ 1.00000e-01, 1.40520e-01, 1.10691e-01, 0.05,  0.02000, 0.05,  5.27930e-01,  0.02000, 7.93503e-02, 0.02000, 1.29546e-01};

double DATAAfb[11]   ={ 3.6193e-01, 5.8762e-02, 5.53442e-02, 0.02, 1.0016e-01, 0.02, -2.8262e-02, 8.84645e-02, 3.9640e-02, 1.18472e-01, 6.47723e-02};  
double DATAAfberr[11]={ 2.19e-01,   3.66e-01,   1.93124e-06, 0.02, 1.07e-01,   0.02,  1.32e-01,   9.54468e-05, 2.60e-04,   2.34817e-01, 4.94694e-02};  
double DATAFh [11]   ={ 1.0000,     1.4123e-01, 1.10705e-01, 0.05, 2.1774e-01, 0.05,  5.2793e-01, 1.77042e-01, 7.9390e-02, 2.36965e-01, 1.29546e-01};
double DATAFherr [11]={ 1.68e-01,   5.05e-01,   9.78095e-07, 0.05, 2.50e-01,   0.05,  3.56e-01,   7.42411e-05, 5.29e-04,   2.09032e-07, 9.06398e-02};

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
//            printf("INFO: readParam, matched %15s!\n",valBuff);
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
//		 cout<<"INFO: readParam, get "<<parName<<"["<<iColumn<<"] = "<<output.at(iColumn)<<endl;
      //  printf("INFO: readParam, get %s[%d]= %f \n",parName,iColumn, output.at(iColumn));
        return output.at(iColumn);
    }else{
//        printf("ERROR: readParam, empty column! Return 0.\n");
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
double readOutput(int iBin, int Index, const char parName[], int iColumn)
{//{{{
    std::vector<double> output;
    char lineBuff[512];
    char *valBuff;
    memset(lineBuff,' ',512*sizeof(char));
    FILE *fp = NULL;
	 fp = fopen(TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt",iBin,iBin,Index),"r");
	 if (fp == NULL) { return 0; }
    while(fgets(lineBuff,512,fp) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
//            printf("INFO: readParam, matched %15s!\n",valBuff);
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
//		 cout<<"INFO: readParam, get "<<parName<<"["<<iColumn<<"] = "<<output.at(iColumn)<<endl;
      //  printf("INFO: readParam, get %s[%d]= %f \n",parName,iColumn, output.at(iColumn));
        return output.at(iColumn);
    }else{
//        printf("ERROR: readParam, empty column! Return 0.\n");
        return 0.;
    }
}//}}}
void writeOutput(int iBin, int Index, const char parName[], double *val, int nVal=3, bool overwrite=true)
{//{{{
    struct stat fiBuff;
    FILE *fi = 0;
    if (stat(TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt",iBin,iBin,Index),&fiBuff) == 0){
        rename(TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt",iBin,iBin,Index),TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt.temp",iBin,iBin,Index));
        fi = fopen(TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt.temp",iBin,iBin,Index),"r");
    }else{
        fi = fopen(TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt.temp",iBin,iBin,Index),"w");
    }
    
    bool parExist = false;
    char lineBuff[512];
    char *valBuff = 0;
    memset(lineBuff,' ',512*sizeof(char));
    FILE *fp = fopen(TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt",iBin,iBin,Index),"w");
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
    remove(TString::Format("./OutputValues/bin%d/OutputValues%d_%d.txt.temp",iBin,iBin,Index));
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
//	RooRealVar fh("fh", "F_{H}", iniFh[iBin], 0., 1.);
//	RooRealVar afb("afb", "A_{FB}", iniAfb[iBin], -1., 1.);
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
		printf("genAfb[%d] =%12.6f +- %12.6f     ",ibin,yafb[i],yerrafb[i]);
		printf("genFh[%d]  =%12.6f +- %12.6f\n",ibin,yfh[i],yerrfh[i]);
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
	frame->SetAxisRange(-0.1,0.5,"Y");  
//	frame->SetAxisRange(-0.02,0.08,"Y");
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
//	frame->SetAxisRange(-0.02,0.02,"Y");
	frame->SetAxisRange(-0.2,0.2,"Y");  
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
std::vector<double> angular_gen_R_bin(int iBin, const char outfile[] = "angular_gen_R")
{//{{{
	RooRealVar genCosThetaL("genCosThetaL", "gencos#theta_{L}", -1., 1.);
	RooRealVar genQ2("genQ2","q^{2}",1.0,22.);
//	RooRealVar fh("fh", "F_{H}", iniFh[iBin], 0., 1.);
//	RooRealVar afb("afb", "A_{FB}", iniAfb[iBin], -1., 1.);
	RooRealVar fh("fh", "F_{H}", genFh[iBin], 0., 1.);
	RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -1., 1.);
	
	RooRealVar nsig("nsig","nsig",1E6,1E2,1E9);
//	RooRealVar nbkg("nbkg","nbkg",10,0.1,1E4);
//	Acceptance
	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));
	accP0.setConstant(kTRUE);
	accP1.setConstant(kTRUE);
	accP2.setConstant(kTRUE);
	accP3.setConstant(kTRUE);
//	accP0.setError(readParam(iBin,"accErr", 0));
//	accP1.setError(readParam(iBin,"accErr", 1));
//	accP2.setError(readParam(iBin,"accErr", 2));
//	accP3.setError(readParam(iBin,"accErr", 3));  
	RooArgSet f_sigA_argset(genCosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
	TString f_sigA_format;
	TString f_ang_format = "( 0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL )";
	TString f_acc_format = "( accP0 + accP1 *exp(-0.5*(((genCosThetaL-accP2)/accP3)**2)) ) ";
	f_sigA_argset.add(RooArgSet(f_acc_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_acc_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sig("f_sig", f_sigA_format,f_sigA_argset);
	RooExtendPdf  f("f","", f_sig, nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaL,genQ2),genQ2range[iBin],0);
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(1));
	f_fitresult->Print();
	
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
	writeParam(iBin, "genfh_R", val);
	val[0] = afb.getVal();val[1] = afb.getError();
	writeParam(iBin, "genafb_R",val);
	
	printf("genAfb_R[%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"genafb_R",0), fabs(readParam(iBin,"genafb_R",1)));
	printf("genFh_R [%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"genfh_R",0),  fabs(readParam(iBin,"genfh_R",1)));

//	write output
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;

}//}}}

void angular_gen_R(const char outfile[] = "angular_gen_R")
{//{{{
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double ygenfh[9], ygenuerrfh[9], ygenderrfh[9], ygenafb[9], ygenerrafb[9]; 
	double yfh[9], yuerrfh[9], yderrfh[9], yafb[9], yerrafb[9]; 
	double yrecofh[9], yrecouerrfh[9], yrecoderrfh[9], yrecoafb[9], yrecoerrafb[9]; 
//	Check input 
	for(int i = 0, ibin = 0; i < 9, ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
	//	reco
		yrecoafb[i]      = readParam(ibin,"recoafb",0);
		yrecoerrafb[i]   = fabs(readParam(ibin,"recoafb",1));
		yrecofh[i]       = readParam(ibin,"recofh",0);
		yrecouerrfh[i]   = fabs(readParam(ibin,"recofh",1));
		if (yrecouerrfh[i] > fabs(yrecofh[i])) { yrecoderrfh[i] = fabs(yrecofh[i]);}
		else { yrecoderrfh[i] = yrecouerrfh[i]; }
		printf("recoAfb[%d]=%6.4f +- %6.4f\n",ibin,yrecoafb[i],yrecoerrafb[i]);
		printf("recoFh [%d]=%6.4f +- %6.4f\n",ibin,yrecofh[i],yrecouerrfh[i]);
	//	gen_R
		yafb[i]      = readParam(ibin,"genafb_R",0);
		yerrafb[i]   = fabs(readParam(ibin,"genafb_R",1));
		yfh[i]       = readParam(ibin,"genfh_R",0);
		yuerrfh[i]   = fabs(readParam(ibin,"genfh_R",1));
		if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		else { yderrfh[i] = yuerrfh[i]; }
		printf("genAfb_R[%d]=%6.4f +- %6.4f\n",ibin,yafb[i],yerrafb[i]);
		printf("genFh_R [%d]=%6.4f +- %6.4f\n",ibin,yfh[i],yuerrfh[i]);
	//	gen
		ygenafb[i]      = readParam(ibin,"genafb",0);
		ygenerrafb[i]   = fabs(readParam(ibin,"genafb",1));
		ygenfh[i]       = readParam(ibin,"genfh",0);
		ygenuerrfh[i]   = fabs(readParam(ibin,"genfh",1));
		if (ygenuerrfh[i] > fabs(ygenfh[i])) { ygenderrfh[i] = fabs(ygenfh[i]);}
		else { ygenderrfh[i] = ygenuerrfh[i]; }
		printf("genAfb[%d]=%6.4f +- %6.4f\n",ibin,ygenafb[i],ygenerrafb[i]);
		printf("genFh [%d]=%6.4f +- %6.4f\n",ibin,ygenfh[i],ygenuerrfh[i]);
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
//	gen_R
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,0.5,"Y");  
//	frame->SetAxisRange(-0.1,1.,"Y");  
//	frame->SetAxisRange(-0.02,0.08,"Y");  
	TGraphAsymmErrors *r_fh  = new TGraphAsymmErrors(7,x,yrecofh,xerr,xerr,yrecoderrfh,yrecouerrfh);
	r_fh->SetMarkerColor(4);
	r_fh->SetMarkerStyle(20);	
	r_fh->SetFillColor(4);
	r_fh->SetFillStyle(3001);
	r_fh->Draw("2");
	r_fh->Draw("P");
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	g_fh->SetMarkerColor(2);
	g_fh->SetMarkerStyle(24);	
	g_fh->SetFillColor(2);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
// gen_R and gen
	TGraphAsymmErrors *gen_fh  = new TGraphAsymmErrors(7,x,ygenfh,xerr,xerr,ygenderrfh,ygenuerrfh);
	gen_fh->SetMarkerColor(1);
	gen_fh->SetMarkerStyle(24);
	gen_fh->SetFillColor(1);
	gen_fh->SetFillStyle(3001);
//	gen_fh->Draw("2");
//	gen_fh->Draw("P");
	
	TLegend *leg =new TLegend(0.74,0.70,0.89,0.89,NULL,"brNDC");
	leg->AddEntry("r_fh"," reco "," P ");
	leg->AddEntry("g_fh"," gen_R "," P ");
	leg->AddEntry("gen_fh"," gen "," P ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
//	leg->Draw();
	
	c->Print(TString::Format("./plots/%s_fh.png",outfile));
	c->Clear();
	
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
//	frame->SetAxisRange(-0.02,0.02,"Y"); 
	frame->SetAxisRange(-0.2,0.2,"Y"); 
//	frame->SetAxisRange(-1.0,1.0,"Y"); 
	frame->Draw();
//	gen_R
	TGraphAsymmErrors *r_afb = new TGraphAsymmErrors(7,x,yrecoafb,xerr,xerr,yrecoerrafb,yrecoerrafb);
	r_afb->SetMarkerColor(4);
	r_afb->SetMarkerStyle(20);
	r_afb->SetFillColor(4);
	r_afb->SetFillStyle(3001);
	r_afb->Draw("2");
	r_afb->Draw("P");
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yerrafb,yerrafb);
	g_afb->SetMarkerColor(2);
	g_afb->SetMarkerStyle(20);
	g_afb->SetFillColor(2);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
	TGraphAsymmErrors *gen_afb = new TGraphAsymmErrors(7,x,ygenafb,xerr,xerr,ygenerrafb,ygenerrafb);
	gen_afb->SetMarkerColor(1);
	gen_afb->SetMarkerStyle(24);
	gen_afb->SetFillColor(1);
	gen_afb->SetFillStyle(3001);
//	gen_afb->Draw("2");
//	gen_afb->Draw("P");
//	leg->Draw();
	c->Print(TString::Format("./plots/%s_afb.png",outfile));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
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
	TString f1_model_format_0 = "[0] + [1]*exp(-0.5* ( ((x-[2])/[3])**2)) ";
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
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_1 = "[0]+[1]*x+[2]*(3*x**2-1)/2.+[3]*(5*x**3-3*x)/2.+[4]*(35*x**4-30*x**2+3)/8.+[5]*(63*x**5-70*x**3+15*x)/8. + [6]*(231*x**6-315*x**4+105*x**2-5)/16."; 
		f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
	//	f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_2 = "( [0]*exp(-0.5*((x-[1])/[2])**2) ) * ( [3]+[4]*x+[5]*(3*x**2-1)/2.+[6]*(5*x**3-3*x)/2.+[7]*(35*x**4-30*x**2+3)/8.+[8]*(63*x**5-70*x**3+15*x)/8. +[9]*(231*x**6-315*x**4+105*x**2-5)/16. )";
	//	f1_model_format_2 = "( ([0]*exp(-0.5*((x-[1])/[2])**2) ) + ([3]*exp(-0.5*((x-[4])/[5])**2)) + ( [6]*exp(-0.5*((x-[7])/[8])**2)) ) + [9]";
	}else { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
		f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]+[8]+[9]"; 
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
	writeParam(iBin,"reco",   arrPar_r,   nPar_r);  // f1_model_format_1
	writeParam(iBin,"recoErr",arrParErr_r,nPar_r);  // f1_model_format_1
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
//	const int nPar = 7;
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, dn, up);
	f1_model->SetParameter(0, readParam(iBin,"acc", 1));
	f1_model->SetParameter(1, readParam(iBin,"acc", 2));
	f1_model->SetParameter(2, readParam(iBin,"acc", 3));
	f1_model->SetParameter(3, readParam(iBin,"reco", 0));
	f1_model->SetParameter(4, readParam(iBin,"reco", 1));
	f1_model->SetParameter(5, readParam(iBin,"reco", 2));
	f1_model->SetParameter(6, readParam(iBin,"reco", 3));
//	f1_model->FixParameter(7, 0.0);
//	f1_model->FixParameter(8, 0.0);
//	f1_model->FixParameter(9, 0.0);
	if ( iBin == 0 || iBin == 1 || iBin == 9 ) {
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
/*		f1_model->SetParError(0, readParam(iBin,"accErr", 1));
		f1_model->SetParError(1, readParam(iBin,"accErr", 2));
		f1_model->SetParError(2, readParam(iBin,"accErr", 3));
		f1_model->SetParError(3, readParam(iBin,"recoErr", 0));
		f1_model->SetParError(4, readParam(iBin,"recoErr", 1));
		f1_model->SetParError(5, readParam(iBin,"recoErr", 2));
		f1_model->SetParError(6, readParam(iBin,"recoErr", 3));
		f1_model->SetParError(7, readParam(iBin,"recoErr", 4));
		f1_model->SetParError(8, readParam(iBin,"recoErr", 5));
		f1_model->SetParError(9, readParam(iBin,"recoErr", 6));
*/	} else {
		f1_model->FixParameter(7, 0.0);
		f1_model->FixParameter(8, 0.0);
		f1_model->FixParameter(9, 0.0);
	}

	h2_eff_fine.SetStats(0);
//	h2_eff_fine.SetMinimum(-0.00005);
	h2_eff_fine.SetMinimum(0.);
	h2_eff_fine.SetTitleOffset(1.3,"XY");
	h2_eff_fine.SetXTitle("CosThetaL");
	h2_eff_fine.SetYTitle("Efficiency");
	h2_eff_fine.SetMaximum(effUpperBound[iBin]); 
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
	f1_model->SetMaximum(effUpperBound[iBin]); 
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> angular_reco_bin(int iBin, const char outfile[] = "angular_reco")
{//{{{
	double up , dn ;
//	if (iBin == 0) {up = 0.80; dn = -0.80;} // 19-09
//	else if (iBin == 1) { up = 0.85; dn = -0.85;}
//	else if (iBin == 9) { up = 0.95; dn = -0.95;}
//	else { up = 1.; dn = -1.;}
//	cout<<genFh[iBin]<<endl;
	up = 1., dn = -1.;
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", dn, up);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
//	RooRealVar fh("fh",   "F_{H}",  readParam(iBin,"genfh",0),   0., 1.);
//	RooRealVar afb("afb", "A_{FB}", readParam(iBin,"genafb",0), -1., 1.);
	RooRealVar fh("fh",   "F_{H}",  iniFh[iBin],   0., 1.);
	RooRealVar afb("afb", "A_{FB}", iniAfb[iBin], -1., 1.);
	
	RooRealVar nsig("nsig","nsig",1E6,1E1,1E9);
//	RooRealVar nbkg("nbkg","nbkg",10,0.1,1E4);
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	Acceptance
	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));
	accP0.setConstant(kTRUE);
	accP1.setConstant(kTRUE);
	accP2.setConstant(kTRUE);
	accP3.setConstant(kTRUE);
//	accP0.setError(readParam(iBin,"accErr", 0));
//	accP1.setError(readParam(iBin,"accErr", 1));
//	accP2.setError(readParam(iBin,"accErr", 2));
//	accP3.setError(readParam(iBin,"accErr", 3));  
//	reco Efficiency
	RooRealVar recoP0("recoP0","recoP0",readParam(iBin,"reco", 0));
	RooRealVar recoP1("recoP1","recoP1",readParam(iBin,"reco", 1));
	RooRealVar recoP2("recoP2","recoP2",readParam(iBin,"reco", 2));
	RooRealVar recoP3("recoP3","recoP3",readParam(iBin,"reco", 3));   
	RooRealVar recoP4("recoP4","recoP4",readParam(iBin,"reco", 4));  
	RooRealVar recoP5("recoP5","recoP5",readParam(iBin,"reco", 5));   
	RooRealVar recoP6("recoP6","recoP6",readParam(iBin,"reco", 6)); 
	recoP0.setConstant(kTRUE);
	recoP1.setConstant(kTRUE);
	recoP2.setConstant(kTRUE);
	recoP3.setConstant(kTRUE);
	recoP4.setConstant(kTRUE);
	recoP5.setConstant(kTRUE);
	recoP6.setConstant(kTRUE);
//	recoP0.setError(readParam(iBin,"recoErr", 0));
//	recoP1.setError(readParam(iBin,"recoErr", 1));
//	recoP2.setError(readParam(iBin,"recoErr", 2));
//	recoP3.setError(readParam(iBin,"recoErr", 3));  
//	recoP4.setError(readParam(iBin,"recoErr", 4)); 
//	recoP5.setError(readParam(iBin,"recoErr", 5)); 
//	recoP6.setError(readParam(iBin,"recoErr", 6)); 
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////

/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
//	Total Efficiency
	RooRealVar effP0("effP0","effP0",readParam(iBin,"accXrecoEff", 0));
	RooRealVar effP1("effP1","effP1",readParam(iBin,"accXrecoEff", 1));
	RooRealVar effP2("effP2","effP2",readParam(iBin,"accXrecoEff", 2));
	RooRealVar effP3("effP3","effP3",readParam(iBin,"accXrecoEff", 3));   
	RooRealVar effP4("effP4","effP4",readParam(iBin,"accXrecoEff", 4));  
	RooRealVar effP5("effP5","effP5",readParam(iBin,"accXrecoEff", 5)); 
	RooRealVar effP6("effP6","effP6",readParam(iBin,"accXrecoEff", 6)); 
	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff", 7)); 
	RooRealVar effP8("effP8","effP8",readParam(iBin,"accXrecoEff", 8)); 
	RooRealVar effP9("effP9","effP9",readParam(iBin,"accXrecoEff", 9)); 
	effP0.setConstant(kTRUE);
	effP1.setConstant(kTRUE);
	effP2.setConstant(kTRUE);
	effP3.setConstant(kTRUE);
	effP4.setConstant(kTRUE);
	effP5.setConstant(kTRUE);
	effP6.setConstant(kTRUE);
	effP7.setConstant(kTRUE);
	effP8.setConstant(kTRUE);
	effP9.setConstant(kTRUE);
//	RooRealVar effP10("effP10","effP10",readParam(iBin,"accXrecoEff", 10));   
//	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
//	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
//	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
//	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
//	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
//	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
//	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
//	effP7.setError(readParam(iBin,"accXrecoEffErr", 7)); 
//	effP8.setError(readParam(iBin,"accXrecoEffErr", 8)); 
//	effP9.setError(readParam(iBin,"accXrecoEffErr", 9)); 
//	effP10.setError(readParam(iBin,"accXrecoEffErr", 10));

//	RooArgSet f_effA_argset(CosThetaL);
//	f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));
//	f_effA_argset.add(RooArgSet(effP8, effP9));      
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	RooArgSet f_accA_argset(CosThetaL);
//	f_accA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
//	RooArgSet f_recoA_argset(CosThetaL);
//	f_recoA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
	
//	RooGenericPdf f_acc("f_acc", "accP0 + accP1 *exp(-0.5*((CosThetaL-accP2)/accP3)^2) ", f_accA_argset);     
//	RooGenericPdf f_reco("f_reco", "recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6", f_recoA_argset);     
//	RooGenericPdf f_sigA("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb));
//	//RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb),TString::Format("fabs(%s) <= (%s)/2.",afb,fh);
	
//	RooProdPdf    f_eff("f_eff","", f_acc, f_reco);
//	//RooProdPdf    f_effXsig("f_effXsig","", f_acc, f_reco, f_sig);
//	RooProdPdf    f_sig("f_effXsig","", f_eff, f_sigA); 
//	RooExtendPdf  f("f","", f_sig, nsig);
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////

/////////////	For bins, except 0, 1, 9 ////////////////////////////////////////////////////////////////////////
	RooArgSet f_sigA_argset(CosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	TString f_sigA_format;
	TString f_rec_format;
	TString f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
	if (iBin != 0 && iBin !=1 && iBin != 9 ) { 
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
	}
	f_sigA_argset.add(RooArgSet(f_rec_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sig("f_sig", f_sigA_format,f_sigA_argset);
	RooExtendPdf  f("f","", f_sig, nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);    // 12-08
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Hesse(0),Minos(1),Strategy(2));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(2));
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(1));
	f_fitresult->Print();
//	Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = -0.;
	double UpperBound[11] ={350,800,1400,0.,1000,0.,800,1000,1200,1600,6000};
		
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,Binning(100)); 
	f.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetMaximum(UpperBound[iBin]);  
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
	writeParam(iBin, "recoafb", val);
	printf("recoAfb[%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"recoafb",0), fabs(readParam(iBin,"recoafb",1)));
	printf("recoFh [%d]=%6.4f +- %6.4f\n", iBin, readParam(iBin,"recofh",0),  fabs(readParam(iBin,"recofh",1)));
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
}//}}}

void angular_reco(const char outfile[] = "angular_reco")
{//{{{
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
		yerrafb[i]   = fabs(readParam(ibin,"recoafb",1));
		yfh[i]       = readParam(ibin,"recofh",0);
		yuerrfh[i]   = fabs(readParam(ibin,"recofh",1));
		if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		else { yderrfh[i] = yuerrfh[i]; }
		printf("recoAfb[%d]=%6.4f +- %6.4f\n",ibin,yafb[i],yerrafb[i]);
		printf("recoFh [%d]=%6.4f +- %6.4f\n",ibin,yfh[i],yuerrfh[i]);
	//	gen
		ygenafb[i]      = readParam(ibin,"genafb",0);
		ygenerrafb[i]   = fabs(readParam(ibin,"genafb",1));
		ygenfh[i]       = readParam(ibin,"genfh",0);
		ygenuerrfh[i]   = fabs(readParam(ibin,"genfh",1));
		if (ygenuerrfh[i] > fabs(ygenfh[i])) { ygenderrfh[i] = fabs(ygenfh[i]);}
		else { ygenderrfh[i] = ygenuerrfh[i]; }
		printf("genAfb[%d]=%6.4f +- %6.4f\n",ibin,ygenafb[i],ygenerrafb[i]);
		printf("genFh [%d]=%6.4f +- %6.4f\n",ibin,ygenfh[i],ygenuerrfh[i]);
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
	frame->SetAxisRange(-0.1,0.5,"Y");  
//	frame->SetAxisRange(-0.1,1.,"Y");  
//	frame->SetAxisRange(-0.02,0.08,"Y");  
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
//	frame->SetAxisRange(-0.02,0.02,"Y"); 
	frame->SetAxisRange(-0.2,0.2,"Y"); 
//	frame->SetAxisRange(-1.0,1.0,"Y"); 
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
   
	RooRealVar nsig("nsig","nsig",0,1E8);
	RooAddPdf f("f", "f",RooArgList(f_sigM),RooArgList(nsig));
	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),Q2range[iBin],0);
	RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(2));
//	RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Extended(kTRUE));  // 12-10-2014
	f_fitresult->Print();

	// Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	RooPlot* frame = Bmass.frame(); 
	data->plotOn(frame,Binning(20)); 
	f.plotOn(frame,LineColor(1)); 
	f.plotOn(frame,Components(f_sigM),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
	f.plotOn(frame,Components(f_sigM),LineStyle(2),LineColor(4),LineWidth(2));

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
//	paveText->AddText(Form(" sigmean   = %.3f #pm %.3f ", sigGauss_mean.getVal(), sigGauss_mean.getError())); 
//	paveText->AddText(Form("sig_sigma1 = %.3f #pm %.3f ", sigGauss1_sigma.getVal(), sigGauss1_sigma.getError())); 
//	paveText->AddText(Form("sig_sigma2 = %.3f #pm %.3f ", sigGauss2_sigma.getVal(), sigGauss2_sigma.getError())); 
//	paveText->AddText(Form("    frac   = %.3f #pm %.3f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->AddText(Form("#Chi^{2}   = %.2f  ", frame->chiSquare())); 
	paveText->Draw(); 
	
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

void angular2D_1b_YpPm(int iBin, const char outfile[] = "angular2D_1b_YpPm", bool keepParam = true)
{//{{{
	bool mctype = false;
	TString func = outfile;
	if (func == "angular2D_1b_YpPm_Jpsi") {
		if (iBin != 10 && iBin != 4 && iBin != 2){    
			mctype = true;
		}
	}else if (func == "angular2D_1b_YpPm_Psi") {
		if (iBin != 10 && iBin != 4 && iBin != 6){
			mctype = true;
		}
	}
	if (keepParam && mctype){
		double val[3]={1,0,0};
		writeParam(iBin, TString::Format("bkgGauss1_mean1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgGauss1_mean2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgGauss1_sigma1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgGauss1_sigma2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgM_frac1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgGauss2_mean1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgGauss2_mean2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgGauss2_sigma1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgGauss2_sigma2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgM_frac2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgM_frac12_%s",outfile), val);
		val[0]=0;
		writeParam(iBin, TString::Format("nbkgPeak_%s",outfile), val);
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
	
	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",1E2,1,1E5);
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
		      f = new RooExtendPdf("f","f",f_bkgPeakM2,nbkgPeak);
				break;
		case 10:
		//2 double guassian ,4+4 deg. ploy
				f = new RooExtendPdf("f","f",f_bkgPeakM12,nbkgPeak);
				break;
		default:
		      break;
	}	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),Q2range[iBin],0);
	RooFitResult *f_fitresult = f->fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Extended());
	f_fitresult->Print();
	
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
		writeParam(iBin, TString::Format("bkgGauss1_mean1_%s",outfile), val);
		val[0] = bkgGauss1_mean2.getVal();val[1] = bkgGauss1_mean2.getError();
		writeParam(iBin, TString::Format("bkgGauss1_mean2_%s",outfile), val);
		val[0] = bkgGauss1_sigma1.getVal();val[1] = bkgGauss1_sigma1.getError();
		writeParam(iBin, TString::Format("bkgGauss1_sigma1_%s",outfile), val);
		val[0] = bkgGauss1_sigma2.getVal();val[1] = bkgGauss1_sigma2.getError();
		writeParam(iBin, TString::Format("bkgGauss1_sigma2_%s",outfile), val);
		val[0] = bkgM_frac1.getVal();val[1] = bkgM_frac1.getError();
		writeParam(iBin, TString::Format("bkgM_frac1_%s",outfile), val);
		val[0] = bkgGauss2_mean1.getVal();val[1] = bkgGauss2_mean1.getError();
		writeParam(iBin, TString::Format("bkgGauss2_mean1_%s",outfile), val);
		val[0] = bkgGauss2_mean2.getVal();val[1] = bkgGauss2_mean2.getError();
		writeParam(iBin, TString::Format("bkgGauss2_mean2_%s",outfile), val);
		val[0] = bkgGauss2_sigma1.getVal();val[1] = bkgGauss2_sigma1.getError();
		writeParam(iBin, TString::Format("bkgGauss2_sigma1_%s",outfile), val);
		val[0] = bkgGauss2_sigma2.getVal();val[1] = bkgGauss2_sigma2.getError();
		writeParam(iBin, TString::Format("bkgGauss2_sigma2_%s",outfile), val);
		val[0] = bkgM_frac2.getVal();val[1] = bkgM_frac2.getError();
		writeParam(iBin, TString::Format("bkgM_frac2_%s",outfile), val);
		val[0] = bkgM_frac12.getVal();val[1] = bkgM_frac12.getError();
		writeParam(iBin, TString::Format("bkgM_frac12_%s",outfile), val);
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
		writeParam(iBin, TString::Format("nbkgPeak_%s",outfile), val);
	}
}//}}}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D_2a_PkPl(int iBin, const char outfile[] = "angular2D_2a_PkPl", bool keepParam = true)
{//{{{
	// Gaussian constraint on yields and mass is needed.
	bool mctype = false;
	TString func = outfile;
	const char read1b[] = "angular2D_1b_YpPm_Jpsi";      /////////////////////////////           CHANGE BY HAND!!!  ///////////////////////
//	const char read1b[] = "angular2D_1b_YpPm_Psi";
	if (func == "angular2D_2a_PkPl_Jpsi") {
		if (iBin != 10 && iBin != 4 && iBin != 2){
			mctype = true;
		}
	}else if (func == "angular2D_2a_PkPl_Psi") {
		if (iBin != 10 && iBin != 4 && iBin != 6){
			mctype = true;
		}
	}else return;
	if (keepParam && mctype){
		double val[3]={0,0,0};
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
//		writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
		return;
	}
	
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
	RooArgSet f_bkgPeakL_argset;
	RooRealVar bkgPeakL_c0("bkgPeakL_c0","c0",0,-5,5);
	RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",0,-5,5);
	RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",0,-5,5);
	RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",0,-5,5);
//	RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",0,-5,5);
	switch (iBin) {
		default:
			f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c0,bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3));
		//	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c4));
		   break;
	}
	RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",50,0.,1E3);
	RooExtendPdf f_bkgPeakL_ext("f_bkgPeakL_ext","f_bkgPeakL_ext",f_bkgPeakL,nbkgPeak);
	
	// Gaussian Constraint
	RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b), 0)),RooConst(readParam(iBin, TString::Format("nbkgPeak_%s",read1b), 1) ) );
	// Get data
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);
	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*data,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));
	f_fitresult->Print();
	
	// Draw CosThetaL
	TCanvas *c = new TCanvas();
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
		val[0] = bkgPeakL_c0.getVal();val[1] = bkgPeakL_c0.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		val[0] = bkgPeakL_c1.getVal();val[1] = bkgPeakL_c1.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		val[0] = bkgPeakL_c2.getVal();val[1] = bkgPeakL_c2.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
		val[0] = bkgPeakL_c3.getVal();val[1] = bkgPeakL_c3.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
	//	val[0] = bkgPeakL_c4.getVal();val[1] = bkgPeakL_c4.getError();
	//	writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
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
//	RooRealVar bkgCombL_c4("bkgCombL_c4","c4",0.,-3.,3.);
//	RooRealVar bkgCombL_c5("bkgCombL_c5","c5",0.,-3.,3.);
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
		case 1:
		case 2:
		case 9:
		case 10:
		   f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
		   f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
			break;
	   default:
		   f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
	   // f_bkgCombL_argset.add(RooArgSet(bkgCombL_c5,bkgCombL_c4));
		   break;
	}
	RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.384 || Bmass < 5.174)",Q2range[iBin]),0);
	RooFitResult *f_fitresult = f_bkgCombL.fitTo(*data,Save(kTRUE),Minimizer("Minuit"));
	f_fitresult->Print();
	
	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
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
		writeParam(iBin, "bkgCombL_c3", val);
		if (iBin == 1 || iBin == 2 || iBin ==9 || iBin == 10) {
		val[0] = bkgCombL_c3.getVal();val[1] = bkgCombL_c3.getError();
		writeParam(iBin, "bkgCombL_c3", val);
	//	val[0] = bkgCombL_c4.getVal();val[1] = bkgCombL_c4.getError();
	//	writeParam(iBin, "bkgCombL_c4", val);
	//	val[0] = bkgCombL_c5.getVal();val[1] = bkgCombL_c5.getError();
	//	writeParam(iBin, "bkgCombL_c5", val);
	   }
	}
}//}}}

//void angular2D_bin(int iBin, const char outfile[] = "angular2D")
std::vector<double> angular2D_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular2D")
//std::vector<double> angular_reco_bin(int iBin, const char outfile[] = "angular_reco")
{//{{{
	// Remark: You must use RooFit!! It's better in unbinned fit.
	//         Extended ML fit is adopted by Mauro, just follow!!
	//         Need some modification for accXrecoEff.
	const char read1b_J[] = "angular2D_1b_YpPm_Jpsi";      ////////
	const char read1b_P[] = "angular2D_1b_YpPm_Psi";
	const char read2a_J[] = "angular2D_2a_PkPl_Jpsi";      ////////
	const char read2a_P[] = "angular2D_2a_PkPl_Psi";
	cout<<endl<<"iBin = "<<iBin<<endl<<endl; 
	// Create parameters and PDFs
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
	RooRealVar Bmass("Bmass","M_{K^{+/-}#Mu#Mu}",5.,5.56);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	// // Angular parameters
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
	RooRealVar fh(  "fh",  "F_{H}", Ifh,   0., 1.);
//	RooRealVar afb("afb", "A_{FB}", dataAfb[iBin], -1., 1.);
//	RooRealVar fh(  "fh",  "F_{H}", dataFh[iBin],   0., 1.);
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	Acceptance
	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));   
	accP0.setConstant(kTRUE);
	accP1.setConstant(kTRUE);
	accP2.setConstant(kTRUE);
	accP3.setConstant(kTRUE);
//	accP0.setError(readParam(iBin,"accErr", 0));
//	accP1.setError(readParam(iBin,"accErr", 1));
//	accP2.setError(readParam(iBin,"accErr", 2));
//	accP3.setError(readParam(iBin,"accErr", 3));  
//	reco Efficiency
	RooRealVar recoP0("recoP0","recoP0",readParam(iBin,"reco", 0));
	RooRealVar recoP1("recoP1","recoP1",readParam(iBin,"reco", 1));
	RooRealVar recoP2("recoP2","recoP2",readParam(iBin,"reco", 2));
	RooRealVar recoP3("recoP3","recoP3",readParam(iBin,"reco", 3));   
	RooRealVar recoP4("recoP4","recoP4",readParam(iBin,"reco", 4));  
	RooRealVar recoP5("recoP5","recoP5",readParam(iBin,"reco", 5));   
	RooRealVar recoP6("recoP6","recoP6",readParam(iBin,"reco", 6)); 
	recoP0.setConstant(kTRUE);
	recoP1.setConstant(kTRUE);
	recoP2.setConstant(kTRUE);
	recoP3.setConstant(kTRUE);
	recoP4.setConstant(kTRUE);
	recoP5.setConstant(kTRUE);
	recoP6.setConstant(kTRUE);
//	recoP0.setError(readParam(iBin,"recoErr", 0));
//	recoP1.setError(readParam(iBin,"recoErr", 1));
//	recoP2.setError(readParam(iBin,"recoErr", 2));
//	recoP3.setError(readParam(iBin,"recoErr", 3));  
//	recoP4.setError(readParam(iBin,"recoErr", 4)); 
//	recoP5.setError(readParam(iBin,"recoErr", 5)); 
//	recoP6.setError(readParam(iBin,"recoErr", 6)); 
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
//	Total Efficiency
	RooRealVar effP0("effP0","effP0",readParam(iBin,"accXrecoEff", 0));
	RooRealVar effP1("effP1","effP1",readParam(iBin,"accXrecoEff", 1));
	RooRealVar effP2("effP2","effP2",readParam(iBin,"accXrecoEff", 2));
	RooRealVar effP3("effP3","effP3",readParam(iBin,"accXrecoEff", 3));   
	RooRealVar effP4("effP4","effP4",readParam(iBin,"accXrecoEff", 4));  
	RooRealVar effP5("effP5","effP5",readParam(iBin,"accXrecoEff", 5)); 
	RooRealVar effP6("effP6","effP6",readParam(iBin,"accXrecoEff", 6)); 
	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff", 7)); 
	RooRealVar effP8("effP8","effP8",readParam(iBin,"accXrecoEff", 8)); 
	RooRealVar effP9("effP9","effP9",readParam(iBin,"accXrecoEff", 9)); 
//	RooRealVar effP10("effP10","effP10",readParam(iBin,"accXrecoEff", 10));   
	effP0.setConstant(kTRUE);
	effP1.setConstant(kTRUE);
	effP2.setConstant(kTRUE);
	effP3.setConstant(kTRUE);
	effP4.setConstant(kTRUE);
	effP5.setConstant(kTRUE);
	effP6.setConstant(kTRUE);
	effP7.setConstant(kTRUE);
	effP8.setConstant(kTRUE);
	effP9.setConstant(kTRUE);
//	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
//	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
//	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
//	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
//	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
//	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
//	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
//	effP7.setError(readParam(iBin,"accXrecoEffErr", 7)); 
//	effP8.setError(readParam(iBin,"accXrecoEffErr", 8)); 
//	effP9.setError(readParam(iBin,"accXrecoEffErr", 9)); 
//	effP10.setError(readParam(iBin,"accXrecoEffErr", 10));
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
	
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	
	// // Signal double gaussian
	RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{+/-}#Mu#Mu}",5.279,5.26,5.30);
	RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0));
	sigGauss1_sigma.setConstant(kTRUE);
//	sigGauss1_sigma.setError(readParam(iBin,"sigGauss1_sigma",1));
	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0));
	sigGauss2_sigma.setConstant(kTRUE);
//	sigGauss2_sigma.setError(readParam(iBin,"sigGauss2_sigma",1));
	RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0));
	sigM_frac.setConstant(kTRUE);
//	sigM_frac.setError(readParam(iBin,"sigM_frac",1));
	// // mass distro of signal
	RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
	RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
	RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
	
	RooArgSet f_sigA_argset(CosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	TString f_sigA_format;
	TString f_rec_format;
	TString f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
	if (iBin != 0 && iBin != 1 && iBin != 9) {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
	}
//	cout<<"\n \n \n \n \nf_rec_format \n"<<f_rec_format<<"\n \n \n \n \n"<<endl;
	f_sigA_argset.add(RooArgSet(f_rec_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sigA("f_sigA", f_sigA_format,f_sigA_argset);
	
	// Create signal distribution
	RooProdPdf f_sig("f_sig","f_sig",f_sigM,f_sigA);
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_sig prepared. <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	
	// Create combinatorial background distribution
//	Build Chebychev polynomial p.d.f.  
/*	RooRealVar     a0("a0", "constant", 0.5, -1, 1);
	RooRealVar     a1("a1", "linear", 0.6, -1, 1);
	RooRealVar     a2("a2", "quadratic", 0.1, -1, 1);
	RooChebychev   f_bkgCombM("f_bkgCombM", "Background", Bmass, RooArgSet(a0, a1, a2));
*/	RooRealVar bkgCombM_c("bkgCombM_c","c",-0.1,-10.,1.);
	RooRealVar offset("offset","offset",-5.);
	RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
	RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",readParam(iBin,"bkgCombL_c0",0));
	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0));
	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0));
	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0)); 
//	RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0));
//	RooRealVar bkgCombL_c5("bkgCombL_c5","c5",readParam(iBin,"bkgCombL_c5",0),-3.,3.);
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
		default:
		      f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
				bkgCombL_c0.setConstant(kTRUE);
				bkgCombL_c1.setConstant(kTRUE);
				bkgCombL_c2.setConstant(kTRUE);
				if (iBin == 1 || iBin == 2 || iBin ==9 || iBin == 10) {
				f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
				bkgCombL_c3.setConstant(kTRUE);}
			//	bkgCombL_c4.setConstant(kTRUE);
			//	bkgCombL_c5.setConstant(kTRUE);
				break;
	}
	RooPolynomial f_bkgCombL("f_bkgCombL","f_bkgCombL",CosThetaL,f_bkgCombL_argset);
	RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombL,f_bkgCombM);
	cout<<">>>>>>>>>>>>>>>> INFO: f_bkgComb prepared. <<<<<<<<<<<<<<<<<<<<<<"<<endl;
	
	// Create peak background distribution
	// // JPsiK MC Peaking Bkg.
	RooRealVar bkgGauss1_mean1_J("bkgGauss1_mean1_J","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_J),0));
	bkgGauss1_mean1_J.setConstant(kTRUE);
//	bkgGauss1_mean1_J.setError(readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_J),1));
	RooRealVar bkgGauss1_mean2_J("bkgGauss1_mean2_J","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_J),0));
	bkgGauss1_mean2_J.setConstant(kTRUE);
//	bkgGauss1_mean2_J.setError(readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_J),1));
	RooRealVar bkgGauss1_sigma1_J("bkgGauss1_sigma1_J","#sigma_{11}",readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_J),0));
	bkgGauss1_sigma1_J.setConstant(kTRUE);
//	bkgGauss1_sigma1_J.setError(readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_J),1));
	RooRealVar bkgGauss1_sigma2_J("bkgGauss1_sigma2_J","#sigma_{12}",readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_J),0));
	bkgGauss1_sigma2_J.setConstant(kTRUE);
//	bkgGauss1_sigma2_J.setError(readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_J),1));
	RooRealVar bkgM_frac1_J("bkgM_frac1_J","bkgM_frac1",readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_J),0));
	bkgM_frac1_J.setConstant(kTRUE);
//	bkgM_frac1_J.setError(readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_J),1));
	
	RooRealVar bkgGauss2_mean1_J("bkgGauss2_mean1_J","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_J),0));
	bkgGauss2_mean1_J.setConstant(kTRUE);
//	bkgGauss2_mean1_J.setError(readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_J),1));
	RooRealVar bkgGauss2_mean2_J("bkgGauss2_mean2_J","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_J),0));
	bkgGauss2_mean2_J.setConstant(kTRUE);
//	bkgGauss2_mean2_J.setError(readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_J),1));
	RooRealVar bkgGauss2_sigma1_J("bkgGauss2_sigma1_J","#sigma_{21}",readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_J),0));
	bkgGauss2_sigma1_J.setConstant(kTRUE);
//	bkgGauss2_sigma1_J.setError(readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_J),1));
	RooRealVar bkgGauss2_sigma2_J("bkgGauss2_sigma2_J","#sigma_{22}",readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_J),0));
	bkgGauss2_sigma2_J.setConstant(kTRUE);
//	bkgGauss2_sigma2_J.setError(readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_J),1));
	RooRealVar bkgM_frac2_J("bkgM_frac2_J","bkgM_frac2",readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_J),0));
	bkgM_frac2_J.setConstant(kTRUE);
//	bkgM_frac2_J.setError(readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_J),1));
	
	RooRealVar bkgM_frac12_J("bkgM_frac12_J","bkgM_frac12",readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_J),0));
	bkgM_frac12_J.setConstant(kTRUE);
//	bkgM_frac12_J.setError(readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_J),1));
	
	RooGaussian f_bkgPeakMGauss11_J("f_bkgPeakMGauss11_J","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean1_J, bkgGauss1_sigma1_J);
	RooGaussian f_bkgPeakMGauss12_J("f_bkgPeakMGauss12_J","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean2_J, bkgGauss1_sigma2_J);
	RooGaussian f_bkgPeakMGauss21_J("f_bkgPeakMGauss21_J","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean1_J, bkgGauss2_sigma1_J);
	RooGaussian f_bkgPeakMGauss22_J("f_bkgPeakMGauss22_J","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean2_J, bkgGauss2_sigma2_J);
	RooAddPdf f_bkgPeakM1_J("f_bkgPeakM1_J","f_bkgPeakM1", RooArgList(f_bkgPeakMGauss11_J, f_bkgPeakMGauss12_J), bkgM_frac1_J);
	RooAddPdf f_bkgPeakM2_J("f_bkgPeakM2_J","f_bkgPeakM2", RooArgList(f_bkgPeakMGauss21_J, f_bkgPeakMGauss22_J), bkgM_frac2_J);
	RooAddPdf f_bkgPeakM12_J("f_bkgPeakM12_J","f_bkgPeakM12", RooArgList(f_bkgPeakM1_J, f_bkgPeakM2_J), bkgM_frac12_J);
	
	RooRealVar bkgPeakL_c0_J("bkgPeakL_c0_J","c0",readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_J),0));
	bkgPeakL_c0_J.setConstant(kTRUE);
//	bkgPeakL_c0_J.setError(readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_J),1));
	RooRealVar bkgPeakL_c1_J("bkgPeakL_c1_J","c1",readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_J),0));
	bkgPeakL_c1_J.setConstant(kTRUE);
//	bkgPeakL_c1_J.setError(readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_J),1));
	RooRealVar bkgPeakL_c2_J("bkgPeakL_c2_J","c2",readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_J),0));
	bkgPeakL_c2_J.setConstant(kTRUE);
//	bkgPeakL_c2_J.setError(readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_J),1));
	RooRealVar bkgPeakL_c3_J("bkgPeakL_c3_J","c3",readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_J),0));
	bkgPeakL_c3_J.setConstant(kTRUE);
//	bkgPeakL_c3_J.setError(readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_J),1));
//	RooRealVar bkgPeakL_c4_J("bkgPeakL_c4_J","c4",readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_J),0));
//	bkgPeakL_c4_J.setError(readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_J),1));
//	RooArgSet f_bkgPeakL_argset_J(bkgPeakL_c1_J,bkgPeakL_c2_J,bkgPeakL_c3_J,bkgPeakL_c4_J,bkgPeakL_c0_J);
	RooArgSet f_bkgPeakL_argset_J(bkgPeakL_c1_J,bkgPeakL_c2_J,bkgPeakL_c3_J,bkgPeakL_c0_J);

	// // PsiK MC Peaking Bkg.
	RooRealVar bkgGauss1_mean1_P("bkgGauss1_mean1_P","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_P),0));
	bkgGauss1_mean1_P.setConstant(kTRUE);
//	bkgGauss1_mean1_P.setError(readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_P),1));
	RooRealVar bkgGauss1_mean2_P("bkgGauss1_mean2_P","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_P),0));
	bkgGauss1_mean2_P.setConstant(kTRUE);
//	bkgGauss1_mean2_P.setError(readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_P),1));
	RooRealVar bkgGauss1_sigma1_P("bkgGauss1_sigma1_P","#sigma_{11}",readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_P),0));
	bkgGauss1_sigma1_P.setConstant(kTRUE);
//	bkgGauss1_sigma1_P.setError(readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_P),1));
	RooRealVar bkgGauss1_sigma2_P("bkgGauss1_sigma2_P","#sigma_{12}",readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_P),0));
	bkgGauss1_sigma2_P.setConstant(kTRUE);
//	bkgGauss1_sigma2_P.setError(readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_P),1));
	RooRealVar bkgM_frac1_P("bkgM_frac1_P","bkgM_frac1",readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_P),0));
	bkgM_frac1_P.setConstant(kTRUE);
//	bkgM_frac1_P.setError(readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_P),1));
	
	RooRealVar bkgGauss2_mean1_P("bkgGauss2_mean1_P","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_P),0));
	bkgGauss2_mean1_P.setConstant(kTRUE);
//	bkgGauss2_mean1_P.setError(readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_P),1));
	RooRealVar bkgGauss2_mean2_P("bkgGauss2_mean2_P","M_{K+#Mu#Mu}",readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_P),0));
	bkgGauss2_mean2_P.setConstant(kTRUE);
//	bkgGauss2_mean2_P.setError(readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_P),1));
	RooRealVar bkgGauss2_sigma1_P("bkgGauss2_sigma1_P","#sigma_{21}",readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_P),0));
	bkgGauss2_sigma1_P.setConstant(kTRUE);
//	bkgGauss2_sigma1_P.setError(readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_P),1));
	RooRealVar bkgGauss2_sigma2_P("bkgGauss2_sigma2_P","#sigma_{22}",readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_P),0));
	bkgGauss2_sigma2_P.setConstant(kTRUE);
//	bkgGauss2_sigma2_P.setError(readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_P),1));
	RooRealVar bkgM_frac2_P("bkgM_frac2_P","bkgM_frac2",readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_P),0));
	bkgM_frac2_P.setConstant(kTRUE);
//	bkgM_frac2_P.setError(readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_P),1));
	
	RooRealVar bkgM_frac12_P("bkgM_frac12_P","bkgM_frac12",readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_P),0));
	bkgM_frac12_P.setConstant(kTRUE);
//	bkgM_frac12_P.setError(readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_P),1));
	
	RooGaussian f_bkgPeakMGauss11_P("f_bkgPeakMGauss11_P","f_bkgPeakMGauss11", Bmass, bkgGauss1_mean1_P, bkgGauss1_sigma1_P);
	RooGaussian f_bkgPeakMGauss12_P("f_bkgPeakMGauss12_P","f_bkgPeakMGauss12", Bmass, bkgGauss1_mean2_P, bkgGauss1_sigma2_P);
	RooGaussian f_bkgPeakMGauss21_P("f_bkgPeakMGauss21_P","f_bkgPeakMGauss21", Bmass, bkgGauss2_mean1_P, bkgGauss2_sigma1_P);
	RooGaussian f_bkgPeakMGauss22_P("f_bkgPeakMGauss22_P","f_bkgPeakMGauss22", Bmass, bkgGauss2_mean2_P, bkgGauss2_sigma2_P);
	RooAddPdf f_bkgPeakM1_P("f_bkgPeakM1_P","f_bkgPeakM1", RooArgList(f_bkgPeakMGauss11_P, f_bkgPeakMGauss12_P), bkgM_frac1_P);
	RooAddPdf f_bkgPeakM2_P("f_bkgPeakM2_P","f_bkgPeakM2", RooArgList(f_bkgPeakMGauss21_P, f_bkgPeakMGauss22_P), bkgM_frac2_P);
	RooAddPdf f_bkgPeakM12_P("f_bkgPeakM12_P","f_bkgPeakM12", RooArgList(f_bkgPeakM1_P, f_bkgPeakM2_P), bkgM_frac12_P);

	RooRealVar bkgPeakL_c0_P("bkgPeakL_c0_P","c0",readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_P),0));
	bkgPeakL_c0_P.setConstant(kTRUE);
//	bkgPeakL_c0_P.setError(readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_P),1));
	RooRealVar bkgPeakL_c1_P("bkgPeakL_c1_P","c1",readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_P),0));
	bkgPeakL_c1_P.setConstant(kTRUE);
//	bkgPeakL_c1_P.setError(readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_P),1));
	RooRealVar bkgPeakL_c2_P("bkgPeakL_c2_P","c2",readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_P),0));
	bkgPeakL_c2_P.setConstant(kTRUE);
//	bkgPeakL_c2_P.setError(readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_P),1));
	RooRealVar bkgPeakL_c3_P("bkgPeakL_c3_P","c3",readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_P),0));
	bkgPeakL_c3_P.setConstant(kTRUE);
//	bkgPeakL_c3_P.setError(readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_P),1));
//	RooRealVar bkgPeakL_c4_P("bkgPeakL_c4_P","c4",readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_P),0));
//	bkgPeakL_c4_P.setError(readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_P),1));
//	RooArgSet f_bkgPeakL_argset_P(bkgPeakL_c1_P,bkgPeakL_c2_P,bkgPeakL_c3_P,bkgPeakL_c4_P,bkgPeakL_c0_P);
	RooArgSet f_bkgPeakL_argset_P(bkgPeakL_c1_P,bkgPeakL_c2_P,bkgPeakL_c3_P,bkgPeakL_c0_P);
	
	switch (iBin) {// Should be fixed constants already.
		default:
		//		bkgPeakL_c0_J.setConstant(kTRUE);
		//		bkgPeakL_c1_J.setConstant(kTRUE);
		//		bkgPeakL_c2_J.setConstant(kTRUE);
		//		bkgPeakL_c3_J.setConstant(kTRUE);
			//	bkgPeakL_c4_J.setConstant(kTRUE);
		//		bkgPeakL_c0_P.setConstant(kTRUE);
		//		bkgPeakL_c1_P.setConstant(kTRUE);
		//		bkgPeakL_c2_P.setConstant(kTRUE);
		//		bkgPeakL_c3_P.setConstant(kTRUE);
			//	bkgPeakL_c4_P.setConstant(kTRUE);
		case 2:
            //1 double guassian ,4+4 deg. ploy
				bkgM_frac12_J.setVal(1.);
				bkgM_frac12_J.setConstant(kTRUE);
				bkgM_frac12_P.setVal(1.);
				bkgM_frac12_P.setConstant(kTRUE);
				break;
		case 4:
            //2 double guassian ,4+4 deg. ploy
				break;
		case 6:
            //1 guassian ,2+2 deg. ploy
				bkgM_frac12_J.setVal(1.);
				bkgM_frac12_J.setConstant(kTRUE);
				bkgM_frac12_P.setVal(1.);
				bkgM_frac12_P.setConstant(kTRUE);
				break;
		case 10:
            //2 double guassian ,4+4 deg. ploy
				break;
	}
	RooPolynomial f_bkgPeakL_J("f_bkgPeakL_J","f_bkgPeakL_J",CosThetaL,f_bkgPeakL_argset_J);
	RooPolynomial f_bkgPeakL_P("f_bkgPeakL_P","f_bkgPeakL_P",CosThetaL,f_bkgPeakL_argset_P);
	RooProdPdf f_bkgPeak_J("f_bkgPeak_J", "f_bkgPeak_J",f_bkgPeakL_J,f_bkgPeakM12_J);
	RooProdPdf f_bkgPeak_P("f_bkgPeak_P", "f_bkgPeak_P",f_bkgPeakL_P,f_bkgPeakM12_P);
	cout<<">>>>>>>>>>>>>>>>> INFO: f_bkgPeak prepared. <<<<<<<<<<<<<<<<<<<"<<endl;
	
	// Observed spectrum = model*fullEfficiency
	RooRealVar nsig("nsig","nsig",20,0,6E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",20,0,1E4);
//	nbkgComb.setConstant(kTRUE);
	RooRealVar nbkgPeak_J("nbkgPeak_J","nbkgPeak_J", (readParam(iBin,TString::Format("nbkgPeak_%s",read1b_J),0) * Lumi_Data / Lumi_Jpsi / Lumi_Scale) );
//	nbkgPeak_J.setError(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_J),1) * Lumi_Data / Lumi_Jpsi / Lumi_Scale);
	nbkgPeak_J.setConstant(kTRUE);
	RooRealVar nbkgPeak_P("nbkgPeak_P","nbkgPeak_P", (readParam(iBin,TString::Format("nbkgPeak_%s",read1b_P),0) * Lumi_Data / Lumi_Psi / Lumi_Scale) );
//	nbkgPeak_P.setError(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_P),1) * Lumi_Data / Lumi_Psi / Lumi_Scale);
	nbkgPeak_P.setConstant(kTRUE);
	if (iBin !=10 && ( iBin == 0 || iBin %2 == 1 || iBin/2 > 3) ){
		nbkgPeak_J.setMin(NULL,0.);
		nbkgPeak_J.setVal(0.);
		nbkgPeak_J.setConstant(kTRUE);
		nbkgPeak_P.setMin(NULL,0.);
		nbkgPeak_P.setVal(0.);
		nbkgPeak_P.setConstant(kTRUE);
	} else if (iBin == 2) {
		nbkgPeak_P.setMin(NULL,0.);
		nbkgPeak_P.setVal(0.);
		nbkgPeak_P.setConstant(kTRUE);
	} else if (iBin == 6) {
		nbkgPeak_J.setMin(NULL,0.);
		nbkgPeak_J.setVal(0.);
		nbkgPeak_J.setConstant(kTRUE);
	}
	
	RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_bkgPeak_J,f_bkgPeak_P,f_sig),RooArgList(nbkgComb,nbkgPeak_J,nbkgPeak_P,nsig));// no penalty term
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_penalty NOT prepared. <<<<<<<<<<<<<<<<<<<<"<<endl;
	
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	

///////////////////////////////////////////////////////////// Gaussian constraints ///////////////////////////////////////////////////	
	// Gaussian constraints
	RooGaussian gaus_sigGauss1_sigma("gaus_sigGauss1_sigma","gaus_sigGauss1_sigma",sigGauss1_sigma,RooConst(readParam(iBin,"sigGauss1_sigma",0)),RooConst(readParam(iBin,"sigGauss1_sigma",1)));
	RooGaussian gaus_sigGauss2_sigma("gaus_sigGauss2_sigma","gaus_sigGauss2_sigma",sigGauss2_sigma,RooConst(readParam(iBin,"sigGauss2_sigma",0)),RooConst(readParam(iBin,"sigGauss2_sigma",1)));
	RooGaussian gaus_sigM_frac("gaus_sigM_frac","gaus_sigM_frac",sigM_frac,RooConst(readParam(iBin,"sigM_frac",0)),RooConst(readParam(iBin,"sigM_frac",1)));
	
	// bkgPeak_J    //////  Lumi. Scale
//	RooGaussian gaus_nbkgPeak_J("gaus_nbkgPeak_J","gaus_nbkgPeak_J",nbkgPeak_J,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_J),0) ),RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_J),1) ));
	RooGaussian gaus_nbkgPeak_J("gaus_nbkgPeak_J","gaus_nbkgPeak_J",nbkgPeak_J,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_J),0) * Lumi_Data / Lumi_Jpsi /Lumi_Scale ),RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_J),1) *  Lumi_Data / Lumi_Jpsi /Lumi_Scale  ));
	
	RooGaussian gaus_bkgGauss1_mean1_J("gaus_bkgGauss1_mean1_J","gaus_bkgGauss1_mean1_J",bkgGauss1_mean1_J,RooConst(readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_J),1)));
	RooGaussian gaus_bkgGauss1_mean2_J("gaus_bkgGauss1_mean2_J","gaus_bkgGauss1_mean2_J",bkgGauss1_mean2_J,RooConst(readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_J),1)));
	RooGaussian gaus_bkgGauss1_sigma1_J("gaus_bkgGauss1_sigma1_J","gaus_bkgGauss1_sigma1_J",bkgGauss1_sigma1_J,RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_J),1)));
	RooGaussian gaus_bkgGauss1_sigma2_J("gaus_bkgGauss1_sigma2_J","gaus_bkgGauss1_sigma2_J",bkgGauss1_sigma2_J,RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_J),1)));
	RooGaussian gaus_bkgM_frac1_J("gaus_bkgM_frac1_J","gaus_bkgM_frac1_J",bkgM_frac1_J,RooConst(readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_J),1)));
	RooGaussian gaus_bkgGauss2_mean1_J("gaus_bkgGauss2_mean1_J","gaus_bkgGauss2_mean1_J",bkgGauss2_mean1_J,RooConst(readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_J),1)));
	RooGaussian gaus_bkgGauss2_mean2_J("gaus_bkgGauss2_mean2_J","gaus_bkgGauss2_mean2_J",bkgGauss2_mean2_J,RooConst(readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_J),1)));
	RooGaussian gaus_bkgGauss2_sigma1_J("gaus_bkgGauss2_sigma1_J","gaus_bkgGauss2_sigma1_J",bkgGauss2_sigma1_J,RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_J),1)));
	RooGaussian gaus_bkgGauss2_sigma2_J("gaus_bkgGauss2_sigma2_J","gaus_bkgGauss2_sigma2_J",bkgGauss2_sigma2_J,RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_J),1)));
	RooGaussian gaus_bkgM_frac2_J("gaus_bkgM_frac2_J","gaus_bkgM_frac2_J",bkgM_frac2_J,RooConst(readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_J),1)));
	RooGaussian gaus_bkgM_frac12_J("gaus_bkgM_frac12_J","gaus_bkgM_frac12_J",bkgM_frac12_J,RooConst(readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_J),0)),RooConst(readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_J),1)));

	RooGaussian gaus_bkgPeakL_c0_J("gaus_bkgPeakL_c0_J","gaus_bkgPeakL_c0_J",bkgPeakL_c0_J,RooConst(readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_J),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_J),1)));
	RooGaussian gaus_bkgPeakL_c1_J("gaus_bkgPeakL_c1_J","gaus_bkgPeakL_c1_J",bkgPeakL_c1_J,RooConst(readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_J),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_J),1)));
	RooGaussian gaus_bkgPeakL_c2_J("gaus_bkgPeakL_c2_J","gaus_bkgPeakL_c2_J",bkgPeakL_c2_J,RooConst(readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_J),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_J),1)));
	RooGaussian gaus_bkgPeakL_c3_J("gaus_bkgPeakL_c3_J","gaus_bkgPeakL_c3_J",bkgPeakL_c3_J,RooConst(readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_J),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_J),1)));
//	RooGaussian gaus_bkgPeakL_c4_J("gaus_bkgPeakL_c4_J","gaus_bkgPeakL_c4_J",bkgPeakL_c4_J,RooConst(readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_J),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_J),1)));
	
	// bkgPeak_P    //////  Lumi. Scale
//	RooGaussian gaus_nbkgPeak_P("gaus_nbkgPeak_P","gaus_nbkgPeak_P",nbkgPeak_P,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_P),0) ),RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_P),1) ));
	RooGaussian gaus_nbkgPeak_P("gaus_nbkgPeak_P","gaus_nbkgPeak_P",nbkgPeak_P,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_P),0) * Lumi_Data / Lumi_Psi /Lumi_Scale ),RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b_P),1) * Lumi_Data / Lumi_Psi /Lumi_Scale  ));
	
	RooGaussian gaus_bkgGauss1_mean1_P("gaus_bkgGauss1_mean1_P","gaus_bkgGauss1_mean1_P",bkgGauss1_mean1_P,RooConst(readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_mean1_%s",read1b_P),1)));
	RooGaussian gaus_bkgGauss1_mean2_P("gaus_bkgGauss1_mean2_P","gaus_bkgGauss1_mean2_P",bkgGauss1_mean2_P,RooConst(readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_mean2_%s",read1b_P),1)));
	RooGaussian gaus_bkgGauss1_sigma1_P("gaus_bkgGauss1_sigma1_P","gaus_bkgGauss1_sigma1_P",bkgGauss1_sigma1_P,RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma1_%s",read1b_P),1)));
	RooGaussian gaus_bkgGauss1_sigma2_P("gaus_bkgGauss1_sigma2_P","gaus_bkgGauss1_sigma2_P",bkgGauss1_sigma2_P,RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss1_sigma2_%s",read1b_P),1)));
	RooGaussian gaus_bkgM_frac1_P("gaus_bkgM_frac1_P","gaus_bkgM_frac1_P",bkgM_frac1_P,RooConst(readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgM_frac1_%s",read1b_P),1)));
	RooGaussian gaus_bkgGauss2_mean1_P("gaus_bkgGauss2_mean1_P","gaus_bkgGauss2_mean1_P",bkgGauss2_mean1_P,RooConst(readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_mean1_%s",read1b_P),1)));
	RooGaussian gaus_bkgGauss2_mean2_P("gaus_bkgGauss2_mean2_P","gaus_bkgGauss2_mean2_P",bkgGauss2_mean2_P,RooConst(readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_mean2_%s",read1b_P),1)));
	RooGaussian gaus_bkgGauss2_sigma1_P("gaus_bkgGauss2_sigma1_P","gaus_bkgGauss2_sigma1_P",bkgGauss2_sigma1_P,RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma1_%s",read1b_P),1)));
	RooGaussian gaus_bkgGauss2_sigma2_P("gaus_bkgGauss2_sigma2_P","gaus_bkgGauss2_sigma2_P",bkgGauss2_sigma2_P,RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgGauss2_sigma2_%s",read1b_P),1)));
	RooGaussian gaus_bkgM_frac2_P("gaus_bkgM_frac2_P","gaus_bkgM_frac2_P",bkgM_frac2_P,RooConst(readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgM_frac2_%s",read1b_P),1)));
	RooGaussian gaus_bkgM_frac12_P("gaus_bkgM_frac12_P","gaus_bkgM_frac12_P",bkgM_frac12_P,RooConst(readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_P),0)),RooConst(readParam(iBin,TString::Format("bkgM_frac12_%s",read1b_P),1)));

	RooGaussian gaus_bkgPeakL_c0_P("gaus_bkgPeakL_c0_P","gaus_bkgPeakL_c0_P",bkgPeakL_c0_P,RooConst(readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_P),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c0_%s",read2a_P),1)));
	RooGaussian gaus_bkgPeakL_c1_P("gaus_bkgPeakL_c1_P","gaus_bkgPeakL_c1_P",bkgPeakL_c1_P,RooConst(readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_P),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c1_%s",read2a_P),1)));
	RooGaussian gaus_bkgPeakL_c2_P("gaus_bkgPeakL_c2_P","gaus_bkgPeakL_c2_P",bkgPeakL_c2_P,RooConst(readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_P),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c2_%s",read2a_P),1)));
	RooGaussian gaus_bkgPeakL_c3_P("gaus_bkgPeakL_c3_P","gaus_bkgPeakL_c3_P",bkgPeakL_c3_P,RooConst(readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_P),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c3_%s",read2a_P),1)));
//	RooGaussian gaus_bkgPeakL_c4_P("gaus_bkgPeakL_c4_P","gaus_bkgPeakL_c4_P",bkgPeakL_c4_P,RooConst(readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_P),0)),RooConst(readParam(iBin,TString::Format("bkgPeakL_c4_%s",read2a_P),1)));

	RooArgSet gausConstraints(gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac);
	switch (iBin) {
		case 2:
            //1 double guassian ,4+4 deg. ploy
			//	gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_J,gaus_bkgPeakL_c2_J,gaus_bkgPeakL_c3_J,gaus_bkgPeakL_c4_J,gaus_bkgPeakL_c0_J));
				gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_J,gaus_bkgPeakL_c2_J,gaus_bkgPeakL_c3_J,gaus_bkgPeakL_c0_J));
				gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1_J,gaus_bkgGauss1_mean2_J,gaus_bkgGauss1_sigma1_J,gaus_bkgGauss1_sigma2_J,gaus_bkgM_frac1_J));
				gausConstraints.add(gaus_nbkgPeak_J);
			//	gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_P,gaus_bkgPeakL_c2_P,gaus_bkgPeakL_c3_P,gaus_bkgPeakL_c4_P,gaus_bkgPeakL_c0_P));
				gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_P,gaus_bkgPeakL_c2_P,gaus_bkgPeakL_c3_P,gaus_bkgPeakL_c0_P));
				gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1_P,gaus_bkgGauss1_mean2_P,gaus_bkgGauss1_sigma1_P,gaus_bkgGauss1_sigma2_P,gaus_bkgM_frac1_P));
				gausConstraints.add(gaus_nbkgPeak_P);
				break;
		
		case 10:
		case 4:
            //2 double guassian ,4+4 deg. ploy
			//	gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_J,gaus_bkgPeakL_c2_J,gaus_bkgPeakL_c3_J,gaus_bkgPeakL_c4_J,gaus_bkgPeakL_c0_J));
				gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_J,gaus_bkgPeakL_c2_J,gaus_bkgPeakL_c3_J,gaus_bkgPeakL_c0_J));
				gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1_J,gaus_bkgGauss1_mean2_J,gaus_bkgGauss1_sigma1_J,gaus_bkgGauss1_sigma2_J,gaus_bkgM_frac1_J));
				gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1_J,gaus_bkgGauss2_mean2_J,gaus_bkgGauss2_sigma1_J,gaus_bkgGauss2_sigma2_J,gaus_bkgM_frac2_J));
				gausConstraints.add(gaus_bkgM_frac12_J);
				gausConstraints.add(gaus_nbkgPeak_J);
			//	gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_P,gaus_bkgPeakL_c2_P,gaus_bkgPeakL_c3_P,gaus_bkgPeakL_c4_P,gaus_bkgPeakL_c0_P));
				gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_P,gaus_bkgPeakL_c2_P,gaus_bkgPeakL_c3_P,gaus_bkgPeakL_c0_P));
				gausConstraints.add(RooArgSet(gaus_bkgGauss1_mean1_P,gaus_bkgGauss1_mean2_P,gaus_bkgGauss1_sigma1_P,gaus_bkgGauss1_sigma2_P,gaus_bkgM_frac1_P));
				gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1_P,gaus_bkgGauss2_mean2_P,gaus_bkgGauss2_sigma1_P,gaus_bkgGauss2_sigma2_P,gaus_bkgM_frac2_P));
				gausConstraints.add(gaus_bkgM_frac12_P);
				gausConstraints.add(gaus_nbkgPeak_P);
				break;
		case 6:
            //1 guassian ,2+2 deg. ploy
			//	gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_J,gaus_bkgPeakL_c2_J,gaus_bkgPeakL_c3_J,gaus_bkgPeakL_c4_J,gaus_bkgPeakL_c0_J));
				gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_J,gaus_bkgPeakL_c2_J,gaus_bkgPeakL_c3_J,gaus_bkgPeakL_c0_J));
				gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1_J,gaus_bkgGauss2_mean2_J,gaus_bkgGauss2_sigma1_J,gaus_bkgGauss2_sigma2_J,gaus_bkgM_frac2_J));
				gausConstraints.add(gaus_nbkgPeak_J);
			//	gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_P,gaus_bkgPeakL_c2_P,gaus_bkgPeakL_c3_P,gaus_bkgPeakL_c4_P,gaus_bkgPeakL_c0_P));
				gausConstraints.add(RooArgSet(gaus_bkgPeakL_c1_P,gaus_bkgPeakL_c2_P,gaus_bkgPeakL_c3_P,gaus_bkgPeakL_c0_P));
				gausConstraints.add(RooArgSet(gaus_bkgGauss2_mean1_P,gaus_bkgGauss2_mean2_P,gaus_bkgGauss2_sigma1_P,gaus_bkgGauss2_sigma2_P,gaus_bkgM_frac2_P));
				gausConstraints.add(gaus_nbkgPeak_P);
				break;
	}
///////////////////////////////////////////////////////////// Gaussian constraints ///////////////////////////////////////////////////	
	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaL,Q2),Q2range[iBin],0);
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"), Minos(kTRUE));
/*
/////////////////////////////////////////////////////////////////// strategy /////////////////////////////
	RooAbsReal* nll = f.createNLL(*data,Extended(kTRUE)) ;
	RooMinuit m(*nll) ;
	m.migrad() ;
	RooFitResult* r = m.save() ;
	RooFitResult* f_fitresult;
	if (r->status() == 0) {
		r->Print();
		delete nll;
		double par1 = afb.getVal();
		double par2 = fh.getVal();
		cout <<">>>>> afb = "<<par1<<">>>>> fh = "<<par2<<endl;
		cout <<"*** STATUS ** REFITTING 1 ***"<<endl;
		afb.setVal(par1);
		fh.setVal(par2);
		RooFitResult *r_1 = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Minos(kTRUE),Strategy(1),Warnings(-1), PrintEvalErrors(-1));
	//	RooAbsReal* nll_1 = f.createNLL(*data,Extended(kTRUE)) ;
	//	RooMinuit m_1(*nll_1);
	//	m_1.migrad() ;
	//	m_1.minos() ;
	//	RooFitResult* r_1 = m_1.save();
		if (r_1->status() == 0) {
			r_1->Print();
		//	delete nll_1;
			double par1_1 = afb.getVal();
			double par2_1 = fh.getVal();
			cout <<">>>>> afb_1 = "<<par1_1<<">>>>> fh_1 = "<<par2_1<<endl;
			cout <<"*** STATUS ** REFITTING 2 ***"<<endl;
			afb.setVal(par1_1);
			fh.setVal(par2_1);
			RooFitResult *r_2 = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Minos(kTRUE),Strategy(1),Warnings(-1), PrintEvalErrors(-1));
		//	RooAbsReal* nll_2 = f.createNLL(*data,Extended(kTRUE)) ;
		//	RooMinuit m_2(*nll_2);
		//	m_2.migrad() ;
		//	m_2.minos() ;
		//	RooFitResult* r_2 = m_2.save();
			if (r_2->status() == 0) {
				r_2->Print();
			//	delete nll_2;
				double par1_2 = afb.getVal();
				double par2_2 = fh.getVal();
				cout <<">>>>> afb_2 = "<<par1_2<<">>>>> fh_2 = "<<par2_2<<endl;
				cout <<"*** STATUS ** REFITTING 3 ***"<<endl;
				afb.setVal(par1_2);
				fh.setVal(par2_2);
				RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Minos(kTRUE),Strategy(1),Warnings(-1), PrintEvalErrors(-1));
			//	RooAbsReal* nll_3 = f.createNLL(*data,Extended(kTRUE)) ;
			//	RooMinuit m_3(*nll_3);
			//	m_3.migrad() ;
			//	m_3.minos() ;
			//	RooFitResult* f_fitresult = m_3.save();
				if (f_fitresult->status() != 0) {
					std::vector<double> output;
					output.push_back(fh.getVal());
					output.push_back(fh.getError());
					return output;
				} else { 
					f_fitresult->Print();
				}				
			} else {
				std::vector<double> output;
				output.push_back(fh.getVal());
				output.push_back(fh.getError());
				return output;
			}
		} else {
			std::vector<double> output;
			output.push_back(fh.getVal());
			output.push_back(fh.getError());
			return output;
		}
	} else {
		std::vector<double> output;
		output.push_back(fh.getVal());
		output.push_back(fh.getError());
		return output;
	}
//	RooFitResult* f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(istra));
//	RooFitResult* f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(istra));
/////////////////////////////////////////////////////////////////// strategy ////////////////////////////////////////
*/
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Strategy(2),Warnings(-1), PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Minos(kTRUE),Strategy(2),Warnings(-1), PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Minos(kTRUE),Strategy(1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Minos(kTRUE),Strategy(2));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) {
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	}
	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
	RooPlot* framemass = Bmass.frame();
	data->plotOn(framemass,RooFit::Name("data"), Binning(20)); 
	f.plotOn(framemass,RooFit::Name("pdf"), LineColor(1)); 
	//f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
	f.plotOn(framemass, Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
	f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(4),LineWidth(2));
	f.plotOn(framemass,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
	//f.plotOn(framemass,RooFit::Name("bkgPeak"), Components(f_bkgPeak),FillStyle(3004),FillColor(3),VLines(), DrawOption("F"));
	f.plotOn(framemass, Components(f_bkgPeak_J),FillStyle(3004),FillColor(3),VLines(), DrawOption("F"));
	f.plotOn(framemass,RooFit::Name("bkgPeak_J"), Components(f_bkgPeak_J),LineColor(3),LineWidth(2),LineStyle(2));
	f.plotOn(framemass, Components(f_bkgPeak_P),FillStyle(3001),FillColor(6),VLines(), DrawOption("F"));
	f.plotOn(framemass,RooFit::Name("bkgPeak_P"), Components(f_bkgPeak_P),LineColor(6),LineWidth(2),LineStyle(2));
	
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
	leg->AddEntry("bkgPeak_J"," Jpsi Peak. bkg. ","L ");
	leg->AddEntry("bkgPeak_P"," Psi Peak. bkg. ","L ");
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
	paveText->AddText(Form(" nbkgPeak_P   = %.0f #pm %.0f ", nbkgPeak_P.getVal(), nbkgPeak_P.getError())); 
	paveText->AddText(Form(" nbkgPeak_J   = %.0f #pm %.0f ", nbkgPeak_J.getVal(), nbkgPeak_J.getError())); 
//	paveText->AddText(Form(" sigmean   = %.3f #pm %.3f ", sigGauss_mean.getVal(), sigGauss_mean.getError())); 
	paveText->AddText(Form("#Chi^{2}_{Bmass}   = %.2f  ", framemass->chiSquare())); 
	paveText->AddText(Form("Fit Status   = %d  ", f_fitresult->status())); 
	paveText->Draw(); 
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	double fixNDC = 0.05;
	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
	c->Update();
	//c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_bin%d_Index_%d.png",outfile,iBin,Index));
	
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
	f.plotOn(framecosl, Components(f_bkgPeak_J),FillStyle(3004),FillColor(3),VLines(), DrawOption("F"));
	f.plotOn(framecosl,RooFit::Name("bkgPeak_J"), Components(f_bkgPeak_J),LineColor(3),LineWidth(2),LineStyle(2));
	f.plotOn(framecosl, Components(f_bkgPeak_P),FillStyle(3001),FillColor(6),VLines(), DrawOption("F"));
	f.plotOn(framecosl,RooFit::Name("bkgPeak_P"), Components(f_bkgPeak_P),LineColor(6),LineWidth(2),LineStyle(2));
	
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
	leg->Draw();
	
	if (iBin == 10) { 
		t1->DrawLatex(.35,.86+fixNDC,TString::Format(" Total Signal Region"));
	} else t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",Q2range[iBin]));
//	t1->DrawLatex(.10,.81+fixNDC,TString::Format("F_{H}  =%9.5f #pm%9.5f",fh.getVal(),fh.getError()));
//	t1->DrawLatex(.10,.76+fixNDC,TString::Format("A_{FB} =%9.5f #pm%9.5f",afb.getVal(),afb.getError()));
	paveText->AddText(Form("#Chi^{2}_{cos#theta}   = %.2f  ", framecosl->chiSquare())); 
	paveText->Draw(); 
	c->Update();
	//c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_cosl_bin%d_Index_%d.png",outfile,iBin,Index));
	delete c;
	delete t1;
	delete data;
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
		val[0] = fh.getVal();val[1] = fh.getError();
		writeParam(iBin, "fh", val);
		val[1]=0; val[2]=0;
		val[0] = afb.getVal();val[1] = afb.getError();
		writeParam(iBin, "afb",val);
	} else if (Index == -2) {
		cout<<endl<<endl<<"This is for a data fitting test!!!"<<endl;
		cout<<"Fit Status  = "<<f_fitresult->status()<<endl;
		cout<<"Afb = "<<afb.getVal()<<" +- "<<afb.getError()<<endl;
		cout<<"Fh  = "<<fh.getVal()<<" +- "<<fh.getError()<<endl;
	} else {
		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeOutput(iBin, Index, "afb", val);
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeOutput(iBin, Index, "fh", val);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeOutput(iBin, Index, "FCN", val);
	}
	
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	
}//}}}

void PlotFCN( int iBin, const char outfile[] = "FCN")
{
	setTDRStyle();
	TGraph *gr_afb = new TGraph();
	TGraph *gr_fh  = new TGraph();
	double FCN = 999, fcn = 999;
	int Index = 0, index = 0, n = -1, NIndex = 0;
	double Iafb, Ifh, afb, afberr, fh, fherr;
	if (iBin == 0 || iBin == 6) NIndex = 200;
	else if (iBin == 8 || iBin == 10) NIndex = 600;
	else if (iBin == 7) NIndex = 1000;
	else NIndex = 400;
	do {
		index+=1;
		while ( ! (fcn = readOutput(iBin, index, "FCN", 0)) && index < NIndex) index+=1;
//		cout<<"fcn = "<<fcn<<endl;
		n+=1;
		afb = readOutput(iBin, index, "afb", 1);
		afberr = readOutput(iBin, index, "afb", 2);
		fh  = readOutput(iBin, index, "fh", 1);
		fherr = readOutput(iBin, index, "fh", 2);
		if (fcn != 0) {
		//	cout<<n<<endl;
			gr_afb->SetPoint(n, afb, fcn);
			gr_fh->SetPoint(n, fh, fcn);
		}
		if ( (FCN > fcn) && (afberr > 0.001) && (fherr > 0.001) ) { FCN = fcn; Index = index; }
	} while ( index < NIndex);
	Iafb   = readOutput(iBin, Index, "afb", 0);
	afb    = readOutput(iBin, Index, "afb", 1);
	afberr = readOutput(iBin, Index, "afb", 2);
	Ifh   = readOutput(iBin, Index, "fh", 0);
	fh    = readOutput(iBin, Index, "fh", 1);
	fherr = readOutput(iBin, Index, "fh", 2);
	TCanvas *c = new TCanvas("c","c",800,600);
//	c->SetTitle("FCN distribution of A_{FB}");
	gr_afb->GetXaxis()->SetTitle("A_{FB}");
	gr_afb->GetYaxis()->SetTitle("NLL");
	gr_afb->Draw("AP");
	c->Print(TString::Format("./plots/%s_afb_bin%d.png",outfile,iBin));
	c->Clear();
//	c->SetTitle("FCN distribution of F_{H}");
	gr_fh->GetXaxis()->SetTitle("F_{H}");
	gr_fh->GetYaxis()->SetTitle("NLL");
	gr_fh->Draw("AP");
	c->Print(TString::Format("./plots/%s_fh_bin%d.png",outfile,iBin));
	c->Clear();
	cout<<"Index = "<<Index<<"   FCN = "<<FCN<<endl;
	cout<<"Iafb  = "<<Iafb<<"   Ifh = "<<Ifh<<endl;
	cout<<"afb = "<<afb<<" +- "<<afberr<<endl;
	cout<<"fh  = "<<fh<<" +- "<<fherr<<endl;
	double val[3]={0,0,0};
	val[0] = Iafb;
	writeParam(iBin, "Iafb",val);
	val[1]=0; val[2]=0;
	val[0] = Ifh;
	writeParam(iBin, "Ifh", val);
	val[1]=0; val[2]=0;
	val[0] = afb;val[1] = afberr;
	writeParam(iBin, "afb",val);
	val[1]=0; val[2]=0;
	val[0] = fh;val[1] = fherr;
	writeParam(iBin, "fh", val);

}

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
	//	yafb[i]      = DATAAfb[ibin];
		yafb[i]      = readParam(ibin,"afb",0);
	//	yuerrafb[i]  = DATAAfberr[ibin];
		yuerrafb[i]  = fabs(readParam(ibin,"afb",1));
		yderrafb[i]  = yuerrafb[i];
	//	yfh[i]       = DATAFh[ibin];
		yfh[i]       = readParam(ibin,"fh",0);
	//	yuerrfh[i]   = DATAFherr[ibin];
		yuerrfh[i]   = fabs(readParam(ibin,"fh",1));
		if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		else { yderrfh[i] = yuerrfh[i]; }
		printf("Afb[%d]=%6.4f +- %6.4f\n",ibin,yafb[i],yuerrafb[i]);
		printf("Fh [%d]=%6.4f +- %6.4f\n",ibin,yfh[i],yuerrfh[i]);
	}
//	Draw
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
	
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
//	frame->SetAxisRange(-0.01,1.,"Y");
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
//	frame->SetAxisRange(-0.2,0.2,"Y");
	frame->SetAxisRange(-1.,1.,"Y");
	frame->Draw();
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
		cout<<"Usage       : ./fit Function infile binID"<<endl;
		cout<<"Functions   :"<<endl;
		cout<<"    0. bmass               Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg."<<endl;
		cout<<"    1. angular_gen         Derive F_{H} and A_{FB} from cosThetaL distribution at GEN level."<<endl;
		cout<<"       acceptance          Get acceptance map from unfiltered signal GEN, |Mu pT| > 2.8 GeV, |Mu eta| < 2.3."<<endl;
		cout<<"       recoEff             Get reconstruction efficiency map from signal simulation."<<endl;
		cout<<"    2. accXrecoEff         Get efficiency map from signal simulation."<<endl;
		cout<<"    3. angular_reco        Derive F_{H} and A_{FB} from cosThetaL distribution at RECO level."<<endl;
		cout<<"    4. angular_gen_R       Derive F_{H} and A_{FB} from gencosThetaL distribution at GEN level with official MC."<<endl;
		cout<<"    5. angular2D_1a_Sm     Leading step1 to angular2D, determine signal shape from simulation."<<endl;
		cout<<"    6. angular2D_1b_YpPm   Leading step2 to angular2D, determine mass spectrum of peaking bkg from simulation."<<endl;
		cout<<"    7. angular2D_2a_PkPl   Leading step3 to angular2D, determine angular dist. of peaking bkg from simulation."<<endl;
		cout<<"    8. angular2D_prior     Leading step4 to angular2D, fit to data sideband to get initial values of combinatorial bkg."<<endl;
		cout<<"    9. angular2D           Derive F_{H} and A_{FB} by fitting to mass and angular distribution."<<endl;
		cout<<"       For data fitting test:        ./fit <function> <input.root> <iBin> <afb> <fh> test "<<endl;
		cout<<"       For initial values scanning:          ./fit <function> <input.root> <iBin>  "<<endl;
		cout<<"   10. PlotFCN             Plot FCN distribution for each q^{2} bin, find the samllest FCN value"<<endl;
		cout<<"                              and it's initial value, save the fitting results."<<endl;
		cout<<"       For refit with final initial values:  ./fit <function> <input.root> <iBin> refit  "<<endl;
		cout<<"Remark      :"<<endl;
		cout<<"    1. Outputs will be stored in ./plots, please keep the directory."<<endl;
		cout<<"    2. Fitted parameters will be stored in ./fitParameters/*.txt, please keep the directory."<<endl;
		cout<<"    3. Scaned fitted parameters will be stored in ./OutputValues/bin*/*.txt, please keep these directories."<<endl;
	//	cout<<"    3. Wildcard is allowed for infile. But you must quote infile like \"inputData_Run*.root\"!"<<endl;
		return 0;
	}
//	main
	if (argc != 4 && argc != 7 && argc != 5){
		cout<<"./fit func infile binID"<<endl;
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
	}else if (func == "angular_gen_R"){
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			std::vector<double> vbin;
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			vbin = angular_gen_R_bin(iBin);
		}else if (iBin == 999) {
			const char outfile[]="angular_gen_R";
			angular_gen_R(outfile);
		}else { 
			cout<<"Refit official gen level,iBin counts from 0 to 10, or 999!"<<endl;
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
	}else if (func == "angular2D_1a_Sm" || func == "angular2D_prior"){
		void (*fx)(int, const char*, bool);
		if ( func == "angular2D_1a_Sm" ){
			fx = angular2D_1a_Sm;
		}else{
			fx = angular2D_prior;
		}	
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		fx(iBin,func,true);// By default overwrite exist parameters.
	}else if ( func == "angular2D_1b_YpPm_Jpsi" || func == "angular2D_1b_YpPm_Psi" || func == "angular2D_2a_PkPl_Jpsi" || func == "angular2D_2a_PkPl_Psi" ) {
		void (*fx)(int, const char*, bool);
		if (func == "angular2D_1b_YpPm_Jpsi" || func == "angular2D_1b_YpPm_Psi") {
			fx = angular2D_1b_YpPm;
		}else if (func == "angular2D_2a_PkPl_Jpsi" || func == "angular2D_2a_PkPl_Psi") {
			fx = angular2D_2a_PkPl;
		}
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		fx(iBin,func,true);// By default overwrite exist parameters.
	/*	for (int iBin = 0; iBin < 11; iBin++) {
			if (iBin == 3 || iBin == 5) continue;     
			fx(iBin,func,true);// By default overwrite exist parameters.
		}
*/	}else if (func == "angular2D"){
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			std::vector<double> vbin;
			
			float Iafb, Ifh;
			int Index = -999;
			if (argc == 7) {
				Iafb  = atof(argv[4]);
				Ifh   = atof(argv[5]);
				TString type  = argv[6];
				if (type == "test") {
					cout<<endl<<"This is for a data fitting test!!!"<<endl;
					Index  = -2;
					cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>>>> Index = "<<Index<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
					vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
				} else {
					cout<<"A test of data fitting with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> <afb> <fh> test "<<endl;
					return 0;
				}
			} else if (argc == 4) {				
				TF1 *f1 = new TF1("f1","gaus",-1,1);
				f1->SetParameters(1,0,0.1);
				int NIndex = 0;
				if (iBin == 0 || iBin == 6) NIndex = 200;
				else if (iBin == 8 || iBin == 10) NIndex = 600;
				else if (iBin == 7) NIndex = 1000;
				else NIndex = 400;
				for (Index = 0; Index < NIndex; Index+=1) {    
					Ifh = fabs(f1->GetRandom());
					Iafb = f1->GetRandom();
					while ( fabs(Iafb) >= 0.5 * Ifh) {
						Iafb = f1->GetRandom();
					}
					cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>>>> Index = "<<Index<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
					vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
				}
			} else if (argc == 5) {
				TString type  = argv[4];
				if ( type == "refit") {
					Iafb = readParam(iBin, "Iafb", 0);
					Ifh  = readParam(iBin, "Ifh", 0);
					Index = -1;
					cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
					vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
				} else {
					cout<<"Refit data with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> refit "<<endl;
					return 0;
				}
			}
		}else if (iBin == 999) {
			const char outfile[]="angular2D";
			angular2D(outfile);
		}else { 
			cout<<"Refit data, iBin counts from 0 to 10; or 999 to plot the results!"<<endl;
			cout<<" Please check the Usage : ./fit"<<endl;
			return 0; 
		}
	} else if (func == "PlotFCN") {
		const char outfile[]="FCN";
		PlotFCN(iBin, outfile);
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
