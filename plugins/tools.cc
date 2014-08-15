// -----------------------------------------------
// Author: Xin Shi <Xin.Shi@cern.ch>
// Created: [2013-08-15 Thu 15:56]
// -----------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
//#include "tools.h"

using namespace std;

void set_root_style(int stat=1110, int grid=0) {
	gROOT->Reset();
	
	gStyle->SetTitleFillColor(0);
	gStyle->SetTitleBorderSize(0);
	
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetCanvasDefX(0);
	gStyle->SetCanvasDefY(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(1);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameFillStyle(0);
	gStyle->SetFrameLineColor(1);
	gStyle->SetFrameLineStyle(1);
	gStyle->SetFrameLineWidth(1);
	
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(grid);
	gStyle->SetPadGridY(grid);
	
	gStyle->SetOptStat(stat);
	gStyle->SetStatColor(0);
	gStyle->SetStatBorderSize(1);
}

TChain* add_chain(TString datatype, TString label, TString cut, int verbose=0) {
	TChain *globChain = new TChain("tree");
//	TString base = Form("%s/sel", getenv("dat"));
	TString fNameList = Form("./sel_%s_%s_%s_rootfiles.txt", datatype.Data(), label.Data(), cut.Data());
//	if (verbose >0 )
	cout << ">>> Load Chain from file: " << fNameList << endl;
	
	ifstream fList(fNameList.Data());
	if (!fList) {
		cerr << "!!! Can't open file " << fNameList << endl;
		return NULL;
	}
	
	char lineFromFile[255];
	while(fList.getline(lineFromFile, 250)) {
		TString fileName = lineFromFile;
		fileName = Form("%s", fileName.Data());
//		fileName = Form("%s/sel/%s/%s", getenv("dat"), datatype.Data(), fileName.Data());
		
		if(globChain->Add(fileName)) {
			if (verbose >0 )
				cout << ">> File '" << fileName << "' has been loaded" << endl;
		}
		else
			cerr << ">> Can't load file '" << fileName << "'" << endl;
	}
	
	if (verbose >0 )
		cout << ">> Total number of entries: " << globChain->GetEntries() << endl;
	fList.close();
	return globChain;
}

char* get_option(char ** begin, char ** end, const std::string & option) {
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end) return *itr;
	return 0;
}

bool option_exists(char** begin, char** end, const std::string& option) {
	return std::find(begin, end, option) != end;
}

