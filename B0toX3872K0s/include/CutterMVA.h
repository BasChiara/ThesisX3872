#ifndef CutterMVA_h
#define CutterMVA_h

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

//#include "RootPlots.h"
 
#include "TFile.h"
#include "TTree.h"
#include "TText.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include <TStyle.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

class CutterMVA{

public:

	// constructor & destructor
	CutterMVA();
	~CutterMVA();

	void SetUpInTreeSGN();
	void SetUpInTreeBKG();
	void LoadFitRanges();

	// getters
	TFile* GetFile(TString file); 
	void GetCorrelation_BDT_MR();
	int  GetSignalPreCut();
	double GetSGNfactor() {return SGNfactor;};
	double GetAvMultiplicity(const double Xcut, const double MRcut, TString dataset );
	void GetControlPlots(const double Xcut, const double MRcut);

	// analysis methods 
	int ApplyBDT();
	double BKG_NevExtraction(const double Xcut, const double MRcut);
	double SGN_NevExtraction(const double Xcut, const double MRcut);
	double PunziSign(const double S, const double B, double* PSerr);
	int    makeSGNvsBKGplot(const double Xcut, const double MRcut);
	int    makeMASSplot2D(const double Xcut, const double MRcut);
	void   Scan_X3872vsB0(const double Xcut, const double MRcut);

	// Fitters
	void  SGNfitter(const double Xcut, const double MRcut);
	void  CMSxxx(TCanvas* c);

private:

	int SGNprecut;
	double SGNfactor;
// I/O file path

	TString outFile_path;

	TString inFileSGN_path, inFileBKG_path;
	TString SGNtreeName, BKGtreeName;
	TString outFileSGN_path, outFileBKG_path;
	TString outPlot_path;
	std::string inFileFitParB0_path, inFileFitParK0s_path, inFileFitParX_path;

// TTree
	float run, LumiBlock, event;
	float pTM_B0, SVprob, LxySign_B0, CosAlpha_B0, DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho;
	float M_Rho, M_B0, M_MuMu, M_X3872, M_K0s, M_K0sPi1, M_K0sPi2;
	float X;
	TFile* InFileSGN;
	TTree* TreeSGN;
	TFile* InFileBKG;
   TTree* TreeBKG;

// BDT model	
	TString BDTweightPath_;

// signal and background edges
	double MX_nearLeft, MX_nearRight, MK0s_nearLeft, MK0s_nearRight;
	double SGNregL, SGNregR;
	double BKGregFarL, BKGregNearL, BKGregFarR, BKGregNearR;
	double Mlow_B0 = 4.9, Mhigh_B0 = 5.6;
		

};



#endif
