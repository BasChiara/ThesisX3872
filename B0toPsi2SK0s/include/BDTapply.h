#ifndef BDTapply_h
#define BDTapply_h

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>
 
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

class BDTapply{

public:

	// constructor & destructor
	BDTapply(const TString& infile,const TString& intree, const TString& outfile, const TString& outtree);
	~BDTapply(){};

	int    Apply();
	// getters
	// setters

private:

//...inputs 
	
	TString inFilePath_;
	TFile* InFile_;
	TString inTreeName_;
	TTree* InTree_;
	// BDT model	
	TString BDTweightPath_;

//...outputs
	TString outFilePath_;
	TString outTreeName_;


// TTree
	float run, LumiBlock, event;
	float pTM_B0, SVprob, LxySign_B0, CosAlpha_B0, DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho;
	float M_Rho, M_B0, M_MuMu, M_X3872, M_K0s;
	float X;

};



#endif
