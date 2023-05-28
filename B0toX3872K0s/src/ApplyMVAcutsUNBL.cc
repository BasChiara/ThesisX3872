#include "../include/CutterMVA.h"
#include "../include/BDTapply.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

int main (int argc, char** argv ){

	// INPUT params	
	if(argc < 3){
		std::cout << "USAGE : ./ApplyMVAcutsUNBL Xcut MRHOcut "<< std::endl; 
		exit(-1);
	}
	const double X_cut = std::stod(argv[1]);
	const double MRho_cut = std::stod(argv[2]);
	// MC simulation
	const TString inFileMC    = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGN_MC.root";
	const TString inTreeMC    = "mcSIGNAL";
	const TString outFileMC   = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonSGN.root";
	const TString outTreeMC   = "TreeBDTx_S";
	// DATA 
	const TString inFileDATA  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data/merged/UNBL_MergData17_v2.root";
	const TString inTreeDATA  = "B0_X3872signal";
	const TString outFileDATA = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonDATAunbl_v2.root";
	const TString outTreeDATA   = "TreeBDTx_B";

	BDTapply* BDTonMC   = new BDTapply(inFileMC, inTreeMC, outFileMC, outTreeMC); 
	BDTonMC->Apply();
	BDTapply* BDTonDATA = new BDTapply(inFileDATA, inTreeDATA, outFileDATA, outTreeDATA); 
	BDTonDATA->Apply();

	delete BDTonMC;
	delete BDTonDATA;

	CutterMVA* cutter = new CutterMVA();
	// pre-cuts

	double Bs_cut = cutter->BKG_NevExtraction(X_cut, MRho_cut);
	double Ss_cut = cutter->SGN_NevExtraction(X_cut, MRho_cut);

	double PsigErr;
	double Psig = cutter->PunziSign(Ss_cut, Bs_cut, &PsigErr);
	cutter->makeSGNvsBKGplot(X_cut, MRho_cut);	
	//cutter->makeMASSplot2D(X_cut, MRho_cut);	
	std::cout << " RESULTS.... " << std::endl;
	std::cout << "  Post-cut (B) (S) "  << Bs_cut   << "\t" << Ss_cut   <<std::endl;
	//std::cout << "  BKGrej    SGNeff Punzi-Sign "  << 1. - Bs_cut/Bs_NOcut << "\t" << Ss_cut/Ss_NOcut << "\t" << Psig << " +/- "<< PsigErr<< std::endl;
	
	return 0;
}
