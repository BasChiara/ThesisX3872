#include "../include/OptimizerMVA.h"
#include "../include/BDTapply.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

int main (int argc, char** argv ){

	// INPUT params	
	if(argc < 3){
		std::cout << "USAGE : ./OptimizeMVA Xcut MRHOcut "<< std::endl; 
		exit(-1);
	}
	const double X_cut = std::stod(argv[1]);
	const double MRho_cut = std::stod(argv[2]);

	const TString inFileMC    = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/SGN_UL17_MC.root";
	const TString inTreeMC    = "mcSIGNAL";
	const TString outFileMC   = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/BDTonSGN_UL17.root";
	const TString outTreeMC   = "TreeBDTx_S";

	BDTapply* BDTonMC   = new BDTapply(inFileMC, inTreeMC, outFileMC, outTreeMC); 
	BDTonMC->Apply();
	delete BDTonMC;

	OptimizerMVA* cutter = new OptimizerMVA();
/*	// pre-cuts
	double Bs_NOcut = cutter->BKG_NevExtraction(-1.1, 0.);
	double Ss_NOcut = cutter->GetSignalPreCut() * cutter->GetSGNfactor();
	cutter->makeSGNvsBKGplot(-1.1, 0);	
	//cutter->makeMASSplot2D(-1.1 , 0);	

	double Bs_cut = cutter->BKG_NevExtraction(X_cut, MRho_cut);
*/
	double Ss_cut = cutter->SGN_NevExtraction(X_cut, MRho_cut);
	std::cout << " #SIGNAL = " << Ss_cut/cutter->GetSGNfactor() << std::endl;
/*
	double PsigErr;
	double Psig = cutter->PunziSign(Ss_cut, Bs_cut, &PsigErr);
	cutter->makeSGNvsBKGplot(X_cut, MRho_cut);	
	//cutter->makeMASSplot2D(X_cut, MRho_cut);	
	//cutter->GetCorrelation_BDT_MR();
	//cutter->GetControlPlots(X_cut, MRho_cut);


	//cutter->FinalFit(X_cut, MRho_cut);

	std::cout << " RESULTS.... " << std::endl;
	std::cout << "  Pre-cut  (B) (S) "  << Bs_NOcut << "\t" << Ss_NOcut <<std::endl;
	std::cout << "  Post-cut (B) (S) "  << Bs_cut   << "\t" << Ss_cut   <<std::endl;
	std::cout << "  BKGrej    SGNeff Punzi-Sign "  << 1. - Bs_cut/Bs_NOcut << "\t" << Ss_cut/Ss_NOcut << "\t" << Psig << " +/- "<< PsigErr<< std::endl;
	
	 
	std::cout << " AVERAGE MULTIPLICITY "<<std::endl; 
	std::cout << "  SGN pre-cut  " << cutter->GetAvMultiplicity(-1.1, 0., "SGN") <<std::endl; 
	std::cout << "  BKG pre-cut  " << cutter->GetAvMultiplicity(-1.1, 0., "BKG") <<std::endl; 
	std::cout << "  SGN post-cut " << cutter->GetAvMultiplicity(X_cut, MRho_cut, "SGN") <<std::endl; 
	std::cout << "  BKG post-cut " << cutter->GetAvMultiplicity(X_cut, MRho_cut, "BKG") <<std::endl; 
	
	cutter->BKGsubtraction(X_cut, MRho_cut);
*/
	delete cutter;
	return 0;
}
