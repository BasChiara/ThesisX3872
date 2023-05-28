#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
 
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TH1F.h"
#include "TH2F.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
 
using namespace TMVA;
 
TString outFile_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAcutApp.root"; 
TString inFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAinputs.root";
TString SGNtreeName = "inputSIGNAL";
TString inFileBKG_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data/merged/BKG_MergData17.root";
TString BKGtreeName = "B0sidebands";


void TMVAClassAppBDT(  )
{
 
   // This loads the library
   TMVA::Tools::Instance();
   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
 
 
   // Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
 
   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t pTM_B0, SVprob, LxySign_B0, CosAlpha_B0;  
   Float_t DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho; 
   reader->AddVariable( "pTM_B0", &pTM_B0);
   reader->AddVariable( "SVprob", &SVprob);
   reader->AddVariable( "LxySign_B0", &LxySign_B0);
   reader->AddVariable( "CosAlpha_B0", &CosAlpha_B0);
   reader->AddVariable( "DR_Pi1B0", &DR_Pi1B0);
   reader->AddVariable( "pT_Pi1", &pT_Pi1);
   reader->AddVariable( "pT_Rho", &pT_Rho);
   reader->AddVariable( "D0_Rho", &D0_Rho);
 
   // Spectator variables declared in the training have to be added to the reader, too
   Float_t M_Rho; 
   reader->AddSpectator( "M_Rho", &M_Rho);
 
   // Book the MVA methods
 
   TString dir    = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/datasetBDT/weights/";
   TString prefix = "TMVAClassification";
 
   // Book BDT method
   TString BDTname = "BDT_nT200_S25_D4_b03_nC35";
	TString methodName = BDTname + TString(" method"); 
	TString weightfile = dir + prefix + TString("_") + BDTname + TString(".weights.xml");
	reader->BookMVA( methodName, weightfile );


   // Book output histograms
   UInt_t nbin = 100;
 	float Xlow = -0.8, Xhigh = 0.4;
   TH1F* h_BDT_SGN       = new TH1F( "MVA_BDT_signal", "", nbin, Xlow, Xhigh);
   TH1F* h_BDT_BKG       = new TH1F( "MVA_BDT_background", "", nbin, Xlow, Xhigh);
	nbin = 75;
	double M_low = 5.0 , M_high = 5.6;
   TH1F* h_SGN_MB0       = new TH1F( "SGN_MB0", "", nbin, M_low, M_high );
   TH1F* h_BKG_MB0       = new TH1F( "BKG_MB0", "", nbin, M_low, M_high);

	TH2F* h_SGN_MB0vsBDTx = new TH2F( "SGN_MB0vsBDTx", "", nbin, Xlow, Xhigh, nbin, M_low, M_high); 
	TH2F* h_BKG_MB0vsBDTx = new TH2F( "BKG_MB0vsBDTx", "", nbin, Xlow, Xhigh, nbin, M_low, M_high); 

	M_low = .4 , M_high = 1.;
   TH1F* h_SGN_MRho      = new TH1F( "SGN_MRho", "", nbin, M_low, M_high );
   TH1F* h_BKG_MRho      = new TH1F( "BKG_MRho", "", nbin, M_low, M_high);
 
   TString fname = inFileSGN_path;
	// Open SIGNAL file
   TFile *input_signal(0);
   if (!gSystem->AccessPathName( fname )) {
      input_signal = TFile::Open( fname ); // check if file in local directory exists
		std::cout << "--- TMVAClassificationApp    : Using input file (SIGNAL): " << input_signal->GetName() << std::endl;
   }
   if (!input_signal) {
      std::cout << "ERROR: could not open signal-data file" << std::endl;
      exit(1);
   }
	// Open BACKGROUND file
   TFile *input_background(0);
	fname = inFileBKG_path;
   if (!gSystem->AccessPathName( fname )) {
      input_background = TFile::Open( fname ); // check if file in local directory exists
		std::cout << "--- TMVAClassificationApp    : Using input file (BACKGROUND): " << input_background->GetName() << std::endl;
   }
   if (!input_background) {
      std::cout << "ERROR: could not open background-data file" << std::endl;
      exit(1);
   }
 
   // Event loop
 
	float M_B0, M_MuMu, M_X3872, M_K0s;
   // ===== SIGNAL TREE =====
   std::cout << "--- Select signal sample" << std::endl;
   TTree* TreeSGN = (TTree*)input_signal->Get(SGNtreeName);

   TreeSGN->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN->SetBranchAddress( "pT_Rho", &pT_Rho);
	TreeSGN->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN->SetBranchAddress( "M_B0", &M_B0);
   TreeSGN->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeSGN->SetBranchAddress( "M_X3872", &M_X3872);
   TreeSGN->SetBranchAddress( "M_K0s", &M_K0s);
	

	// ===== BACKGROUND TREE =====

   std::cout << "--- Select background sample" << std::endl;
   TTree* TreeBKG = (TTree*)input_background->Get(BKGtreeName);

   TreeBKG->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeBKG->SetBranchAddress( "SVprob", &SVprob);
   TreeBKG->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeBKG->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeBKG->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeBKG->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeBKG->SetBranchAddress( "pT_Rho", &pT_Rho);
   TreeBKG->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeBKG->SetBranchAddress( "M_Rho", &M_Rho);
   TreeBKG->SetBranchAddress( "M_B0", &M_B0);
   TreeBKG->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeBKG->SetBranchAddress( "M_X3872", &M_X3872);
   TreeBKG->SetBranchAddress( "M_K0s", &M_K0s);

	
	const double Xcut = -0.271178;
	Double_t X, BKGrej = 0., SGNeff = 0.;

	// ===== OUTPUT TREE ===== // 
	Double_t SGN_x, SGN_MB0, SGN_MRho;
	Double_t BKG_x, BKG_MB0, BKG_MRho;
	TTree *OutTreeSGN = new TTree("OutTreeSGN", "Tree with SIGNAL variables post BDT cut");	 
	OutTreeSGN->Branch("SGN_x", &X);
	OutTreeSGN->Branch("SGN_MB0", &M_B0);
	OutTreeSGN->Branch("SGN_MRho", &M_Rho);
	TTree *OutTreeBKG = new TTree("OutTreeBKG", "Tree with BACKGROUND variables post BDT cut");	 
	OutTreeBKG->Branch("BKG_x", &X);
	OutTreeBKG->Branch("BKG_MB0", &M_B0);
	OutTreeBKG->Branch("BKG_MRho", &M_Rho);

	// ---> SIGNAL evaluation 
   std::cout << "--- Processing: " << TreeSGN->GetEntries() << "signal events" << std::endl;
   TStopwatch sw;
   sw.Start();

	for (Long64_t ievt=0; ievt<TreeSGN->GetEntries();ievt++) {
 
		if (ievt%100 == 0) std::cout << "--- ... Processing signal-event: " << ievt << std::endl;
		TreeSGN->GetEntry(ievt);
		X = reader->EvaluateMVA(methodName);
		h_SGN_MB0vsBDTx->Fill(X,M_B0);
		h_BDT_SGN->Fill( X );
		h_SGN_MRho->Fill(M_Rho);
		if (X > Xcut){
			SGNeff += 1.;
			h_SGN_MB0->Fill(M_B0);

			OutTreeSGN->Fill();
		}
 
   }
	SGNeff /= TreeSGN->GetEntries(); 
	// ---> BACKGROUND evaluation 
   std::cout << "--- Processing: " << TreeBKG->GetEntries() << "background events" << std::endl;

	for (Long64_t ievt=0; ievt<TreeBKG->GetEntries();ievt++) {
 
		if (ievt%1000 == 0) std::cout << "--- ... Processing background-event: " << ievt << std::endl;
		TreeBKG->GetEntry(ievt);
		X = reader->EvaluateMVA(methodName);
		h_BKG_MB0vsBDTx->Fill(X,M_B0);
		h_BDT_BKG->Fill( X );
		h_BKG_MRho->Fill(M_Rho);
		if (X > Xcut){	
			BKGrej+=1.;
			h_BKG_MB0->Fill(M_B0);

			OutTreeBKG->Fill();
		}
   }
	BKGrej = 1. - BKGrej/TreeBKG->GetEntries(); 


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
 
   // Write histograms
 
   TFile *target  = new TFile(outFile_path,"RECREATE" );
   h_BDT_SGN->Write();
   h_BDT_BKG->Write();
	h_SGN_MB0vsBDTx->Write();
	h_BKG_MB0vsBDTx->Write();
	h_SGN_MB0->Write();
	h_BKG_MB0->Write();
	h_BKG_MRho->Write();
	h_SGN_MRho->Write();
	OutTreeSGN->Write();
	OutTreeBKG->Write();
 
   std::cout << "--- Created root file: "<< target->GetName() << " containing the MVA output histograms & TTree" << std::endl;
	target->Close();
 
   delete reader;
 
   std::cout << "==> TMVAClassAppBDT is done!" << std::endl << std::endl;
	std::cout << "SIGNAL EFFICIENCY = " << SGNeff << std::endl;
	std::cout << "BACKGROUND REJECTION = "<< BKGrej << std::endl;
}
 
int main( int argc, char** argv )
{
   TMVAClassAppBDT();
	



   return 0;
}
