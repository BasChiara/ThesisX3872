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
 
TString outFile_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/TMVAcompareMC.root"; 

TString inFileSGN_X3872_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGN_MC.root";
TString inFileSGN_PSI2S_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/SGN_MC.root";
TString SGNtreeName = "mcSIGNAL";
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
   UInt_t nbin = 50;
 	float Xlow = -0.8, Xhigh = 0.4;
   TH1F* h_BDT_SGN_X       = new TH1F( "BDT_X", "", nbin, Xlow, Xhigh);
   TH1F* h_BDT_SGN_Psi     = new TH1F( "BDT_Psi", "", nbin, Xlow, Xhigh);
	nbin = 75;
	double M_low = 5.0 , M_high = 5.6;
   TH1F* h_SGN_MB0_X       = new TH1F( "SGN_MB0_X", "", nbin, M_low, M_high );
   TH1F* h_SGN_MB0_Psi     = new TH1F( "SGN_MB0_Psi", "", nbin, M_low, M_high);

	TH2F* h_SGN_MB0vsBDTx_X = new TH2F( "SGN_MB0vsBDTx_X", "", nbin, Xlow, Xhigh, nbin, M_low, M_high); 
	TH2F* h_SGN_MB0vsBDTx_Psi  = new TH2F( "SGN_MB0vsBDTx_Psi", "", nbin, Xlow, Xhigh, nbin, M_low, M_high); 

	nbin = 40, M_low = .2, M_high = 1.;
   TH1F* h_SGN_MRho_X      = new TH1F( "SGN_MRho_X", "", nbin, M_low, M_high );
   TH1F* h_SGN_MRho_Psi      = new TH1F( "SGN_MRho_Psi", "", nbin, M_low, M_high);
 
   TString fname = inFileSGN_X3872_path;
	// Open SIGNAL file
   TFile *input_signal_X3872(0);
   if (!gSystem->AccessPathName( fname )) {
      input_signal_X3872 = TFile::Open( fname ); // check if file in local directory exists
		std::cout << "--- TMVAClassificationApp    : Using input file SIGNAL X(3872): " << input_signal_X3872->GetName() << std::endl;
   }
   if (!input_signal_X3872) {
      std::cout << "ERROR: could not open signal-data file" << std::endl;
      exit(1);
   }
	// Open BACKGROUND file
   TFile *input_signal_Psi2S(0);
	fname = inFileSGN_PSI2S_path;
   if (!gSystem->AccessPathName( fname )) {
      input_signal_Psi2S = TFile::Open( fname ); // check if file in local directory exists
		std::cout << "--- TMVAClassificationApp    : Using input file SIGNAL Psi(2S): " << input_signal_Psi2S->GetName() << std::endl;
   }
   if (!input_signal_Psi2S) {
      std::cout << "ERROR: could not open background-data file" << std::endl;
      exit(1);
   }
 
   // Event loop
 
	float M_B0, M_MuMu, M_X3872, M_K0s;
   // ===== SIGNAL TREE =====
   std::cout << "--- Select signal sample" << std::endl;
   TTree* TreeSGN_X = (TTree*)input_signal_X3872->Get(SGNtreeName);

   TreeSGN_X->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN_X->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN_X->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN_X->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN_X->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN_X->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN_X->SetBranchAddress( "pT_Rho", &pT_Rho);
	TreeSGN_X->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN_X->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN_X->SetBranchAddress( "M_B0", &M_B0);
   TreeSGN_X->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeSGN_X->SetBranchAddress( "M_X3872", &M_X3872);
   TreeSGN_X->SetBranchAddress( "M_K0s", &M_K0s);
	

	// ===== BACKGROUND TREE =====

   std::cout << "--- Select background sample" << std::endl;
   TTree* TreeSGN_Psi = (TTree*)input_signal_Psi2S->Get(SGNtreeName);

   TreeSGN_Psi->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN_Psi->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN_Psi->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN_Psi->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN_Psi->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN_Psi->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN_Psi->SetBranchAddress( "pT_Rho", &pT_Rho);
   TreeSGN_Psi->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN_Psi->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN_Psi->SetBranchAddress( "M_B0", &M_B0);
   TreeSGN_Psi->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeSGN_Psi->SetBranchAddress( "M_X3872", &M_X3872);
   TreeSGN_Psi->SetBranchAddress( "M_K0s", &M_K0s);

	
	const double Xcut = -0.08;
	Double_t X, BKGrej = 0., X_SGNeff = 0., Psi_SGNeff = 0.;
	bool inSignalRegB0 = false;
	const double MB_sgnR	= 5.35438;
	const double MB_sgnL	= 5.20443;
	// ===== OUTPUT TREE ===== // 
	Double_t SGN_x, SGN_MB0, SGN_MRho;
	Double_t BKG_x, BKG_MB0, BKG_MRho;
	TTree *OutTreeSGN_X = new TTree("BDToutSGN_X3872", "Tree with SIGNAL X(3872) variables post BDT cut");	 
	OutTreeSGN_X->Branch("SGN_x", &X);
	OutTreeSGN_X->Branch("SGN_MB0", &M_B0);
	OutTreeSGN_X->Branch("SGN_MRho", &M_Rho);
	TTree *OutTreeSGN_Psi2S= new TTree("BDToutSGN_Psi2S", "Tree with SIGNAL Psi(2S) variables post BDT cut");	 
	OutTreeSGN_Psi2S->Branch("BKG_x", &X);
	OutTreeSGN_Psi2S->Branch("BKG_MB0", &M_B0);
	OutTreeSGN_Psi2S->Branch("BKG_MRho", &M_Rho);

	// ---> X(3872) SIGNAL evaluation 
   std::cout << "--- Processing: " << TreeSGN_X->GetEntries() << " X(3872) signal events" << std::endl;
   TStopwatch sw;
   sw.Start();

	for (Long64_t ievt=0; ievt<TreeSGN_X->GetEntries();ievt++) {
 
		if (ievt%100 == 0) std::cout << "--- ... Processing X(3872) signal-event: " << ievt << std::endl;
		TreeSGN_X->GetEntry(ievt);

		inSignalRegB0 = (M_B0 > MB_sgnL) &&  (M_B0 < MB_sgnR);
		if (!inSignalRegB0) continue;
		X = reader->EvaluateMVA(methodName);
		h_SGN_MB0vsBDTx_X->Fill(X,M_B0);
		h_BDT_SGN_X->Fill( X );
		h_SGN_MRho_X->Fill(M_Rho);
		if (X > Xcut){
			h_SGN_MB0_X->Fill(M_B0);

			OutTreeSGN_X->Fill();
		}
 
   }
	X_SGNeff /= TreeSGN_X->GetEntries(); 
	// ---> BACKGROUND evaluation 
   std::cout << "--- Processing: " << TreeSGN_Psi->GetEntries() << " Psi(2S) signal-events" << std::endl;

	for (Long64_t ievt=0; ievt<TreeSGN_Psi->GetEntries();ievt++) {
 
		if (ievt%1000 == 0) std::cout << "--- ... Processing Psi(2S) signal event: " << ievt << std::endl;
		TreeSGN_Psi->GetEntry(ievt);

		inSignalRegB0 = (M_B0 > MB_sgnL) &&  (M_B0 < MB_sgnR);
		if (!inSignalRegB0) continue;
		X = reader->EvaluateMVA(methodName);
		h_SGN_MB0vsBDTx_Psi->Fill(X,M_B0);
		h_BDT_SGN_Psi->Fill( X );
		h_SGN_MRho_Psi->Fill(M_Rho);
		if (X > Xcut){	
			Psi_SGNeff +=1.;
			h_SGN_MB0_Psi->Fill(M_B0);

			OutTreeSGN_Psi2S->Fill();
		}
   }


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
 
   // Write histograms
 
   TFile *target  = new TFile(outFile_path,"RECREATE" );
   h_BDT_SGN_X->Write();
   h_BDT_SGN_Psi->Write();
	h_SGN_MB0vsBDTx_X->Write();
	h_SGN_MB0vsBDTx_Psi->Write();
	h_SGN_MB0_X->Write();
	h_SGN_MB0_Psi->Write();
	h_SGN_MRho_X->Write();
	h_SGN_MRho_Psi->Write();
	OutTreeSGN_X->Write();
	OutTreeSGN_Psi2S->Write();
 
   std::cout << "--- Created root file: "<< target->GetName() << " containing the MVA output histograms & TTree" << std::endl;
	target->Close();
 
   delete reader;
 
   std::cout << "==> TMVAClassAppBDT is done!" << std::endl << std::endl;
}
 
int main( int argc, char** argv )
{
   TMVAClassAppBDT();
	



   return 0;
}
