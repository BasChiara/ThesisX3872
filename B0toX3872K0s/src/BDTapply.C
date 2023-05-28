#include "../include/BDTapply.h"


BDTapply::BDTapply(const TString& infile,const TString& intree, const TString& outfile, const TString& outtree){

//...inputs 
	
	inFilePath_= infile;
	inTreeName_= intree;

//...outputs
	outFilePath_ = outfile;
	outTreeName_= outtree;
	
	// BDT model
	BDTweightPath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/datasetBDT/weights/TMVAClassification_BDT_nT200_S25_D4_b03_nC35.weights.xml";


}//BDTapply()



int BDTapply::Apply(){

// --> TMVA SET UP
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
 
   // Book the BDT method
   TString methodName = TString("BDT method"); 
	reader->BookMVA( methodName, BDTweightPath_);

// --> SIGNAL & BACKGROUND files reading
   TFile *InFile(0);
   if (!gSystem->AccessPathName(inFilePath_)) {
      InFile = TFile::Open(inFilePath_); // check if file in local directory exists
   }
   if (!InFile) {
      std::cout << "ERROR: could not open the data file : [" << inFilePath_ << "]" << std::endl;
      exit(1);
	}
	std::cout << "     TMVAClassificationApp    : Using input file : " << InFile->GetName() << std::endl;
 
// --> EVENT LOOP 
 

	float M_B0, M_MuMu, M_X3872, M_K0s;
   // ===== SIGNAL TREE =====
   std::cout << "--- Select signal sample" << std::endl;
   TTree* InTree = (TTree*)InFile->Get(inTreeName_);

   InTree->SetBranchAddress( "pTM_B0", &pTM_B0);
   InTree->SetBranchAddress( "SVprob", &SVprob);
   InTree->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   InTree->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   InTree->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   InTree->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   InTree->SetBranchAddress( "pT_Rho", &pT_Rho);
	InTree->SetBranchAddress( "D0_Rho", &D0_Rho);
   InTree->SetBranchAddress( "M_Rho", &M_Rho);
   InTree->SetBranchAddress( "M_B0", &M_B0);
   InTree->SetBranchAddress( "M_mumu", &M_MuMu);
   InTree->SetBranchAddress( "M_X3872", &M_X3872);
   InTree->SetBranchAddress( "M_K0s", &M_K0s);
	 // friend-tree with BDT response
	Float_t Xs;
	TFile* output = new TFile(outFilePath_, "RECREATE"); 
	TTree* TreeBDTx = new TTree(outTreeName_, "friend with BDT response"); 
	TreeBDTx->Branch("BDTx", &Xs, "BDTx/F");

	// ---> BDT evaluation 
   std::cout << "--- Processing: " << InTree->GetEntries() << " signal events" << std::endl;
   TStopwatch sw;
   sw.Start();

	for (Long64_t ievt=0; ievt<InTree->GetEntries();ievt++) {
 
		if (ievt%100 == 0) std::cout << "    ... Processing event: " << ievt << std::endl;
		InTree->GetEntry(ievt);
		Xs = reader->EvaluateMVA(methodName);

		TreeBDTx->Fill(); 
	}

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
	
	output->cd();
	TreeBDTx->Write();

	// close files
	output->Close();
	InFile->Close();

	return 0;

} //ApplyBDT()



