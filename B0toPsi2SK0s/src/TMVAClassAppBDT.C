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
 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
 
using namespace TMVA;
 
TString outFile_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/TMVAapponMC.root"; 

TString BDTweight_dir = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/";
TString inFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/TMVAinputs.root";
TString SGNtreeName = "inputSIGNAL";
//TString inFileBKG_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BKG_data17.root";
TString BKGtreeName = "inputBACKGROUND";


void TMVAClassAppBDT( TString myMethodList = "" )
{
 
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();
 
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
 
   // Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
 
   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
 
   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
 
      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);
 
         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }
 
   // --------------------------------------------------------------------------------------------------
 
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
 
   TString dir    = BDTweight_dir + "datasetBDT/weights/";
   TString prefix = "TMVAClassification";
 
   // Book method
   TString BDTname = "BDT_nT200_S25_D4_b03_nC35";
	TString methodName = BDTname + TString(" method"); 
	TString weightfile = dir + prefix + TString("_") + BDTname + TString(".weights.xml");
	reader->BookMVA( methodName, weightfile );
  // for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
  //    if (it->second) {
  //       TString methodName = TString(it->first) + TString(" method");
  //       TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
  //       reader->BookMVA( methodName, weightfile );
  //    }
  // }
 
   // Book output histograms
   UInt_t nbin = 100;
   TH1F *h_BDT_SGN(0);
   TH1F *h_BDT_BKG(0);
   TH1F *h_EV_SGN(0);
   TH1F *h_EV_BKG(0);
 
   if (Use["BDT"])           h_BDT_SGN     = new TH1F( "MVA_BDT_signal",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDT"])           h_BDT_BKG     = new TH1F( "MVA_BDT_background",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDT"])           h_EV_SGN     = new TH1F( "MVA_EV_signal",           "",           3913, 0., 3913 );
   if (Use["BDT"])           h_EV_BKG     = new TH1F( "MVA_EV_background",       "",           3913, 0., 3913 );
 
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input_signal(0);
   TString fname = inFileSGN_path;
   if (!gSystem->AccessPathName( fname )) {
      input_signal = TFile::Open( fname ); // check if file in local directory exists
		std::cout << "--- TMVAClassificationApp    : Using input file (SIGNAL): " << input_signal->GetName() << std::endl;
   }
   if (!input_signal) {
      std::cout << "ERROR: could not open signal-data file" << std::endl;
      exit(1);
   }
   TFile *input_background(0);
	fname = inFileSGN_path;
   if (!gSystem->AccessPathName( fname )) {
      input_background = TFile::Open( fname ); // check if file in local directory exists
		std::cout << "--- TMVAClassificationApp    : Using input file (BACKGROUND): " << input_background->GetName() << std::endl;
   }
   if (!input_background) {
      std::cout << "ERROR: could not open background-data file" << std::endl;
      exit(1);
   }
 
   // Event loop
 
   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   // ===== SIGNAL TREE =====
   //
	Long64_t iEv;
   std::cout << "--- Select signal sample" << std::endl;
   TTree* TreeSGN = (TTree*)input_signal->Get(SGNtreeName);

   //Float_t pTM_B0, SVprob, LxySign_B0, CosAlpha_B0;  
   //Float_t DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho; 
   TreeSGN->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN->SetBranchAddress( "pT_Rho", &pT_Rho);
   TreeSGN->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN->SetBranchAddress( "iEv", &iEv);

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
   TreeBKG->SetBranchAddress( "iEv", &iEv);

	float X;
	// ---> SIGNAL evaluation 
   std::cout << "--- Processing: " << TreeSGN->GetEntries() << "signal events" << std::endl;
   TStopwatch sw;
   sw.Start();

	for (Long64_t ievt=0; ievt<TreeSGN->GetEntries();ievt++) {
 
		if (ievt%100 == 0) std::cout << "--- ... Processing signal-event: " << ievt << std::endl;
		TreeSGN->GetEntry(ievt);
		if (Use["BDT"]) h_BDT_SGN->Fill( reader->EvaluateMVA(methodName ) );
		X = reader->EvaluateMVA(methodName ) ;	
		if((X > -0.01) && (M_Rho > 0.695))h_EV_SGN->Fill(iEv);
 
   }
 
	// ---> BACKGROUND evaluation 
   std::cout << "--- Processing: " << TreeBKG->GetEntries() << "background events" << std::endl;

	for (Long64_t ievt=0; ievt<TreeBKG->GetEntries();ievt++) {
 
		if (ievt%1000 == 0) std::cout << "--- ... Processing background-event: " << ievt << std::endl;
		TreeBKG->GetEntry(ievt);
		h_BDT_BKG->Fill( reader->EvaluateMVA(methodName ) );
		X = reader->EvaluateMVA(methodName ) ;	
		if((X > -0.01) && (M_Rho > 0.695)) {
			h_EV_BKG->Fill(iEv);
			std::cout << "  #EVbkg " << iEv + 1<< std::endl; 
		}
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
 
   // Write histograms
 
   TFile *target  = new TFile(outFile_path,"RECREATE" );
   h_BDT_SGN->Write();
   h_BDT_BKG->Write();
   h_EV_SGN->Write();
   h_EV_BKG->Write();
 
	target->Close();
   std::cout << "--- Created root file: \"TMVAprovaApp.root\" containing the MVA output histograms" << std::endl;
 
   delete reader;
 
   std::cout << "==> TMVAClassAppBDT is done!" << std::endl << std::endl;
}
 
int main( int argc, char** argv )
{
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVAClassAppBDT(methodList);
   return 0;
}
