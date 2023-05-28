#include "../include/CutterMVA.h"

// RooFit
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "TLatex.h"
using namespace RooFit;

CutterMVA::CutterMVA(){

	// SGN scale factor
	SGNfactor = 0.107/3.;

	// I/O files

	inFileSGN_path  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGN_v2_MC.root";
	SGNtreeName     = "mcSIGNAL";
	outFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonSGN.root";

	inFileBKG_path  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data/merged/UNBL_MergData17_v2.root";
	BKGtreeName     = "B0_X3872signal";
	outFileBKG_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonDATAunbl_v2.root";

	inFileFitParB0_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/B0params.txt";
	inFileFitParK0s_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/K0sparams.txt";
	inFileFitParX_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/X3872params.txt";
	outPlot_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOnUnblind/";

	LoadFitRanges();

	SetUpInTreeSGN();
	SetUpInTreeBKG();

	// BDT model
	BDTweightPath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/datasetBDT/weights/TMVAClassification_BDT_nT200_S25_D4_b03_nC35.weights.xml";


}//CutterMVA()

CutterMVA::~CutterMVA(){ 

	InFileSGN->Close();
	InFileBKG->Close();
}

//    ==== UTILITY ====     //
TFile* CutterMVA::GetFile(TString file){

   TFile *File(0);
	TString Fpath("");
	if (file == "InS") Fpath = inFileSGN_path; 
	if (file == "InB") Fpath = inFileBKG_path;
	if (file == "xS")  Fpath = outFileSGN_path;
	if (file == "xB")  Fpath = outFileBKG_path;

	// Open SIGNAL file
   if (!gSystem->AccessPathName(Fpath)) {
      File = TFile::Open(Fpath); // check if file in local directory exists
   }
   if (!File) {
      std::cout << "ERROR: could not open signal-data file" << std::endl;
      exit(1);
   }

	return File; 

}//GetFile()


void CutterMVA::SetUpInTreeSGN(){

   InFileSGN = GetFile("InS");
   // ===== SIGNAL TREE =====
	TreeSGN = (TTree*)InFileSGN->Get(SGNtreeName);

   TreeSGN->SetBranchAddress( "run", &run);
   TreeSGN->SetBranchAddress( "LumiBlock", &LumiBlock);
   TreeSGN->SetBranchAddress( "event", &event);

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
   TreeSGN->SetBranchAddress( "M_K0sPi1", &M_K0sPi1);
   TreeSGN->SetBranchAddress( "M_K0sPi2", &M_K0sPi2);

	 // friend-tree with BDT response
   TFile* outputBDT = GetFile("xS");
   TTree* Tree_xS= (TTree*)outputBDT->Get("TreeBDTx_S");
	Tree_xS->SetBranchAddress( "BDTx", &X);
	TreeSGN->AddFriend("TreeBDTx_S");


}//SetUpInTreeSGN()

void CutterMVA::SetUpInTreeBKG(){

	// ===== BACKGROUND TREE =====

   InFileBKG = GetFile("InB");
   TreeBKG = (TTree*)InFileBKG->Get(BKGtreeName);

   TreeBKG->SetBranchAddress( "run", &run);
   TreeBKG->SetBranchAddress( "LumiBlock", &LumiBlock);
   TreeBKG->SetBranchAddress( "event", &event);

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
   TreeBKG->SetBranchAddress( "M_K0sPi1", &M_K0sPi1);
   TreeBKG->SetBranchAddress( "M_K0sPi2", &M_K0sPi2);
	
   TFile* outputBDT = GetFile("xB");
   TTree* Tree_xB= (TTree*)outputBDT->Get("TreeBDTx_B");
	Tree_xB->SetBranchAddress( "BDTx", &X);
	TreeBKG->AddFriend("TreeBDTx_B");



}// SetUpInTreeBKG()



void CutterMVA::LoadFitRanges(){
	
	std::string line;
	int Nline = 0;

	char ParName[30];
	double err;
	// --- B0 FIT 
	std::ifstream inFileParB0(inFileFitParB0_path);	
	if(!inFileParB0.is_open()) std::cout << "ERROR cannot open " << inFileFitParB0_path << std::endl;
	while(!inFileParB0.eof()){

		getline(inFileParB0, line); Nline++;
		//std::cout << Nline << "\t" << line << std::endl;
		if(line.c_str()[0] == '#') continue;
		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &BKGregNearL);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &BKGregFarL);
		if(Nline == 5) sscanf(line.c_str(), "%s %lf", ParName, &BKGregNearR);
		if(Nline == 6) sscanf(line.c_str(), "%s %lf", ParName, &BKGregFarR);
		if(Nline == 9) sscanf(line.c_str(), "%s %lf", ParName, &SGNregL);
		if(Nline ==10) sscanf(line.c_str(), "%s %lf", ParName, &SGNregR);
	}
	inFileParB0.close();

	// --- K0s FIT 
	Nline = 0;
	std::ifstream inFileParK0s(inFileFitParK0s_path);	
	if(!inFileParK0s.is_open()) std::cout << "ERROR cannot open " << inFileFitParK0s_path << std::endl;
	while(!inFileParK0s.eof()){

		getline(inFileParK0s, line); Nline++;
		if(line.c_str()[0] == '#') continue;

		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &MK0s_nearLeft);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &MK0s_nearRight);
	}
	inFileParK0s.close();

	
	// --- X FIT 
	Nline = 0;
	std::ifstream inFileParX(inFileFitParX_path);	
	if(!inFileParX.is_open()) std::cout << "ERROR cannot open " << inFileFitParX_path << std::endl;
	while(!inFileParX.eof()){

		getline(inFileParX, line); Nline++;
		if(line.c_str()[0] == '#') continue;

		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &MX_nearLeft);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &MX_nearRight);
	}
	inFileParX.close();

	std::cout << " ---> MASS FIT RESULTS " << std::endl;
	std::cout << "    B0 sidebands   [" << BKGregFarL    << "," << BKGregNearL    << "]" << " + [" << BKGregNearR << "," << BKGregFarR << "]"  <<std::endl; 
	std::cout << "    B0 SGNregion   [" << SGNregL       << "," << SGNregR        << "]" << std::endl;
	std::cout << " X(3872) SGNregion [" << MX_nearLeft   << "," << MX_nearRight   << "]" << std::endl; 
	std::cout << "   K0s SGNregion   [" << MK0s_nearLeft << "," << MK0s_nearRight << "]" << std::endl; 

	

}//LoadFitRanges()

int CutterMVA::GetSignalPreCut(){
	
	int S_precut = 0;
	int Nentries = TreeSGN->GetEntries();
	bool inSGNregion;

	for(int ien = 0; ien < Nentries; ien++){
		TreeSGN->GetEntry(ien);
		inSGNregion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if(inSGNregion)S_precut++;
	}
	
	return S_precut;
	
}//GetSignalPreCut()

void CutterMVA::GetCorrelation_BDT_MR(){

	double CorrF;
	const int NbinsX = 100; 
	const double Xlow = -0.8, Xhigh = 0.5;
	const int NbinsMR = 100; 
	const double MRlow = 0.4, MRhigh = .9;
	TH2F* h_CorrS = new TH2F("CorrS_XvsMrho", "", NbinsMR, MRlow, MRhigh, NbinsX, Xlow, Xhigh);	
	TH2F* h_CorrB = new TH2F("CorrB_XvsMrho", "", NbinsMR, MRlow, MRhigh, NbinsX, Xlow, Xhigh);	
	
	TreeSGN->Draw("TreeBDTx_S.BDTx:M_Rho>>CorrS_XvsMrho");
	h_CorrS->GetXaxis()->SetTitle("M(\\rho(770))\\ [GeV]");
	h_CorrS->GetYaxis()->SetTitle("BDT response [x]");
	h_CorrS->GetYaxis()->SetTitleOffset(1.05);
	h_CorrS->SetLineColor(kAzure +1);
	h_CorrS->SetFillColorAlpha(kAzure +1, 0.60);
	TreeBKG->Draw("TreeBDTx_B.BDTx:M_Rho>>CorrB_XvsMrho");
	h_CorrB->SetLineColor(kRed);
	h_CorrB->SetFillColorAlpha(kRed, 0.50);


	auto legend= new TLegend(0.525, 0.8,.89,.89);
	legend->SetTextSize(0.025);
	legend->SetBorderSize(0);
	
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);	
	gStyle->SetOptStat(0);
	h_CorrS->Draw("BOX");
	legend->AddEntry(h_CorrS, Form("SIGNAL (MC) Corr = %.3f", h_CorrS->GetCorrelationFactor()));
	h_CorrB->Draw("BOX SAME");
	legend->AddEntry(h_CorrB, Form("DATA 2017   Corr = %.3f", h_CorrB->GetCorrelationFactor()));
	legend->Draw();

	TString outPath = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/";
	c1->SaveAs(outPath + "CorrelationBDTxMRho.png");

}

double CutterMVA::GetAvMultiplicity(const double Xcut, const double MRcut, TString dataset = "SGN"){
	
	TTree* inTree;
	if (dataset == "SGN") inTree = TreeSGN;
	if (dataset == "BKG") inTree = TreeBKG;

	double AvMultiplicity;
	int Nentries = inTree->GetEntries(), Ncand = 0, Nevents = 0;
	float prevRun = -1, prevLumiBlock = -1, prevEvent = -1;
	bool isNewEvent, passCut;

	for(int ientry = 0; ientry < Nentries; ientry++){
		inTree->GetEntry(ientry);
		passCut = (M_Rho > MRcut) && (X > Xcut );
		if (!passCut) continue;
		Ncand++;
		isNewEvent = !((run == prevRun) && (LumiBlock == prevLumiBlock ) && (event == prevEvent));
		prevRun = run; prevLumiBlock = LumiBlock; prevEvent = event;
		if ( isNewEvent ) Nevents++;

	}
	std::cout << "#cand vs #events\t" << Ncand << " " << Nevents << std::endl;

	return ((double)Ncand/Nevents);

}//GetAvMultiplicitySGN()



// ==== ANALYSIS TOOLS ==== //

int CutterMVA::ApplyBDT(){

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
	TFile* input_signal = GetFile("InS");
	std::cout << "     TMVAClassificationApp    : Using input file (SIGNAL): " << input_signal->GetName() << std::endl;
 
   TFile *input_background = GetFile("InB");
	std::cout << "     TMVAClassificationApp    : Using input file (BACKGROUND): " << input_background->GetName() << std::endl;

// --> EVENT LOOP 
 

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
	 // friend-tree with BDT response
	Float_t Xs;
	TFile* output_signal = new TFile(outFileSGN_path, "RECREATE"); 
	TTree* TreeBDTx_S = new TTree("TreeBDTx_S", "friend with BDT response"); 
	TreeBDTx_S->Branch("BDTx", &Xs, "BDTx/F");


	// ===== BACKGROUND TREE =====

   std::cout << "--- Select background sample" << std::endl;
	TreeBKG = (TTree*)input_background->Get(BKGtreeName);

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
	 // friend-tree with BDT response
	Float_t Xb;
	TFile* output_background = new TFile(outFileBKG_path, "RECREATE"); 
	TTree* TreeBDTx_B = new TTree("TreeBDTx_B", "friend with BDT response"); 
	TreeBDTx_B->Branch("BDTx", &Xb, "BDTx/F");
	


	// ---> SIGNAL evaluation 
   std::cout << "--- Processing: " << TreeSGN->GetEntries() << " signal events" << std::endl;
   TStopwatch sw;
   sw.Start();

	for (Long64_t ievt=0; ievt<TreeSGN->GetEntries();ievt++) {
 
		if (ievt%100 == 0) std::cout << "    ... Processing signal-event: " << ievt << std::endl;
		TreeSGN->GetEntry(ievt);
		Xs = reader->EvaluateMVA(methodName);

		TreeBDTx_S->Fill(); 
   }

	// ---> BACKGROUND evaluation 
   std::cout << "--- Processing: " << TreeBKG->GetEntries() << " background events" << std::endl;

	for (Long64_t ievt=0; ievt<TreeBKG->GetEntries();ievt++) {
 
		if (ievt%1000 == 0) std::cout << "    ... Processing background-event: " << ievt << std::endl;
		TreeBKG->GetEntry(ievt);
		Xb = reader->EvaluateMVA(methodName);

		TreeBDTx_B->Fill();
   }


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
	
	output_signal->cd();
	TreeBDTx_S->Write();

	output_background->cd();
	TreeBDTx_B->Write();

	// close files
	output_signal->Close();
	output_background->Close();
	input_signal->Close();
	input_background->Close();	



	return 0;

} //ApplyBDT()



double CutterMVA::BKG_NevExtraction(const double Xcut, const double MRcut){

	double Bs;

	// ==== ROOFIT SET UP ==== //
	
	RooRealVar B0m("B0m", "M(B_0)\\ [GeV]", Mlow_B0, Mhigh_B0);

	B0m.setRange(   "SB1"      , Mlow_B0, BKGregNearL);       // SIDEBAND 1
	B0m.setRange(   "SB2"      , BKGregNearR, Mhigh_B0);        // SIDEBAND 2
	B0m.setRange("BLINDregion" , BKGregNearL, BKGregNearR);        // SIDEBAND 2
	B0m.setRange("FULLregion"  , Mlow_B0    , Mhigh_B0);        // FULLregion RANGE 
	B0m.setRange("SGNregion"   , SGNregL    , SGNregR);           // SGNregion RANGE 

	RooDataSet B0m_ds("B0m_ds","B0m_ds", RooArgSet(B0m));

	// Expo + C	
	RooRealVar tau("tau", "", -15.0, -30.0 , -5.);
	RooExponential expo("expo", "", B0m, tau);
	RooProdPdf BKGdesc("BKGdesc", "BKGdesc", expo); 

	RooPolynomial Const("Const", "", B0m);
	RooProdPdf BKGbase("BKGbase", "BKGbase", Const); 

	
	// Fermi 
	RooRealVar Slope("Slope", "", 10.0, 1.0 , 100.);
	RooRealVar Flex("Flex", "", 4., 1. , 5.);
	RooRealVar C("C", "", 0.5, 0. , 1.);
	
	RooGenericPdf Fermi("Fermi", "", "1./(1. + exp(( @0 - @1)*@2)) + @3", RooArgList(B0m, Flex, Slope, C));

	// Build PDF
	RooRealVar bkg_yield("nBKG", "", 1600, 100, 10000);
	RooRealVar f("f", "", 0., 1.);
	RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Fermi, bkg_yield);
	//RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", RooArgList(expo, Const), f);


	// ==== LOOP ON EVENTS ==== // 
	int N = TreeBKG->GetEntries();
	bool isOutRangeB0, isOutSgnRangeX3872, isOutSgnRangeK0s;
	for (int i = 0; i < N; i++){

		TreeBKG->GetEntry(i);
		isOutRangeB0 = (M_B0 < Mlow_B0) || (M_B0 > Mhigh_B0) || ((M_B0 > BKGregNearL ) && (M_B0 < BKGregNearR));
		isOutSgnRangeX3872 = (M_X3872 < MX_nearLeft) || (M_X3872 > MX_nearRight);
		isOutSgnRangeK0s = (M_K0s < MK0s_nearLeft) || (M_K0s > MK0s_nearRight);
		if( ( X < Xcut) || ( M_Rho < MRcut) || isOutRangeB0 || isOutSgnRangeX3872 || isOutSgnRangeK0s) continue;
		B0m = M_B0;
		B0m_ds.add(RooArgSet(B0m));

	}
	B0m_ds.Print("v");
	// ==== FIT PERFORM ==== //
	
	RooFitResult *ResBKGmodelSB2 = BKGmodel.fitTo(B0m_ds, Range("SB2"), Save());
	//ResBKGmodelSB2->Print();
	RooFitResult *ResBKGmodelSB1 = BKGmodel.fitTo(B0m_ds, Range("SB1"), Save());
	//ResBKGmodelSB1->Print();
	RooFitResult *ResBKGmodel= BKGmodel.fitTo(B0m_ds, Range("SB1,SB2"), Save());
	//ResBKGmodel->Print("v");

	// ==== SAVE RESULTS ==== //
	RooPlot* B0m_fr = B0m.frame(Bins(35));
	B0m_ds.plotOn(B0m_fr);
	BKGmodel.plotOn(B0m_fr, Range("FULLregion"), RooFit::NormRange("SB1,SB2"));

	
	TText *TxtChi2= new TText(5.2, B0m_fr->GetMaximum()*0.85, Form("Chi^2 = %.3f", B0m_fr->chiSquare()));
   TxtChi2->SetTextSize(0.035);
   TxtChi2->SetTextColor(kRed);
   B0m_fr->addObject(TxtChi2);
	std::cout << " ---> Chi^2 = " << B0m_fr->chiSquare() << std::endl;
	BKGmodel.paramOn(B0m_fr, Layout(0.60));
	B0m_fr->SetTitle(Form("X-cut = %.2f Mrho-cut = %.3f", Xcut, MRcut));

	TH2 *h_Corr = ResBKGmodel->correlationHist();
	
	// ==== Bs ==== //
	RooAbsReal* IntSreg = BKGmodel.createIntegral(B0m, NormSet(B0m), Range("SGNregion")); //integrate sgn region
	double Is = IntSreg->getVal();
	RooAbsReal* IntBreg = BKGmodel.createIntegral(B0m, NormSet(B0m), Range("SB1,SB2"));   //integrate bkg region
	double Ib = IntBreg->getVal();

	Bs = B0m_ds.sumEntries() * Is/Ib; 

	// ==== SAVE ON FILE ==== //
	TString outFileName(outPlot_path + "FitRes/FitBKG_X");
	if (Xcut < 0.) outFileName.Append(Form("m"));
	outFileName.Append(Form("%.0f_MR%.0f.root", fabs(Xcut)*100., MRcut*1000.));

	TFile* outFitFile = new TFile(outFileName, "RECREATE");
	B0m_fr->Write("FitPlot");
	ResBKGmodel->Write("FitResults");
	h_Corr->Write();

	outFitFile->Close();
	

	return Bs;

}//BKG_NevExtraction()




double CutterMVA::SGN_NevExtraction(const double Xcut, const double MRcut){

	int iSs = 0;
	double Ss;
	bool inSGNregion = false;
	int Nentries = TreeSGN->GetEntries();
	bool isOutSgnRangeK0s,  isOutSgnRangeX3872;
	for(int ien = 0; ien < Nentries; ien++){
		TreeSGN->GetEntry(ien);
		isOutSgnRangeK0s = (M_K0s < MK0s_nearLeft) || (M_K0s > MK0s_nearRight);
		isOutSgnRangeX3872 = (M_X3872 < MX_nearLeft) || (M_X3872 > MX_nearRight);
		inSGNregion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if( (X < Xcut) || (M_Rho < MRcut)  || (!inSGNregion) || isOutSgnRangeX3872 || isOutSgnRangeK0s) continue;
		iSs ++;
	}
	Ss = iSs * SGNfactor;
	return Ss;

}//SGN_NevExtraction()
	



double CutterMVA::PunziSign(const double S, const double B, double* PSerr){

	double PunziSign, Sign_denom, Sign_denom_error;
	double SGNeff, SGNeff_error;
	const double b = 5.0;    // #sigmas corresp. to 5sigma significance level
	const double a = 2.0; // #sigmas corresp. to CL (90%--> 1.2816) (95% --> 1.6448) (CMStwiki --> 2.)

	SGNeff       = S / (GetSignalPreCut() * SGNfactor);
	SGNeff_error = SGNfactor / S; // error square
	
	Sign_denom        = b*b + 2.*a*sqrt(B) + b*sqrt(b*b + 4.*a*sqrt(B) + 4.*B );
	Sign_denom_error  = a/sqrt(B) + b * (a / sqrt(B) + 2. ) / sqrt(b*b + 4.*a*sqrt(B) + 4.*B );	
	Sign_denom_error /= Sign_denom;	

	PunziSign = SGNeff/Sign_denom;
	*PSerr = PunziSign*sqrt(SGNeff_error + Sign_denom_error*Sign_denom_error);


	return PunziSign;

}



void  CutterMVA::SGNfitter(const double Xcut, const double MRcut){

	int Nbins = 70;
	double Mlow = 5.0, Mhigh = 5.6;
	TH1F* h_SGN_B0  = new TH1F("SGN_B0" , "", Nbins, Mlow, Mhigh); 
	RooRealVar B0m("B0m", "M(B_0)\\ [GeV]", BKGregFarL, BKGregFarR);
	B0m.setRange(   "SB1"      , BKGregFarL , BKGregNearL);       // SIDEBAND 1
	B0m.setRange(   "SB2"      , BKGregNearR, BKGregFarR);        // SIDEBAND 2
	B0m.setRange("BLINDregion" , BKGregNearL, BKGregNearR);        // SIDEBAND 2
	B0m.setRange("SGNregion"   , SGNregL    , SGNregR);           // SGNregion RANGE 
	B0m.setRange("FULLregion"  , BKGregFarL , BKGregFarR);        // FULLregion RANGE 

	RooDataSet SB_B0m_ds("SB_B0m_ds","SB_B0m_ds", RooArgSet(B0m));		
	// LOOP ON SIGNAL EVENTS
	int Nsignal = TreeSGN->GetEntries();
	bool inSGNregion;
	for (int i = 0; i < Nsignal; i++){
		TreeSGN->GetEntry(i);
		inSGNregion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if( ( X < Xcut) || ( M_Rho < MRcut) || !inSGNregion) continue;
		h_SGN_B0->Fill(M_B0);
	}
	h_SGN_B0->Scale(SGNfactor);
	RooDataHist h_Signal("h_Signal", "h_Signal", B0m, Import(*h_SGN_B0));

	// LOOP ON BACKGROUND EVENTS
	int Nbkg = TreeBKG->GetEntries();
	for (int i = 0; i < Nbkg ; i++){
		TreeBKG->GetEntry(i);
		if( ( X < Xcut) || ( M_Rho < MRcut) ) continue;
		B0m = M_B0;
		SB_B0m_ds.add(RooArgSet(B0m));
	}
	SB_B0m_ds.Print("v");
	// --> BKG model
	RooRealVar Slope("Slope", "", 9.3, 1.0 , 100.);
	RooRealVar Flex("Flex", "", 4.39, 1. , 5.2);
	RooRealVar C("C", "", 0.0005, 0. , 1.);
	
	RooGenericPdf Fermi("Fermi", "", "1./(1. + exp(( @0 - @1)*@2)) + @3", RooArgList(B0m, Flex, Slope, C));
	RooDataSet* FR_B0m_ds = Fermi.generate(B0m, 190); 
	RooDataSet* BR_B0m_ds = (RooDataSet*)FR_B0m_ds->reduce(Form("(B0m < %f) && (B0m > %f)", BKGregNearR, BKGregNearL)); 
	//SB_B0m_ds.append(*BR_B0m_ds);	
	
	// --> SGN models
	// Breit-Wigner 
	RooRealVar B0mass ("B0mass" , "", 5.279, 5.275, 5.285);
	RooRealVar B0width("B0width", "", 0.01, 0.  , 0.02);
	RooBreitWigner BWsignal("BWsignal", "", B0m, B0mass, B0width);
	RooProdPdf SGNpdfBW("SGNpdfBW", "SGNpdfBW", BWsignal ); 

	// CrystalBall + Gauss
	RooRealVar alpha("alpha", "", 2., 0.,  10.);
	RooRealVar     N(    "N", "",25., 10. ,  50 );
	RooCBShape  CBsignal("CBsignal", "", B0m, B0mass, B0width, alpha, N);
	RooProdPdf SGNpdfCB("SGNpdfCB", "SGNpdfCB", CBsignal ); 
	RooRealVar B0shift("B0shift", "", 0);
	RooRealVar B0sigma("B0sigma", "", 0.02, 0.01  , 2.);
	RooGaussian Gsignal("Gsignal", "", B0m, B0shift, B0sigma);
	RooProdPdf SGNpdfG("SGNpdfG", "SGNpdfG", Gsignal ); 
	
	RooFFTConvPdf BWxG("BWxG", "BreitWig (X) Gauss", B0m, BWsignal, Gsignal);

	// PDF
	RooRealVar sgn_yield("nSGN", "", 79, 1, 100);
	RooRealVar f("f", "", 0., 1.);
	RooRealVar BKGfrac("BKGfrac", "", 0.5, 0., 1.);
	RooAddPdf SGNmodel("SGNmodel", "SGNmodel", RooArgList(CBsignal, Gsignal), f); 
	RooAddPdf FULLmodel("FULLmodel", "FULLmodel", RooArgList(Fermi, BWxG), BKGfrac); 

	// FIT the signal only
	//RooFitResult *ResSGNmodel = SGNmodel.fitTo(h_Signal, Range("SGNregion"), Save());
	//RooFitResult *ResSGNmodel = BWxG.fitTo(h_Signal, Range("SGNregion"), Save());
	//h_Signal.add(SB_B0m_ds);
	//RooFitResult *ResFULLmodel = FULLmodel.fitTo(h_Signal, Range("FULLregion"), Save());
	//ResSGNmodel->Print("v");
	//ResFULLmodel->Print("v");

	// SAVE RESULTS //
	RooPlot* SgnB0m_fr = B0m.frame(Bins(47));
	//h_Signal.plotOn(SgnB0m_fr,LineColor(kBlue), MarkerColor(kBlue)); 
	//h_Signal.plotOn(SgnB0m_fr);
	BR_B0m_ds->plotOn(SgnB0m_fr, LineColor(kGreen), MarkerColor(kGreen));
	SB_B0m_ds.plotOn(SgnB0m_fr);
	//FULLmodel.plotOn(SgnB0m_fr);
	//TText *TxtChi2= new TText(5.31, SgnB0m_fr->GetMaximum()*0.85, Form("ChiSq = %.3f", SgnB0m_fr->chiSquare()));
	//TxtChi2->SetTextSize(0.035);   TxtChi2->SetTextColor(kRed);
   //SgnB0m_fr->addObject(TxtChi2);
	//std::cout << " ----> Chi2 SIGNAL FIT \t" << SgnB0m_fr->chiSquare() << std::endl;	
	//SB_B0m_ds.plotOn(SgnB0m_fr, LineColor(kRed), MarkerColor(kRed));
	//FULLmodel.plotOn(SgnB0m_fr, Components(Fermi), LineColor(kRed), LineStyle(kDashed), NormRange("FULLregion"));
	//FULLmodel.plotOn(SgnB0m_fr, Components(BWxG), LineColor(kGreen), LineStyle(kDashed), NormRange("FULLregion"));

	//RooAbsReal* IntSreg = FULLmodel.createIntegral(B0m, NormSet(B0m), Range("SGNregion"));
	//std::cout << "\n ---> INEGRAL full model sgn region " << IntSreg->getVal() << std::endl;

	// DRAW RESULTS	
	h_SGN_B0->GetXaxis()->SetTitle("M(B_0)\\ [GeV]");
	h_SGN_B0->GetYaxis()->SetTitle(Form("Events/%.0f MeV", h_SGN_B0->GetXaxis()->GetBinWidth(10)* 1000.));
	h_SGN_B0->SetMaximum(30);
	h_SGN_B0->SetLineColor(kAzure +1 ); h_SGN_B0->SetFillColorAlpha(kAzure +1, 0.3);
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	//h_SGN_B0->Draw("HIST");
	SgnB0m_fr->Draw("SAME");
	c1->SaveAs("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/FULLfit_dataset.png");
	gPad->SetLogy();
	c1->SaveAs("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/FULLfit_datasetLOG.png");

}//SGNfitter()



int CutterMVA::makeSGNvsBKGplot(const double Xcut, const double MRcut){
	int Nbins = 35;
	TH1F* h_Data_B0 = new TH1F("Data_B0", "", Nbins, Mlow_B0, Mhigh_B0);
	TH1F* h_SGN_B0  = new TH1F("SGN_B0" , "", Nbins, Mlow_B0, Mhigh_B0); 
	Nbins = 15;
	double Mlow = 3.8, Mhigh = 3.95;
	TH1F* h_Data_X  = new TH1F("Data_X", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_X   = new TH1F("SGN_X" , "", Nbins, Mlow, Mhigh); 
	Nbins = 20;
	Mlow = 0.45, Mhigh = 0.55;
	TH1F* h_Data_K0s  = new TH1F("Data_K0s", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_K0s   = new TH1F("SGN_K0s" , "", Nbins, Mlow, Mhigh); 
	Nbins = 35;
	Mlow = 0.3, Mhigh = 1.;
	TH1F* h_Data_Rho  = new TH1F("Data_Rho", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_Rho   = new TH1F("SGN_Rho" , "", Nbins, Mlow, Mhigh); 
	TH2F* h_Data_XvsRho = new TH2F("Data_XvsRho", "", 20, 0.6, Mhigh, 30, 3.8, 4.1);
	TH2F* h_SGN_XvsRho = new TH2F("SGN_XvsRho", "", 20, 0.6, Mhigh, 30, 3.8, 4.1);
	Nbins = 30; 	
	Mlow = 0.5, Mhigh = 2.;
	TH1F* h_Data_K0sPi  = new TH1F("Data_K0sPi", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_K0sPi  = new TH1F("SGN_K0sPi", "", Nbins, Mlow, Mhigh);

	//MX_nearRight = 5.;
	//MX_nearLeft = 3.88;
	//MK0s_nearRight = 0.6;
	//MK0s_nearLeft = 0.4;
	// ==== LOOP ON EVENTS ==== // 
	// ...BKG
	int Nb = TreeBKG->GetEntries();
	bool isOutRangeB0, isOutSgnRangeX3872, isOutSgnRangeK0s, isB0sgnRegion;
	for (int i = 0; i < Nb; i++){
		TreeBKG->GetEntry(i);
		isOutRangeB0 = (M_B0 < Mlow_B0) || (M_B0 > Mhigh_B0);
		isOutSgnRangeX3872 = (M_X3872 < MX_nearLeft) || (M_X3872 > MX_nearRight);
		isOutSgnRangeK0s = (M_K0s < MK0s_nearLeft) || (M_K0s > MK0s_nearRight);
		isB0sgnRegion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if( ( X < Xcut) || ( M_Rho < MRcut) || isOutRangeB0 || isOutSgnRangeX3872 || isOutSgnRangeK0s) continue;
		h_Data_B0->Fill(M_B0);
		if(!isB0sgnRegion) continue;
		h_Data_X->Fill(M_X3872);
		h_Data_K0s->Fill(M_K0s);
		h_Data_K0sPi->Fill(M_K0sPi1);
		h_Data_K0sPi->Fill(M_K0sPi2);
		h_Data_Rho->Fill(M_Rho);
		h_Data_XvsRho->Fill(M_Rho, M_X3872);
	}
	int Ns = TreeSGN->GetEntries();
	for (int i = 0; i < Ns; i++){
		TreeSGN->GetEntry(i);
		isOutRangeB0 = (M_B0 < Mlow_B0) || (M_B0 > Mhigh_B0);
		isOutSgnRangeX3872 = (M_X3872 < MX_nearLeft) || (M_X3872 > MX_nearRight);
		isOutSgnRangeK0s = (M_K0s < MK0s_nearLeft) || (M_K0s > MK0s_nearRight);
		isB0sgnRegion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if( ( X < Xcut) || ( M_Rho < MRcut) || isOutRangeB0 || isOutSgnRangeX3872 || isOutSgnRangeK0s) continue;
		h_SGN_B0->Fill(M_B0);
		if(!isB0sgnRegion) continue;
		h_SGN_X->Fill(M_X3872);
		h_SGN_K0s->Fill(M_K0s);
		h_SGN_K0sPi->Fill(M_K0sPi1);
		h_SGN_K0sPi->Fill(M_K0sPi2);
		h_SGN_Rho->Fill(M_Rho);
		h_SGN_XvsRho->Fill(M_Rho, M_X3872);
	}
	// Histo set-up
	// ---> B0
	//h_SGN_B0->SetTitle(Form("B0 %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_B0->GetXaxis()->SetTitle("M(B^{0})[GeV]");
	h_SGN_B0->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_B0->GetXaxis()->GetBinWidth(1)));
	h_SGN_B0->GetXaxis()->SetTitleOffset(1.1); h_SGN_B0->GetXaxis()->SetTitleSize(0.04);
	h_SGN_B0->GetYaxis()->SetTitleSize(0.04);
	h_SGN_B0->SetLineWidth(3);
	h_SGN_B0->SetLineColor(kAzure + 1); h_SGN_B0->SetFillColorAlpha(kAzure + 1, 0.30);
	h_Data_B0->SetLineColor(kBlack);
	h_Data_B0->SetMarkerStyle(20);
	h_Data_B0->SetLineWidth(2);
	// ---> X(3872)
	//h_SGN_X->SetTitle(Form("X %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_X->GetXaxis()->SetTitle("M(J/#Psi #pi^{+}#pi^{-})[GeV]");
	h_SGN_X->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_X->GetXaxis()->GetBinWidth(1)));
	h_SGN_X->GetXaxis()->SetTitleOffset(1.1); h_SGN_X->GetXaxis()->SetTitleSize(0.04);
	h_SGN_X->GetYaxis()->SetTitleSize(0.04);
	h_SGN_X->SetLineWidth(3);
	h_SGN_X->SetLineColor(kViolet + 1); h_SGN_X->SetFillColorAlpha(kViolet + 1, 0.30);
	h_Data_X->SetLineColor(kBlack);
	h_Data_X->SetMarkerStyle(20);
	h_Data_X->SetLineWidth(2);
	// ---> K0short 
	//h_SGN_K0s->SetTitle(Form("K0s %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_K0s->GetXaxis()->SetTitle("M(K^{0}_{s})[GeV]");
	h_SGN_K0s->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_K0s->GetXaxis()->GetBinWidth(1)));
	h_SGN_K0s->GetXaxis()->SetTitleOffset(1.1); h_SGN_K0s->GetXaxis()->SetTitleSize(0.04);
	h_SGN_K0s->GetYaxis()->SetTitleSize(0.04);
	h_SGN_K0s->SetLineWidth(3);
	h_SGN_K0s->SetLineColor(kGreen + 1); h_SGN_K0s->SetFillColorAlpha(kGreen + 1, 0.30);
	h_Data_K0s->SetLineColor(kBlack);
	h_Data_K0s->SetMarkerStyle(20);
	h_Data_K0s->SetLineWidth(2);
	// ---> K0short + Pi1 
	//h_SGN_K0s->SetTitle(Form("K0s %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_K0sPi->GetXaxis()->SetTitle("M(K^{0}_{s}#pi_{1/2})[GeV]");
	h_SGN_K0sPi->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_K0sPi->GetXaxis()->GetBinWidth(1)));
	h_SGN_K0sPi->GetXaxis()->SetTitleOffset(1.1); h_SGN_K0sPi->GetXaxis()->SetTitleSize(0.04);
	h_SGN_K0sPi->GetYaxis()->SetTitleSize(0.04);
	h_SGN_K0sPi->SetLineWidth(3);
	h_SGN_K0sPi->SetLineColor(kGreen + 1); h_SGN_K0sPi->SetFillColorAlpha(kGreen + 1, 0.30);
	h_Data_K0sPi->SetLineColor(kBlack);
	h_Data_K0sPi->SetMarkerStyle(20);
	h_Data_K0sPi->SetLineWidth(2);
	// ---> Rho
	//h_SGN_Rho->SetTitle(Form("Rho %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_Rho->GetXaxis()->SetTitle("M(\\rho)\\ [GeV]");
	h_SGN_Rho->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_Rho->GetXaxis()->GetBinWidth(1)));
	h_SGN_Rho->SetLineWidth(3);
	h_SGN_Rho->SetLineColor(kOrange + 1); h_SGN_Rho->SetFillColorAlpha(kOrange + 1, 0.30);
	h_Data_Rho->SetLineColor(kBlack);
	h_Data_Rho->SetMarkerStyle(20);
	// ---> X vs Rho
	h_SGN_XvsRho->GetXaxis()->SetTitle("M(\\rho)\\ [GeV]");
	h_SGN_XvsRho->GetYaxis()->SetTitle("M(X)\\ [GeV]");
	h_SGN_XvsRho->SetLineWidth(3);
	h_SGN_XvsRho->SetLineColor(kOrange + 1); h_SGN_XvsRho->SetFillColor(kOrange + 1);
	h_Data_XvsRho->SetLineColor(kBlack);
	h_Data_XvsRho->SetLineWidth(3);

	// get FIT
	TString inFileName = outPlot_path + "FitRes/FitBKG_X";
	if (Xcut < 0.) inFileName.Append(Form("m"));
	inFileName.Append(Form("%.0f_MR%.0f.root", fabs(Xcut)*100., MRcut*1000.));
	TFile* inFile = new TFile(inFileName);
	RooPlot* FitPlot = (RooPlot*)inFile->Get("FitPlot");
	RooCurve* FitCurve = (RooCurve*)FitPlot->getObject(1);
	FitCurve->SetLineColor(kRed);

	// Legend 
	auto legendB0 = new TLegend(0.15, 0.70,.45,.83);
	legendB0->SetTextSize(0.035);
	legendB0->SetBorderSize(0);
	auto legendX= new TLegend(0.15, 0.70,.45,.80);
	legendX->SetTextSize(0.035);
	legendX->SetBorderSize(0);
	auto legendK0s= new TLegend(0.15, 0.70,.45,.80);
	legendK0s->SetTextSize(0.035);
	legendK0s->SetBorderSize(0);
	auto legendRho= new TLegend(0.60, 0.80,.89,.89);
	legendRho->SetTextSize(0.025);
	legendRho->SetBorderSize(0);

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
	gStyle->SetLineWidth(3);
	h_SGN_B0->Scale(SGNfactor);
	h_SGN_B0->GetYaxis()->SetRangeUser(0., 1.2 * std::max(h_SGN_B0->GetMaximum(), h_Data_B0->GetMaximum()));
	h_SGN_B0->Draw("HISTE");
	legendB0->AddEntry(h_SGN_B0, "SIGNAL (MC)");
	h_Data_B0->Draw("PE0 SAME");
	legendB0->AddEntry(h_Data_B0, "DATA 2017");
	//FitCurve->Draw("SAME");
	//legendB0->AddEntry(FitCurve, "SIDEBANDS-FIT");

	legendB0->Draw();
	TString outPath = outPlot_path + "CutOptPlots/V2_";
	gPad->SetLeftMargin(0.12), gPad->SetBottomMargin(0.12);
	gPad->Update();
	CMSxxx(c1);
	c1->SaveAs(outPath + "B0massCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	c1->SaveAs(outPath + "B0massCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".pdf");

	h_SGN_X->Scale(SGNfactor);
	h_SGN_X->GetYaxis()->SetRangeUser(0.,1.4 * std::max(h_SGN_X->GetMaximum(), h_Data_X->GetMaximum()));
	h_SGN_X->Draw("HISTE");
	legendX->AddEntry(h_SGN_X, "SIGNAL (MC)");
	h_Data_X->Draw("PE0 SAME");
	legendX->AddEntry(h_Data_X, "DATA 2017");
	legendX->Draw();
	gPad->SetLeftMargin(0.12), gPad->SetBottomMargin(0.12);
	gPad->Update();
	CMSxxx(c1);
	c1->SaveAs(outPath + "X3872massCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	c1->SaveAs(outPath + "X3872massCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".pdf");

	h_SGN_K0s->Scale(SGNfactor);
	h_SGN_K0s->GetYaxis()->SetRangeUser(0.,1.2 * std::max(h_SGN_K0s->GetMaximum(), h_Data_K0s->GetMaximum()));
	h_SGN_K0s->Draw("HISTE");
	legendK0s->AddEntry(h_SGN_K0s, "SIGNAL (MC)");
	h_Data_K0s->Draw("PE0 SAME");
	legendK0s->AddEntry(h_Data_K0s, "DATA 2017");
	legendK0s->Draw();
	gPad->SetLeftMargin(0.12), gPad->SetBottomMargin(0.12);
	gPad->Update();
	CMSxxx(c1);
	c1->SaveAs(outPath + "K0smassCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	c1->SaveAs(outPath + "K0smassCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".pdf");

	h_SGN_K0sPi->Scale(SGNfactor);
	h_SGN_K0sPi->GetYaxis()->SetRangeUser(0.,1.2 * std::max(h_SGN_K0sPi->GetMaximum(), h_Data_K0sPi->GetMaximum()));
	h_SGN_K0sPi->Draw("HISTE");
	h_Data_K0sPi->Draw("PE0 SAME");
	legendK0s->Draw();
	gPad->SetLeftMargin(0.12), gPad->SetBottomMargin(0.12);
	gPad->Update();
	CMSxxx(c1);
	c1->SaveAs(outPath + "K0sPi12massCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	h_SGN_Rho->Scale(SGNfactor);
	h_SGN_Rho->GetYaxis()->SetRangeUser(0.,1.2 * std::max(h_SGN_Rho->GetMaximum(), h_Data_Rho->GetMaximum()));
	h_SGN_Rho->Draw("HISTE");
	legendRho->AddEntry(h_SGN_Rho, "SIGNAL (MC)");
	h_Data_Rho->Draw("PE0 SAME");
	legendRho->AddEntry(h_Data_Rho, "DATA 2017");
	legendRho->Draw();
	CMSxxx(c1);
	c1->SaveAs(outPath + "RhomassCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	c1->SaveAs(outPath + "RhomassCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".pdf");

	//h_SGN_XvsRho->Scale(SGNfactor);
	//h_SGN_XvsRho->Draw("BOX");
	h_Data_XvsRho->Draw("COLZ0");
	//legendRho->Draw();
	//c1->SaveAs(outPath + "XvsRhomassCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	

	inFile->Close();
	delete h_Data_B0;
	delete h_SGN_B0;
	delete h_Data_X;
	delete h_SGN_X;
	delete h_Data_K0s;
	delete h_SGN_K0s;
	delete h_Data_Rho;
	delete h_SGN_Rho;

	return 0;
}


int CutterMVA::makeMASSplot2D(const double Xcut, const double MRcut){

	int Nbins = 35;
	double Mlow_B0 = 4.90 , Mhigh_B0 = 5.60;
	double Mlow_X  = 3.80 , Mhigh_X  = 3.96;
	double Mlow_K0s= 0.45 , Mhigh_K0s= 0.55; 
	TH2F* h_XvsB0_m       = new TH2F( "XvsB0_m"     , "", Nbins, Mlow_B0 , Mhigh_B0,  Nbins, Mlow_X  , Mhigh_X  );
	TH2F* h_XvsB0_m_bkg   = new TH2F( "XvsB0_m_bkg" , "", Nbins, Mlow_B0 , Mhigh_B0,  Nbins, Mlow_X  , Mhigh_X  );
	TH2F* h_K0svsB0_m     = new TH2F("K0svsB0_m"    , "", Nbins, Mlow_B0 , Mhigh_B0,   Nbins, Mlow_K0s, Mhigh_K0s); 
	TH2F* h_K0svsB0_m_bkg = new TH2F("K0svsB0_m_bkg", "", Nbins, Mlow_B0 , Mhigh_B0,   Nbins, Mlow_K0s, Mhigh_K0s); 
	TH2F* h_XvsK0s_m      = new TH2F("XvsK0s_m"     , "", Nbins, Mlow_K0s, Mhigh_K0s , Nbins, Mlow_X, Mhigh_X); 
	TH2F* h_XvsK0s_m_bkg  = new TH2F("XvsK0s_m_bkg" , "", Nbins, Mlow_K0s, Mhigh_K0s , Nbins, Mlow_X, Mhigh_X); 

	TString CutOpt = Form("TreeBDTx_S.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	TString CutOptB = Form("TreeBDTx_B.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	TreeSGN->Draw("M_X3872:M_B0>>XvsB0_m",  CutOpt); 
	TreeBKG->Draw("M_X3872:M_B0>>XvsB0_m_bkg",  CutOptB); 
	h_XvsB0_m->SetFillColor(kAzure+1);  h_XvsB0_m_bkg->SetFillColor(kBlack);
	h_XvsB0_m_bkg->GetXaxis()->SetTitle("M(B_0)\\ [GeV]");
	h_XvsB0_m_bkg->GetYaxis()->SetTitle("M(X)\\ [GeV]");
	TreeSGN->Draw("M_K0s:M_B0>>K0svsB0_m", CutOpt); 
	TreeBKG->Draw("M_K0s:M_B0>>K0svsB0_m_bkg", CutOptB); 
	h_K0svsB0_m_bkg->SetFillColor(kGreen); h_K0svsB0_m_bkg->SetFillColor(kBlack);
	h_K0svsB0_m_bkg->GetXaxis()->SetTitle("M(B_0)\\ [GeV]");
	h_K0svsB0_m_bkg->GetYaxis()->SetTitle("M(K_0^s)\\ [GeV]");
	TreeSGN->Draw("M_X3872:M_K0s>>XvsK0s_m", CutOpt); 
	TreeBKG->Draw("M_X3872:M_K0s>>XvsK0s_m_bkg", CutOptB); 
	h_XvsK0s_m->SetFillColor(kViolet);  h_XvsK0s_m_bkg->SetFillColor(kBlack);
	h_XvsK0s_m_bkg->GetXaxis()->SetTitle("M(K_0^s)\\ [GeV]");
	h_XvsK0s_m_bkg->GetYaxis()->SetTitle("M(X)\\ [GeV]");

	
	TString outPath = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/";
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	TPad* pad = new TPad("pad", "", 0.,0.,1., 1.);
	pad->SetLeftMargin(0.15); pad->SetBottomMargin(0.15); 
	pad->Draw();
	pad->cd();
	//h_XvsB0_m->Draw("BOX"); 
	h_XvsB0_m_bkg->Draw("COLZ0");
	pad->Update();
	c1->SaveAs(outPath + "M_X3872vsM_B0" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	//h_K0svsB0_m->Draw("BOX"); 
	gStyle->SetPalette(kViridis);
	h_K0svsB0_m_bkg->Draw("COLZ0");
	pad->Update();
	c1->SaveAs(outPath + "M_K0svsM_B0" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	//h_XvsK0s_m->Draw("BOX"); 
	gStyle->SetPalette(kBird);
	h_XvsK0s_m_bkg->Draw("COLZ0");
	pad->Update();
	c1->SaveAs(outPath + "M_X3872vsM_K0s" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	delete c1;
	return 0;
}//makeMASSplot2D



void CutterMVA::Scan_X3872vsB0(const double Xcut, const double MRcut){

	bool maybeKstar = false;

	TString VarName[] = {"D0_Rho", "pT_Rho", "pT_Pi1", "DR_Pi1B0", "CosAlpha_B0", "LxySign_B0", "SVprob", "pTM_B0", "M_K0s", "M_Rho", "M_B0", "M_X3872"};
	TString VarLabel[] = {"D0(#pi_{1})/#sigma", "p_{T}(#pi#pi)/p_{T}(B_{0})", "p_{T}(#pi_{1})/p_{T}(B_{0})", "#Delta R(#pi_{1}, B_{0})", "cos(#alpha)", "L_{xy}/#sigma", "P(SV)", "p_{T}(B_{0})/M_{B_{0}}", "M(K_{s}^{0}) [GeV]", "M(#rho (770)) [GeV]", "M(B^{0}) [GeV]", "M(X_{3872}) [GeV]"};
	int    VarNbins[] = { 15,  10,  10,  10,   15,  20, 10, 14,   10,  10,   15,   15};
	double VarXlow[]  = { 0., .05,  0.,  0.,  .98,  0., 0., 0., 0.45, 0.6, 5.05,  3.8};
	double VarXhigh[] = {30., .30, .25,  .5, 1.01, 200, 1., 14, 0.55, 1.0, 5.20, 3.95};
	bool   VarNorm[]  = {true, true, true, true, true, true, true, true, false, false, false, false};
	const int Nvars = sizeof(VarName)/sizeof(TString); 
	std::cout << " --- # variables = " << Nvars << std::endl;

	TH1F* h_InVar_S = new TH1F("InVar_S", "", 100, 0., 500.);
	TH1F* h_InVar_B = new TH1F("InVar_B", "", 100, 0., 500.);
	double M_X3872_low = MX_nearLeft, M_B0_high = 5.2;
	TString Mcuts = Form("M_X3872>%f && M_B0<%f",M_X3872_low, M_B0_high);
	if (maybeKstar) Mcuts = Form("M_X3872>%f && M_X3872<%f && M_B0>5.08 && M_B0<5.18", MX_nearLeft, MX_nearRight); 
	TString CutOptS = Form("TreeBDTx_S.BDTx>%f && M_Rho>%f && ", Xcut, MRcut) + Mcuts;
	TString CutOptB = Form("TreeBDTx_B.BDTx>%f && M_Rho>%f && ", Xcut, MRcut) + Mcuts;
	TCanvas* c = new TCanvas("c","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	TString eosPath = "/eos/user/c/cbasile/www/B0toX3872K0s/FitResults/Unblind/XvsB_InVarsScan/";	
	if (maybeKstar) eosPath = "/eos/user/c/cbasile/www/B0toX3872K0s/FitResults/Unblind/B0around5_1GeV/";

	for(int i = 0; i < Nvars; i++){

		h_InVar_S->SetBins(VarNbins[i], VarXlow[i], VarXhigh[i]);	
		h_InVar_B->SetBins(VarNbins[i], VarXlow[i], VarXhigh[i]);
		TreeSGN->Draw(VarName[i]+">>InVar_S",  CutOptS); 
		TreeBKG->Draw(VarName[i]+">>InVar_B",  CutOptB); 
		if (VarNorm[i]){
			h_InVar_S->Scale(1./h_InVar_S->Integral());
			h_InVar_B->Scale(1./h_InVar_B->Integral());
		}
		if (!maybeKstar) h_InVar_S->SetTitle(Form("M_{X(3872)} > %.2f M_{B^{0}} < %.2f", M_X3872_low, M_B0_high));	
		h_InVar_S->GetXaxis()->SetTitle(VarLabel[i]);
		h_InVar_S->GetYaxis()->SetTitle(Form("Events/%.3f", h_InVar_S->GetXaxis()->GetBinWidth(1)));
		std::cout << Form("MaxS %.2f @ %d \t MaxB %.2f ",h_InVar_S->GetMaximum(), h_InVar_S->GetMaximumBin(), h_InVar_B->GetMaximum()) << std::endl;
		if (i < 9) h_InVar_S->SetMaximum(1.3 * std::max(h_InVar_S->GetMaximum(), h_InVar_B->GetMaximum()));
		else h_InVar_S->SetMaximum( 1.3 * h_InVar_B->GetMaximum());
		std::cout << Form("MaxFinale %.2f",h_InVar_S->GetMaximum()) << std::endl;
		h_InVar_S->SetLineWidth(3);	h_InVar_B->SetLineWidth(3);	
		h_InVar_S->SetLineColor(kAzure +1); h_InVar_S->SetFillColorAlpha(kAzure +1, 0.30);
		h_InVar_B->SetLineColor(kBlack); h_InVar_B->SetMarkerStyle(20);

		h_InVar_S->Draw("HIST");
		h_InVar_B->Draw("PE0 same");

		c->SaveAs(eosPath + VarName[i] +"_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
		c->SaveAs(eosPath + VarName[i] +"_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".pdf");
	}
	delete c;
}//Scan_X3872vsB0()


/*void CutterMVa::M_K0sPi(){


	TH1F* h_K0sPi1_dat= new TH1F("K0sPi1_dat", "", 30, 0.5, 2.0);
	Mcuts = Form("M_X3872>%f && M_X3872<%f && M_B0>5.08 && M_B0<5.18", MX_nearLeft, MX_nearRight); 
	TString CutOpt = Form("TreeBDTx_B.BDTx>%f && M_Rho>%f && ", Xcut, MRcut) + Mcuts;
	TCanvas* c = new TCanvas("c","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	TString eosPath = "/eos/user/c/cbasile/www/B0toX3872K0s/FitResults/Unblind/B0around5_1GeV/";	
	TreeBKG->Draw("M_K0sPi1>>K0sPi1_dat",  CutOpt); 

	h_K0sPi1_dat->GetYaxis()->SetTitle(Form("Events/%.3f", h_K0sPi1_dat->GetXaxis()->GetBinWidth(1)));
	std::cout << Form("MaxS %.2f @ %d \t MaxB %.2f ",h_K0sPi1_dat->GetMaximum(), h_K0sPi1_dat->GetMaximumBin(), h_InVar_B->GetMaximum()) << std::endl;
	if (i < 9) h_K0sPi1_dat->SetMaximum(1.3 * std::max(h_K0sPi1_dat->GetMaximum(), h_InVar_B->GetMaximum()));
	else h_K0sPi1_dat->SetMaximum( 1.3 * h_InVar_B->GetMaximum());
	std::cout << Form("MaxFinale %.2f",h_K0sPi1_dat->GetMaximum()) << std::endl;
	h_K0sPi1_dat->SetLineWidth(3);	h_InVar_B->SetLineWidth(3);	
	h_K0sPi1_dat->SetLineColor(kAzure +1); h_K0sPi1_dat->SetFillColorAlpha(kAzure +1, 0.30);
	h_InVar_B->SetLineColor(kBlack); h_InVar_B->SetMarkerStyle(20);

	h_K0sPi1_dat->Draw("HIST");






}//M_K0sPi()
*/
void CutterMVA::GetControlPlots(const double Xcut, const double MRcut){

	int Nbins = 100;
	double Mlow_Rho = 0.4 , Mhigh_Rho= 0.9;
	double low_X = -0.80 , high_X  = 0.2;

	TH1F* h_MrhoCutS= new TH1F("MrhoCutS", "", Nbins, Mlow_Rho, Mhigh_Rho);
	TH1F* h_MrhoCutB= new TH1F("MrhoCutB", "", Nbins, Mlow_Rho, Mhigh_Rho);
	TH1F* h_BDTxCutS= new TH1F("BDTxCutS", "", Nbins, low_X, high_X);
	TH1F* h_BDTxCutB= new TH1F("BDTxCutB", "", Nbins, low_X, high_X);

	TString CutOptS = Form("TreeBDTx_S.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	TString CutOptB = Form("TreeBDTx_B.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	// Mrho
	TreeSGN->Draw("M_Rho>>MrhoCutS",  CutOptS); 
	TreeBKG->Draw("M_Rho>>MrhoCutB",  CutOptB); 
	h_MrhoCutS->GetXaxis()->SetTitle("M(\\rho(770))\\ GeV");
	h_MrhoCutS->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_MrhoCutS->GetXaxis()->GetBinWidth(1)));
	h_MrhoCutS->SetLineWidth(3);	
	h_MrhoCutS->SetLineColor(kAzure +1); h_MrhoCutS->SetFillColorAlpha(kAzure +1, 0.30);
	h_MrhoCutB->SetLineColor(kRed); h_MrhoCutB->SetFillColorAlpha(kRed, 0.30);
	// X BDT
	TreeSGN->Draw("TreeBDTx_S.BDTx>>BDTxCutS",  CutOptS); 
	TreeBKG->Draw("TreeBDTx_B.BDTx>>BDTxCutB",  CutOptB); 
	h_BDTxCutS->GetXaxis()->SetTitle("BDT response [x]");
	h_BDTxCutS->GetYaxis()->SetTitle(Form("Events/%.3f ", h_BDTxCutS->GetXaxis()->GetBinWidth(1)));
	h_BDTxCutS->SetLineWidth(3);	
	h_BDTxCutS->SetLineColor(kAzure +1); h_BDTxCutS->SetFillColorAlpha(kAzure +1, 0.30);
	h_BDTxCutB->SetLineColor(kRed); h_BDTxCutB->SetFillColorAlpha(kRed, 0.30);


	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	h_MrhoCutS->Draw();
	h_MrhoCutB->Draw("same");

	c1->SaveAs(outPlot_path + "M_RhoCUT" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	h_BDTxCutS->Draw();
	h_BDTxCutB->Draw("same");
	c1->SaveAs( outPlot_path + "BDTxCUT" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
}//GetControlPlots()

void CutterMVA::CMSxxx(TCanvas* c){
	c->cd();
	TLatex RunDetails; RunDetails.SetNDC(); 
	RunDetails.SetTextFont(61);
	RunDetails.SetTextAlign(10);
	RunDetails.SetTextSize(0.035);
	RunDetails.DrawLatex(.12, .92, "CMS");
	RunDetails.SetTextFont(52);
	RunDetails.SetTextSize(0.035);
	RunDetails.DrawLatex(.20, .91, "Work in progress");
	RunDetails.SetTextFont(42);
	RunDetails.SetTextSize(0.030);
	RunDetails.DrawLatex(.70, .91, "41 fb^{-1} (13 TeV)");

}
