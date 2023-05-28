#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooCBShape.h"

#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooGenericPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "RooStats/SPlot.h"

#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"

#include "TSystem.h"

#include <map>
#include <iostream>


using namespace RooFit;
using namespace RooStats;
using namespace std;


void AddModel(RooWorkspace*, float, float);
void AddData(RooWorkspace* , float , float);
void FitToData(RooWorkspace* , float , float);

void Psi2S_fitter(){

    Float_t lowM  = 4.9; 
    Float_t highM = 5.6;     
    RooWorkspace* wspace = new RooWorkspace("myWS");

    AddModel(wspace, lowM, highM);
    AddData(wspace, lowM, highM);

    wspace->Print();

    FitToData(wspace, lowM, highM);

}


void AddModel(RooWorkspace* ws, float Mlow, float Mhigh){

	RooRealVar B0m("B0m", "B0 reco mass", Mlow, Mhigh, "GeV");
    cout << " ... signali+background model for M(B0) [B0 -> Psi(2S) K0s]" << endl;

	// --> BKG model
	// Poisson
	RooRealVar SlopeP("SlopeP", "", -19.78, -25.0 , -10.);
	RooRealVar Min("Min", "", 4.8, 4. , 5.);
	RooRealVar C("C", "", 0.001, 0. , 0.005);
	RooRealVar EXP("EXP", "", 4.);
	
	RooGenericPdf Pois("Pois", "", "pow((@0 - @1), @4) * exp(( @0 - @1)*@2) + @3", RooArgList(B0m, Min,  SlopeP, C, EXP));
	
	// --> SGN model

	// CrystalBall + Gauss
	RooRealVar B0mass ("B0mass" , "", 5.28, 5.279, 5.289);

	RooRealVar B0width("B0width", "", 0.005, 0.001  , 0.01);
	RooRealVar alpha("alpha", "", 2., 0.,  10.);
	RooRealVar     N(    "N", "",100., 50. ,  200);
	RooCBShape CBsignal("CBsignal", "", B0m, B0mass, B0width, alpha, N);
	RooProdPdf SGNpdfCB("SGNpdfCB", "SGNpdfCB", CBsignal ); 

	RooRealVar B0sigma("B0sigma", "", 0.02, 0.0  , 2.);
	RooGaussian Gsignal("Gsignal", "", B0m, B0mass, B0sigma);
	RooProdPdf SGNpdfG("SGNpdfG", "SGNpdfG", Gsignal ); 

	// PDF
	RooRealVar sgn_yield("nSGN", "", 1000, 500, 10000);
	RooRealVar bkg_yield("nBKG", "", 1000, 800, 10000);
	RooRealVar f("f", "", 0., 1.);
	RooRealVar BKGfrac("BKGfrac", "", 0.5, 0., 1.);
	RooAddPdf SGNmodel("SGNmodel", "SGNmodel", RooArgList(CBsignal, Gsignal), f); 
	RooAddPdf FULLmodel("FULLmodel", "FULLmodel", RooArgList(Pois, SGNmodel), RooArgList(bkg_yield,sgn_yield)); 
	
    // ... save on workspace
    cout << " --> import FULLmodel model" << endl;
    ws->import(FULLmodel);


}//AddModel()


void AddData(RooWorkspace* ws,  float Mlow, float Mhigh){

    // ===== DATA TREE =====
    // tree variables & path
    Float_t run, LumiBlock, event;
    Float_t pTM_B0, SVprob, LxySign_B0, CosAlpha_B0, DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho; 
    Float_t M_Rho, M_B0, M_MuMu, M_X3872, M_K0s;
    Float_t X;
    
    TString MainTreePath = "./SGNPsi2S_MergData17.root";
    TString TreeName = "B0_Psi2Ssignal"; 
    TString BDTTreePath = "./BDTonDATA.root";
    // open file 
    TFile* InFile(0);
    if (!gSystem->AccessPathName(MainTreePath)) {
        InFile = TFile::Open(MainTreePath); // check if file in local directory exists
    }else {
        cout << "ERROR: could not open data file " << MainTreePath  << endl;
        exit(1);
    }
    TTree* DataTree= (TTree*)InFile->Get(TreeName);

    DataTree->SetBranchAddress( "run", &run);
    DataTree->SetBranchAddress( "LumiBlock", &LumiBlock);
    DataTree->SetBranchAddress( "event", &event);

    DataTree->SetBranchAddress( "pTM_B0", &pTM_B0);
    DataTree->SetBranchAddress( "SVprob", &SVprob);
    DataTree->SetBranchAddress( "LxySign_B0", &LxySign_B0);
    DataTree->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
    DataTree->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
    DataTree->SetBranchAddress( "pT_Pi1", &pT_Pi1);
    DataTree->SetBranchAddress( "pT_Rho", &pT_Rho);
    DataTree->SetBranchAddress( "D0_Rho", &D0_Rho);
    DataTree->SetBranchAddress( "M_Rho", &M_Rho);
    DataTree->SetBranchAddress( "M_B0", &M_B0);
    DataTree->SetBranchAddress( "M_mumu", &M_MuMu);
    DataTree->SetBranchAddress( "M_X3872", &M_X3872);
    DataTree->SetBranchAddress( "M_K0s", &M_K0s);

    TFile* outputBDT = TFile::Open(BDTTreePath); 
    if(!outputBDT){
        cout << "ERROR: could not open BDT file " << BDTTreePath  << endl;
        exit(1);
    }
    TTree* TreeBDT= (TTree*)outputBDT->Get("TreeBDTx_B");
	TreeBDT->SetBranchAddress( "BDTx", &X);
	DataTree->AddFriend("TreeBDTx_B");

    // build the RooDataSet
    RooRealVar *B0m = ws->var("B0m");

    RooDataSet data("data","data", RooArgSet(*B0m));		

	// LOOP ON BACKGROUND EVENTS
	float Xcut = -0.08;                 // cut on BDT output
    float MRcut = 0.5;                  // cut on PiPi mass
	bool isOutRange = true;
	int Nbkg = DataTree->GetEntries();
	for (int i = 0; i < Nbkg ; i++){
		DataTree->GetEntry(i);
		isOutRange = (M_B0 < Mlow) ||(M_B0 > Mhigh);
		if( ( X < Xcut) || ( M_Rho < MRcut) || isOutRange) continue;
		//h_Data_B0->Fill(M_B0);
		*B0m = M_B0;
		data.add(RooArgSet(*B0m));
	}
	data.Print("v");

    // import in workspace
    ws->import(data);

}//AddData()


void FitToData(RooWorkspace* ws, float Mlow, float Mhigh){

    cout << " ... start fitting data ..." << std::endl;
    
    RooAbsPdf* model = ws->pdf("FULLmodel");
    RooAbsPdf* SGNmodel = ws->pdf("SGNmodel");
    RooAbsPdf* BKGmodel = ws->pdf("Pois");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    RooFitResult* results = model->fitTo(*data, Save());
	results->Print("v");

    // save parameters values after the fit
    RooRealVar* B0m     = ws->var("B0m");
    RooRealVar* SlopeP  = ws->var("SlopeP");
    RooRealVar* Min     = ws->var("Min");
    RooRealVar* C       = ws->var("C");
    RooRealVar* B0mass  = ws->var("B0mass");
    RooRealVar* B0sigma = ws->var("B0sigma");
    RooRealVar* B0width = ws->var("B0width");
    RooRealVar* alpha   = ws->var("alpha");
    RooRealVar* N       = ws->var("N");

    // create the plot
    int Nbins = 70;
    RooPlot *frame = B0m->frame(Title("B0 -> Psi(2S) K0s"));
    data->plotOn(frame, Binning(Nbins));
    model->plotOn(frame, LineColor(kAzure +1));
    double Chi2 = frame->chiSquare(8);
    data->plotOn(frame, Binning(Nbins)); 

	TText *TxtChi2= new TText(5.31, frame->GetMaximum()*0.85, Form("ChiSq = %.3f", Chi2)); 
	TxtChi2->SetTextSize(0.035);   TxtChi2->SetTextColor(kRed);
    frame->addObject(TxtChi2);
    frame->Draw();
    cout << "  == Chi2 : " << Chi2 << endl;

    // save results on file
    TString outFileParams = "./Psi2S_fitterOUT.root";
    TFile* outFile = new TFile(outFileParams, "RECREATE");
    frame->Write("FitPlot");
    results->Write("FitParams");
    outFile->Close();
    cout << " ... fit results saved in " << outFileParams << endl;

    // compute number of events
	B0m->setRange("SGNregion"   , 5.22, 5.35);           // SGNregion RANGE 
	std::cout << " TOTAL ENTRIES = " << data->sumEntries() << std::endl;
	
	double I_SGN_T  = SGNmodel->createIntegral(*B0m, NormSet(*B0m))->getVal(); //integrate sgn region
	double I_SGN_Sr = SGNmodel->createIntegral(*B0m, NormSet(*B0m), Range("SGNregion"))->getVal(); //integrate sgn region
	double nSGN = ws->var("nSGN")->getVal();
	double Ss = nSGN*I_SGN_Sr/I_SGN_T;
	std::cout << " nSGN = " << nSGN << " I(s | TOT ) = " << I_SGN_T << " I(s | SR ) = " << I_SGN_Sr << " --> #S in sgn region " << Ss << std::endl; 
	double I_BKG_T  = BKGmodel->createIntegral(*B0m, NormSet(*B0m))->getVal(); //integrate sgn region
	double I_BKG_Sr = BKGmodel->createIntegral(*B0m, NormSet(*B0m), Range("SGNregion"))->getVal(); //integrate sgn region
	double nBKG =  ws->var("nBKG")->getVal();
	double Bs = nBKG*I_BKG_Sr/I_BKG_T;
	std::cout << " nBKG = " << nBKG << "  I(b | TOT ) = " << I_BKG_T << " I(b | SR ) = " << I_BKG_Sr << " --> #B in sgn region " << Bs << std::endl; 
    
    

}//FitToData()

