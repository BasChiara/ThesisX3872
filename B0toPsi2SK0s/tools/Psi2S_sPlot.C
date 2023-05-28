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
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "RooStats/SPlot.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"

#include "TSystem.h"

#include <iostream>


// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*, float, float);

void AddData(RooWorkspace*, float, float);
void DoSPlot(RooWorkspace*);
void MakePlots(RooWorkspace*);
void MakeHistos(RooWorkspace*);


void Psi2S_sPlot(){

// set range of observable
  Float_t lowMB  = 4.9; 
  Float_t highMB = 5.6;  
  
  // Create a new workspace to manage the project.
  RooWorkspace* wspace = new RooWorkspace("myWS");
  
  // add the signal and background models to the workspace.
  // Inside this function you will find a description our model.
  AddModel(wspace, lowMB, highMB);
  
  // add dataset from converted root tree
  //getDataSet("/eos/cms/store/user/crovelli/LowPtEle/TnpData/Sept/Jan16/Formatted_Parking_Run2018ALL_probeLowPt__tagIdCutsAt0.root", wspace, lowRange, highRange);
  AddData(wspace, lowMB, highMB);
  // inspect the workspace if you wish
  wspace->Print();
  
  // make a new dataset with sWeights added for every event.
  DoSPlot(wspace);

  // Make some plots showing the discriminating variable and
  // the control variable after unfolding.
  MakePlots(wspace);
  
  // Save variables in histos
  MakeHistos(wspace);

  // cleanup
  delete wspace;


}// Psi2S_sPlot()

// Signal and background fit models
void AddModel(RooWorkspace* ws, float low, float high){

    // make a RooRealVar for the observables
    // discriminant observable --> M(B0)
    RooRealVar B0m("B0m", "B0 reco mass", low, high, "GeV");
    
    // 
    // signal model for M(B0)
    cout << " ... signal model for M(B0)" << endl;
    // CrystalBall + Gauss with the same mean value
    RooRealVar B0mass ("B0mass" , "", 5.279, 5.275, 5.285);
    // Crystal Ball
    RooRealVar B0width("B0width", "", 0.01, 0.  , 0.05);
	RooRealVar alpha("alpha", "", 2., 0.,  10.);
	RooRealVar     N(    "N", "",100., 50. ,  200);
	
    RooCBShape CBsignal("CBsignal", "", B0m, B0mass, B0width, alpha, N);
	RooProdPdf SGNpdfCB("SGNpdfCB", "SGNpdfCB", CBsignal ); 
	// Gaussian
	RooRealVar B0sigma("B0sigma", "", 0.02, 0.0  , 2.);
	
    RooGaussian Gsignal("Gsignal", "", B0m, B0mass, B0sigma);
	
    RooRealVar f("f", "", 0., 1.);
    RooAddPdf SGNmodel("SGNmodel", "SGNmodel", RooArgList(CBsignal, Gsignal), f);

    // background model for M(B0)
    cout << " ... background model for M(B0)" << endl;
    RooRealVar SlopeP("SlopeP", "", -19.78, -25.0 , -10.);
	RooRealVar Min("Min", "", 4.8, 4. , 5.);
	RooRealVar C("C", "", 0.5, 0. , 1.);
	RooRealVar EXP("EXP", "", 4.);
	
	RooGenericPdf Pois("Pois", "", "pow((@0 - @1), @4) * exp(( @0 - @1)*@2) + @3", RooArgList(B0m, Min,  SlopeP, C, EXP));

    // combined signal + background model
    cout << " ... full signal + background model" << endl;
    RooRealVar sgn_yield("sgn_yield", "", 1000, 500, 10000);
	RooRealVar bkg_yield("bkg_yield", "", 1000, 800, 10000);

    RooAddPdf FULLmodel("FULLmodel", "FULLmodel", 
                        RooArgList(Pois, SGNmodel), 
                        RooArgList(bkg_yield,sgn_yield)); 
	
    // import model in the workspace
    cout << " --> import the full model" << endl;
    ws->import(FULLmodel);


}// AddModel()

void AddData(RooWorkspace* ws, float Mlow, float Mhigh){

    Float_t X_cut = -1.1;
    Float_t Mpipi_cut = 0.01;
    float run, LumiBlock, event;
    float pTM_B0, SVprob, LxySign_B0, CosAlpha_B0, DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho;
    float M_Rho, M_B0, M_MuMu, M_X3872, M_K0s;
    float X;

    // ===== DATA TREE =====

    TFile* InFileData = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/data/merged/SGNPsi2S_MergData17.root");
    TTree* TreeData = (TTree*)InFileData->Get("B0_Psi2Ssignal");

    TreeData->SetBranchAddress( "run", &run);
    TreeData->SetBranchAddress( "LumiBlock", &LumiBlock);
    TreeData->SetBranchAddress( "event", &event);

    TreeData->SetBranchAddress( "pTM_B0", &pTM_B0);
    TreeData->SetBranchAddress( "SVprob", &SVprob);
    TreeData->SetBranchAddress( "LxySign_B0", &LxySign_B0);
    TreeData->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
    TreeData->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
    TreeData->SetBranchAddress( "pT_Pi1", &pT_Pi1);
    TreeData->SetBranchAddress( "pT_Rho", &pT_Rho);
    TreeData->SetBranchAddress( "D0_Rho", &D0_Rho);
    TreeData->SetBranchAddress( "M_Rho", &M_Rho);
    TreeData->SetBranchAddress( "M_B0", &M_B0);
    TreeData->SetBranchAddress( "M_mumu", &M_MuMu);
    TreeData->SetBranchAddress( "M_X3872", &M_X3872);
    TreeData->SetBranchAddress( "M_K0s", &M_K0s);

    TFile* outputBDT = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/BDTonDATA.root");
    TTree* Tree_xB= (TTree*)outputBDT->Get("TreeBDTx_B");
    Tree_xB->SetBranchAddress( "BDTx", &X);
    TreeData->AddFriend("TreeBDTx_B");


    // get the what is needed
    RooRealVar *B0m = ws->var("B0m");
    // control observable BDT output
    RooRealVar BDTout("BDTout", "BDT output", -1.0, 1.0, "");

    RooDataSet data("data", "data", RooArgSet(*B0m, BDTout));

    bool isOutRange = true;
	int N = TreeData->GetEntriesFast();
	for (int i = 0; i < N ; i++){
		TreeData->GetEntry(i);
		isOutRange = (M_B0 < Mlow) ||(M_B0 > Mhigh);
		if( ( X < X_cut) || ( M_Rho < Mpipi_cut) || isOutRange) continue;
		*B0m = M_B0;
        BDTout = X;
		data.add(RooArgSet(*B0m, BDTout));
	}
	data.Print("v");
    // import in the workspace
    ws->import(data);

    InFileData->Close();
    outputBDT->Close();

}// AddData()



void DoSPlot(RooWorkspace* ws){
  
  std::cout << "Calculate sWeights" << std::endl;
  
  // get what we need out of the workspace to do the fit
  RooAbsPdf* model = ws->pdf("FULLmodel");
  RooRealVar* sgn_yield = ws->var("sgn_yield");
  RooRealVar* bkg_yield = ws->var("bkg_yield");
  RooDataSet* data = (RooDataSet*) ws->data("data");

  // fit the model to the data.
  model->fitTo(*data, Extended());

  //TCanvas* cdata = new TCanvas("cdata","data fit", 1);
  //RooRealVar* pair_mass = ws->var("pair_mass");
  //RooPlot* frame = pair_mass->frame() ;
  //data->plotOn(frame) ;
  //model->plotOn(frame) ;
  //model->plotOn(frame,Components(*jpsiModel),LineStyle(kDashed), LineColor(kRed)) ;
  //model->plotOn(frame4,Components(*bkgModel),LineStyle(kDashed),LineColor(kGreen)) ;
  //frame->SetTitle("Fit of model to discriminating variable");
  //frame->Draw() ;
  //cdata->SaveAs("Fit.png");

  // The sPlot technique requires that we fix the parameters
  // of the model that are not yields after doing the fit.
  // This *could* be done with the lines below, however this is taken care of
  // by the RooStats::SPlot constructor (or more precisely the AddSWeight method).
  RooRealVar* B0mass  = ws->var("B0mass");
  RooRealVar* B0width = ws->var("B0width");
  RooRealVar* alpha = ws->var("alpha");
  RooRealVar* N = ws->var("N");
  RooRealVar* B0sigma = ws->var("B0sigma");
  RooRealVar* SlopeP = ws->var("SlopeP");
  RooRealVar* Min  = ws->var("Min");
  RooRealVar* C  = ws->var("C");
  B0mass->setConstant();   
  B0width->setConstant(); 
  alpha->setConstant();   
  N->setConstant();   
  B0sigma->setConstant(); 
  SlopeP->setConstant();   
  Min->setConstant();   
  C->setConstant();   

  RooMsgService::instance().setSilentMode(true);

  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",
					       *data, model, RooArgList(*sgn_yield,*bkg_yield) );


  // Check that our weights have the desired properties
  std::cout << "Check SWeights:" << std::endl;

  std::cout << std::endl <<  "Yield of B0 is "
	    << sgn_yield->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("sgn_yield") << std::endl;
  
  std::cout << "Yield of background is "
	    << bkg_yield->getVal() << ".  From sWeights it is "
	    << sData->GetYieldFromSWeight("bkg_yield") << std::endl
	    << std::endl;

  std::cout << "sWeights for the first 10 events \n" << std::endl;
  for(Int_t i=0; i < 10; i++) {
      std::cout << "B0 Weight = " << sData->GetSWeight(i,"sgn_yield")
		<< ", bkg Weight = " << sData->GetSWeight(i,"bkg_yield")
		<< ", Total Weight = " << sData->GetSumOfEventSWeight(i)
		<< std::endl;
    }
  
  std::cout << std::endl;

  // import this new dataset with sWeights
  std::cout << "import new dataset with sWeights" << std::endl;
  ws->import(*data, Rename("dataWithSWeights"));
}

// Control plots
void MakePlots(RooWorkspace* ws){
  
  // Here we make plots of the discriminating variable (B0m) after the fit
  // and of the control variable (BDToutput) after unfolding with sPlot.
  std::cout << std::endl;
  std::cout << "make plots" << std::endl;

  // make our canvas
  TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 800, 1200);
  cdata->Divide(1,3);

  // get what we need out of the workspace
  RooAbsPdf* model = ws->pdf("FULLmodel");
  RooAbsPdf* SGNmodel = ws->pdf("SGNmodel");
  RooAbsPdf* BKGmodel = ws->pdf("Pois");
  
  RooRealVar* BDTout = ws->var("BDTout");
  RooRealVar* B0m = ws->var("B0m");

  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");
  
  // this shouldn't be necessary, need to fix something with workspace
  // do this to set parameters back to their fitted values.
  model->fitTo(*data, Extended() );

  //plot B0 mass for data with full model and individual components overlaid
  cdata->cd(1);
  RooPlot* frame = B0m->frame() ;
  data->plotOn(frame ) ;
  model->plotOn(frame) ;
  model->plotOn(frame,Components(*SGNmodel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame,Components(*BKGmodel),LineStyle(kDashed),LineColor(kGreen)) ;
  frame->SetTitle("Fit of model to discriminating variable");
  frame->Draw() ;


  // Now use the sWeights to show our variable distribution for B0 and background.
  //
  // Plot our variable for B0 component.
  // Do this by plotting all events weighted by the sWeight for the JPsi component.
  // The SPlot class adds a new variable that has the name of the corresponding
  // yield + "_sw".
   cdata->cd(2);
   
  // create weighted data set
  RooDataSet * dataw_B0 = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"sgn_yield_sw") ;
  
  RooPlot* frame2 = BDTout->frame();
  // RooPlot* frame2 = probePfmvaId->frame() ;
  dataw_B0->plotOn(frame2, DataError(RooAbsData::SumW2) ) ;
  frame2->SetTitle("M(B0) distribution with sWeights B0");
  frame2->Draw() ;

  // Plot interesting variables for background
  cdata->cd(3);
  RooDataSet * dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkg_yield_sw") ;
  RooPlot* frame3 = BDTout->frame() ;
  // RooPlot* frame3 = probePfmvaId->frame() ;
  dataw_bkg->plotOn(frame3,DataError(RooAbsData::SumW2) ) ;
  frame3->SetTitle("M(B0) distribution for background");
  frame3->Draw() ;
  
  cdata->SaveAs("/eos/user/c/cbasile/www/B0toX3872K0s/Psi2S/MVAcuts/BEST_CUT/BDTout_SPlot.png");


  // Fit variable
  TCanvas* cdata2 = new TCanvas("cdata2","data fit", 1);
  RooPlot* frame4 = B0m->frame() ;
  data->plotOn(frame4 ) ;
  model->plotOn(frame4) ;
  model->plotOn(frame4,Components(*SGNmodel),LineStyle(kDashed), LineColor(kRed)) ;
  model->plotOn(frame4,Components(*BKGmodel),LineStyle(kDashed),LineColor(kGreen)) ;
  frame4->SetTitle("Fit of model to discriminating variable");
  frame4->Draw() ;
  cdata2->SaveAs("/eos/user/c/cbasile/www/B0toX3872K0s/Psi2S/MVAcuts/BEST_CUT/B0massPsi2S_Fit.png");
}//MakePlots()

void MakeHistos(RooWorkspace* ws){

	gStyle->SetOptStat(0);
    //Int_t status = gSystem->Load("../lib/RootPlots.o");
    //std::cout << "status: " << status << std::endl;
	std::cout << std::endl;
	std::cout << "save histos" << std::endl;
    Float_t X_cut = -1.1;
    Float_t Mpipi_cut = 0.01;

	RooRealVar* BDTout = ws->var("BDTout");
	// note, we get the dataset with sWeights
	RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

	RooDataSet * dataw_B0 = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"sgn_yield_sw") ;
	RooDataSet * dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkg_yield_sw") ;

	TH1 *h_BDTout_sgnData = dataw_B0->createHistogram("BDTout_sgn",*BDTout,Binning(40));

    TFile* InFileMC = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/SGN_MC.root");
    TTree* TreeMC = (TTree*)InFileMC->Get("mcSIGNAL");
    TFile* outputBDT = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/BDTonSGN.root");
    TTree* Tree_xS= (TTree*)outputBDT->Get("TreeBDTx_S");
    TreeMC->AddFriend("TreeBDTx_S");

	TH1* h_BDTout_sgnMC  = new TH1F("BDTout_sgnMC" , "", 40, -1., 1.);	
	TreeMC->Draw("BDTx>>BDTout_sgnMC",Form("(M_B0 > 4.9)&&(M_B0 < 5.6)&&(BDTx>%f)&&(M_Rho>%f)",X_cut, Mpipi_cut));

		
	h_BDTout_sgnMC->Scale(3.786);	
	std::cout << " MC integral prima " << h_BDTout_sgnMC->Integral() << " vs " << h_BDTout_sgnData->Integral()  << std::endl;
    //myRootLib::histoSetUp(h_BDTout_sgnMC, "BDT output", Form("Events/%f", h_BDTout_sgnMC->GetBinWidth(2)), kRed, true, false);
	h_BDTout_sgnMC->SetLineColor(kRed); h_BDTout_sgnMC->SetFillColorAlpha(kRed, 0.3); h_BDTout_sgnMC->SetLineWidth(3);
	h_BDTout_sgnData->SetLineColor(kBlack); h_BDTout_sgnData->SetMarkerStyle(20); h_BDTout_sgnData->SetLineWidth(3);
	
    h_BDTout_sgnMC->Scale(1./h_BDTout_sgnMC->Integral());
    h_BDTout_sgnData->Scale(1./h_BDTout_sgnData->Integral());	
	std::cout << " MC integral dopo " << h_BDTout_sgnMC->Integral() << " vs " << h_BDTout_sgnData->Integral()  << std::endl;

    h_BDTout_sgnMC->SetMaximum(1.3*std::max(h_BDTout_sgnMC->GetMaximum(), h_BDTout_sgnData->GetMaximum()));
    h_BDTout_sgnMC->GetXaxis()->SetTitle("BDT output"); 
    h_BDTout_sgnMC->GetYaxis()->SetTitle(Form("Events/%.2f", h_BDTout_sgnMC->GetBinWidth(2)));


	TCanvas* c = new TCanvas("c","", 1024, 1024);
	h_BDTout_sgnMC->Draw("HIST");
	h_BDTout_sgnData->Draw("PE SAME");
	c->SaveAs("/eos/user/c/cbasile/www/B0toX3872K0s/Psi2S/MVAcuts/BEST_CUT/MCvsSPlot_BDTout.png");
    
    TFile* outPlotFile = new TFile("./outRoot/sPlot_MCpreUL.root", "RECREATE");
    outPlotFile->cd();
    h_BDTout_sgnMC->Write();
    h_BDTout_sgnData->Write();
    std::cout << " [OUTPUT] output saved in file " << outPlotFile->GetName() << std::endl;
    outPlotFile->Close();


    InFileMC->Close();
	outputBDT->Close();

}//MakeHistos()



