#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TLine.h>
#include <TMarker.h>
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <THStack.h>

TString FilePath_       = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/";
TString RootFileNameBDT = "TMVAoutDataBDTbest.root";
TString RootFileNameDNN = "TMVAoutDataDNNbest.root";
TString OutPlotFile     = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/plots/TMVA/";
TString DirBDT_0        = "datasetBDT";
TString subDirBDT_1     = "Method_BDT";
int N_BDT = 7;
TString subDirBDT_2[]   =  {"BDT_nT200_S25_D4_b03_nC35", "BDT_nT200_S25_D5_b03_nC30", "BDT_nT250_S25_D4_b03_nC35", "BDT_nT250_S25_D5_b03_nC30", "BDT_nT300_S25_D4_b03_nC35","BDT_nT300_S25_D3_b03_nC30", "BDT_nT200_S25_D3_b03_nC30"};
TString BDTlegend[]   =  {"N_{T}=200 #Phi=2.5 D_{M}=4 #gamma=0.3 #delta_{x}=35", "N_{T}=200 #Phi=2.5 D_{M}=5 #gamma=0.3 #delta_{x}=30", "N_{T}=250 #Phi=2.5 D_{M}=4 #gamma=0.3 #delta_{x}=35", "N_{T}=250 #Phi=2.5 D_{M}=5 #gamma=0.3 #delta_{x}=30", "N_{T}=300 #Phi=2.5 D_{M}=4 #gamma=0.3 #delta_{x}=35","N_{T}=300 #Phi=2.5 D_{M}=3 #gamma=0.3 #delta_{x}=30", "N_{T}=200 #Phi=2.5 D_{M}=3 #gamma=0.3 #delta_{x}=30"};
TString subDirBDT_2_1   = "BDT_nT100_b08_nC30"; TString subDirBDT_2_2   = "BDT_nT100_b08_nC30";


TString DirDNN_0        = "datasetDNN";
TString subDirDNN_1     = "Method_DL";
int N_DNN = 5;
TString subDirDNN_2[]   =  {"DNN_CPU_L2_LR1e4_E1000_P09_BS30", "DNN_CPU_L2_LR1e4_E1000_P07_BS30", "DNN_CPU_L0_LR1e4_E1000_P08_BS30", "DNN_CPU_L0_LR1e4_E1000_P07_BS30", "DNN_CPU_L1_LR1e4_E1000_P08_BS30"};
TString DNNlegend[]   =  {"Lyt.3|#eta=10^{-4} #varepsilon_{M}=1000 #beta_{1}=0.9 #Theta=30", "Lyt.3|#eta=10^{-4} #varepsilon_{M}=1000 #beta_{1}=0.7 #Theta=30", "Lyt.1|#eta=10^{-4} #varepsilon_{M}=1000 #beta_{1}=0.8 #Theta=30", "Lyt.1|#eta=10^{-4} #varepsilon_{M} =1000 #beta_{1}=0.7 #Theta=30", "Lyt.2|#eta=10^{-4} #varepsilon_{M}=1000 #beta_{1}=0.8 #Theta=30"};
TString subDirDNN_2_1   = "DNN_CPU";
TString subDirANN_1     = "Method_MLP";
TString subDirANN_2_1   = "MLP";

TFile* openTMVArootFile(TString method = "BDT"){

	TString FileToOpen = FilePath_ ;
	if (method == "BDT") FileToOpen.Append(RootFileNameBDT);
	if (method == "DNN") FileToOpen.Append(RootFileNameDNN);
	TFile* RootFile = new TFile(FileToOpen); 
	if ( !RootFile->IsOpen() ) {
		std::cout << "ERROR IN OPENING FILE "<< RootFile<< std::endl;
		exit(-1);
	}

	return RootFile;

}//openTMVArootFile()

TDirectoryFile* accesTMVAmethodDir(TFile* RootFile, TString method, int Ndir = 0){

	TDirectoryFile *subDir1   = new TDirectoryFile();
	TDirectoryFile *subDir2_1 = new TDirectoryFile();
	TDirectoryFile *dir = new TDirectoryFile(); 
	if (method == "BDT"){
		dir = (TDirectoryFile*)RootFile->Get(DirBDT_0);
		subDir1 = (TDirectoryFile*)dir->Get(subDirBDT_1);
		subDir2_1 = (TDirectoryFile*)subDir1->Get(subDirBDT_2[Ndir]);
	}
	if (method == "DNN"){
		dir = (TDirectoryFile*)RootFile->Get(DirDNN_0);
		subDir1 = (TDirectoryFile*)dir->Get(subDirDNN_1);
		subDir2_1 = (TDirectoryFile*)subDir1->Get(subDirDNN_2[Ndir]);
	}
	if (method == "MLP"){
		subDir1 = (TDirectoryFile*)dir->Get(subDirANN_1);
		subDir2_1 = (TDirectoryFile*)subDir1->Get(subDirANN_2_1);
	}

	return subDir2_1;

}//accesTMVAmethodDir()


Color_t CategoryColor(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["BDT"]           = kGreen; 
  Color["MLP"]           = kBlue;
  Color["DNN"]           = kOrange;
  Color["SGN_B0"]        = kAzure  + 1;
  Color["SGN_B0_pstcut"] = kAzure  - 9;
  Color["BKG_B0"]        = kRed;
  Color["BKG_B0_pstcut"] = kRed    - 9;
  Color["BKG_mc_B0"]     = kOrange + 1;

  Color["BKG_X3872"]     = kBlue   - 9;
  Color["BKG_K0s"]       = kTeal   + 2;
  Color["BKG_JPsi"]      = kRed    - 7;
  Color["SGN_Rho"]		 = kGreen  + 1;
  Color["BKG_Rho"]		 = kYellow - 7;

  return Color[category];
}//CategoryLegend()


TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["BDT"] = "BDT";
  Leg_entry["MLP"] = "ANN";
  Leg_entry["DNN"] = "DNN";
  Leg_entry["SGN_B0"] = "SIGNAL(MC)";
  Leg_entry["BKG_B0"] = "SIDE-BANDS(DATA)";
  Leg_entry["SGN_B0_precut"] = "SIGNAL-MC (PRE CUT)";
  Leg_entry["BKG_B0_precut"] = "DATA (PRE CUT)";
  Leg_entry["SGN_B0_pstcut"] = "SIGNAL-MC (POST CUT)";
  Leg_entry["BKG_B0_pstcut"] = "DATA (POST CUT)";
  Leg_entry["BKG_mc_B0"] = "MC non matching";

  return Leg_entry[category];
}//CategoryLegend()


void CMSxxx(TCanvas* c){
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
void histo_SetUp(TH1* histo, const TString& category, const TString& x_name, const TString& y_name, bool fill = true , bool norm = true){
  //AXIS LABEL
  histo->SetTitle("");
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.03);
  histo->GetXaxis()->SetTitleOffset(1.25);
  histo->GetYaxis()->SetTitle(y_name);
  histo->GetYaxis()->SetLabelSize(0.03);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.25);


  //WIDTH & COLOR                                                                                                                                                               
  histo->SetLineWidth(4);
  histo->SetLineColor(CategoryColor(category));
  if (fill)histo->SetFillColorAlpha(CategoryColor(category), 0.1);
  //NORMALIZATION
  if(norm) histo->Scale(1./histo->Integral());
}


int MVA_SignalEff(){ 

	TFile* RootFile = openTMVArootFile();
	TDirectoryFile* souceDir = accesTMVAmethodDir(RootFile, "BDT");

	TH1F* h_SGNeff = (TH1F*)souceDir->Get("MVA_" + subDirBDT_2[0] + "_effS");
	histo_SetUp(h_SGNeff, "SGN_B0","BDT-output (Y_{BDT})", "#varepsilon(SGN)", false, false); 
	h_SGNeff->GetYaxis()->SetNdivisions(11, 2, 0);

	// WORKING POINT
	int binEff = -1;
	float Eff = 0.9, INprecision = 0.0001, OUTgap = 10.;
	double OUTeff = 0., OUTx = 0.; 
	std::cout << "(efficiency required is " << Eff << " )" << std::endl;
	while(OUTgap > 1){
		std::cout << "precision = " << INprecision  << std::endl;
		OUTgap = h_SGNeff->GetBinWithContent(Eff, binEff, 0, h_SGNeff->GetNbinsX() + 1, INprecision);
		OUTeff = h_SGNeff->GetBinContent(binEff);
		OUTx = h_SGNeff->GetBinLowEdge(0) + (binEff + 0.5) * h_SGNeff->GetBinWidth(0);
		
		INprecision *= 10.;
	}

	std::cout << "	 ==> GAP " << OUTgap << " at x " << OUTx << " --> efficiency = " << h_SGNeff->GetBinContent(binEff) << std::endl; 	
	// OPTIMIZATION
	double Wx, X2, Deff2, Deff = Eff - OUTeff;
	if (Deff > 0){
		X2 = h_SGNeff->GetBinLowEdge(0) + (binEff - 0.5) * h_SGNeff->GetBinWidth(0); // previous bin-center
		Deff2 = Eff - h_SGNeff->GetBinContent(binEff - 1);
	} else {
		X2 = h_SGNeff->GetBinLowEdge(0) + (binEff + 1.5) * h_SGNeff->GetBinWidth(0); // previous bin-center
		Deff2 = Eff - h_SGNeff->GetBinContent(binEff + 1);
	}

	Wx = ( fabs(1./Deff) * OUTx + fabs(1./Deff2) * OUTx)/ (fabs(1./Deff) + fabs(1./Deff2) );
	//std::cout << " Next efficiency " << Deff2 << " @ X2 = " << X2  << " FINAL X --> " << Wx << std::endl;
	std::cout <<" FINAL X --> " << Wx << std::endl;

	TMarker EffPoint(OUTx,Eff,34);
	EffPoint.SetMarkerSize(3);
	TLine* lineX = new TLine( OUTx, 0, OUTx, Eff); 
	lineX->SetLineColor(kBlack); lineX->SetLineStyle(10);
	lineX->SetLineWidth(3);
	TLine* lineX2 = new TLine( Wx, 0, Wx, Eff); 
	lineX2->SetLineColor(kRed);
	TLine* lineY = new TLine( h_SGNeff->GetBinLowEdge(0) , Eff, OUTx, Eff); 
	lineY->SetLineColor(kBlack); lineY->SetLineStyle(10);
	lineY->SetLineWidth(3);

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	c1->SetGrid();
	gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
	h_SGNeff->Draw();
	lineX->Draw("same");	
	//lineX2->Draw("same");	
	lineY->Draw("same");	
	EffPoint.Draw();

	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(3);
	c1->SaveAs(OutPlotFile + subDirBDT_2[0] + "_SGNeffREDO.png");
	c1->SaveAs(OutPlotFile + subDirBDT_2[0] + "_SGNeffREDO.pdf");

	return 0;

}//MVA_SignalEff()


int MVAcutMB0_MRho(){
	
	int Nbins = 60;
	double M_low = 5.0 , M_high = 5.6;
	// PRE - CUT 	
	TFile* RootFileMC   = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAinputs.root");
	TTree* TreeMC = (TTree*)RootFileMC->Get("inputSIGNAL");
	TH1F* h_SGN_B0 = new TH1F("h_SGN_B0", "", Nbins, M_low, M_high); 
	TreeMC->Draw("M_B0>>h_SGN_B0");

	TFile* RootFileDATA = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data/merged/BKG_MergData17.root");
	TTree* Tree17 = (TTree*)RootFileDATA->Get("B0sidebands");
	Int_t N = Tree17->GetEntriesFast();
	TH1F* h_BKG_B0 = new TH1F("h_BKG_B0", "", Nbins, M_low,M_high );
	Tree17->Draw("M_B0>>h_BKG_B0");

	// POST - CUT
	TFile* RootFileCUT = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAcutApp.root");
	//TTree* TreeMC = (TTree*)RootFileMC->Get("");
	TH1F* h_SGN_B0cut = (TH1F*)RootFileCUT->Get("SGN_MB0"); 
	TH1F* h_BKG_B0cut = (TH1F*)RootFileCUT->Get("BKG_MB0"); 
	TH1F* h_SGN_Rhocut = (TH1F*)RootFileCUT->Get("SGN_MRho"); 
	TH1F* h_BKG_Rhocut = (TH1F*)RootFileCUT->Get("BKG_MRho"); 


	histo_SetUp(h_SGN_B0,    "SGN_B0", 			"M(B_{0}) [GeV]", Form("%s %4.3f", "dN/" ,h_SGN_B0->GetXaxis()->GetBinWidth(1)), true, false);
	h_SGN_B0->SetFillColor(CategoryColor("SGN_B0"));
	histo_SetUp(h_BKG_B0,    "BKG_B0", 			"M(B_{0}) [GeV]", Form("%s %4.3f %s", "dN/" ,h_SGN_B0->GetXaxis()->GetBinWidth(1), " [GeV]"), true, false);
	h_BKG_B0->SetFillColor(CategoryColor("BKG_B0"));
	histo_SetUp(h_SGN_B0cut, "SGN_B0_pstcut", "M(B_{0}) [GeV]", Form("%s %4.3f", "dN/" ,h_SGN_B0cut->GetXaxis()->GetBinWidth(1)), true, false);
	h_SGN_B0cut->SetFillColor(CategoryColor("SGN_B0_pstcut"));
	histo_SetUp(h_BKG_B0cut, "BKG_B0_pstcut", "M(B_{0}) [GeV]", Form("%s %4.3f %s", "dN/" ,h_SGN_B0cut->GetXaxis()->GetBinWidth(1), " [GeV]"), true, false);
	h_BKG_B0cut->SetFillColor(CategoryColor("BKG_B0_pstcut"));
	
	histo_SetUp(h_SGN_Rhocut, "SGN_Rho", "M[#rho (770)]# [GeV]", Form("%s %4.3f %s", "dN/" ,h_SGN_Rhocut->GetXaxis()->GetBinWidth(1), " [GeV]"), true, false);
	histo_SetUp(h_BKG_Rhocut, "BKG_Rho", "M[#rho (770)]\\ [GeV]", Form("%s %4.3f %s", "dN/" ,h_SGN_Rhocut->GetXaxis()->GetBinWidth(1), " [GeV]"), true, false);


	// Legend
	auto legend = new TLegend(0.45, 0.65,.85,.85);
	legend->SetTextSize(0.03);
	legend->SetBorderSize(0);

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	h_BKG_B0->Draw();
	h_BKG_B0cut->Draw("SAME");
	h_SGN_B0->Scale(0.173);
	h_SGN_B0->Draw("HIST SAME");
	h_SGN_B0cut->Scale(0.173);
	h_SGN_B0cut->Draw("HIST SAME");

	legend->AddEntry(h_SGN_B0,CategoryLegend("SGN_B0_precut"), "f"); 
	legend->AddEntry(h_SGN_B0cut,CategoryLegend("SGN_B0_pstcut"), "f"); 
	legend->AddEntry(h_BKG_B0,CategoryLegend("BKG_B0_precut"), "f"); 
	legend->AddEntry(h_BKG_B0cut,CategoryLegend("BKG_B0_pstcut"), "f"); 
	legend->Draw();

	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(3);
	gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
	gPad->RedrawAxis();
	gPad->Update();

	c1->SaveAs(OutPlotFile + "MB0_BDTcut.png");
	c1->SaveAs(OutPlotFile + "MB0_BDTcut.pdf");

	std::cout << "EFF --> " << (double)h_SGN_B0cut->Integral()/h_SGN_B0->Integral() << std::endl;
	std::cout << "REJ --> " << 1. - (double)h_BKG_B0cut->Integral()/h_BKG_B0->Integral() << std::endl;


	TCanvas* c2 = new TCanvas("c2","canvas", 1024, 1024);
	h_SGN_Rhocut->Draw();
	h_BKG_Rhocut->Draw("SAME");

	c2->SaveAs(OutPlotFile + "MRho_BDTcut.png");
	

	RootFileMC->Close();
	RootFileDATA->Close();
	RootFileCUT->Close();

	return 0;

}//MVAcutMB0_MRho()


int SGNvsBKGblinded(){

	TFile* RootFileMC   = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/plots/SB_variables.root");
	TH1F* h_SGN_B0 = (TH1F*)RootFileMC->Get("SGN_B0_M");

	TFile* RootFileDATA = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data/merged/BKG_MergData17.root");
	TTree* Tree17 = (TTree*)RootFileDATA->Get("B0sidebands");
	Int_t N = Tree17->GetEntriesFast();
	TH1F* h_BKG_B0 = new TH1F("h_BKG_B0", "", 60, 5., 5.6);
	Tree17->Draw("M_B0>>h_BKG_B0");

	histo_SetUp(h_SGN_B0, "SGN_B0", "M(B^{0}) GeV", Form("%s %4.3f", "dN/" ,h_SGN_B0->GetXaxis()->GetBinWidth(1)), true, false);
	h_SGN_B0->Scale(0.17345103);
	histo_SetUp(h_BKG_B0, "BKG_B0", "M(B^{0}) GeV", Form("%s %4.3f %s", "dN/" ,h_SGN_B0->GetXaxis()->GetBinWidth(1), " [GeV]"), true, false);
	// Legend
	auto legend = new TLegend(0.50, 0.70,.79,.85);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	h_BKG_B0->Draw();
	//h_BKG_B0->SetLineColor(kBlack);
	//h_BKG_B0->SetMarkerStyle(20);
	//h_BKG_B0->Draw("E1 SAME");
	legend->AddEntry(h_BKG_B0,CategoryLegend("BKG_B0")); 
	h_SGN_B0->Draw("SAME HIST");
	legend->AddEntry(h_SGN_B0,CategoryLegend("SGN_B0")); 
	legend->Draw();
	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(3);
	gPad->RedrawAxis();
	gPad->Update();

	CMSxxx(c1);
	c1->SaveAs(OutPlotFile + "SvsBblind.png");
	c1->SaveAs(OutPlotFile + "SvsBblind.pdf");
	c1->SaveAs(OutPlotFile + "SvsBblind.eps");

	RootFileMC->Close();
	RootFileDATA->Close();

	return 0;

}//SGNvsBKGblinded()

int MassBlinded(){
		
		TFile* RootFileDATA = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BKG_data17.root");
		TTree* Tree17 = (TTree*)RootFileDATA->Get("B0sidebands");
		int nbins = 75;
		TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);


		TH1F* h_MuMu_M= new TH1F("h_MuMu_M", "", nbins, 2.9, 3.25);
		Tree17->Draw("M_mumu>>h_MuMu_M");
		histo_SetUp(h_MuMu_M, "BKG_JPsi", "M(\\mu^+\\mu^-)", Form("%s %4.3f", "dN/" ,h_MuMu_M->GetXaxis()->GetBinWidth(1)), true, false);
		h_MuMu_M->Draw();	
		c1->SaveAs(OutPlotFile + "Mmumu.png");

		TH1F* h_Rho_M= new TH1F("h_Rho_M", "", nbins, 0.2, 1.);
		Tree17->Draw("M_Rho>>h_Rho_M");
		histo_SetUp(h_Rho_M, "BKG_Rho", "M(\\pi^+\\pi^-)", Form("%s %4.3f", "dN/" ,h_Rho_M->GetXaxis()->GetBinWidth(1)), true, false);
		h_Rho_M->Draw();
		c1->SaveAs(OutPlotFile + "Mrho.png");

		TH1F* h_X3872_M= new TH1F("h_X3872_M", "", nbins, 3.7, 4.05);
		Tree17->Draw("M_X3872>>h_X3872_M");
		histo_SetUp(h_X3872_M, "BKG_X3872", "M(\\mu^+\\mu^-\\pi^+\\pi^-)", Form("%s %4.3f", "dN/" ,h_X3872_M->GetXaxis()->GetBinWidth(1)), true, false);
		h_X3872_M->Draw();
		c1->SaveAs(OutPlotFile + "MX3872.png");

		TH1F* h_K0s_M= new TH1F("h_K0s_M", "", nbins, 0.4, 0.6);
		Tree17->Draw("M_K0s>>h_K0s_M");
		histo_SetUp(h_K0s_M, "BKG_K0s", "M(K_s^0)\\ ", Form("%s %4.3f", "dN/" ,h_K0s_M->GetXaxis()->GetBinWidth(1)), true, false);
		h_K0s_M->Draw();	
		c1->SaveAs(OutPlotFile + "MK0s.png");


		RootFileDATA->Close();

		return 0;
	}



	int SGNvsBKGinput(){
		
		// FILE	
		TFile* RootFile = openTMVArootFile();
		TDirectoryFile *dir0 = (TDirectoryFile*)RootFile->Get(DirBDT_0);
		TDirectoryFile *dir = (TDirectoryFile*)dir0->Get("InputVariables_Id");
		
		// Var structure
		TString SGNid = "__Signal_Id";
		TString BKGid = "__Background_Id";
		int NVar = 8;
		TString VarName[] = {"D0_Rho", "pT_Rho", "pT_Pi1", "DR_Pi1B0", "LxySign_B0", "SVprob", "pTM_B0", "CosAlpha_B0"};
		std::string VarLabel[] = {"D0(#pi_{1})/#sigma_{xy}", "p_{T}(#rho)/p_{T}(B_{0})", "p_{T}(#pi_{1})/p_{T}(B_{0})", "#Delta R(#pi_{1}, B_{0})", "L_{xy}/#sigma_{xy}", "P(SV)", "p_{T}(B_{0})/M_{B}", "cos(#alpha)"};
		TString Ylabel("");
		double MaxS, MaxB;

		// canva
		TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
		gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		
		// Legend
		auto legend = new TLegend(0.45,0.75,.80,.85);
		legend->SetTextSize(0.035);
		legend->SetBorderSize(0);

		// HISTOS
		TH1F* hS(0);
		TH1F* hB(0);
		for (int i = 0; i < NVar; i++){

			hS = (TH1F*)dir->Get(VarName[i] + SGNid);
			Ylabel = Form("1/N dN/%.3f" ,hS->GetXaxis()->GetBinWidth(1));
			histo_SetUp(hS, "SGN_B0", VarLabel[i], "", false);
			hS->GetYaxis()->SetTitle(Ylabel);
			if(i == 0) legend->AddEntry(hS,CategoryLegend("SGN_B0"), "l");

			hB = (TH1F*)dir->Get(VarName[i] + BKGid);
			histo_SetUp(hB, "BKG_B0", VarLabel[i], "dN/", false);
			if(i == 0) legend->AddEntry(hB,CategoryLegend("BKG_B0"), "l");
			hS->Draw("HIST");
			hB->Draw("HIST same");
			legend->Draw();
		
			MaxS = hS->GetBinContent(hS->GetMaximumBin());
			MaxB = hB->GetBinContent(hB->GetMaximumBin());
			hS->SetMaximum(1.2 * std::max(MaxS, MaxB));
			if(i==7)gPad->SetLogy();
			gStyle->SetOptStat(0);
			gStyle->SetLineWidth(3);
			c1->SaveAs(OutPlotFile + "input_" + VarName[i] + "REDO.png");
			c1->SaveAs(OutPlotFile + "input_" + VarName[i] + "REDO.pdf");

		}

		RootFile->Close();

		return 0;

	}//SGNvsBKGinput()


	int BKGinputMCvsDATA(){
		
		// FILE	
		TFile* RootFileDATA = openTMVArootFile();
		TDirectoryFile *dir0 = (TDirectoryFile*)RootFileDATA->Get(DirBDT_0);
		TDirectoryFile *dir = (TDirectoryFile*)dir0->Get("InputVariables_Id");

		TFile* RootFileMC = new TFile("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/plots/SB_variables.root");
		
		// Var structure
		TString SGNid = "__Signal_Id";
		TString BKGid = "__Background_Id";
		int NVar = 8;
		TString VarName[] = {"D0_Rho", "pT_Rho", "pT_Pi1", "DR_Pi1B0", "CosAlpha_B0", "LxySign_B0", "SVprob", "pTM_B0"};
		TString HistoName[] = {"BKGb_Rho_D0", "BKGb_Rho_pT", "BKGb_Pi1_pT", "BKGb_DR_Pi1B0_Rho", "BKG_B0_cosA", "BKG_B0_LxySign", "BKG_B0_SVp", "BKG_B0_pT"};
		std::string VarLabel[] = {"D0(\\pi_1)/\\sigma", "p_T(\\rho)/p_T(B_0)", "p_T(\\pi_1)/p_T(B_0)", "\\Delta R(\\pi_1, B_0)", "cos(\\alpha)", "L_{xy}/\\sigma", "P(SV)", "p_T(B_0)/M_B \\ "};
		TString Ylabel("");
		double MaxD, MaxMC;

		// canva
		TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
		
		// Legend
		auto legend = new TLegend(0.55,0.80,.89,.89);
		legend->SetTextSize(0.025);
		legend->SetBorderSize(0);

		// HISTOS
		TH1F* hS(0);
		TH1F* hMC(0);
		TH1F* hD(0);
		for (int i = 0; i < NVar; i++){
			hMC = (TH1F*)RootFileMC->Get(HistoName[i]);
			histo_SetUp(hMC, "BKG_mc_B0", VarLabel[i], Form("%s %4.3f", "dN/" ,hMC->GetXaxis()->GetBinWidth(10)), false);
			if(i == 0) legend->AddEntry(hMC,CategoryLegend("BKG_mc_B0"), "l");

			hD = (TH1F*)dir->Get(VarName[i] + BKGid);
			histo_SetUp(hD, "BKG_B0", VarLabel[i], Form("%s %4.3f", "dN/" ,hD->GetXaxis()->GetBinWidth(10)), false);
			if(i == 0) legend->AddEntry(hD,CategoryLegend("BKG_B0"), "l");

			hS = (TH1F*)dir->Get(VarName[i] + SGNid);
			histo_SetUp(hS, "SGN_B0", VarLabel[i], "", false);
			if(i == 0) legend->AddEntry(hS,CategoryLegend("SGN_B0"), "l");

			hS->Draw("HIST");
			hMC->Draw("HIST same");
			hD->Draw("HIST same");
			legend->Draw();
		
			MaxD = hD->GetBinContent(hD->GetMaximumBin());
			MaxMC = hMC->GetBinContent(hMC->GetMaximumBin());
			hD->SetMaximum(1.2 * std::max(MaxD, MaxMC));

			gStyle->SetOptStat(0);
			gPad->RedrawAxis();
			gPad->Update();
			c1->SaveAs(OutPlotFile + "MCvsDATAinput_" + VarName[i] + ".png");

		}

		RootFileMC->Close();
		RootFileDATA->Close();

		return 0;

	}//BKGinputMCvsDATA()




	int ROC_methods(){
		
		TFile* RootFile = openTMVArootFile();
		TDirectoryFile* souceDir = new TDirectoryFile(); 
		// Legend
		auto legend = new TLegend(0.15,0.15,.70,.50);
		legend->SetBorderSize(0);
		//legend->SetTextSize(0.027);
		float ROCintegral;
		
		TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
		c1->SetGrid();
		TCanvas* c2 = new TCanvas("c2","canvas", 1248, 1024);
		c2->SetGrid();

		// HISTO init
		TH1F* hBDT(0);
		TH1F* hDNN(0);

		// BDT
		c1->cd();
		for (int i = 0; i < N_BDT; i++){
			souceDir = accesTMVAmethodDir(RootFile, "BDT", i);
			hBDT = (TH1F*)souceDir->Get("MVA_" + subDirBDT_2[i] + "_rejBvsS");
			histo_SetUp(hBDT, "BDT", "#varepsilon(SGN)", "1- #varepsilon(BKG)", false, false);
			//hBDT->GetYaxis()->SetRangeUser(0., 1.05);
			ROCintegral = hBDT->Integral();
			legend->AddEntry(hBDT, BDTlegend[i] + "|AUC=" + Form("%4.3f",ROCintegral / 100.),"l");
			hBDT->SetLineColor(2+ i);
			hBDT->Draw("same");
		}
		// Draw ...
		legend->Draw();
		
		gStyle->SetOptStat(0);
		gStyle->SetLineWidth(3);
		c1->SaveAs(OutPlotFile+"ROCmethodsBDT.png");
		c1->SaveAs(OutPlotFile+"ROCmethodsBDT.pdf");
		legend->Clear();



		//DNN
		TFile* RootFile2 = openTMVArootFile("DNN");
		c2->cd();
		for (int i = 0; i < N_DNN; i++){
			souceDir = accesTMVAmethodDir(RootFile2, "DNN", i);
			hDNN= (TH1F*)souceDir->Get("MVA_" + subDirDNN_2[i] + "_rejBvsS");
			histo_SetUp(hDNN, "DNN", "#varepsilon(SGN)", "1-#varepsilon(BKG)", false, false);
			legend->AddEntry(hDNN, DNNlegend[i] + "|AUC = " + Form("%4.3f", hDNN->Integral() / 100.), "l");
			hDNN->SetLineColor(51 + 10*i);
			hDNN->Draw("same");
		}
		// Draw ...
		legend->Draw();
		
		gStyle->SetOptStat(0);
		gStyle->SetLineWidth(3);
		c2->SaveAs(OutPlotFile+"ROCmethodsDNN.png");
		c2->SaveAs(OutPlotFile+"ROCmethodsDNN.pdf");

		RootFile2->Close();
		RootFile->Close();
		return 0;	

	}//ROC_methods()




	int Covariance(){
		
		TFile* RootFile = openTMVArootFile();
		TDirectoryFile *dir = (TDirectoryFile*)RootFile->Get(DirBDT_0);
		TH2F* CovMatrix = (TH2F*)dir->Get("CorrelationMatrixS"); 

		// Label 
		int NVar = 8;
		std::string VarLabel[] = {"D0(#pi_{1})/#sigma_{xy}", "p_{T}(#rho)/p_{T}(B_{0})", "p_{T}(#pi_{1})/p_{T}(B_{0})", "#Delta R(#pi_{1}, B_{0})", "cos(#alpha)", "L_{xy}/#sigma_{xy}", "P(SV)", "p_{T}(B_{0})/M_{B}"};
		//std::string VarLabel[] = {"D0(\\pi_1)", "p_T(\\rho)", "p_T(\\pi_1)", "\\Delta R(\\pi_1, B_0)", "cos(\\alpha)", "L_{xy}/\\sigma", "P(SV)", "p_T(B_0) \\ "};
		for (int i = 0; i < NVar; i++){
			CovMatrix->GetXaxis()->SetBinLabel(i+1, (VarLabel[i].c_str()));
			CovMatrix->GetYaxis()->SetBinLabel(i+1, (VarLabel[NVar - 1 - i].c_str()));
		}
		CovMatrix->GetXaxis()->SetLabelSize(0.05);
		CovMatrix->GetYaxis()->SetLabelSize(0.045);	
		CovMatrix->SetTitle("");

		// Draw ...
		TCanvas* c1 = new TCanvas("c1","canvas", 1124, 1024);
		gStyle->SetOptStat(0);
		gStyle->SetLineWidth(3);
		//gStyle->SetPalette(kRainBow);
		float margin = 0.17;
		TPad* pad = new TPad("pad", "", 0.,0.,1., 1.);
		pad->SetRightMargin(margin);	pad->SetLeftMargin(margin);
		pad->Draw();
		pad->cd();
		CovMatrix->SetMarkerColor(kBlack);
		CovMatrix->Draw("COLZ0 TEXT");
		pad->Update();


		c1->SaveAs(OutPlotFile+"CorrMatrixSGN.png");
		c1->SaveAs(OutPlotFile+"CorrMatrixSGN.pdf");
		c1->Delete();
		RootFile->Close();
		return 0;	

	}//Covariance()


	int SGNvsBKGdistribution(){

		TString Tree = subDirBDT_2[0];
		TFile* RootFile = openTMVArootFile();
		// BDT
		TDirectoryFile* souceDir = accesTMVAmethodDir(RootFile, "BDT");
		TH1F* hTestSGN = (TH1F*)souceDir->Get("MVA_" + Tree + "_S");
		TH1F* hTrainSGN = (TH1F*)souceDir->Get("MVA_" + Tree + "_Train_S");
		TH1F* hTestBKG = (TH1F*)souceDir->Get("MVA_" + Tree + "_B");
		TH1F* hTrainBKG = (TH1F*)souceDir->Get("MVA_" + Tree + "_Train_B");

		// SetUP
		//histo_SetUp(TH1* histo, const TString& category, const TString& x_name, const TString& y_name, bool fill = true , bool norm = true){
		histo_SetUp(hTestSGN,  "SGN_B0", "BDT-output (Y_{BDT})", "(1/N)dN/dx", true, false);	
		hTestSGN->SetMaximum(4.5);
		histo_SetUp(hTrainSGN, "SGN_B0", "BDT-output (Y_{BDT})", "(1/N)dN/dx", false, false);	
		hTrainSGN->SetMarkerStyle(20); hTrainSGN->SetMarkerColor(CategoryColor("SGN_B0"));
		histo_SetUp(hTestBKG,  "BKG_B0", "BDT-output (Y_{BDT})", "(1/N)dN/dx", true, false);	
		histo_SetUp(hTrainBKG, "BKG_B0", "BDT-output (Y_{BDT})", "(1/N)dN/dx", false, false);	
		hTrainBKG->SetMarkerStyle(20); hTrainBKG->SetMarkerColor(CategoryColor("BKG_B0"));

		// Legend
		auto legend = new TLegend(0.15 ,0.69 ,0.45,0.89);
		legend->SetBorderSize(0);	
		legend->AddEntry(hTestSGN, CategoryLegend("SGN_B0") + " TEST-set", "f");
		legend->AddEntry(hTrainSGN, CategoryLegend("SGN_B0")+ " TRAIN-set", "lep");
		legend->AddEntry(hTestBKG, CategoryLegend("BKG_B0") + " TEST-set", "f");
		legend->AddEntry(hTrainBKG, CategoryLegend("BKG_B0")+ " TRAIN-set", "lep");

		
		TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
		hTestSGN->Draw("hist");	
		hTrainSGN->Draw("P same");	
		hTestBKG->Draw("hist same");	
		hTrainBKG->Draw("E1 same");	
		legend->Draw();

		gStyle->SetLineWidth(3);
		gStyle->SetOptStat(0);
		gStyle->SetLegendTextSize(0.03);

		c1->SaveAs(OutPlotFile+"SGNvsBKGdistBDT1.png");
		c1->SaveAs(OutPlotFile+"SGNvsBKGdistBDT1.pdf");
		c1->Delete();
		RootFile->Close();
		return 0;	

	}//SGNvsBKGdistribution()

	int ROC_SplitSeed(){

		TString Seeds[] = {"10", "30", "50", "100", "150", "200"};
		TString RootFileName = "TMVAoutputsBDT_S";
		TString FileToOpen("");
		TFile* RootFile = new TFile(); 

		TH1F* ROC_BDT_nT100, *ROC_BDT_nT50;
		TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
		TCanvas* c2 = new TCanvas("c2","canvas", 1024, 1024);
		//LEGEND
		TString leg_entry = "SEED ";
		auto legend = new TLegend(0.12,0.12,0.24,.36);

		for (int f = 0; f < 6; f++){
			
			FileToOpen = FilePath_+RootFileName+Seeds[f]+".root";
			RootFile = new TFile(FileToOpen); 

			if ( !RootFile->IsOpen() ) {
				std::cout << "ERROR IN OPENING FILE "<< RootFile<< std::endl;
				exit(-1);
			}else std::cout << "... reading " << FileToOpen << std::endl;

			TDirectoryFile *dir = (TDirectoryFile*)RootFile->Get(DirBDT_0);
			TDirectoryFile *subDir1 = (TDirectoryFile*)dir->Get(subDirBDT_1);
			TDirectoryFile *subDir2_1 = (TDirectoryFile*)subDir1->Get(subDirBDT_2_1);
			TDirectoryFile *subDir2_2 = (TDirectoryFile*)subDir1->Get(subDirBDT_2_2);

			ROC_BDT_nT100 = (TH1F*)subDir2_1->Get("MVA_BDT_nT100_b08_nC30_rejBvsS");
			ROC_BDT_nT50  = (TH1F*)subDir2_2->Get("MVA_BDT_nT50_b05_nC30_rejBvsS");


			if (f == 0){
				ROC_BDT_nT100->SetTitle("BDT NTrees = 100 beta = 0.8");
				ROC_BDT_nT50->SetTitle("BDT NTrees = 50 beta = 0.5");
				ROC_BDT_nT100->GetXaxis()->SetTitle("signal efficiency");
				ROC_BDT_nT50->GetXaxis()->SetTitle("signal efficiency");
				ROC_BDT_nT100->GetYaxis()->SetTitle("background rejection");
				ROC_BDT_nT50->GetYaxis()->SetTitle("background rejection");

			}

			ROC_BDT_nT100->SetLineWidth(2);
			ROC_BDT_nT50->SetLineWidth(2);
			ROC_BDT_nT100->SetLineColor(f+2);
			ROC_BDT_nT50->SetLineColor(f+2);

			c1->cd();
			ROC_BDT_nT100->Draw("same");			
			c2->cd();
			ROC_BDT_nT50->Draw("same");			

			leg_entry += Seeds[f];
			legend->AddEntry(ROC_BDT_nT50,leg_entry  ,"l");
			leg_entry = "SEED ";
		}	

		gStyle->SetOptStat(0);	
		c1->cd();
		legend->Draw();
		c2->cd();
		legend->Draw();
		c1->SaveAs(FilePath_+ DirBDT_0 + "/plots/SeedROC_nT100.png");
		c2->SaveAs(FilePath_+ DirBDT_0 + "/plots/SeedROC_nT50.png");
		RootFile->Close();

		return 0;
	}



int CountCandidates_B0(TString filePath, TString treeName){

	TFile* inFile = new TFile(filePath); 
	if(!inFile->IsOpen()){
		std::cout << " ERROR : cannot open file " << filePath << std::endl;
		exit(-1);
	}

	TTree* inTree = (TTree*)inFile->Get(treeName);	
	int Nentries = inTree->GetEntries(), Nevents = 0;	
	std::cout<< " # CANDIDATES " << Nentries << std::endl; 
	float prevRun = -1, prevLumiBlock = -1, prevEvent = -1;
	float Run, LumiBlock, Event;
	inTree->SetBranchAddress("run", &Run );
	inTree->SetBranchAddress("LumiBlock", &LumiBlock);
	inTree->SetBranchAddress("event", &Event);

	bool isNewEvent;
	for (int ientry = 0; ientry < Nentries; ientry++){
		inTree->GetEntry(ientry);
		isNewEvent = !((prevRun == Run) && (LumiBlock == prevLumiBlock ) && (Event == prevEvent));
		prevRun = Run; prevLumiBlock = LumiBlock; prevEvent = Event;
		if ( isNewEvent )	Nevents++;
	}
	std::cout << " # EVENTS " << Nevents << " => average multiplicity " << (float)Nentries / Nevents<< std::endl;
	

	return 0;

}
