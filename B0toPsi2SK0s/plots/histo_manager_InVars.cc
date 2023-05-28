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
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <THStack.h>
#include <TLatex.h>


TString SGNregionOptPsi = "(M_X3872 > 3.65455) && (M_X3872 < 3.71843)";
TString SGNregionOptX   = "(M_X3872 > 3.81) && (M_X3872 < 3.93)";
TString SGNregionOptB0  = "(M_B0 > 5.20443) && (M_B0 < 5.35438)";
TString SGNregionOptK0s = "(M_K0s > 0.454494) && (M_K0s < 0.540385)";


TString FileNameBDT_MC     = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/TMVAcompareMC.root";
TString FileNameX3872_MC   = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGN_MC.root";
TString TreeNameX3872_MC   = "mcSIGNAL";
TString FileNamePsi2S_MC   = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/SGN_MC.root";
TString TreeNamePsi2S_MC   = "mcSIGNAL";
TString FileNamePsi2S_DATA = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/data/merged/SGNPsi2S_MergData17.root";
TString TreeNamePsi2S_DAT  = "B0_Psi2Ssignal";

TString OutPlotPath        = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/plots/InVarsCompare/";
TString VarName[] = {"D0_Rho", "pT_Rho", "pT_Pi1", "DR_Pi1B0", "CosAlpha_B0", "LxySign_B0", "SVprob", "pTM_B0"};
std::string VarLabel[] = {"D0(#pi_{1})/#sigma", "p_{T}(#pi#pi)/p_{T}(B_{0})", "p_{T}(#pi_{1})/p_{T}(B_{0})", "#Delta R(#pi_{1}, B_{0})", "cos(#alpha)", "L_{xy}/#sigma", "P(SV)", "p_{T}(B_{0})/M_{B_{0}}"};
int    VarNbins[] = { 25,  35,  25,  25, 30,  50, 25, 20};
double VarXlow[]  = { 0.,  0.,  0.,  0.,.97,  0., 0., 0.};
double VarXhigh[] = {50., .35, .25, .75, 1., 450, 1., 20};


TFile* openFile(TString method = "BDT"){

	TString FileToOpen(""); 
	if (method == "X_MC")     FileToOpen = FileNameX3872_MC;
	if (method == "PSI_MC")   FileToOpen = FileNamePsi2S_MC;
	if (method == "PSI_DATA") FileToOpen = FileNamePsi2S_DATA;
	if (method == "BDT")      FileToOpen = FileNameBDT_MC;

	TFile* RootFile = new TFile(FileToOpen); 
	if ( !RootFile->IsOpen() ) {
		std::cout << "ERROR IN OPENING FILE "<< RootFile<< std::endl;
		exit(-1);
	}

	return RootFile;

}//openFile()

Color_t CategoryColor(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["X_MC"]           = kAzure + 1; 
  Color["PSI_MC"]         = kPink - 4; 
  Color["PSI_DAT"]        = kBlack; 

  return Color[category];
}//CategoryLegend()


TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["X_MC"] = "X(3872) SIGNAL (MC)";
  Leg_entry["PSI_MC"] = "#Psi(2S) SIGNAL (MC)";
  Leg_entry["PSI_DAT"] = TString("#Psi(2S) DATA2017");

  return Leg_entry[category];
}//CategoryLegend()


void histo_SetUp(TH1* histo, const TString& category, const TString& x_name, const TString& y_name, bool fill = true , bool norm = true){
  //AXIS LABEL
  histo->SetTitle("");
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.03);
  histo->GetYaxis()->SetTitle(y_name);
  histo->GetYaxis()->SetLabelSize(0.035);
  histo->GetYaxis()->SetTitleSize(0.04);


  //WIDTH & COLOR 
  gStyle->SetLineWidth(3);
  histo->SetLineWidth(4);
  histo->SetLineColor(CategoryColor(category));
  if (fill)histo->SetFillColorAlpha(CategoryColor(category), 0.1);
  //NORMALIZATION
  if(norm) histo->Scale(1./histo->Integral());
}

int InVarsMCvsDAT_Psi(){

	
	// FILE	
	TFile* RootFileMC_Psi = openFile("PSI_MC");
	TTree* TreeMC_Psi = (TTree*)RootFileMC_Psi->Get(TreeNamePsi2S_MC);
	std::cout << TreeMC_Psi->GetEntries() << std::endl;
	
	TFile* RootFileDAT_Psi = openFile("PSI_DATA");
	TTree* TreeDAT_Psi = (TTree*)RootFileDAT_Psi->Get(TreeNamePsi2S_DAT);
	std::cout << TreeDAT_Psi->GetEntries() << std::endl;
	

	TH1F* hMC = new TH1F("hMC", "", 100, 0, 500 );  
	TH1F* hDAT = new TH1F("hDAT", "", 100, 0, 500); 

	// canva
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

	// Legend
	auto legend = new TLegend(0.55,0.80,.89,.89);
	legend->SetTextSize(0.025);
	legend->SetBorderSize(0);
	const int NVars = 8;
	TString DrawMC = ">>hMC"; 
	TString DrawDAT= ">>hDAT"; 
	for(int i = 0; i < NVars; i++){

		hMC->SetBins(VarNbins[i], VarXlow[i], VarXhigh[i]);
		hDAT->SetBins(VarNbins[i], VarXlow[i], VarXhigh[i]);
		TreeMC_Psi->Draw (VarName[i] + DrawMC,  SGNregionOptPsi + " && " + SGNregionOptB0 + " && " + SGNregionOptK0s);
		TreeDAT_Psi->Draw(VarName[i] + DrawDAT, SGNregionOptPsi + " && " + SGNregionOptB0 + " && " + SGNregionOptK0s);

		histo_SetUp(hMC, "PSI_MC", VarLabel[i], Form("%s %4.3f", "dN/" ,hMC->GetXaxis()->GetBinWidth(10)));
		histo_SetUp(hDAT, "PSI_DAT", VarLabel[i], Form("%s %4.3f", "dN/" ,hDAT->GetXaxis()->GetBinWidth(10)), false);
		hDAT->SetLineWidth(1);
		hMC->Draw("HIST");
		hDAT->Draw("E1 SAME");
		hMC->SetMaximum(1.2 * std::max(hDAT->GetBinContent(hDAT->GetMaximumBin()), hMC->GetBinContent(hMC->GetMaximumBin()) ));
		hDAT->SetMarkerStyle(kFullCircle);
		if(i == 0) legend->AddEntry(hDAT, CategoryLegend("PSI_DAT"));
		if(i == 0) legend->AddEntry(hMC, CategoryLegend("PSI_MC"));
		legend->Draw();	
		gStyle->SetOptStat(0);
		gPad->RedrawAxis();
		gPad->Update();
		c1->SaveAs(OutPlotPath + "MCvsDAT_Psi2S" + VarName[i] + ".png");
		c1->SaveAs(OutPlotPath + "MCvsDAT_Psi2S" + VarName[i] + ".pdf");


	}
	std::cout << "MC " << TreeMC_Psi->GetEntries() << " ---> " << hMC->GetEntries() << std::endl; 
	std::cout << "DAT " << TreeDAT_Psi->GetEntries() << " ---> " << hDAT->GetEntries() << std::endl; 

	return 0;

}//InVarsMCvsDAT_Psi()

int OutVarMCvsMC_PsiX(){

	// FILE	
	TFile* RootFile = openFile("BDT");
	TH1F* hBDTx_X = (TH1F*)RootFile->Get("BDT_X");
	TH1F* hMRho_X = (TH1F*)RootFile->Get("SGN_MRho_X");
	TH1F* hBDTx_Psi = (TH1F*)RootFile->Get("BDT_Psi");
	TH1F* hMRho_Psi = (TH1F*)RootFile->Get("SGN_MRho_Psi");

	histo_SetUp(hBDTx_X, "X_MC", "BDT output (Y_{BDT})", Form("%s %4.3f", "dN/" ,hBDTx_X->GetXaxis()->GetBinWidth(10)));
	histo_SetUp(hMRho_X, "X_MC", "M(#pi^{+}#pi^{-})", Form("dN/%.0f MeV", hMRho_X->GetXaxis()->GetBinWidth(10) * 1000.));
	histo_SetUp(hBDTx_Psi, "PSI_MC", "BDT output (Y_{BDT})", Form("%s %4.3f", "dN/" ,hBDTx_Psi->GetXaxis()->GetBinWidth(10)));
	histo_SetUp(hMRho_Psi, "PSI_MC", "M(#pi^{+}#pi^{-})", Form("%s %4.3f", "dN/" ,hMRho_Psi->GetXaxis()->GetBinWidth(10)));
	hBDTx_X->SetMaximum(1.2 * std::max(hBDTx_X->GetBinContent(hBDTx_X->GetMaximumBin()), hBDTx_Psi->GetBinContent(hBDTx_X->GetMaximumBin()) ));
	hMRho_X->SetMaximum(1.2 * std::max(hMRho_X->GetBinContent(hMRho_X->GetMaximumBin()), hMRho_Psi->GetBinContent(hMRho_X->GetMaximumBin()) ));

	// canva
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

	auto legend = new TLegend(0.15,0.73,.45,.82);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);
	legend->AddEntry(hBDTx_X, CategoryLegend("X_MC"));
	legend->AddEntry(hBDTx_Psi, CategoryLegend("PSI_MC"));

	hBDTx_X->Draw("HIST");
	hBDTx_Psi->Draw("HIST SAME");

	gStyle->SetOptStat(0);
	gPad->RedrawAxis();
	gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
	gPad->Update();
	legend->Draw();
	c1->SaveAs(OutPlotPath +  "BDTx_SGN.png");
	c1->SaveAs(OutPlotPath +  "BDTx_SGN.pdf");

	hMRho_X->Draw("HIST");
	hMRho_Psi->Draw("HIST SAME");

	gStyle->SetOptStat(0);
	gPad->RedrawAxis();
	gPad->Update();
	legend->Draw();
	c1->SaveAs(OutPlotPath +  "MRho_SGN.png");
	c1->SaveAs(OutPlotPath +  "MRho_SGN.pdf");

	return 0;

}//OutVarMCvsMC_PsiX()

int InVarsMCvsMC_PsiX(){


	// FILE	
	TFile* RootFileMC_Psi = openFile("PSI_MC");
	TTree* TreeMC_Psi = (TTree*)RootFileMC_Psi->Get(TreeNamePsi2S_MC);
	
	TFile* RootFileMC_X = openFile("X_MC");
	TTree* TreeMC_X = (TTree*)RootFileMC_X->Get(TreeNameX3872_MC);

	// Var structure
	TString SGNid = "__Signal_Id";
	const int NVar = 8;
	TH1F* hMC_X = new TH1F("hMC_X", "", 100, 0, 500 );  
	TH1F* hMC_Psi = new TH1F("hMC_Psi", "", 100, 0, 500); 
	TString DrawMCX = ">>hMC_X"; 
	TString DrawMCP= ">>hMC_Psi"; 

	// canva
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

	// Legend
	auto legend = new TLegend(0.45,0.75,.84,.89);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);

	// HISTOS

	for(int i =0;i < NVar; i++){

		hMC_X->SetBins(VarNbins[i], VarXlow[i], VarXhigh[i]);
		hMC_Psi->SetBins(VarNbins[i], VarXlow[i], VarXhigh[i]);
		TreeMC_X->Draw (VarName[i] + DrawMCX,  SGNregionOptX+ " && " + SGNregionOptB0 + " && " + SGNregionOptK0s);
		TreeMC_Psi->Draw(VarName[i] + DrawMCP, SGNregionOptPsi + " && " + SGNregionOptB0 + " && " + SGNregionOptK0s);

		histo_SetUp(hMC_X, "X_MC", VarLabel[i], Form("%s %4.3f", "dN/" ,hMC_X->GetXaxis()->GetBinWidth(10)));
		histo_SetUp(hMC_Psi, "PSI_MC", VarLabel[i], Form("%s %4.3f", "dN/" ,hMC_Psi->GetXaxis()->GetBinWidth(10)));
		hMC_X->Draw("HIST");
		hMC_Psi->Draw("HIST SAME");
		hMC_X->SetMaximum(1.3 * std::max(hMC_Psi->GetBinContent(hMC_Psi->GetMaximumBin()), hMC_X->GetBinContent(hMC_X->GetMaximumBin()) ));
		if(i == 0) legend->AddEntry(hMC_Psi, CategoryLegend("PSI_MC"));
		if(i == 0) legend->AddEntry(hMC_X, CategoryLegend("X_MC"));
		legend->Draw();	
		gStyle->SetOptStat(0);
		gStyle->SetLineWidth(3);
		if( VarName[i] == "CosAlpha_B0")gPad->SetLogy();
		else gPad->SetLogy(0);
		gPad->RedrawAxis();
		gPad->Update();
		gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
		c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[i] + ".png");
		c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[i] + ".pdf");

	}

	/*TH1F hMC_D0_Rho_X("hMC_D0_Rho_X", "", Nbins, low, high); // D0_Rho
	TH1F hMC_D0_Rho_Psi("hMC_D0_Rho_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("D0_Rho>>hMC_D0_Rho_X", SGNregionOptB0);	
	TreeMC_Psi->Draw("D0_Rho>>hMC_D0_Rho_Psi", SGNregionOptB0);	
	histo_SetUp(&hMC_D0_Rho_X, "X_MC", VarLabel[0], Form("%s %4.3f", "dN/" ,hMC_D0_Rho_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_D0_Rho_Psi, "PSI_MC", VarLabel[0], Form("%s %4.3f", "dN/" ,hMC_D0_Rho_Psi.GetXaxis()->GetBinWidth(10)), false);
	legend->AddEntry(&hMC_D0_Rho_X, CategoryLegend("X_MC"));
	legend->AddEntry(&hMC_D0_Rho_Psi, CategoryLegend("PSI_MC"));
	
	hMC_D0_Rho_X.Draw("HIST");
	hMC_D0_Rho_Psi.Draw("HIST SAME");
	MaxX = hMC_D0_Rho_X.GetBinContent(hMC_D0_Rho_X.GetMaximumBin());
	MaxPsi = hMC_D0_Rho_Psi.GetBinContent(hMC_D0_Rho_Psi.GetMaximumBin());
	hMC_D0_Rho_X.SetMaximum(1.2 * std::max(MaxX, MaxPsi));	
	legend->Draw();

	gStyle->SetOptStat(0);
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[0] + ".png");

	Nbins = 35 ,low = 0. , high = .35;
	TH1F hMC_pT_Rho_X("hMC_pT_Rho_X", "", Nbins, low, high); // pT_Rho
	TH1F hMC_pT_Rho_Psi("hMC_pT_Rho_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("pT_Rho>>hMC_pT_Rho_X",SGNregionOptB0);	
	TreeMC_Psi->Draw("pT_Rho>>hMC_pT_Rho_Psi",SGNregionOptB0);	
	histo_SetUp(&hMC_pT_Rho_X, "X_MC", VarLabel[1], Form("%s %4.3f", "dN/" ,hMC_pT_Rho_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_pT_Rho_Psi, "PSI_MC", VarLabel[1], Form("%s %4.3f", "dN/" ,hMC_pT_Rho_Psi.GetXaxis()->GetBinWidth(10)), false);

	hMC_pT_Rho_X.Draw("HIST");
	hMC_pT_Rho_Psi.Draw("HIST SAME");
	MaxX = hMC_pT_Rho_X.GetBinContent(hMC_pT_Rho_X.GetMaximumBin());
	MaxPsi = hMC_pT_Rho_Psi.GetBinContent(hMC_pT_Rho_Psi.GetMaximumBin());
	hMC_pT_Rho_X.SetMaximum(1.2 * std::max(MaxX, MaxPsi));	
	legend->Draw();
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[1] + ".png");

	Nbins = 25 ,low = 0. , high = .25;
	TH1F hMC_pT_Pi1_X("hMC_pT_Pi1_X", "", Nbins, low, high); // pT_Pi1
	TH1F hMC_pT_Pi1_Psi("hMC_pT_Pi1_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("pT_Pi1>>hMC_pT_Pi1_X",SGNregionOptB0);	
	TreeMC_Psi->Draw("pT_Pi1>>hMC_pT_Pi1_Psi",SGNregionOptB0);	
	histo_SetUp(&hMC_pT_Pi1_X, "X_MC", VarLabel[2], Form("%s %4.3f", "dN/" ,hMC_pT_Pi1_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_pT_Pi1_Psi, "PSI_MC", VarLabel[2], Form("%s %4.3f", "dN/" ,hMC_pT_Pi1_Psi.GetXaxis()->GetBinWidth(10)), false);

	hMC_pT_Pi1_X.Draw("HIST");
	hMC_pT_Pi1_Psi.Draw("HIST SAME");
	MaxX = hMC_pT_Pi1_X.GetBinContent(hMC_pT_Pi1_X.GetMaximumBin());
	MaxPsi = hMC_pT_Pi1_Psi.GetBinContent(hMC_pT_Pi1_Psi.GetMaximumBin());
	hMC_pT_Pi1_X.SetMaximum(1.2 * std::max(MaxX, MaxPsi));	
	legend->Draw();
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[2] + ".png");

	Nbins = 25 ,low = 0. , high = .75;
	TH1F hMC_DR_Pi1B0_X("hMC_DR_Pi1B0_X", "", Nbins, low, high); // DR_Pi1B0
	TH1F hMC_DR_Pi1B0_Psi("hMC_DR_Pi1B0_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("DR_Pi1B0>>hMC_DR_Pi1B0_X",SGNregionOptB0);	
	TreeMC_Psi->Draw("DR_Pi1B0>>hMC_DR_Pi1B0_Psi",SGNregionOptB0);	
	histo_SetUp(&hMC_DR_Pi1B0_X, "X_MC", VarLabel[3], Form("%s %4.3f", "dN/" ,hMC_DR_Pi1B0_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_DR_Pi1B0_Psi, "PSI_MC", VarLabel[3], Form("%s %4.3f", "dN/" ,hMC_DR_Pi1B0_Psi.GetXaxis()->GetBinWidth(10)), false);

	hMC_DR_Pi1B0_X.Draw("HIST");
	hMC_DR_Pi1B0_Psi.Draw("HIST SAME");
	MaxX = hMC_DR_Pi1B0_X.GetBinContent(hMC_DR_Pi1B0_X.GetMaximumBin());
	MaxPsi = hMC_DR_Pi1B0_Psi.GetBinContent(hMC_DR_Pi1B0_Psi.GetMaximumBin());
	hMC_DR_Pi1B0_X.SetMaximum(1.2 * std::max(MaxX, MaxPsi));	
	legend->Draw();
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[3] + ".png");

	Nbins = 30,low = 0.97 , high = 1.;
	TH1F hMC_CosAlpha_B0_X("hMC_CosAlpha_B0_X", "", Nbins, low, high); // CosAlpha_B0
	TH1F hMC_CosAlpha_B0_Psi("hMC_CosAlpha_B0_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("CosAlpha_B0>>hMC_CosAlpha_B0_X", SGNregionOptB0);	
	TreeMC_Psi->Draw("CosAlpha_B0>>hMC_CosAlpha_B0_Psi", SGNregionOptB0);	
	histo_SetUp(&hMC_CosAlpha_B0_X, "X_MC", VarLabel[4], Form("%s %4.3f", "dN/" ,hMC_CosAlpha_B0_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_CosAlpha_B0_Psi, "PSI_MC", VarLabel[4], Form("%s %4.3f", "dN/" ,hMC_CosAlpha_B0_Psi.GetXaxis()->GetBinWidth(10)), false);

	hMC_CosAlpha_B0_X.Draw("HIST");
	hMC_CosAlpha_B0_Psi.Draw("HIST SAME");
	MaxX = hMC_CosAlpha_B0_X.GetBinContent(hMC_CosAlpha_B0_X.GetMaximumBin());
	MaxPsi = hMC_CosAlpha_B0_Psi.GetBinContent(hMC_CosAlpha_B0_Psi.GetMaximumBin());
	hMC_CosAlpha_B0_X.SetMaximum(3.);	
	legend->Draw();
	gPad->SetLogy();
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[4] + ".png");

	Nbins = 50 ,low = 0. , high = 450;
	TH1F hMC_LxySign_B0_X("hMC_LxySign_B0_X", "", Nbins, low, high); // LxySign_B0
	TH1F hMC_LxySign_B0_Psi("hMC_LxySign_B0_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("LxySign_B0>>hMC_LxySign_B0_X", SGNregionOptB0);	
	TreeMC_Psi->Draw("LxySign_B0>>hMC_LxySign_B0_Psi", SGNregionOptB0);	
	histo_SetUp(&hMC_LxySign_B0_X, "X_MC", VarLabel[5], Form("%s %4.3f", "dN/" ,hMC_LxySign_B0_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_LxySign_B0_Psi, "PSI_MC", VarLabel[5], Form("%s %4.3f", "dN/" ,hMC_LxySign_B0_Psi.GetXaxis()->GetBinWidth(10)), false);

	hMC_LxySign_B0_X.Draw("HIST");
	hMC_LxySign_B0_Psi.Draw("HIST SAME");
	MaxX = hMC_LxySign_B0_X.GetBinContent(hMC_LxySign_B0_X.GetMaximumBin());
	MaxPsi = hMC_LxySign_B0_Psi.GetBinContent(hMC_LxySign_B0_Psi.GetMaximumBin());
	hMC_LxySign_B0_X.SetMaximum(1.2 * std::max(MaxX, MaxPsi));	
	legend->Draw();
	gPad->SetLogy(0);
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[5] + ".png");

	Nbins = 40 ,low = 0. , high = 1.;
	TH1F hMC_SVprob_X("hMC_SVprob_X", "", Nbins, low, high); // SVprob
	TH1F hMC_SVprob_Psi("hMC_SVprob_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("SVprob>>hMC_SVprob_X", SGNregionOptB0);	
	TreeMC_Psi->Draw("SVprob>>hMC_SVprob_Psi", SGNregionOptB0);	
	histo_SetUp(&hMC_SVprob_X, "X_MC", VarLabel[6], Form("%s %4.3f", "dN/" ,hMC_SVprob_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_SVprob_Psi, "PSI_MC", VarLabel[6], Form("%s %4.3f", "dN/" ,hMC_SVprob_Psi.GetXaxis()->GetBinWidth(10)), false);

	hMC_SVprob_X.Draw("HIST");
	hMC_SVprob_Psi.Draw("HIST SAME");
	MaxX = hMC_SVprob_X.GetBinContent(hMC_SVprob_X.GetMaximumBin());
	MaxPsi = hMC_SVprob_Psi.GetBinContent(hMC_SVprob_Psi.GetMaximumBin());
	hMC_SVprob_X.SetMaximum(1.2 * std::max(MaxX, MaxPsi));	
	legend->Draw();
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[6] + ".png");

	Nbins = 20,low = 0. , high = 20;
	TH1F hMC_pTM_B0_X("hMC_pTM_B0_X", "", Nbins, low, high); // pTM_B0
	TH1F hMC_pTM_B0_Psi("hMC_pTM_B0_Psi", "", Nbins, low, high);
	TreeMC_X->Draw("pTM_B0>>hMC_pTM_B0_X",SGNregionOptB0);	
	TreeMC_Psi->Draw("pTM_B0>>hMC_pTM_B0_Psi",SGNregionOptB0);	
	histo_SetUp(&hMC_pTM_B0_X, "X_MC", VarLabel[7], Form("%s %4.3f", "dN/" ,hMC_pTM_B0_X.GetXaxis()->GetBinWidth(10)), false);
	histo_SetUp(&hMC_pTM_B0_Psi, "PSI_MC", VarLabel[7], Form("%s %4.3f", "dN/" ,hMC_pTM_B0_Psi.GetXaxis()->GetBinWidth(10)), false);

	hMC_pTM_B0_X.Draw("HIST");
	hMC_pTM_B0_Psi.Draw("HIST SAME");
	MaxX = hMC_pTM_B0_X.GetBinContent(hMC_pTM_B0_X.GetMaximumBin());
	MaxPsi = hMC_pTM_B0_Psi.GetBinContent(hMC_pTM_B0_Psi.GetMaximumBin());
	hMC_pTM_B0_X.SetMaximum(1.2 * std::max(MaxX, MaxPsi));	
	legend->Draw();
	gPad->RedrawAxis();
	gPad->Update();
	c1->SaveAs(OutPlotPath + "MCvsMC_" + VarName[7] + ".png");

*/

	/*for (int i = 0; i < NVar; i++){
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

	}*/

	RootFileMC_Psi->Close();
	RootFileMC_X->Close();




	return 0;

}
