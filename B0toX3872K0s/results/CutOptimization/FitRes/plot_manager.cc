#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TLine.h>
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "RooPlot.h"
using namespace RooFit;


TFile* open_file(double& Xcut, double& MRcut){
	
	TString root_file = Form("FitBKG_X");
	if (Xcut < -0.001) root_file.Append("m");
	root_file.Append(Form("%.0f_MR%.0f.root", fabs(Xcut)*100., MRcut*1000.));
	TFile* input_file = new TFile(root_file);

	if ( !input_file->IsOpen() ) {
		std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
		exit(-1);
	}

	return input_file;
}

TString pngName(double& Xcut, double& MRcut, std::string category){
	TString png_name = Form("./fitPlots/%sX", category.c_str());
	if (Xcut < -0.001) png_name.Append("m");
	png_name.Append(Form("%.0f_MR%.0f.png", fabs(Xcut)*100., MRcut*1000.)); 
	return png_name;
}


int FitPlot(double Xcut, double MRcut){


	TFile* inFile = open_file(Xcut, MRcut);
	RooPlot* hFit = (RooPlot*)inFile->Get("FitPlot");
	
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	hFit->Draw();
	c1->SaveAs(pngName(Xcut, MRcut, "Fit"));

	inFile->Close();
	return 0;

}

int FitPlotBEST(double Xcut, double MRcut){


	TFile* inFile = open_file(Xcut, MRcut);
	RooPlot* hFit = (RooPlot*)inFile->Get("FitPlot");
	TGraph* h_Dati = (TGraph*)hFit->getObject(0);	
	//std::cout << h_Dati->GetXaxis()->GetNbins() << std::endl;

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	h_Dati->Draw();
	c1->SaveAs(pngName(Xcut, MRcut, "FIT"));

	inFile->Close();
	return 0;

}

int CorrMtx(double Xcut, double MRcut){

	TFile* inFile = open_file(Xcut, MRcut);
	TH2* hCorrMtx = (TH2*)inFile->Get("correlation_matrix");
	hCorrMtx->SetTitle(Form("X-cut = %.2f Mrho-cut = %.3f", Xcut, MRcut));

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	float margin = 0.17;
	TPad* pad = new TPad("pad", "", 0.,0.,1., 1.);
	pad->SetRightMargin(margin);	pad->SetLeftMargin(margin);
	pad->Draw();
	pad->cd();
	hCorrMtx->Draw("COLZ TEXT");
	pad->Update();
	gStyle->SetOptStat(0);
	c1->SaveAs(pngName(Xcut, MRcut, "CorrM"));

	inFile->Close();


	return 0;
}


