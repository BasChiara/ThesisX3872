#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TLine.h>
#include <TBox.h>
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include "RooPlot.h"
using namespace RooFit;


TFile* open_file(TString mode = "FINE"){
	
	TString root_file = "../CutOptimization";
	if (mode == "FINE") root_file.Append("Plus");
	root_file.Append(".root");

	TFile* input_file = new TFile(root_file);
	if ( !input_file->IsOpen() ) {
		std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
		exit(-1);
	}

	return input_file;
}

TString pngName(TString name){
	TString png_name = "./" + name +".png"; 
	return png_name;
}
TString pdfName(TString name){
	TString pdf_name = "./" + name +".pdf"; 
	return pdf_name;
}


int PunziPlot(TString mode = "FINE"){
	 
	TFile* inFile = open_file(mode);
	TTree* inTree = (TTree*)inFile->Get("CutOpt");
	double PS, PSerr, Xcut, Mcut, Ss, Bs;
	inTree->SetBranchAddress("cut_X", &Xcut);
	inTree->SetBranchAddress("cut_MRho", &Mcut);
	inTree->SetBranchAddress("SignPunzi", &PS);
	inTree->SetBranchAddress("SignPunzi_error", &PSerr);
	inTree->SetBranchAddress("nS_Sreg", &Ss);
	inTree->SetBranchAddress("nB_Sreg", &Bs);
	float minX, maxX, minM, maxM;
	int Nx, Nm;
	if(mode == "ROUGH"){ 
		std::cout << " ==> Large Scan" << std::endl;
		minX = -0.525, maxX = 0.175;
		Nx = 14;
		minM =  0.475, maxM = 0.725;
		Nm = 5;
	}
	if(mode == "FINE"){
		std::cout << " ==> Fine Scan" << std::endl;
		minX = -0.205, maxX = 0.005;
		Nx = 21;
		minM =  0.5975, maxM = 0.7025;
		Nm = 21;
	}
	// Punzi Significance
	TH2F* h_PunziScan = new TH2F("PunziScan", "", Nx, minX, maxX, Nm, minM, maxM);
	inTree->Draw("cut_X:cut_MRho>>PunziScan");
	h_PunziScan->GetXaxis()->SetTitle("cut on BDT-output (Y_{BDT})");
	h_PunziScan->GetYaxis()->SetTitle("cut on M(#pi^{+} #pi^{-}) [GeV]");
	h_PunziScan->GetZaxis()->SetTitle("Punzi-Sign.");
	h_PunziScan->GetZaxis()->SetTitleOffset(1.8);
	//h_PunziScan->SetTitle("Punzi Significance");
	TBox* FineBOX = new TBox(-0.225, 0.575, 0.025, 0.725);
	FineBOX->SetLineColor(kRed);
	FineBOX->SetLineWidth(4);
	FineBOX->SetFillStyle(0);

	// Signal events 
	TH2F* h_NevSGN= new TH2F("NevSGN", "", Nx, minX, maxX, Nm, minM, maxM);
	inTree->Draw("cut_X:cut_MRho>>NevSGN");
	h_NevSGN->GetXaxis()->SetTitle("cut on BDT-output (Y_{BDT})");
	h_NevSGN->GetYaxis()->SetTitle("cut on M(#pi^{+} #pi^{-}) [GeV]");
	h_NevSGN->GetZaxis()->SetTitle("#SGN");
	h_NevSGN->SetTitle("# signal events");

	// Background events 
	TH2F* h_NevBKG= new TH2F("NevBKG", "", Nx, minX, maxX, Nm, minM, maxM);
	inTree->Draw("cut_X:cut_MRho>>NevBKG");
	h_NevBKG->GetXaxis()->SetTitle("cut on BDT-output (Y_{BDT})");
	h_NevBKG->GetYaxis()->SetTitle("cut on M(#pi^{+} #pi^{-}) [GeV]");
	h_NevBKG->GetZaxis()->SetTitle("#BKG");
	h_NevBKG->SetTitle("# background events");
		
	int Nent = inTree->GetEntries();	
	int binX, binM;
	std::cout << "X-cut \t M-cut \t PunziSign  +/-  error \t Signal" << std::endl;
	for(int i = 0; i < Nent; i++){
		inTree->GetEntry(i);
		binX = h_PunziScan->GetXaxis()->FindBin(Xcut );
		binM = h_PunziScan->GetYaxis()->FindBin(Mcut );
		h_PunziScan->SetBinContent(binX,binM, PS);
		h_NevSGN->SetBinContent(binX,binM, Ss);
		h_NevBKG->SetBinContent(binX,binM, Bs);

		if((PS > 0.00295) && (Ss > 78)){
			std::cout << Xcut << "\t" << Mcut << "\t" << PS  << "  +/-  " << PSerr  << "\t" << Ss << std::endl;
		//std::cout << binX << "\t" << binM << "\t" << h_PunziScan->GetXaxis()->GetBinCenter(binX)<< "\t" << h_PunziScan->GetYaxis()->GetBinCenter(binM) << std::endl;
		}
	}
	
	TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("4.0f");
	float margin = 0.15;
	TPad* pad = new TPad("pad", "", 0.,0.,1., 1.);
	pad->SetRightMargin(margin+.05);	//pad->SetLeftMargin(margin); pad->SetBottomMargin(margin);
	pad->Draw();
	pad->cd();

	gStyle->SetPalette(kBird);
	gStyle->SetLineWidth(3);
	h_PunziScan->Draw("COLZ");
	if(mode == "ROUGH") FineBOX->Draw();
	pad->Update();
	c1->SaveAs(pngName("PunziScan" + mode));
	c1->SaveAs(pdfName("PunziScan" + mode));

	h_NevSGN->SetBarOffset(0.2);
	h_NevSGN->Draw("SAME TEXT");
	h_NevBKG->SetBarOffset(-0.2);
	h_NevBKG->SetMarkerColor(kRed);
	h_NevBKG->Draw("SAME TEXT");
	c1->SaveAs(pngName("PunziScanSB" + mode));
	c1->SaveAs(pdfName("PunziScanSB" + mode));

	gStyle->SetPalette(kLake);
	gStyle->SetLineWidth(3);
	h_NevSGN->Draw("COLZ");
	pad->Update();
	c1->SaveAs(pngName("SGNevents" + mode));
	c1->SaveAs(pdfName("SGNevents" + mode));

	gStyle->SetPalette(kCherry);
	gStyle->SetLineWidth(3);
	h_NevBKG->Draw("COLZ");
	pad->Update();
	c1->SaveAs(pngName("BKGevents" + mode));
	c1->SaveAs(pdfName("BKGevents" + mode));

	inFile->Close();
	return 0;

}




