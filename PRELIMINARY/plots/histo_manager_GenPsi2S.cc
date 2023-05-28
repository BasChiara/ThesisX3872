#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


const TString eosPath = "/eos/user/c/cbasile/www/B0toX3872K0s/GEN_LEVEL/GenUL17Psi2S/";

TFile* open_file(){
    TString root_file = "../outRoot/CheckGenLev_UL17_Psi2S";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}


Color_t PtlColorMap(const TString& particle){

  std::map <TString , Color_t> PtlColor{};
  PtlColor["mu"] = kPink + 5;
  PtlColor["muL"] = kAzure + 5;
  PtlColor["muSL"] = kPink + 5;
  PtlColor["pi"] = kBlue - 4;
  PtlColor["PiPi"] = kOrange + 6;
  PtlColor["JPsi"] = kRed;
  PtlColor["Rho"] = kOrange + 8;
  PtlColor["Ks"] = kGreen;
  PtlColor["X"] = kBlack + 3;
  PtlColor["Psi2S"] = kBlack + 3;
  PtlColor["B0"] = kViolet + 8;

  return PtlColor[particle];
}


TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["mu"] = "#mu";
  Leg_entry["muL"] = "Leading #mu";
  Leg_entry["muSL"] = "Sub-leading #mu";
  Leg_entry["pi"] = "pions";
  Leg_entry["PiPi"] = "#pi^{+} #pi^{-}";
  Leg_entry["JPsi"] = "J/#psi";
  Leg_entry["Rho"] = "#rho(770)";
  Leg_entry["Ks"] = "K_{s}^{0}";
  Leg_entry["X"] = "X(3872)";
  Leg_entry["Psi2S"] = "#psi (2S)";
  Leg_entry["B0"] = "B_{0} ";

  return Leg_entry[category];
}

void histoSetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true , bool norm = true){

  //AXIS LABEL 
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitle(Form("Events/ %.1f", histo->GetXaxis()->GetBinWidth(1)));
  histo->GetYaxis()->SetTitleOffset(1.6);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.04);

  //WIDTH & COLOR
  histo->SetLineWidth(3);
  histo->SetLineColor(PtlColorMap(category));
  if (fill) histo->SetFillColorAlpha(PtlColorMap(category), 0.3);

  //NORMALIZATION
  if(norm )histo->Scale(1./histo->Integral());



}

TString pngName(TString histo_name){

	TString pngName = eosPath + histo_name + ".png";
	return pngName;
}
TString pdfName(TString histo_name){

	TString pdfName = eosPath + histo_name + ".pdf";
	return pdfName;
}

void CMSxxx(TCanvas* c){
	c->cd();
	TLatex RunDetails; RunDetails.SetNDC(); 
	RunDetails.SetTextFont(61);
	RunDetails.SetTextAlign(10);
	RunDetails.SetTextSize(0.03);
	RunDetails.DrawLatex(.10, .91, "CMS");
	RunDetails.SetTextFont(52);
	RunDetails.DrawLatex(.17, .91, "Simulation");
	RunDetails.SetTextFont(42);
	RunDetails.SetTextSize(0.025);
	RunDetails.DrawLatex(.70, .91, "41 fb^{-1} (13 TeV)");

}

int draw_pT_histo(const TString histo_name, const TString particle){
    
    TFile* input_file = open_file();

    TH1F* h = (TH1F*)input_file->Get(histo_name);
    
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
	  histoSetUp(h, particle, "p_{T} [GeV]");
	  h->GetYaxis()->SetTitle(Form("1/N dN/%.2f GeV", h->GetXaxis()->GetBinWidth(1)));

    auto legend = new TLegend(0.60,0.75,.80,.80);
	  legend->SetBorderSize(0);
	  legend->SetTextSize(0.035);
    legend->AddEntry(h,CategoryLegend(particle),"f");

    //STATISTICS
    gStyle->SetOptStat(0);

	 //TEXT

    TString png_name = pngName(histo_name); 
    TString pdf_name = pdfName(histo_name); 
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
	  c1->DrawFrame(0,0,1,1);
    h->Draw("HIST");
	  legend->Draw();
	  gPad->SetLeftMargin(0.13);
	  gPad->SetBottomMargin(0.13);
    c1->Update(); 
    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

    input_file->Close();
    return 0;

}//draw_pT()


int draw_two_histograms(const TString histo1,const TString& category1, const TString histo2, const TString& category2, const TString& x_name, bool fill = true){

    TFile* input_file = open_file();
    TH1F* h1 = (TH1F*)input_file->Get(histo1);
    TH1F* h2 = (TH1F*)input_file->Get(histo2);
    if ( !h1 ){
      std::cout<< "null pointer for histogram named " << histo1 << std::endl;
      exit(-1);
    }    
    if ( !h2 ){
      std::cout<< "null pointer for histogram named " << histo2 << std::endl;
      exit(-1);
    } 
    
    histoSetUp(h1, category1, x_name);
    histoSetUp(h2, category2, x_name);
	 h1->GetYaxis()->SetTitle(Form("1/N dN/%.2f GeV", h1->GetXaxis()->GetBinWidth(1)));

    //STATISTICS
    gStyle->SetOptStat(0);
    
    //SETMAXIMUM                                                                                                                                                                  
    double M1 = h1->GetBinContent(h1->GetMaximumBin());
    double M2 = h2->GetBinContent(h2->GetMaximumBin());
    if (M1 > M2){ h1->SetMaximum(1.2*M1);
    }else {h1->SetMaximum(1.2*M2);}
    
    //LEGEND
    auto legend = new TLegend(0.50,0.70,.80,.80);
	 legend->SetBorderSize(0);
	 legend->SetTextSize(0.035);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    legend->AddEntry(h2, CategoryLegend(category2) ,"f");
    
    TString png_name = pngName(histo1);
    TString pdf_name = pdfName(histo1);
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.13);
    gPad->RedrawAxis();
    legend->Draw();
	 //c1->SetLogy();
    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

    input_file->Close();

    return 0;
}

int draw_mass_histo(const TString histo_name, const TString particle){
    
    TFile* input_file = open_file();

    TH1F* h = (TH1F*)input_file->Get(histo_name);
    
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
	 histoSetUp(h, particle, "M [GeV]");
    
    //LEGEND
    auto legend = new TLegend(0.75,0.82,.89,.89);
	 legend->SetBorderSize(0);
	 legend->SetTextSize(0.04);
    legend->AddEntry(h,CategoryLegend(particle),"l");

    //STATISTICS
    gStyle->SetOptStat(0);
    TString png_name = pngName(histo_name);
    TString pdf_name = pdfName(histo_name);
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
	
	 h->Draw("HIST");
    legend->Draw();
    
    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

    input_file->Close();
    return 0;

}//draw_pT()


int draw_Eta_histo(const TString histo_name, const TString particle){

  TFile* input_file = open_file();

  TH1F* h = (TH1F*)input_file->Get(histo_name);

  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }
  histoSetUp(h, particle, "\\eta");

  //LEGEND                                                                                                               
  auto legend = new TLegend(0.75,0.82,.89,.89);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(h,CategoryLegend(particle),"f");

  //STATISTICS                                                                                                           
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name); 
  TString pdf_name = pdfName(histo_name); 
  TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.13);
  h->Draw("HIST");
  legend->Draw();

  c1->SaveAs(png_name);
  c1->SaveAs(pdf_name);

  input_file->Close();

  return 0;

}//drae_Eta()

int draw_Mul(const TString histo_name, const TString particle){

  TString root_file = "analysis_Multiplicity_UL17.root";
  TFile* input_file = new TFile(root_file);

  TH1I* h = (TH1I*)input_file->Get(histo_name);

  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }

  histoSetUp(h, particle, "# candidates", true, false);
  h->GetYaxis()->SetTitle("counts");
  if(particle == "JPsi"){
	  int Nbins = h->GetNbinsX();
	  for (int i = 0; i < Nbins; i++) h->GetXaxis()->SetBinLabel(i+1,std::to_string(i).c_str());
	  h->GetXaxis()->SetLabelSize(0.05);
  }
if(particle == "PiPi") h->GetXaxis()->SetNdivisions(5,kTRUE); 
	
  //LEGEND                                                                                                                                                                    
  auto legend = new TLegend(0.65,0.75,.80,.80);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.035);
  legend->AddEntry(h,CategoryLegend(particle),"f");

  //STATISTICS                                                                                                                                                                
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name); 
  TString pdf_name = pdfName(histo_name); 
  TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.13);
  h->Draw("HIST");
  legend->Draw();

  c1->SaveAs(png_name);
  c1->SaveAs(pdf_name);

  input_file->Close();
  return 0;

}//draw_mul()   
