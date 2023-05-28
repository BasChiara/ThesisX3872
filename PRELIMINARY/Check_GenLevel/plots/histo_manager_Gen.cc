#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


TFile* open_file(){
    //TString root_file = "CheckGenLev_Chiara.root";
    TString root_file = "CheckGenLev_Livia.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}


TString pngName(const TString& histo ){
  //TString png_name = "./CheckGenLev_C/" + histo ;
  TString png_name = "./CheckGenLev_L/" + histo;
  
  return png_name + ".png";
}

Color_t PtlColorMap(const TString& particle){

  std::map <TString , Color_t> PtlColor{};
  PtlColor["mu"] = kPink + 5;
  PtlColor["pi"] = kBlue - 4;
  PtlColor["PiPi"] = kOrange + 6;
  PtlColor["JPsi"] = kRed;
  PtlColor["Rho"] = kOrange + 8;
  PtlColor["Ks"] = kGreen;
  PtlColor["X"] = kBlack + 3;
  PtlColor["B0"] = kViolet + 8;

  return PtlColor[particle];
}


TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["mu"] = "\\mu";
  Leg_entry["pi"] = "\\pi";
  Leg_entry["PiPi"] = "\\pi^+ \\pi^-";
  Leg_entry["JPsi"] = "J\\Psi";
  Leg_entry["Rho"] = "\\rho";
  Leg_entry["Ks"] = "K_s^0 \\";
  Leg_entry["X"] = "X(3872)";
  Leg_entry["B0"] = "B_0 \\";

  return Leg_entry[category];
}


int draw_pT_histo(const TString histo_name, const TString particle){
    
    TFile* input_file = open_file();

    TH1F* h = (TH1F*)input_file->Get(histo_name);
    
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }

    //AXIS & TITLE
    h->GetXaxis()->SetTitle("\\ p_T [GeV]");
    h->GetYaxis()->SetTitle("counts"); 

    //LINE COLOR MAP & LINE WIDTH
    Color_t color = PtlColorMap(particle);
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->SetFillColor(color);
    
    //LEGEND
    auto legend = new TLegend(0.75,0.8,.89,.89);
    legend->AddEntry(h,CategoryLegend(particle),"l");

    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = pngName(histo_name); 
    TCanvas* c1 = new TCanvas("c1","canvas", 1248,1024);
    h->Draw();
    legend->Draw();
    
    c1->SaveAs(png_name);

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

  //AXIS & TITLE                                                                                                         
  h->GetXaxis()->SetTitle("\\eta");
  h->GetYaxis()->SetTitle("counts");

  //LINE COLOR MAP & LINE WIDTH                                                                                          
  Color_t color = PtlColorMap(particle);
  h->SetLineWidth(2);
  h->SetLineColor(color);
  h->SetFillColor(color);

  //LEGEND                                                                                                               
  auto legend = new TLegend(0.75,0.8,.89,.89);
  legend->AddEntry(h,CategoryLegend(particle),"l");

  //STATISTICS                                                                                                           
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name); 
  TCanvas* c1 = new TCanvas("c1","canvas", 1248,1024);
  h->Draw();
  legend->Draw();

  c1->SaveAs(png_name);

  input_file->Close();

  return 0;

}//drae_Eta()

int draw_Mul(const TString histo_name, const TString particle){

  TString root_file = "analysis_Multiplicity.root";
  TFile* input_file = new TFile(root_file);

  TH1I* h = (TH1I*)input_file->Get(histo_name);

  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }

  //AXIS & TITLE                                                                                                                                                              
  h->GetXaxis()->SetTitle("n");
  h->GetYaxis()->SetTitle("counts");
  if(particle == "JPsi"){
	  int Nbins = h->GetNbinsX();
	  for (int i = 0; i < Nbins; i++) h->GetXaxis()->SetBinLabel(i+1,std::to_string(i).c_str());

  }
  h->GetXaxis()->SetLabelSize(0.05);
  //LINE COLOR MAP & LINE WIDTH                                                                                                                                               
  Color_t color = PtlColorMap(particle);
  h->SetLineWidth(2);
  h->SetLineColor(color);
  h->SetFillColor(color);

  //LEGEND                                                                                                                                                                    
  auto legend = new TLegend(0.75,0.8,.89,.89);
  legend->AddEntry(h,CategoryLegend(particle),"l");

  //STATISTICS                                                                                                                                                                
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name); 
  TCanvas* c1 = new TCanvas("c1","canvas", 1248,1024);
  h->Draw();
  legend->Draw();

  c1->SaveAs(png_name);

  input_file->Close();
  return 0;

}//draw_mul()   
