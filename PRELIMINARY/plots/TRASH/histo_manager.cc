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


int draw_single_histogram(const TString histo_name, const TString x_name){
    
    
    TFile* input_file = new TFile("analysis.root");

    if ( !input_file->IsOpen() ) {
       std::cout << "ERRORE ... " << std::endl;
       exit(-1);
    }

    TH1F* h = (TH1F*)input_file->Get(histo_name);
    
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }

    //AXIS
    h->GetXaxis()->SetTitle(x_name);
    h->GetYaxis()->SetTitle("counts"); 

    //LINE COLOR AND WIDTH
    h->SetLineWidth(2);
    h->SetLineColor(4);
    
    //LEGEND
    auto legend = new TLegend(0.7,0.8,.9,.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header                                   
    legend->AddEntry(h,"ALL MUONS","l");
    //legend->AddEntry(h2,"MINIMUM\\ \\Delta R","l")

    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = histo_name + ".png";
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
    h->Draw();
    legend->Draw();
    
    c1->SaveAs(png_name);

    input_file->Close();
    return 0;

}

int draw_two_histograms(const TString histo1, const TString histo2, const TString x_name){

    TFile* input_file = new TFile("analysis.root");

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR OPENING FILE" << std::endl;
       exit(-1);
    }

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
    
    //AXIS
    //TString x_name = "\\Delta R \\mu" ;
    h1->GetXaxis()->SetTitle(x_name);
    h1->GetYaxis()->SetTitle("counts"); 
    h1->GetXaxis()->SetTitle(x_name);
    h1->GetYaxis()->SetTitle("counts"); 
    
    //LINE COLOR AND WIDTH
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kRed);

    //STATISTICS
    gStyle->SetOptStat(0);

    //LEGEND
    auto legend = new TLegend(0.7,0.8,.9,.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(h1,"ALL MUONS","l");
    legend->AddEntry(h2,"MINIMUM\\ \\Delta R","l");
    

    TString png_name = histo1 + ".png";
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
    h1->Draw();
    h2->Draw("SAME");
    legend->Draw();

    c1->SaveAs(png_name);

    input_file->Close();

    return 0;
}


int draw_lead_pt(const TString h_leading, const TString h_subleading){

  TFile* input_file = new TFile("analysis.root");

  if ( !input_file->IsOpen() ) {
    std::cout << "ERROR OPENING FILE" << std::endl;
    exit(-1);
  }

  TH1F* h_l = (TH1F*)input_file->Get(h_leading);
  TH1F* h_sl = (TH1F*)input_file->Get(h_subleading);

  if ( !h_l ){
    std::cout<< "null pointer for histogram named " << h_leading << std::endl;
    exit(-1);
  }
  if ( !h_sl){
    std::cout<< "null pointer for histogram named " << h_subleading << std::endl;
    exit(-1);
  }

  //AXIS                                                                                                                   
  TString x_name = "p_T [GeV]";
  h_l->GetXaxis()->SetTitle(x_name);
  h_l->GetYaxis()->SetTitle("counts");
  h_sl->GetXaxis()->SetTitle(x_name);
  h_sl->GetYaxis()->SetTitle("counts");

  //LINE COLOR AND WIDTH                                                                                                   
  h_l->SetLineWidth(3);
  h_sl->SetLineWidth(3);
  h_l->SetLineColor(kMagenta + 1);
  h_sl->SetLineColor(kCyan - 6);
  h_l->SetFillColor(kMagenta + 1);
  h_sl->SetFillColor(kCyan - 6);
  //STATISTICS                                                                                                             
  gStyle->SetOptStat(0);

  //LEGEND                                                                                                                 
  auto legend = new TLegend(0.7,0.8,.9,.9);
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry(h_l,"LEADING","l");
  legend->AddEntry(h_sl,"SUBLEADING","l");


  TString png_name = h_leading + ".png";
  TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
  h_sl->Draw();
  h_l->Draw("SAME");
  legend->Draw();

  c1->SaveAs(png_name);

  input_file->Close();

  return 0;
}

int draw_matching(const TString histo_name, const TString x_name){


  TFile* input_file = new TFile("analysis.root");

  if ( !input_file->IsOpen() ) {
    std::cout << "ERRORE ... " << std::endl;
    exit(-1);
  }

  TH1F* h = (TH1F*)input_file->Get(histo_name);

  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }
  
  //AXIS                                                                                                                                                                      
  h->GetXaxis()->SetTitle(x_name);
  h->GetXaxis()->SetRangeUser(-0.5, 2.);
  h->GetYaxis()->SetTitle("counts/events");
  
  h->Scale(1./h->Integral());  
  
  //LINE COLOR AND WIDTH
  h->SetStats(0);
  //h->SetLineWidth(2);
  h->SetBarWidth(1.);
  h->SetLineColor(kOrange);
  h->SetFillColor(kOrange);

  const int nx = 3;
  std::string os_X[nx]   = {"NON MATCHING"," ","MATCHING"};
  for (int i=1; i<=nx; i++) {
    h->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
  }

  //STATISTICS
  gStyle->SetOptStat(0);

  TString png_name = histo_name + ".png";
  TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
  h->Draw("HIST");

  c1->SaveAs(png_name);

  input_file->Close();
  

  return 0;
}
