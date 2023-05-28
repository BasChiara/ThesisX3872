#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

#define NEVENTS_ 3936
using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


TFile* open_file(){
    TString root_file = "analysis_Tracks.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["disc"] = kBlue - 7;
  Color["rec"] = kRed + 1;
  Color["sel2"] = kOrange - 7;
  Color["gen"] = kOrange;
  Color["DR"] = kRed;
  Color["NReco"] = kAzure + 1;
  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["disc"] = "OTHER TRACKS";
  Leg_entry["rec"] = "SIGNAL PIONS ";
  Leg_entry["sel2"] = "RECO DR < 0.1 PIONS ";
  Leg_entry["gen"] = "GENERATED PIONS";
  Leg_entry["NReco"] = "SIGNAL PIONS";

  return Leg_entry[category];
}


TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./Rho_PiPi/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
}


void histo_SetUp(TH1* histo, const TString& category, const TString& x_name,  int norm = 1, bool fill = true ){
  //AXIS LABEL                                                                                                                                                                  
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.03);
  histo->GetXaxis()->SetLabelSize(0.03);
  histo->GetYaxis()->SetTitle("counts");
  histo->GetYaxis()->SetLabelSize(0.03);

  if (category == "DR"){
    histo->GetXaxis()->SetTitleOffset(1.);
    histo->GetYaxis()->SetTitleOffset(1.);
  }
  //WIDTH & COLOR
  histo->SetLineWidth(4);
  histo->SetLineColor(CategoryColorMap(category));
  if(fill) histo->SetFillColorAlpha(CategoryColorMap(category), 0.4);
  
  //NORMALIZATION
  float factor = 1.;
  if (norm == 1) factor/= histo->Integral(); //1 to normalize to 1
  if (norm == 2) factor/= NEVENTS_;//2 to normalize with the totality of the events
  histo->Scale(factor);
  
}

int draw_single_histogram(const TString histo_name, const TString x_name){
    
    TFile* input_file = open_file();
    TH1F * h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    
 
    //AXIS & TITLE
    h->SetTitle("Minimum\\ \\Delta R\\ &\\ \\Delta p_T < 0.5");
    h->GetXaxis()->SetTitle(x_name);
    h->GetYaxis()->SetTitle("counts"); 

    //LINE COLOR AND WIDTH
    Color_t color = kRed;
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->SetFillColor(color);
    
    //LEGEND
    auto legend = new TLegend(0.7,0.8,.9,.9);
    
    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = pngName(histo_name, "");
    TCanvas* c1 = new TCanvas("c1","canvas", 1248,1024);
    h->Draw();
    //legend->Draw();
    
    c1->SaveAs(png_name);

    input_file->Close();
    return 0;

}


int draw_two_histograms(const TString histo1,const TString& category1, const TString histo2, const TString& category2, const TString& x_name, const int& norm = 2 , bool fill = true){

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

    //HISTO SETUP
    histo_SetUp(h1, category1, x_name, norm, fill);
    histo_SetUp(h2, category2, x_name, norm, fill);
    if (category1 == "disc") h2->Scale(100.);  
    
    //SETMAXIMUM
    double M1 = h1->GetBinContent(h1->GetMaximumBin());
    double M2 = h2->GetBinContent(h2->GetMaximumBin());
    if (M1 > M2){ h1->SetMaximum(1.2*M1);
    }else { h1->SetMaximum(1.2*M2);}

    //STATISTICS
    gStyle->SetOptStat(0);

    //LEGEND
    auto legend = new TLegend(0.6,0.8,.89,.89);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    if (category1 == "disc")legend->AddEntry(h2, CategoryLegend(category2)+" x 100" ,"f");    
    else legend->AddEntry(h2, CategoryLegend(category2) ,"f");
    

    TString png_name = pngName(histo1, category2);
    TCanvas* c1 = new TCanvas("c1","canvas", 1248,1024);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    gPad->RedrawAxis();
    legend->Draw();

    c1->SaveAs(png_name);

    input_file->Close();

    return 0;
}



int draw_QC_histo(const TString hReco, const TString hDisc, const TString& title){

  TFile* input_file = open_file();

  TH1F* h1 = (TH1F*)input_file->Get(hReco);
  TH1F* h2 = (TH1F*)input_file->Get(hDisc);

  if ( !h1 ){
    std::cout<< "null pointer for histogram named " << hReco << std::endl;
    exit(-1);
  }
  if ( !h2 ){
    std::cout<< "null pointer for histogram named " << hDisc << std::endl;
    exit(-1);
  }

  TString category1 = "rec", category2 = "disc";
  //SETUP                                                                                                                                                                       
  histo_SetUp(h1, category1, "", true);
  h1->SetTitle(title);
  histo_SetUp(h2, category2, "", true);

  //STATISTICS
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.4f");
  //LEGEND
  auto legend = new TLegend(0.62,0.8,.89,.89);
  legend->AddEntry(h1, CategoryLegend(category1) ,"f");
  legend->AddEntry(h2, CategoryLegend(category2) ,"f");

  TString png_name = pngName(hReco, category2);
  TCanvas* c1 = new TCanvas("c1","canvas", 864, 720);
  gPad->SetLogy();
  h1->Draw("HIST TEXT0");
  h2->Draw("HIST TEXT0 SAME");
  gPad->RedrawAxis();
  legend->Draw();
  c1->SaveAs(png_name);

  input_file->Close();

  return 0;

}



int draw_2Dhisto(const TString histo2D_name){

  TFile* input_file = open_file();
  TH2F * h = (TH2F*)input_file->Get(histo2D_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo2D_name << std::endl;
    exit(-1);
  }

  //AXIS & TITLE
  h->GetXaxis()->SetTitle("\\Delta R");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("\\frac{\\Delta p_T}{p_T^G}");
  h->GetYaxis()->CenterTitle();
  
  //LINE COLOR AND WIDTH                                                                                                                                                      

  TLine* lh = new TLine(0.,0.5,0.03,0.5);
  lh->SetLineColor(CategoryColorMap("DR"));
  lh->SetLineWidth(3);
  TLine* lv = new TLine(0.03,0.,0.03,0.5);
  lv->SetLineColor(CategoryColorMap("DR"));
  lv->SetLineWidth(3);
  //STATISTICS
  gStyle->SetOptStat(0);
  
  TString png_name = pngName(histo2D_name,"");

  TCanvas* c1 = new TCanvas("c1","canvas", 1520,1024);
  c1->SetPhi(200);

  h->Draw("COLZ");
  lh->Draw("same");
  lv->Draw("same");
  c1->SaveAs(png_name);
  
  TCanvas* c2 = new TCanvas("c2","canvas", 1248,1024);
  
  TH1* hx = h->ProjectionX();
  histo_SetUp(hx, "DR" ,"\\Delta R");
  hx->Draw("HIST");
  hx->SetTitle("");
  
  c2->SaveAs(pngName("projX_"+histo2D_name, ""));

  TH1* hy = h->ProjectionY();
  histo_SetUp(hy, "DR" ,"\\Delta p_T");
  hy->SetTitle("");
  hy->Draw("HIST");
  
  c2->SaveAs(pngName("projY_"+histo2D_name, ""));

  input_file->Close();

  return 0;

}//draw_2Dhisto()

int draw_matching(const TString histo_name, const TString title){

  TFile* input_file = open_file();
  TH1F* h = (TH1F*)input_file->Get(histo_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }
  
  //TITLE & AXIS
  h->SetTitle(title);
  h->GetXaxis()->SetRangeUser(-0.5, 2.);
  h->GetYaxis()->SetTitle("counts/total");
  std::cout << " NM " << h->GetBinContent(1) << std::endl;
  h->Scale(1./h->Integral());
  h->GetYaxis()->SetRangeUser(0., 1.);  
  h->SetCanExtend(TH1::kAllAxes);
  
  //PRINT SCORE
  float percent = h->GetBinContent(3);
  TString text = std::to_string(percent); 
  TText *t = new TText(1.,.3, text);
  t->SetTextAlign(22);
  t->SetTextColor(kBlack);
  t->SetTextFont(43);
  t->SetTextSize(25);

  //LINE COLOR AND WIDTH
  Color_t color = kOrange - 3;
  h->SetStats(0);
  h->SetBarWidth(1.);
  h->SetLineColor(color);
  h->SetFillColor(color);

  const int nx = 3;
  std::string os_X[nx]   = {"LOST"," ","RECONSTRUCTED"};
  for (int i=1; i<=nx; i++) {
    h->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
  }

  //STATISTICS
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
  h->Draw("HIST");
  t->Draw("same");
  c1->SaveAs(png_name);

  input_file->Close();
  

  return 0;
}

int draw_NReco(const TString histo_name){

  TFile* input_file = open_file();
  TH1F* h = (TH1F*)input_file->Get(histo_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }

  TString category = "NReco";   
  //TITLE & AXIS                                                                                                                                                                  
  h->SetTitle("Reconstructed tracks per event");
  const int nx = 3;
  std::string os_X[nx]   = {"0","1","2"};
  for (int i=1; i<=nx; i++) {
    h->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
  }
  
  //SETUP histo_SetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true )
  histo_SetUp( h, category, "");
  h->GetXaxis()->SetLabelSize(0.05);
  h->SetMaximum(2.);
  //STATISTICS
  gStyle->SetOptStat(0);
  
  auto legend = new TLegend(0.62,0.85,.89,.89);
  legend->AddEntry(h, CategoryLegend(category) ,"f");

  TString png_name = pngName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
  gPad->SetLogy();
  h->Draw("HIST TEXT0");
  legend->Draw();
  
  c1->SaveAs(png_name);
  input_file->Close();


  return 0;
}
