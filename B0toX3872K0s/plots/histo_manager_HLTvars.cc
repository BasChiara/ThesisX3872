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
#include <THStack.h>


using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


TFile* open_file(){
    TString root_file = "HLT_observables.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["fired"] = kOrange;
  Color["pT"] = kViolet;
  Color["eta"] = kGreen + 1;
  Color["M"] = kBlue - 7;
  Color["DCA"] = kOrange + 1;
  Color["sign"] = kAzure;
  Color["cosA"] = kMagenta;
  Color["SVp"] = kSpring -9;

  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["PRE"] = "PREFIT";
  Leg_entry["VTX"] = "VERTEX FIT";
  Leg_entry["TOT"] = "FINAL FIT";

  return Leg_entry[category];
}

TString PlotTitle(const TString& category){

  std::map <TString , TString> Title{};
  Title["MuTrigger_L0"] = "Selection on single muon";
  Title[""] = "VERTEX FIT";
  Title[""] = "FINAL FIT";

  return Title[category];
}



TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./HLT_observables/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
}


void histo_SetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true ){
  //AXIS LABEL 
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitle("counts");
  histo->GetYaxis()->SetLabelSize(0.03);

  if (category == "??"){
    histo->GetXaxis()->SetTitleOffset(1.);
    histo->GetYaxis()->SetTitleOffset(1.);
  }

  //WIDTH & COLOR                                                                                                                                                               
  histo->SetLineWidth(4);
  histo->SetLineColor(CategoryColorMap(category));
  if (fill) histo->SetFillColor(CategoryColorMap(category));
  histo->SetFillColorAlpha(CategoryColorMap(category), 0.3);
  //NORMALIZATION
  histo->Scale(1./histo->Integral());
}

int draw_multiple_histo(std::vector<TString> histos, std::vector<TString> categories, std::vector<TString> x_name, TString name){

   TFile* input_file = open_file();

   TCanvas* c1 = new TCanvas("c1","canvas", 3072, 2048);
   c1->Divide(2 ,2, 0.00001, 0.00001);
   gStyle->SetOptStat(0);
   int NH = histos.size();

   for (int i = 0; i < NH; i++){

      TH1F * h = (TH1F*)input_file->Get(histos[i]);
      if ( !h ){
	 std::cout<< "null pointer for histogram named " << histos[i] << std::endl;
	 exit(-1);
      }

      histo_SetUp(h, categories[i], x_name[i]);
      c1->cd(i+1);
      h->Draw("HIST");
   }
   c1->SaveAs(pngName(name, ""));
   return 0;
}//draw_multiple_histo()

int draw_single_histogram(const TString& histo_name, const TString& x_name, const TString& category){
    
    TFile* input_file = open_file();
    TH1F * h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    
    //AXIS & TITLE
    h->SetTitle("");
    h->GetXaxis()->SetTitle(x_name);
    h->GetYaxis()->SetTitle("counts"); 

    //LINE COLOR AND WIDTH
    Color_t color = CategoryColorMap(category);
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->SetFillColor(color);
    
    //LEGEND
    
    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = pngName(histo_name, "");
    TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
    h->Draw();
    c1->SaveAs(png_name);

    input_file->Close();
    return 0;

}


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
    
    histo_SetUp(h1, category1, x_name, fill);
    histo_SetUp(h2, category2, x_name, fill);

    //STATISTICS
    gStyle->SetOptStat(0);
    
    //SETMAXIMUM                                                                                                                                                                  
    double M1 = h1->GetBinContent(h1->GetMaximumBin());
    double M2 = h2->GetBinContent(h2->GetMaximumBin());
    if (M1 > M2){ h1->SetMaximum(1.2*M1);
    }else {h1->SetMaximum(1.2*M2);}
    
    //LEGEND
    auto legend = new TLegend(0.62,0.8,.89,.89);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    legend->AddEntry(h2, CategoryLegend(category2) ,"f");
    

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



int draw_2Dhisto(const TString histo2D_name){

  TFile* input_file = open_file();
  TH2F * h = (TH2F*)input_file->Get(histo2D_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo2D_name << std::endl;
    exit(-1);
  }

  //AXIS & TITLE
  TString title = "Minimum\\ \\Delta R\\ vs\\ \\Delta p_T";
  h->SetTitle(title);
  h->GetXaxis()->SetTitle("\\Delta R");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("\\Delta p_T");
  h->GetYaxis()->CenterTitle();

  //LINE COLOR AND WIDTH                                                                                                                                                      
  Color_t color = CategoryColorMap("DR");
  h->SetLineWidth(2);
  h->SetLineColor(color);

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


int draw_Nfired(const TString histo_name, const TString title){

  TFile* input_file = open_file();
  TH1F* h = (TH1F*)input_file->Get(histo_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }
  
  //TITLE & AXIS
  h->SetTitle(title);
  h->GetYaxis()->SetTitle("counts/total");
  h->Scale(1./h->Integral());
  h->GetYaxis()->SetRangeUser(0., 1.1);  
  
  //LINE COLOR AND WIDTH
  Color_t color = CategoryColorMap("fired");
  h->SetStats(0);
  h->SetBarWidth(1.);
  h->SetLineColor(color);
  h->SetFillColorAlpha(color, 0.3);

  const int nx = 2;
  std::string os_X[nx]   = {"\\ K_S^0","\\rho(770)"};
  for (int i=1; i<=nx; i++) {
    h->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
  }

  //STATISTICS
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
  h->Draw("HIST TEXT0");
  c1->SaveAs(png_name);

  input_file->Close();
  

  return 0;
}

