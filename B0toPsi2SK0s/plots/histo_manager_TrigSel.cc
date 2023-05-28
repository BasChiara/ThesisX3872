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
#include <TLatex.h>


using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


TFile* open_file(){
    TString root_file = "TriggerSelection.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["MC"] = kViolet - 9;
  Color["DR"] = kRed;
  Color["JPsi"] = kRed;
  Color["Psi2S"] = kAzure + 1;
  Color["B0"] = kBlue;

  Color["PRE_SGN_JPsi"] = kRed;
  Color["PRE_BKG_JPsi"] = kOrange;
  Color["VTX_SGN_JPsi"] = kRed - 7;
  Color["VTX_BKG_JPsi"] = kGreen -7;

  Color["PRE_SGN_Rho"] = kGreen + 1;
  Color["PRE_BKG_Rho"] = kYellow - 7;

  Color["VTX_SGN_K0s"] = kTeal  + 2;
  Color["VTX_BKG_K0s"] = kPink - 6;

  Color["PRE_SGN_Psi2S"] = kAzure + 1;
  Color["PRE_BKG_Psi2S"] = kOrange + 6;
  Color["TOT_SGN_Psi2S"] = kBlue - 9;
  Color["TOT_BKG_Psi2S"] = kMagenta - 6;

  Color["PRE_SGN_B0"] = kViolet + 1;
  Color["PRE_BKG_B0"] = kRed - 7;
  Color["VTX_SGN_B0"] = kGreen - 6;
  Color["VTX_BKG_B0"] = kOrange + 1;
  Color["TOT_SGN_B0"] = kAzure + 1;
  Color["TOT_BKG_B0"] = kRed;

  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["JPsi"] = "J/#psi CANDIDATES ";
  Leg_entry["B0"] = "B_{0} CANDIDATES ";

  Leg_entry["PRE_SGN_JPsi"] = "MC-matching (pre-fit)";
  Leg_entry["PRE_BKG_JPsi"] = "FAKE (prefit)";
  Leg_entry["VTX_SGN_JPsi"] = "MC-matching J/#psi";
  Leg_entry["VTX_BKG_JPsi"] = "FAKE J/#psi";

  Leg_entry["PRE_SGN_Rho"] = "MC-matching #pi^{+}#pi^{-}";
  Leg_entry["PRE_BKG_Rho"] = "FAKE #pi^{+}#pi^{-}";

  Leg_entry["VTX_SGN_K0s"] = "MC-matching K_{0}^{s}";
  Leg_entry["VTX_BKG_K0s"] = "FAKE K_{0}^{s}";

  Leg_entry["PRE_SGN_Psi2S"] = "MC-matching pre-fit";
  Leg_entry["PRE_BKG_Psi2S"] = "FAKE pre-fit";
  Leg_entry["TOT_SGN_Psi2S"] = "MC-matching #psi(2S)";
  Leg_entry["TOT_BKG_Psi2S"] = "FAKE #psi(2S)";

  Leg_entry["PRE_SGN_B0"] = "MC-matching prefit";
  Leg_entry["PRE_BKG_B0"] = "FAKE prefit";
  Leg_entry["VTX_SGN_B0"] = "MC-matching vtx-fit";
  Leg_entry["VTX_BKG_B0"] = "FAKE vtx-fit";
  Leg_entry["TOT_SGN_B0"] = "MC-matching B^{0}";
  Leg_entry["TOT_BKG_B0"] = "FAKE B^{0}";

  return Leg_entry[category];
}


TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./TriggerSelection/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
}

TString pdfName(const TString& histo1_name, const TString& category2){
  TString pdf_name = "./TriggerSelection/" + histo1_name ;
  if (category2 != "") pdf_name += "_" + category2;
  
  return pdf_name + ".pdf";
}

void histo_SetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true , bool norm = true){
  //AXIS LABEL
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitle("counts");
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.04);
	gStyle->SetLineWidth(3);

  //WIDTH & COLOR                                                                                                                                                               
  histo->SetLineWidth(4);
  histo->SetLineColor(CategoryColorMap(category));
  if (fill)histo->SetFillColorAlpha(CategoryColorMap(category), 0.3);
  //NORMALIZATION
  if(norm) histo->Scale(1./histo->Integral());
}



int draw_single_histogram(const TString& histo_name, const TString& category, const TString& x_name){
    
    TFile* input_file = open_file();
    TH1F * h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    
	 histo_SetUp(h, category, x_name);   
    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = pngName(histo_name, "");
    TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
    h->Draw("HIST");
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

int draw_many_histo(std::vector<TString> histos, std::vector<TString> categories, const TString & x_name){
  TFile* input_file = open_file();
  //FILL THE STACK
  THStack* Stk = new THStack("hStack",";"+x_name+";counts");;
  
  //LEGEND
  UInt_t NH = histos.size();
  auto legend = new TLegend(0.62,1 - 0.08*NH,.89,.89);
	legend->SetTextSize(0.035);
  
  for (UInt_t i = 0; i < NH; i++){
    TH1* h = (TH1*)input_file->Get(histos[i]);
    histo_SetUp(h, categories[i], "", 0);
    
    Stk->Add(h);
    legend->AddEntry(h, CategoryLegend(categories[i]) ,"f");
  }
  
  Stk->SetMaximum(1.2*Stk->GetMaximum());
 
  //DRAW
  TString png_name = pngName("STK"+histos[1], "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
  Stk->Draw("nostack HIST");
  legend->Draw();
  c1->SaveAs(png_name);

  input_file->Close();
  return 0;
}


int draw2StackedHistos(const TString histo1,const TString& category1, const TString histo2, const TString& category2, const TString& x_name, bool fill = true, bool norm = false){

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

  //STACKING
  histo_SetUp(h1, category1, x_name, fill, norm);
  histo_SetUp(h2, category2, x_name, fill, norm);
  TString y_name = Form("Entries / %.3f GeV", h1->GetXaxis()->GetBinWidth(1));
  THStack * hStack = new THStack("hStack",";"+x_name+";"+ y_name);
  hStack->Add(h1);
  hStack->Add(h2);
  hStack->SetMaximum(1.35*hStack->GetMaximum());

  //STATISTICS                                                                                                                                                                
  gStyle->SetOptStat(0);
  
  //LEGEND                                                                                                                                                                    
  auto legend = new TLegend(0.50,0.73,.85,.85);
  legend->SetBorderSize(0);
	legend->SetTextSize(0.035);
  legend->AddEntry(h1, CategoryLegend(category1) ,"f");
  legend->AddEntry(h2, CategoryLegend(category2) ,"f");


  TString png_name = pngName("STK_"+histo1, category2);
  TString pdf_name = pdfName("STK_"+histo1, category2);
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

  hStack->Draw();
  legend->Draw();
	gPad->SetLeftMargin(0.13);  gPad->SetBottomMargin(0.13);
  c1->SaveAs(png_name);
  c1->SaveAs(pdf_name);

  input_file->Close();

  return 0;
}



int draw_2Dhisto(const TString histo2D_name, const TString ptl){

  TFile* input_file = open_file();
  TH2F * h = (TH2F*)input_file->Get(histo2D_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo2D_name << std::endl;
    exit(-1);
  }

  //AXIS & TITLE
  h->GetXaxis()->SetTitle("\\Delta R_{min}");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("\\frac{\\Delta p_T}{p_T}");
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
  if(ptl != "mu") lh->Draw("same");
  lv->Draw("same");
  c1->SaveAs(png_name); 

  input_file->Close();

  return 0;

}//draw_2Dhisto()


int draw_Ncandidates(const TString histo_name, const TString category, const TString title){

  TFile* input_file = open_file();
  TH1F* h = (TH1F*)input_file->Get(histo_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }
  
  //TITLE & AXIS
  h->SetTitle(title);
  histo_SetUp(h,category, "candidates");
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetRangeUser(0., 1.);

  const int Nb = h->GetNbinsX();
  for (int i=0; i < Nb; i++) {
    h->GetXaxis()->SetBinLabel(i+1, (std::to_string(i)).c_str());
  }

  //STATISTICS
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
  h->Draw("HIST");
  c1->SaveAs(png_name);

  input_file->Close();
  

  return 0;
}

int draw_ptlMCmissed(const TString histo_name, const TString category, const TString title){

  TFile* input_file = open_file();
  TH1F* h = (TH1F*)input_file->Get(histo_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }
  
  //TITLE & AXIS
  h->SetTitle(title);
  histo_SetUp(h,category, "", true, false);
  h->GetXaxis()->SetLabelSize(0.05);

  const int Nb = h->GetNbinsX();
  std::string os_X[Nb];
  if (category == "MC"){ 
	  os_X[0] = "J\\Psi";
	  os_X[1] = "\\rho";
	  os_X[2] = "\\ K_s^0";}
  if (category == "B0"){
	  os_X[0] = "J\\Psi";
	  os_X[1] = "\\pi \\ (\\rho)";
	  os_X[2] = "TRIG.\\pi \\ from\\ \\rho";
	  os_X[3] = "\\ K_s^0";
	  os_X[4] = "TRIG.\\pi \\(K_s^0)";}

  for (int i=0; i < Nb; i++) {
    h->GetXaxis()->SetBinLabel(i+1, (os_X[i].c_str()));
  }

  //STATISTICS
  gStyle->SetOptStat(0);

  TString png_name = pngName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
  h->Draw("HIST");
  c1->SaveAs(png_name);

  input_file->Close();
  

  return 0;
}

