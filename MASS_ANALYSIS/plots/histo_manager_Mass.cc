#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLine.h>
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


TFile* open_file(){
    TString root_file = "invMassCutON.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["JPsiSGN"] = kRed;
  Color["JPsiBKG"] = kOrange;
  Color["RhoSGN"] = kGreen + 1;
  Color["PiPiBKG"] = kYellow - 7; 
  Color["Pions"] = kCyan + 2 ;
  Color["XSGN"] = kAzure + 1; //kViolet+ 1;
  Color["XGEN"] = kPink + 1;
  Color["JPsiPi"] = kOrange + 8;
  Color["XBBKG"] = kViolet;
  Color["XBKG"] = kTeal - 8;
  Color["BSGN"] = kAzure + 1;
  Color["BGEN"] = kOrange+ 1;
  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["JPsiSGN"] = "J/#psi MC-matched";
  Leg_entry["JPsiBKG"] = "Combinatory #mu^{+}#mu^{-}";
  Leg_entry["RhoSGN"] = "#rho(770)";
  Leg_entry["PiPiBKG"] = "Combinatory #pi^{+}#pi^{-}";
  Leg_entry["Pions"] = "ALL #pi";
  Leg_entry["XSGN"] = "X(3872) MC-matched"; 
  Leg_entry["XGEN"] = "X(3872) Generated"; 
  Leg_entry["JPsiPi"] = "J/#psi#pi MC-matched"; 
  Leg_entry["XBBKG"] = "Combinatory #mu^{+}#mu^{-} #pi^{+}#pi^{-}"; 
  Leg_entry["XBKG"] = "Combinatory J/#psi #pi^{+}#pi^{-}"; 
  Leg_entry["BSGN"] = "B_{0} MC-matched"; 
  Leg_entry["BGEN"] = "B_{0} Generated";

  return Leg_entry[category];
}


TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./CombinatorialCUTON/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
}


TString pdfName(const TString& histo1_name, const TString& category2){
  TString pdf_name = "./CombinatorialCUTON/" + histo1_name ;
  if (category2 != "") pdf_name += "_" + category2;
  
  return pdf_name + ".pdf";
}

void histo_SetUp(TH1* histo, const TString& category, const TString& x_name, const bool fill = true, const bool norm = true ){
  //AXIS LABEL                                                                                                                                                                  
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.4);
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitle(Form("Entries / %.3f GeV", histo->GetXaxis()->GetBinWidth(1)));
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.04);

  //WIDTH & COLOR                                                                                                                                                               
  histo->SetLineWidth(3);
  histo->SetLineColor(CategoryColorMap(category));
  if (fill) histo->SetFillColor(CategoryColorMap(category));
  else histo->SetFillColorAlpha(CategoryColorMap(category), .3);
  //NORMALIZATION
  if (norm) {histo->Scale(1./histo->Integral());
	  histo->GetYaxis()->SetTitle(Form("1/N dN/ %.3f GeV", histo->GetXaxis()->GetBinWidth(1)));}
  
  //SET MAXIMUM
  double M = histo->GetBinContent(histo->GetMaximumBin());
  histo->SetMaximum(1.2*M);

	gStyle->SetLineWidth(3);
}



int draw_single_histogram(const TString& histo_name, const TString& x_name, const TString& category, const bool norm = true){
    
    TFile* input_file = open_file();
    TH1F * h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    
    //AXIS & TITLE
    h->SetTitle("");
    histo_SetUp(h, category, x_name, true, norm);
    
    //STATISTICS
    gStyle->SetOptStat(0);

    //LEGEND
	 auto legend = new TLegend(0.55,0.75,.85,.80);
	 legend->SetBorderSize(0);
	 legend->SetTextSize(0.035);
	 legend->AddEntry(h, CategoryLegend(category) ,"f");

    TString png_name = pngName(histo_name, "");
    TString pdf_name = pdfName(histo_name, "");
    TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    h->Draw("HIST");
    legend->Draw();
    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

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
  THStack * hStack = new THStack("hStack",";"+x_name+";"+y_name);
  hStack->Add(h1);
  hStack->Add(h2);
  hStack->SetMaximum(1.45*hStack->GetMaximum());
  //hStack->GetXaxis()->SetTitleOffset(1.5);

  //STATISTICS                                                                                                                                                                
  gStyle->SetOptStat(0);
  
  //LEGEND                                                                                                                                                                    
  auto legend = new TLegend(0.20 ,0.73,.50,.85);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.035);
  legend->AddEntry(h1, CategoryLegend(category1) ,"f");
  legend->AddEntry(h2, CategoryLegend(category2) ,"f");


  TString png_name = pngName("STK_"+histo1, category2);
  TString pdf_name = pdfName("STK_"+histo1, category2);
  TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  hStack->Draw();
  legend->Draw();

  c1->SaveAs(png_name);
  c1->SaveAs(pdf_name);

  input_file->Close();

  return 0;
}

int draw_two_histograms(const TString histo1,const TString& category1, const TString histo2, const TString& category2, const TString& x_name,bool norm = false,  bool fill = true, const bool log = false){

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
    
    histo_SetUp(h1, category1, x_name, fill, norm);
    histo_SetUp(h2, category2, x_name, fill, norm);

    //STATISTICS
    gStyle->SetOptStat(0);
    
    //SETMAXIMUM                                                                                                                                                                  
    double M1 = h1->GetBinContent(h1->GetMaximumBin());
    double M2 = h2->GetBinContent(h2->GetMaximumBin());
    if (M1 > M2){ h1->SetMaximum(1.4*M1);
    }else {h1->SetMaximum(1.4*M2);}
    
    //LEGEND
    auto legend = new TLegend(0.50,0.70,.76,.85);
	legend->SetBorderSize(0);
	legend->SetTextSize(0.035);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    legend->AddEntry(h2, CategoryLegend(category2) ,"f");
    

    TString png_name = pngName(histo1, category2);
    TString pdf_name = pdfName(histo1, category2);
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
    if(log) gPad->SetLogy();
    h1->Draw("HIST");
    h2->Draw("HIST SAME");
	gPad->SetBottomMargin(0.13); gPad->SetLeftMargin(0.13);
    gPad->RedrawAxis();
    legend->Draw();

    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

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
  
  TCanvas* c2 = new TCanvas("c2","canvas", 1024,1024);
  
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
  Color_t color = CategoryColorMap("match");
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
  TString pdf_name = pdfName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
  h->Draw("HIST");
  t->Draw("same");
  c1->SaveAs(png_name);
  c1->SaveAs(pdf_name);

  input_file->Close();
  

  return 0;
}

int draw_NReco(const TString histo_name, const TString& category = "Pions"){

  TFile* input_file = open_file();
  TH1F* h = (TH1F*)input_file->Get(histo_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo_name << std::endl;
    exit(-1);
  }

  //TITLE & AXIS
  TString title = "Reconstructed tracks per event";
  if (category == "OtherP") title = "Number of muons from other particles per event";
  h->SetTitle(title);
  const int nx = 3;
  std::string os_X[nx]   = {"0","1","2"};
  for (int i=1; i<=nx; i++) {
    h->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
  }

  histo_SetUp(h, "NReco" , "");
  h->GetXaxis()->SetLabelSize(0.05);
  h->SetMaximum(1.2);
  //STATISTICS                                                                                                                                                                  
  gStyle->SetOptStat(0);

  //LEGEND                                                                                                                                                                      
  auto legend = new TLegend(0.62,0.8,.89,.89);
	legend->SetBorderSize(0);
  legend->AddEntry(h, CategoryLegend(category) ,"f");

  TString png_name = pngName(histo_name, "");
  TString pdf_name = pdfName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

  gPad->SetLogy();
  h->Draw("HIST");
  legend->Draw();
  c1->SaveAs(png_name);
  c1->SaveAs(pdf_name);
  input_file->Close();

  return 0;
}
