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

#define NEVENTS_ 3527
using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


TFile* open_file(){
    TString root_file = "analysis_K0REDO.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["disc"] = kAzure - 9;
  Color["rec"] = kGreen + 1;
  Color["gen"] = kOrange +1;
  Color["prefit"] = kOrange + 7;
  Color["womc"] = kViolet;
  Color["DR"] = kRed;
  
  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["disc"] = "OTHER K^{0}_{s}";
  Leg_entry["rec"] = "MC-MATCHED K^{0}_{s}";
  Leg_entry["gen"] = "GENERATED K^{0}_{s}";
  Leg_entry["prefit"] = " K0s (No fit)";
  Leg_entry["womc"] = "K0s (VTX)fit";  

  return Leg_entry[category];
}


TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./K0short/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
}
TString pdfName(const TString& histo1_name, const TString& category2){
  TString pdf_name = "./K0short/" + histo1_name ;
  if (category2 != "") pdf_name += "_" + category2;
  
  return pdf_name + ".pdf";
}


void histo_SetUp(TH1* histo, const TString& category, const TString& x_name,  int norm = 1, bool fill = true ){
  //AXIS LABEL                                                                                                                                                                  
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitle("counts");
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.04);
  gStyle->SetLineWidth(3);

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


void GetCMS(TCanvas* c){
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
	RunDetails.DrawLatex(.70, .91, "(#sqrt{s} = 13 TeV) 2017");
}

int draw_single_histogram(const TString histo_name, const TString category, const TString x_name){
    
    TFile* input_file = open_file();
    TH1F * h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    
    //AXIS & TITLE
    histo_SetUp(h, category, x_name, 1); 
	 h->SetMaximum(1.3*h->GetMaximum());
   h->GetXaxis()->SetNdivisions(505,kTRUE); 
    //LEGEND
    auto legend = new TLegend(0.53,0.75,.78,.80);
	 legend->SetBorderSize(0);
	 legend->SetTextSize(0.035);
    
    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = pngName(histo_name, "");
    TString pdf_name = pdfName(histo_name, "");
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
    h->Draw("HIST");
    legend->AddEntry(h, CategoryLegend(category), "f");
    legend->Draw();
	 //GetCMS(c1);    
	 gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

    input_file->Close();
    return 0;

}

int draw_many_histo(std::vector<TString> histos, std::vector<TString> categories, const TString & x_name){
	TFile* input_file = open_file();
	//FILL THE STACK

	THStack* Stk = new THStack("hStack",";"+x_name+";counts");

	//LEGEND
	UInt_t NH = histos.size();
	auto legend = new TLegend(0.55,1 - 0.08*NH,.89,.89);
	legend->SetBorderSize(0);

	for (UInt_t i = 0; i < NH; i++){
		TH1* h = (TH1*)input_file->Get(histos[i]);
		histo_SetUp(h, categories[i], "", 0);
		Stk->Add(h);
		legend->AddEntry(h, CategoryLegend(categories[i]) ,"f");
	}

	Stk->SetMaximum(1.1*Stk->GetMaximum());

	//DRAW
	TString png_name = pngName("STK"+histos[1], "");
	TString pdf_name = pdfName("STK"+histos[1], "");
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	Stk->Draw("nostack HIST");
	legend->Draw();
	//GetCMS(c1);    
	c1->SaveAs(png_name);
	c1->SaveAs(pdf_name);

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
    
    //SETMAXIMUM
    double M1 = h1->GetBinContent(h1->GetMaximumBin());
    double M2 = h2->GetBinContent(h2->GetMaximumBin());
    if (M1 > M2){ h1->SetMaximum(1.2*M1);
    }else { h1->SetMaximum(1.2*M2);}

    //STATISTICS
    gStyle->SetOptStat(0);

    //LEGEND
    auto legend = new TLegend(0.50,0.70,.80,.85);
	 legend->SetBorderSize(0);
	 legend->SetTextSize(0.035);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    legend->AddEntry(h2, CategoryLegend(category2) ,"f");

    TString png_name = pngName(histo1, category2);
    TString pdf_name = pdfName(histo1, category2);
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
	 gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    gPad->RedrawAxis();
    legend->Draw();
	 //GetCMS(c1);    

    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

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
  auto legend = new TLegend(0.55,0.8,.89,.89);
	legend->SetBorderSize(0);
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
  h->GetXaxis()->SetTitle("#DeltaR_{min}");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("#frac{#Deltap_{T}}{p_{T}^{G}}");
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

  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
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
  TString pdf_name = pdfName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
  h->Draw("HIST");
  t->Draw("same");
  c1->SaveAs(png_name);
  c1->SaveAs(pdf_name);

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

  TString category = "disc";   
  //TITLE & AXIS                                                                                                                                                                  
  const int nx = 10;
  std::string os_X[nx]   = {"0","1","2"};
  for (int i=1; i<=nx; i++) {
    h->GetXaxis()->SetBinLabel(i,std::to_string(i-1).c_str());
  }
  
  //SETUP histo_SetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true )
  histo_SetUp( h, category, "# candidate K_{s}^{0}");
  h->GetXaxis()->SetLabelSize(0.05);
  //STATISTICS
  gStyle->SetOptStat(0);
  
  auto legend = new TLegend(0.55,0.85,.89,.89);
  legend->SetBorderSize(0);
  legend->AddEntry(h, CategoryLegend(category) ,"f");

  TString png_name = pngName(histo_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
  gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
  h->Draw("HIST");
  //legend->Draw();
  //GetCMS(c1);    
  
  c1->SaveAs(png_name);
  input_file->Close();


  return 0;
}
