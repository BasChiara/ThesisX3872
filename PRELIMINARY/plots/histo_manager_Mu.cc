#include <iostream>
#include <vector>
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

#define NEVENTS_ 3527 //3936

using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms()


TFile* open_file(){
    TString root_file = "analysis_MuonsREDO.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["disc"] = kMagenta - 7;
  Color["rec"] = kOrange + 1;
  Color["gen"] = kPink -7;
  Color["OtherP"] = kSpring - 6;
  Color["DR"] = kRed;
  Color["match"] = kOrange;
  Color["LMu"] = kAzure + 1;
  Color["sLMu"] = kViolet + 1;
  Color["NReco"] = kAzure + 1;
  
  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["disc"] = "OTHER MUONS";
  Leg_entry["rec"] = "MC-MATCHED MUONS ";
  Leg_entry["gen"] = "GENERATED MUONS";
  Leg_entry["OtherP"] = "OTHER PTL MUONS";
  Leg_entry["LMu"] = "LEADING MUON";
  Leg_entry["sLMu"] = "SUB-LEADING MUON";

  return Leg_entry[category];
}


TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./JPsi_MuMuREDO/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
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

void histo_SetUp(TH1* histo, const TString& category, const TString& x_name, int norm = 1, bool fill = true ){
  //AXIS LABEL                                                                                                                                                                  
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetYaxis()->SetTitle("counts");

  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitleSize(0.04);

  if (category == "DR"){
    histo->GetXaxis()->SetTitleOffset(1.);
    histo->GetYaxis()->SetTitleOffset(1.);
  }

  //WIDTH & COLOR
  gStyle->SetLineWidth(3);
  histo->SetLineWidth(4);
  histo->SetLineColor(CategoryColorMap(category));
  if (fill) histo->SetFillColorAlpha(CategoryColorMap(category), 0.3);

  //NORMALIZATION: 0 for noNorm
  float factor = 1.;
  if (norm == 1){
	 factor/= histo->Integral(); //1 to normalize to 1
	}
  if (norm == 2) factor/= NEVENTS_;//2 to normalize with the totality of the events
  histo->Scale(factor);
}



int draw_single_histogram(const TString& histo_name, const TString& x_name, const TString& category){
    
    TFile* input_file = open_file();
    TH1F * h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    
    //AXIS & TITLE
    h->SetTitle("");
    histo_SetUp(h, category, x_name);
	 h->GetYaxis()->SetTitle(Form("1/N_{EV} counts")); 
    
    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = pngName(histo_name, "");
    TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	 gPad->SetBottomMargin(0.13); gPad->SetLeftMargin(0.13);
    h->Draw("HIST");
	 //GetCMS(c1);
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
  auto legend = new TLegend(0.45,.85 - 0.05*NH,.85,.85);
  legend->SetBorderSize(0);
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
  TCanvas* c1 = new TCanvas("c1","canvas", 1064, 1024);

  Stk->Draw("nostack HIST");
  legend->Draw();
  gPad->SetBottomMargin(0.13); gPad->SetLeftMargin(0.13);
  //GetCMS(c1);
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
    
    histo_SetUp(h1, category1, x_name, 1, fill);
    histo_SetUp(h2, category2, x_name, 1, fill);
	 //h1->GetYaxis()->SetTitle(Form("1/N dN/%.2f GeV", h1->GetXaxis()->GetBinWidth(1)));

    //STATISTICS
    gStyle->SetOptStat(0);
    
    //SETMAXIMUM                                                                                                                                                                  
    double M1 = h1->GetBinContent(h1->GetMaximumBin());
    double M2 = h2->GetBinContent(h2->GetMaximumBin());
    if (M1 > M2){ h1->SetMaximum(1.4*M1);
    }else {h1->SetMaximum(1.4*M2);}
    
    //LEGEND
    auto legend = new TLegend(0.50,0.75,.85,.85);
	 legend->SetBorderSize(0);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    legend->AddEntry(h2, CategoryLegend(category2) ,"f");
    

    TString png_name = pngName(histo1, category2);
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);

	 gPad->SetBottomMargin(0.13); gPad->SetLeftMargin(0.13);
    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    gPad->RedrawAxis();
    legend->Draw();
	 //GetCMS(c1);

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
	  h->GetXaxis()->SetTitle("#Delta R_{min}");
	  h->GetXaxis()->CenterTitle(); 
	  h->GetXaxis()->SetTitleOffset(1.2);
	  h->GetXaxis()->SetNdivisions(5,5,kTRUE);
	  h->GetYaxis()->SetTitle("#frac{#Deltap_{T}}{p_{T}^{G}}");
	  h->GetYaxis()->CenterTitle(); 
	  h->GetYaxis()->SetTitleOffset(1.8);
	  h->LabelsDeflate("v");

	  //LINE COLOR AND WIDTH                                                                                                                                                      
	  Color_t color = CategoryColorMap("DR");
	  h->SetLineWidth(2);
	  h->SetLineColor(color);

	  //TLine* lh = new TLine(0.,0.5,0.03,0.5);
	  //lh->SetLineColor(CategoryColorMap("DR"));
	  //lh->SetLineWidth(3);
	  TLine* lv = new TLine(0.03,0.,0.03,1.);
	  lv->SetLineColor(CategoryColorMap("DR"));
	  lv->SetLineWidth(3);

	  //STATISTICS
	  gStyle->SetOptStat(0);
	  
	  TString png_name = pngName(histo2D_name,"");

	  TCanvas* c1 = new TCanvas("c1","canvas", 1024,768);
	  c1->SetPhi(200);
	  gPad->SetLeftMargin(0.15);
	  gPad->SetRightMargin(0.15);

	  h->Draw("COLZ");
	  //lh->Draw("same");
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
	  histo_SetUp(h1, category1, "");
	  h1->SetTitle(title);
	  h1->GetYaxis()->SetRangeUser(0.01, 3.);
	  histo_SetUp(h2, category2, "");

	  //STATISTICS
	  gStyle->SetOptStat(0);
	  gStyle->SetPaintTextFormat("1.4f");
	  

	  //LEGEND
	  auto legend = new TLegend(0.62,0.8,.89,.89);
		legend->SetBorderSize(0);
	  legend->AddEntry(h1, CategoryLegend(category1) ,"f");
	  legend->AddEntry(h2, CategoryLegend(category2) ,"f");

	  TString png_name = pngName(hReco, category2);
	  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	  gPad->SetLogy();
	  h1->Draw("HIST TEXT0");
	  h2->Draw("HIST TEXT0 SAME");

	  gPad->RedrawAxis();
	  legend->Draw();

	  c1->SaveAs(png_name);

	  input_file->Close();

	  return 0;

	}

	int draw_matching(const TString histo_name, const TString title){

	  TFile* input_file = open_file();
	  TH1F* h = (TH1F*)input_file->Get(histo_name);
	  if ( !h ){
		 std::cout<< "null pointer for histogram named " << histo_name << std::endl;
		 exit(-1);
	  }
	  
	  //TITLE & AXIS
	  //h->SetTitle(title);
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
	  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
	  gPad->SetBottomMargin(0.13); gPad->SetLeftMargin(0.13);
	  h->Draw("HIST");
	  t->Draw("same");
	  c1->SaveAs(png_name);

	  input_file->Close();
	  

	  return 0;
	}


	int draw_lead_pt(const TString h_leading, const TString h_subleading){

	  TFile* input_file = open_file();

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
	  //SETUP
	  histo_SetUp(h_l, "LMu" , "\\ p_T [GeV]");
	  histo_SetUp(h_sl, "sLMu" , "\\ p_T [GeV]");

	  //STATISTICS
																																		 
	  gStyle->SetOptStat(0);

	  //LEGEND

	  auto legend = new TLegend(0.55,0.8,.89,.89);
		legend->SetBorderSize(0);
	  legend -> AddEntry(h_sl, CategoryLegend("sLMu") ,"f");
	  legend -> AddEntry(h_l, CategoryLegend("LMu") ,"f");
	  legend->SetTextAlign(22);

	  TString png_name = pngName(h_leading,"");

	  TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);

	  h_sl->Draw("HIST");
	  h_l->Draw("HIST SAME");
	  legend->Draw();
	  GetCMS(c1);

	  c1->SaveAs(png_name);

	  input_file->Close();

	  return 0;
	}

int draw_NReco(const TString histo_name, const TString& category = "rec"){

	TFile* input_file = open_file();
	TH1F* h = (TH1F*)input_file->Get(histo_name);
	if ( !h ){
		std::cout<< "null pointer for histogram named " << histo_name << std::endl;
		exit(-1);
	}

	//TITLE & AXIS
	TString title = "Reconstructed tracks per event";
	if (category == "OtherP") title = "Number of muons from other particles per event"; 
	// h->SetTitle(title);
	const int nx = 3;
	std::string os_X[nx]   = {"0","1","2"};
	for (int i=1; i<=nx; i++) {
		h->GetXaxis()->SetBinLabel(i,os_X[i-1].c_str());
	}

	histo_SetUp(h, "NReco" , "");
	h->GetXaxis()->SetLabelSize(0.05);
	h->SetMaximum(3.);
	//STATISTICS
	gStyle->SetOptStat(0);

	//LEGEND
	auto legend = new TLegend(0.20,0.80,.40,.85);
	legend->AddEntry(h, CategoryLegend(category) ,"f");
	legend->SetBorderSize(0);

	TString png_name = pngName(histo_name, "");
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetPaintTextFormat(".4f"); 
	h->SetBarOffset(0.5);
	h->SetMarkerSize(1.5);
	gPad->SetLogy();
	gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.11);
	h->Draw("HIST TEXT0");
	legend->Draw();
	//GetCMS(c1);
	c1->SaveAs(png_name);
	input_file->Close();

	return 0;
}

