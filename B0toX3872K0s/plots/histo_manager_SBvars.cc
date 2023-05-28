#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TLine.h>
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TLatex.h>
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>

using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// draw_single_histogram(const TString histo_name)
// draw_two_histograms(const TString histo1,const TString& category1, const TString histo2, const TString& category2, const TString& x_name, bool fill = true, bool log = false){
//


TFile* open_file(){
    TString root_file = "SB_variables.root";
    TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};
  Color["SGN_Dphi"] = kGreen - 6;
  Color["BKG_Dphi"] = kViolet;

  Color["SGN_JPsi"] = kRed;
  Color["BKG_JPsi"] = kOrange;

  Color["SGN_Pi"] = kRed;
  Color["BKG_Pi"] = kBlack;
  Color["BKG_B0_Pi"] = kCyan + 1;

  Color["SGN_Rho"] = kRed - 7;
  Color["BKG_Rho"] = kOrange - 2;
  Color["BKG_B0_Rho"] = kTeal + 2;

  Color["SGN_K0trk"] = kViolet - 1;
  Color["BKG_K0trk"] = kAzure  - 9;

  Color["SGN_K0s"] = kGreen + 1;
  Color["BKG_K0s"] = kOrange + 1;
  Color["BKG_B0_K0s"] = kPink- 4;

  Color["SGN_X3872"] = kBlue - 9;
  Color["BKG_X3872"] = kAzure  - 4;

  Color["SGN_B0"] = kAzure + 1;
  Color["BKG_B0"] = kRed;
  Color["BKG_Rho_B0"] = kOrange;
  Color["BKG_K0s_B0"] = kGreen;
  Color["BKG_RK_B0"] = kMagenta - 3;

  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["JPsi"] = "J/#Psi CANDIDATES ";

  Leg_entry["SGN_Dphi"] = "MC-matching ";
  Leg_entry["BKG_Dphi"] = "FAKE";
  Leg_entry["SGN_JPsi"] = "MC-matching J/#Psi"; 
  Leg_entry["BKG_JPsi"] = "FAKE B_{0}";

  Leg_entry["SGN_Pi"] = "MC-matching #pi [#rho]";
  Leg_entry["BKG_Pi"] = "FAKE #pi [#rho]";
  Leg_entry["BKG_B0_Pi"] = "FAKE B_{0}";

  Leg_entry["SGN_Rho"] = "MC-matching #rho(770)";
  Leg_entry["BKG_Rho"] = "FAKE #rho(770)";
  Leg_entry["BKG_B0_Rho"] = "FAKE B_{0}";

  Leg_entry["SGN_K0trk"] = "MC-matching #pi (K_{0}^{s})";
  Leg_entry["BKG_K0trk"] = "FAKE #pi (K_{0}^{s})";

  Leg_entry["SGN_K0s"] = "MC-matching K_{0}^{s}";
  Leg_entry["BKG_K0s"] = "FAKE K_{0}^{s}";
  Leg_entry["BKG_B0_K0s"] = "FAKE B_{0}";

  Leg_entry["SGN_X3872"] = "MC-matching X(3872)"; 
  Leg_entry["BKG_X3872"] = "FAKE #rho(770)";

  Leg_entry["SGN_B0"] = "MC-matching B_{0}";
  Leg_entry["BKG_B0"] = "FAKE B_{0}";
  Leg_entry["BKG_Rho_B0"] = "FAKE #rho(770)";
  Leg_entry["BKG_K0s_B0"] = "FAKE K_{s}^{0}";
  Leg_entry["BKG_RK_B0"] = "FAKE #rho(770) & K_{s}^{0}";

  return Leg_entry[category];
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

TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./SB_variables/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
}

TString pdfName(const TString& histo1_name, const TString& category2){
  TString pdf_name = "./SB_variables/" + histo1_name ;
  if (category2 != "") pdf_name += "_" + category2;
  
  return pdf_name + ".pdf";
}

void histo_SetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true , bool norm = true){
  //AXIS LABEL                                                                                                                                                                  
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitle("counts");
  histo->GetYaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetTitleSize(0.04);

  if (category == "??"){
    histo->GetXaxis()->SetTitleOffset(1.);
    histo->GetYaxis()->SetTitleOffset(1.);
  }

  //WIDTH & COLOR
  gStyle->SetLineWidth(3);
  histo->SetLineWidth(3);
  histo->SetLineColor(CategoryColorMap(category));
  if (fill)histo->SetFillColorAlpha(CategoryColorMap(category), 0.3);
  //NORMALIZATION
  if(norm) histo->Scale(1./histo->Integral());
}



int draw_single_histogram(const TString& histo_name, const TString& category, const TString& x_name, bool fill = true){
    
    TFile* input_file = open_file();
    TH1F * h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    
    //AXIS & TITLE
    h->SetTitle("");
	histo_SetUp(h, category, x_name, fill);

    
    //STATISTICS
    gStyle->SetOptStat(0);

    TString png_name = pngName(histo_name, "");
    TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
    h->Draw("hIST");
    c1->SaveAs(png_name);

    input_file->Close();
    return 0;

}


int draw_two_histograms(const TString histo1,const TString& category1, const TString histo2, const TString& category2, const TString& x_name, bool fill = true, bool log = false){

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
	 if(log) h1->SetMaximum(3.); 
    //LEGEND
    auto legend = new TLegend(0.62,0.8,.89,.89);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    legend->AddEntry(h2, CategoryLegend(category2) ,"f");
    

    TString png_name = pngName(histo1, category2);
    TCanvas* c1 = new TCanvas("c1","canvas", 1248,1024);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    gPad->RedrawAxis();
	 if(log) gPad->SetLogy();
    legend->Draw();

    c1->SaveAs(png_name);

    input_file->Close();

    return 0;
}

int draw_many_histo(std::vector<TString> histos, std::vector<TString> categories, const TString & x_name, const TString MyPngName ,bool stack = true){
  TFile* input_file = open_file();
  //FILL THE STACK
  TH1* h1 = (TH1*)input_file->Get(histos[0]);
  TString y_name = Form("Entries / %.3f GeV", h1->GetXaxis()->GetBinWidth(1));
  THStack* Stk = new THStack("hStack",";" + x_name + ";" + y_name);
  
  //LEGEND
  UInt_t NH = histos.size();
  auto legend = new TLegend(0.53,0.89 - 0.07*NH,.84,.80);
  legend->SetBorderSize(0);
  
  for (UInt_t i = 0; i < NH; i++){
    TH1* h = (TH1*)input_file->Get(histos[i]);
    histo_SetUp(h, categories[i], "", true, false);
    
    Stk->Add(h);
    legend->AddEntry(h, CategoryLegend(categories[i]) ,"f");
  }
  
  Stk->SetMaximum(1.1*Stk->GetMaximum());

  //DRAW
  TString png_name, pdf_name;
  if(MyPngName == ""){ 
	  png_name = pngName("STK_"+histos[1], "");
	  pdf_name = pdfName("STK_"+histos[1], "");
  }else{ 
	  png_name = pngName(MyPngName,"");
	  pdf_name = pdfName(MyPngName,"");
  }

  TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);

  TString DRAWopt = "HIST";
  if(!stack) DRAWopt = "nostack HIST";
  Stk->Draw(DRAWopt);
  legend->Draw();
  //GetCMS(c1);
  Stk->GetXaxis()->SetTitleSize(0.04);
  Stk->GetYaxis()->SetTitleSize(0.04);
  c1->Update();
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
  THStack * hStack = new THStack("hStack",";"+x_name+";counts");
  hStack->Add(h1);
  hStack->Add(h2);
  hStack->SetMaximum(1.3*hStack->GetMaximum());

  //STATISTICS                                                                                                                                                                
  gStyle->SetOptStat(0);
  
  //LEGEND
  auto legend = new TLegend(0.62,0.8,.89,.89);
  legend->AddEntry(h1, CategoryLegend(category1) ,"f");
  legend->AddEntry(h2, CategoryLegend(category2) ,"f");


  TString png_name = pngName("STK_"+histo1, category2);
  TCanvas* c1 = new TCanvas("c1","canvas", 1248,1024);

  hStack->Draw();
  legend->Draw();

  c1->SaveAs(png_name);

  input_file->Close();

  return 0;
}



int draw_2Dhisto(const TString histo2D_name, const TString xname, const TString yname,const TString category,const int option = 0){

  TFile* input_file = open_file();
  TH2F * h = (TH2F*)input_file->Get(histo2D_name);
  if ( !h ){
    std::cout<< "null pointer for histogram named " << histo2D_name << std::endl;
    exit(-1);
  }

  //AXIS & TITLE
  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();

  //LINE COLOR AND WIDTH                                                                                                                                                      
  Color_t color = CategoryColorMap(category);
  h->SetLineWidth(2);
  h->SetLineColor(color);

  //STATISTICS
  gStyle->SetOptStat(0);
  
  TString png_name = pngName(histo2D_name,"");

  TCanvas* c1 = new TCanvas("c1","canvas", 1520,1024);
  h->Draw("COLZ");
  c1->SaveAs(png_name);
  
  TCanvas* c2 = new TCanvas("c2","canvas", 1248,1024);
  
  TH1* hx = h->ProjectionX();
  histo_SetUp(hx, category, xname);
  hx->Draw("HIST");
  hx->SetTitle("");
  
  c2->SaveAs(pngName("projX_"+histo2D_name, ""));

  TH1* hy = h->ProjectionY();
  histo_SetUp(hy, category, yname);
  hy->SetTitle("");
  hy->Draw("HIST");
  
  c2->SaveAs(pngName("projY_"+histo2D_name, ""));

  input_file->Close();

  return 0;

}//draw_2Dhisto()

int draw_two_2Dhisto(const TString h2D_name1,const TString category1, const TString h2D_name2, const TString category2, const TString xname, const TString yname, const int option = 0){

  TFile* input_file = open_file();
  TH2F * h1 = (TH2F*)input_file->Get(h2D_name1);
  if ( !h1 ){
    std::cout<< "null pointer for histogram named " << h2D_name1 << std::endl;
    exit(-1);
  }
  TH2F * h2 = (TH2F*)input_file->Get(h2D_name2);
  if ( !h2 ){
    std::cout<< "null pointer for histogram named " << h2D_name2 << std::endl;
    exit(-1);
  }

  //AXIS & TITLE
  h1->GetXaxis()->SetTitle(xname);
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle(yname);
  h1->GetYaxis()->CenterTitle();

  //LINE COLOR AND WIDTH                                                                                                                                                      
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h1->SetLineColor(CategoryColorMap(category1));
  h2->SetLineColor(CategoryColorMap(category2));

  //STATISTICS
  gStyle->SetOptStat(0);
  
  TString png_name = pngName(h2D_name1,category2);

  TCanvas* c1 = new TCanvas("c1","canvas", 1520,1024);
  h1->Draw("BOX0");
  h2->Draw("BOX0 SAME");
  c1->SaveAs(png_name);
  
  input_file->Close();

	return 0;
}//draw_two_2Dhisto()




int draw_many_2Dhisto(std::vector<TString> histos, std::vector<TString> categories, const TString & x_name, const TString & y_name,const TString MyPngName ){

  TFile* input_file = open_file();
  
  //LEGEND
  UInt_t NH = histos.size();
  auto legend = new TLegend(0.62,1 - 0.08*NH,.89,.89);
  TString draw_opt = "BOX";

	//DRAW
  TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);
  for (UInt_t i = 0; i < NH; i++){
	  TH2* h = (TH2*)input_file->Get(histos[i]);
	  histo_SetUp(h, categories[i], x_name, true, false);
		h->GetYaxis()->SetTitle(y_name);
		h->SetFillColor(CategoryColorMap(categories[i]));

	  if(i > 0) draw_opt = "BOX same";
	  h->Draw(draw_opt);

	  legend->AddEntry(h, CategoryLegend(categories[i]) ,"f");
  }
  
 
  //SAVE
  TString png_name;
  if(MyPngName == "") png_name = pngName(histos[1], "");
  else png_name = pngName(MyPngName,"");


  legend->Draw();
  c1->SaveAs(png_name);

  input_file->Close();
  return 0;


}//draw_two_2Dhisto()


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


int  makeROCcurve(TString sigHist_name, TString bkgHist_name){

	
  TFile* input_file = open_file();
  TH1F* sigHist = (TH1F*)input_file->Get(sigHist_name);
  if ( !sigHist ){
    std::cout<< "null pointer for histogram named " << sigHist_name << std::endl;
    exit(-1);}

  TH1F* bkgHist = (TH1F*)input_file->Get(bkgHist_name);
  if ( !bkgHist ){
    std::cout<< "null pointer for histogram named " << bkgHist_name << std::endl;
    exit(-1);}

  int nbins = sigHist->GetNbinsX();
  float sig_integral = sigHist->Integral(1,nbins);
  float bkg_integral = bkgHist->Integral(1,nbins);
  std::cout<<" total int  sig: "<<sig_integral<<" bkg: "<<bkg_integral<<std::endl;
  std::vector<float> sigPoints(nbins);
  std::vector<float> bkgPoints(nbins);
  for ( int i = nbins; i > 0; i-- ) {
    float sig_slice_integral = sigHist->Integral(i,nbins);
    float bkg_slice_integral = bkgHist->Integral(i,nbins);
    sigPoints.push_back(sig_slice_integral/sig_integral);
    bkgPoints.push_back(bkg_slice_integral/bkg_integral);

    std::cout<<i<<" "<<sig_slice_integral<<" "<<sig_slice_integral/sig_integral<<" "<<bkg_slice_integral<<" "<<bkg_slice_integral/bkg_integral<<std::endl;
  }

  TGraph *g = new TGraph(sigPoints.size(),&bkgPoints[0],&sigPoints[0]);


  TString png_name = pngName("ROC"+sigHist_name, "");
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 768);
  g->Draw();
  c1->SaveAs(png_name);

  input_file->Close();



  return 0;
}// makeROCcurve()




int Fit_mass_X(TString func = "GAUS"){

	// load HISTO
	TString histo_name = "SGN_X3872_M";
	TFile* input_file = open_file();
	TH1F* h = (TH1F*)input_file->Get(histo_name);
	if ( !h ){
		std::cout<< "null pointer for histogram named " << histo_name << std::endl;
		exit(-1);}


	// optimal FUNC
	TF1 * fitFunc(0); 
	TF1 * fitSubFunc_G1(0); 
	TF1 * fitSubFunc_G2(0); 

	double Mlow = 3.8 , Mhigh = 3.95;
	if(func == "GAUS") {
		fitFunc	= new TF1("fitFunc","gaus(0)+gaus(3)",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
		fitFunc->SetParNames("A_G1", "mean_G1", "Sigma_G1", "A_G2", "mean_G2", "Sigma_G2");
	}
	if(func == "CB"){
		fitFunc = new TF1("fitFunc", "crystalball", Mlow, Mhigh);
	}
	if(func == "GAUS+CB"){
		fitFunc = new TF1("fitFunc", "gaus(0) + crystalball(3)", Mlow, Mhigh);
		fitFunc->SetParNames("A_G", "mean_G", "Sigma_G", "A_CB", "mean_CB", "Sigma_CB", "alpha_CB", "N_CB");
		std::cout << "\n FIT function " << fitFunc->GetExpFormula() << std::endl;
	}

	// FIT options
	ROOT::Fit::DataOptions opt;
	opt.fIntegral = true;

	ROOT::Fit::BinData data(opt);
	ROOT::Fit::FillData(data, h);

	ROOT::Math::WrappedMultiTF1 fitFunction(*fitFunc, fitFunc->GetNdim() );

	ROOT::Fit::Fitter fitter;
	fitter.SetFunction( fitFunction, false);

	if(func == "GAUS"){
		double par0[] = {130., 3.87, 0.05, 50, 3.87, 0.1};
		fitter.Config().SetParamsSettings(fitFunc->GetNumberFreeParameters(), par0);
		fitter.Config().ParSettings(0).SetLimits(50, 200);
		fitter.Config().ParSettings(1).SetLimits(3.870,3.875);
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.1);
		fitter.Config().ParSettings(3).SetLimits(0, 100);
		fitter.Config().ParSettings(4).SetLimits(3.870,3.875);
		fitter.Config().ParSettings(5).SetLimits(0.0, 0.08);
	}
	if (func == "CB"){
		double par0[] = {130., 5.27, 0.02, 1, 3};
		fitter.Config().ParSettings(0).SetLimits(0, 150);
		fitter.Config().ParSettings(1).SetLimits(5.26,5.29);
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.3);
		fitter.Config().ParSettings(3).SetLimits(0., 80.);
		fitter.Config().ParSettings(4).SetLimits(0., 80.);
	}
	if (func == "GAUS+CB"){
		double par0[] = {180., 5.27, 0.02, 150., 5.27, 0.02, 1, 3};
		// ... gauss
		fitter.Config().ParSettings(0).SetLimits(20, 150);   // ampl
		fitter.Config().ParSettings(1).SetLimits(5.27,5.28); // mean 
		//fitter.Config().ParSettings(1).Fix();
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.05); // sigma
		// ... crystal ball
		fitter.Config().ParSettings(3).SetLimits(0, 150);   // ampl
		fitter.Config().ParSettings(4).SetLimits(5.27,5.28); // mean
		//fitter.Config().ParSettings(4).Fix();
		fitter.Config().ParSettings(5).SetLimits(0.0, 0.1);  // sigma
		fitter.Config().ParSettings(6).SetLimits(0., 50.);  // N
		fitter.Config().ParSettings(7).SetLimits(0., 50.);  // alpha		
	}

	// FIT
	fitter.Fit(data);
	ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);
	const double * FitParams = result.GetParams();
	const double * FitParErr = result.GetErrors();

	// FUNCTION SET UP
	fitFunc->SetFitResult(result);
	fitFunc->SetLineColor(kRed);
	if (func == "GAUS"){
		fitSubFunc_G1 = new TF1("fitSubFunc_G1","gaus",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
		fitSubFunc_G1->SetParameters(*FitParams, *(FitParams+1), *(FitParams+2));
		fitSubFunc_G1->SetLineColor(kGreen); fitSubFunc_G1->SetLineWidth(2);
		fitSubFunc_G2 = new TF1("fitSubFunc_G2", "gaus", Mlow, Mhigh);
		fitSubFunc_G2->SetParameters(*(FitParams+3), *(FitParams+4), *(FitParams+5));
		fitSubFunc_G2->SetLineColor(kOrange); fitSubFunc_G2->SetLineWidth(2);
	}

	// SIDEBANDS
	float meanXm = (*(FitParams+1) + *(FitParams+4))*0.5;	
	float wG1 = fitSubFunc_G1->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float wG2 = fitSubFunc_G2->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float SIGMA = 2 * sqrt( (*(FitParams+2) * *(FitParams+2) * wG1*wG1  )  + (*(FitParams+5) * *(FitParams+5) * wG2*wG2));
	std::cout << "\nMEAN " << meanXm << "\tSIGMA " << SIGMA << std::endl;
	float nSIG = 4;
	float MXnearLeft = meanXm - nSIG * SIGMA; 
	float MXnearRight = meanXm + nSIG * SIGMA; 
	std::cout << " X Signal Region : [ " << MXnearLeft << ", " << MXnearRight << "]" << std::endl;

	TLine* NearR = new TLine(MXnearRight ,0., MXnearRight ,20);
	TLine* NearL = new TLine(MXnearLeft  ,0., MXnearLeft  ,20);
	NearR->SetLineColor(kBlack); NearL->SetLineColor(kBlack); 
	NearR->SetLineWidth(3); NearL->SetLineWidth(3) ; 

	double SGNefficiency= h->Integral(h->FindBin(MXnearLeft), h->FindBin(MXnearRight));
	SGNefficiency = SGNefficiency/h->Integral();
	std::cout << "SGN EFFICIENCY: " << SGNefficiency << std::endl;

	//LEGEND
	auto legend = new TLegend(0.11,0.75,.41,.89);
	legend->SetBorderSize(0);
	gStyle->SetLegendTextSize(0.025);

	// DRAW results 
	TString png_name = pngName(histo_name + "_Fit" + func, "");
	TString png_name_log = pngName(histo_name + "_FitLog_" + func, "");
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

	histo_SetUp(h, "SGN_X3872", "", false, false);

	h->Draw(); fitSubFunc_G1->Draw("same"); fitSubFunc_G2->Draw("same");
	fitFunc->Draw("same");
	NearR->Draw("same"); NearL->Draw("same");
	legend->AddEntry(h, CategoryLegend("SGN_X3872") ,"f");
	legend->AddEntry("fitFunc", "FIT function = GAUS + GAUS" ,"l");
	legend->AddEntry("fitSubFunc_G1", "GAUSIAN 1" ,"l");
	legend->AddEntry("fitSubFunc_G2", "GAUSSIAN 2" ,"l");
	legend->Draw();

	gStyle->SetOptStat(0);

	c1->SaveAs(png_name);
	c1->SetLogy();
	c1->SaveAs(png_name_log);
	input_file->Close();
	
	// Write on file
	std::string outFileParams = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/X3872params.txt";	
	std::ofstream outFile;
	outFile.open(outFileParams);
	if(!outFile.is_open()) std::cout << "ERROR cannot open file" << outFileParams << std::endl;
	outFile << "# X(3872) fit parameters ..." << std::endl;	
	outFile << func << std::endl;
	outFile << "MXnearLeft \t" << MXnearLeft << "\n";
	outFile << "MXnearRight \t" << MXnearRight << "\n";
	outFile << "SGN EFFICIENCY: " << SGNefficiency << "\n"; 
	for (unsigned int p = 0; p < result.NFreeParameters(); p++){
		//result.GetParameterName(p)
		outFile << fitter.Config().ParSettings(p).Name() << "\t" << *(FitParams + p) << "\t" << *(FitParErr + p)<< "\n"; 
	}	
	outFile << "Chi2/Ndf" << "\t" << result.Chi2()/result.Ndf() << "\n";
	outFile.close();

	return 0;
}//Fit_mass_X()



int Fit_mass_K0s(TString func = "GAUS"){

	// load HISTO
	TString histo_name = "SGN_K0s_M";
	TFile* input_file = open_file();
	TH1F* h = (TH1F*)input_file->Get(histo_name);
	if ( !h ){
		std::cout<< "null pointer for histogram named " << histo_name << std::endl;
		exit(-1);}


	// optimal FUNC
	TF1 * fitFunc(0); 
	TF1 * fitSubFunc_G(0); 
	TF1 * fitSubFunc_CB(0);
	TF1 * fitSubFunc_G1(0);
	TF1 * fitSubFunc_G2(0);
	double Mlow = 0.44 , Mhigh = 0.6;
	if(func == "GAUS") {
		fitFunc	= new TF1("fitFunc","gaus(0)+gaus(3)",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
		fitFunc->SetParNames("A_G1", "mean_G1", "Sigma_G1", "A_G2", "mean_G2", "Sigma_G2");
	}
	if(func == "CB"){
		fitFunc = new TF1("fitFunc", "crystalball", Mlow, Mhigh);
	}
	if(func == "GAUS+CB"){
		fitFunc = new TF1("fitFunc", "gaus(0) + crystalball(3)", Mlow, Mhigh);
		fitFunc->SetParNames("A_G", "mean_G", "Sigma_G", "A_CB", "mean_CB", "Sigma_CB", "alpha_CB", "N_CB");
		std::cout << "\n FIT function " << fitFunc->GetExpFormula() << std::endl;
	}

	// FIT options
	ROOT::Fit::DataOptions opt;
	opt.fIntegral = true;

	ROOT::Fit::BinData data(opt);
	ROOT::Fit::FillData(data, h);

	ROOT::Math::WrappedMultiTF1 fitFunction(*fitFunc, fitFunc->GetNdim() );

	ROOT::Fit::Fitter fitter;
	fitter.SetFunction( fitFunction, false);

	if(func == "GAUS"){
		double par0[] = {180., 0.5, 0.02, 50, 0.5, 0.02};
		fitter.Config().SetParamsSettings(fitFunc->GetNumberFreeParameters(), par0);
		fitter.Config().ParSettings(0).SetLimits(150, 200);
		fitter.Config().ParSettings(1).SetLimits(0.49, 0.51);
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.01);
		fitter.Config().ParSettings(3).SetLimits(0, 50);
		fitter.Config().ParSettings(4).SetLimits(0.49, 0.51);
		fitter.Config().ParSettings(5).SetLimits(0.0, 0.08);
	}
	if (func == "CB"){
		double par0[] = {130., 5.27, 0.02, 1, 3};
		fitter.Config().ParSettings(0).SetLimits(0, 150);
		fitter.Config().ParSettings(1).SetLimits(5.26,5.29);
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.3);
		fitter.Config().ParSettings(3).SetLimits(0., 80.);
		fitter.Config().ParSettings(4).SetLimits(0., 80.);
	}
	if (func == "GAUS+CB"){
		double par0[] = {170, 0.5, 0.02, 20., 0.5, 0.02, 1, 3};
		// ... gauss
		fitter.Config().ParSettings(0).SetLimits(20, 180);   // ampl
		fitter.Config().ParSettings(1).SetLimits(0.49,0.51); // mean 
		//fitter.Config().ParSettings(1).Fix();
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.05); // sigma
		// ... crystal ball
		fitter.Config().ParSettings(3).SetLimits(0, 180);   // ampl
		fitter.Config().ParSettings(4).SetLimits(0.49,0.51); // mean
		//fitter.Config().ParSettings(4).Fix();
		fitter.Config().ParSettings(5).SetLimits(0.0, 0.1);  // sigma
		fitter.Config().ParSettings(6).SetLimits(0., 10.);  // N
		fitter.Config().ParSettings(7).SetLimits(0., 10.);  // alpha
	}

	// FIT
	fitter.Fit(data);
	//fitter.FitFCN(3, *fitFunc , 0, data.Size(), true);
	ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);
	const double * FitParams = result.GetParams();
	const double * FitParErr = result.GetErrors();

	// FUNCTION SET UP
	fitFunc->SetFitResult(result);
	fitFunc->SetLineColor(kRed);
	if (func == "GAUS+CB"){
		fitSubFunc_G = new TF1("fitSubFunc_G","gaus",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
		fitSubFunc_G->SetParameters(*FitParams, *(FitParams+1), *(FitParams+2));
		fitSubFunc_G->SetLineColor(kBlue); fitSubFunc_G->SetLineWidth(2);
		fitSubFunc_CB = new TF1("fitiSubFunc_CB", "crystalball", Mlow, Mhigh);
		fitSubFunc_CB->SetParameters(*(FitParams+3), *(FitParams+4), *(FitParams+5), *(FitParams+6), *(FitParams+7));
		fitSubFunc_CB->SetLineColor(kOrange); fitSubFunc_CB->SetLineWidth(2);
	}
	if (func == "GAUS"){
		fitSubFunc_G1 = new TF1("fitSubFunc_G1","gaus",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
		fitSubFunc_G1->SetParameters(*FitParams, *(FitParams+1), *(FitParams+2));
		fitSubFunc_G1->SetLineWidth(2); fitSubFunc_G1->SetLineColor(kBlue);
		fitSubFunc_G2 = new TF1("fitSubFunc_G2","gaus",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
		fitSubFunc_G2->SetParameters(*(FitParams+3), *(FitParams+4), *(FitParams+5));
		fitSubFunc_G2->SetLineWidth(2); fitSubFunc_G2->SetLineColor(kOrange);
	}
	// SIGNAL REGION  
	float meanKm = (*(FitParams+1) + *(FitParams+4))*0.5;	
	float wG1 = fitSubFunc_G1->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float wG2 = fitSubFunc_G2->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float SIGMA = 2 * sqrt( (*(FitParams+2) * *(FitParams+2) * wG1*wG1  )  + (*(FitParams+5) * *(FitParams+5) * wG2*wG2));
	std::cout << "\nMEAN " << meanKm << "\tSIGMA " << SIGMA << std::endl;
	float nSIG = 4;
	float MK0s_nearRight = meanKm + nSIG * SIGMA; 
	float MK0s_nearLeft  = meanKm - nSIG * SIGMA; 
	std::cout << " K0s Signal Region : [ " << MK0s_nearLeft << ", " << MK0s_nearRight << "]" << std::endl;

	TLine* NearR = new TLine(MK0s_nearRight, 0., MK0s_nearRight, 20);
	TLine* NearL = new TLine(MK0s_nearLeft,  0., MK0s_nearLeft, 20);
	NearR->SetLineColor(kBlack); NearL->SetLineColor(kBlack); 
	NearR->SetLineWidth(3); NearL->SetLineWidth(3); 

	double SGNefficiency= h->Integral(h->FindBin(MK0s_nearLeft), h->FindBin(MK0s_nearRight));
	SGNefficiency /= h->Integral();
	std::cout << "SGN EFFICIENCY: " << SGNefficiency << std::endl;

	//LEGEND
	auto legend = new TLegend(0.11,0.75,.41,.89);
	legend->SetBorderSize(0);
	gStyle->SetLegendTextSize(0.025);

	// DRAW results 
	TString png_name = pngName(histo_name + "_Fit" + func, "");
	TString png_name_log = pngName(histo_name + "_FitLog" + func, "");
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

	histo_SetUp(h, "SGN_K0s", "", false, false);
	h->Draw();
	if(func == "GAUS+CB"){
		fitSubFunc_G->Draw("same");
		fitSubFunc_CB->Draw("same");
	}
	if (func == "GAUS"){
		fitSubFunc_G1->Draw("same");
		fitSubFunc_G2->Draw("same");
	}
	fitFunc->Draw("same");
	NearR->Draw("same"); NearL->Draw("same");
	legend->AddEntry(h, CategoryLegend("SGN_K0s") ,"f");
	legend->AddEntry("fitFunc", "FIT function = GAUS + GAUS" ,"l");
	legend->AddEntry("fitSubFunc_G1", "GAUSIAN 1" ,"l");
	legend->AddEntry("fitSubFunc_G2", "GAUSSIAN 2" ,"l");
	legend->Draw();

	gStyle->SetOptStat(0);

	c1->SaveAs(png_name);
	c1->SetLogy();
	c1->SaveAs(png_name_log);
	input_file->Close();

	// Write on file
	std::string outFileParams = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/K0sparams.txt";	
	std::ofstream outFile;
	outFile.open(outFileParams);
	if(!outFile.is_open()) std::cout << "ERROR cannot open file" << outFileParams << std::endl;
	outFile << "# K0s fit parameters ..." << std::endl;	
	outFile << func << std::endl;
	outFile << "MK0s_nearLeft \t" << MK0s_nearLeft << "\n";
	outFile << "MK0s_nearRight \t" << MK0s_nearRight << "\n";
	outFile << "SGN EFFICIENCY: " << SGNefficiency << "\n"; 
	for (unsigned int p = 0; p < result.NFreeParameters(); p++){
		//result.GetParameterName(p)
		outFile << fitter.Config().ParSettings(p).Name() << "\t" << *(FitParams + p) << "\t" << *(FitParErr + p) <<"\n"; 
	}	
	outFile << "Chi2/Ndf" << "\t" << result.Chi2()/result.Ndf() << "\n";
	outFile.close();

	return 0;
}//Fit_mass_K0s()




int Fit_mass_B0(std::string func = "GAUS+CB"){

	// load HISTO
	TString histo_name = "SGN_B0_M";
	TFile* input_file = open_file();
	TH1F* h = (TH1F*)input_file->Get(histo_name);
	if ( !h ){
		std::cout<< "null pointer for histogram named " << histo_name << std::endl;
		exit(-1);}


	// optimal FUNC
	TF1 * fitFunc(0); 
	TF1 * fitSubFunc_G(0); 
	TF1 * fitSubFunc_CB(0); 

	double Mlow = 5.1 , Mhigh = 5.4;
	fitFunc = new TF1("fitFunc", "gaus(0) + crystalball(3)", Mlow, Mhigh);
	fitFunc->SetParNames("A_G", "mean_G", "Sigma_G", "A_CB", "mean_CB", "Sigma_CB", "alpha_CB", "N_CB");
	std::cout << "\n FIT function " << fitFunc->GetExpFormula() << std::endl;

	// FIT options
	ROOT::Fit::DataOptions opt;
	opt.fIntegral = true;

	ROOT::Fit::BinData data(opt);
	ROOT::Fit::FillData(data, h);

	ROOT::Math::WrappedMultiTF1 fitFunction(*fitFunc, fitFunc->GetNdim() );

	ROOT::Fit::Fitter fitter;
	fitter.SetFunction( fitFunction, false);

	double par0[] = {130., 5.27, 0.02, 130., 5.27, 0.02, 1, 3};
	// ... gauss
	fitter.Config().ParSettings(0).SetLimits(20, 150);   // ampl
	fitter.Config().ParSettings(1).SetLimits(5.27,5.28); // mean 
	//fitter.Config().ParSettings(1).Fix();
	fitter.Config().ParSettings(2).SetLimits(0.0, 0.05); // sigma
	// ... crystal ball
	fitter.Config().ParSettings(3).SetLimits(0, 100);   // ampl
	fitter.Config().ParSettings(4).SetLimits(5.27,5.28); // mean
	//fitter.Config().ParSettings(4).Fix();
	fitter.Config().ParSettings(5).SetLimits(0.0, 0.1);  // sigma
	fitter.Config().ParSettings(6).SetLimits(0., 50.);  // N
	fitter.Config().ParSettings(7).SetLimits(0., 50.);  // alpha		

	// FIT
	fitter.Fit(data);
	//fitter.FitFCN(3, *fitFunc , 0, data.Size(), true);
	ROOT::Fit::FitResult result = fitter.Result();
   result.Print(std::cout);

	const double * FitParams = result.GetParams();
	const double * FitParErr = result.GetErrors();

	// FUNCTION SET UP
	fitFunc->SetFitResult(result);
	fitFunc->SetLineColor(kRed);
	fitSubFunc_G = new TF1("fitSubFunc_G","gaus",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
	fitSubFunc_G->SetParameters(*FitParams, *(FitParams+1), *(FitParams+2));
	fitSubFunc_G->SetLineColor(kGreen); fitSubFunc_G->SetLineWidth(2);
	fitSubFunc_CB = new TF1("fitSubFunc_CB", "crystalball", Mlow, Mhigh);
	fitSubFunc_CB->SetParameters(*(FitParams+3), *(FitParams+4), *(FitParams+5), *(FitParams+6), *(FitParams+7));
	fitSubFunc_CB->SetLineColor(kOrange); fitSubFunc_CB->SetLineWidth(2);

	// SIDEBANDS
	float meanBm = (*(FitParams+1) + *(FitParams+4))*0.5;	
	float wG = fitSubFunc_G->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float wCB= fitSubFunc_CB->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float SIGMA = 2 * sqrt( (*(FitParams+2) * *(FitParams+2) * wG*wG  )  + (*(FitParams+5) * *(FitParams+5) * wCB*wCB));
	std::cout << "SIGMA " << SIGMA << std::endl;

	float nSIGmin = 4., nSIGmax = 8;
	float MB_nearLeft = meanBm - nSIGmin * SIGMA;	
	float MB_nearRight = meanBm + nSIGmin * SIGMA;	
	float MB_farRight = meanBm + nSIGmax * SIGMA; 
	float MB_farLeft = meanBm - nSIGmax * SIGMA; 

	TLine* NearR = new TLine(MB_nearRight, 0., MB_nearRight,20);
	TLine* NearL = new TLine(MB_nearLeft,  0., MB_nearLeft, 20);
	TLine* FarR  = new TLine(MB_farRight,  0., MB_farRight, 20);
	TLine* FarL  = new TLine(MB_farLeft,   0., MB_farLeft,  20);
	std::cout << "B0 interval [" << MB_farLeft << "," << MB_nearLeft << "]" << " + [" << MB_nearRight << "," << MB_farRight << "]"  <<std::endl;

	NearR->SetLineColor(kBlack); NearL->SetLineColor(kBlack); FarR->SetLineColor(kBlack); FarL->SetLineColor(kBlack);
	NearR->SetLineWidth(3); NearL->SetLineWidth(3) , FarR->SetLineWidth(3) , FarL->SetLineWidth(3) ;

	double SGNcontamination = h->Integral(h->FindBin(MB_farLeft), h->FindBin(MB_nearLeft)) + h->Integral(h->FindBin(MB_nearRight), h->FindBin(MB_farRight));
	SGNcontamination /= h->Integral();
	std::cout << "SGN CONTAMINATION : " << SGNcontamination << std::endl;

	//LEGEND
	auto legend = new TLegend(0.11,0.75,.41,.89);
	legend->SetBorderSize(0);

	// DRAW results 
	TString png_name = pngName(histo_name + "_Fit" + func, "");
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	TString png_name_log = pngName(histo_name + "_FitLOG" + func, "");

	histo_SetUp(h, "SGN_B0", "", false, false);
	h->GetListOfFunctions()->Add(fitSubFunc_G);
	h->GetListOfFunctions()->Add(fitSubFunc_CB);
	h->GetListOfFunctions()->Add(fitFunc);
	
	h->Draw(); fitFunc->Draw("same"); fitSubFunc_G->Draw("same"); fitSubFunc_CB->Draw("same");
	NearR->Draw("same"); NearL->Draw("same");
	FarR->Draw("same"); FarL->Draw("same");
	legend->AddEntry(h, CategoryLegend("SGN_B0") ,"f");
	legend->AddEntry("fitFunc", "FIT function = GAUS + CB" ,"l");
	legend->AddEntry("fitSubFunc_G", "GAUSIAN" ,"l");
	legend->AddEntry("fitSubFunc_CB", "CRYSTALBALL" ,"l");
	legend->Draw();

	gStyle->SetOptStat(0);

	c1->SaveAs(png_name);
	c1->SetLogy();
	c1->SaveAs(png_name_log);

	input_file->Close();

	std::string outFileParams = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/B0params.txt";	
	std::ofstream outFile;
	outFile.open(outFileParams);
	if(!outFile.is_open()) std::cout << "ERROR cannot open file" << outFileParams << std::endl;
	outFile << "# B0 fit parameters ..." << std::endl;	
	outFile << func << std::endl;
	outFile << "MB_nearLeft" << "\t" << MB_nearLeft << "\n";
	outFile << "MB_farLeft" << "\t" << MB_farLeft << "\n";
	outFile << "MB_nearRight" << "\t" << MB_nearRight << "\n";
	outFile << "MB_farRight" << "\t" << MB_farRight << "\n";
	outFile << "SGN CONTAMINATION : " << SGNcontamination << "\n"; 
	outFile << "# PAR-NAME \t # VALUE \t # ERROR" << "\n";
	for (unsigned int p = 0; p < result.NFreeParameters(); p++){
		//result.GetParameterName(p)
		outFile << fitter.Config().ParSettings(p).Name() << "\t" << *(FitParams + p) << "\t" << *(FitParErr + p) << "\n"; 
	}	
	outFile << "Chi2/Ndf" << "\t" << result.Chi2()/result.Ndf() << "\n";
	outFile.close();

	return 0;
}//Fit_mass_B0()
