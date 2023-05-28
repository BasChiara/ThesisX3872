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
#include <TTree.h>
#include <TLine.h>
#include <TText.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TLatex.h>

// RooFit
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooPlot.h"
using namespace RooFit;

// ROOT::Fit
#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>

using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html


TFile* open_file(TString file = ""){
	TString root_file = "TriggerSelection.root";
	if (file == "DAT") root_file = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/data/merged/SGNPsi2S_MergData17.root";	
	if (file == "MVA") root_file = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAcutApp.root";
	TFile* input_file = new TFile(root_file);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< root_file << std::endl;
       exit(-1);
    }

    return input_file;
}

Color_t CategoryColorMap(const TString& category){

  std::map <TString , Color_t> Color{};

  Color["SGN_JPsi"] = kRed;
  Color["BKG_JPsi"] = kOrange;

  Color["SGN_Rho"] = kRed - 7;
  Color["BKG_Rho"] = kOrange - 2;
  Color["BKG_B0_Rho"] = kTeal + 2;

  Color["SGN_K0s"] = kGreen + 1;
  Color["BKG_K0s"] = kOrange + 1;
  Color["BKG_B0_K0s"] = kPink- 4;

  Color["SGN_Psi2S"] = kBlue - 9;
  Color["BKG_Psi2S"] = kAzure  - 4;

  Color["SGN_B0"] = kAzure + 1;
  Color["BKG_B0"] = kRed;
  Color["BKG_B0_pstcut"] = kRed    - 9;

  return Color[category];
}

TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["JPsi"] = "J\\Psi \\ CANDIDATES ";

  Leg_entry["SGN_JPsi"] = "SIGNAL J#Psi"; 
  Leg_entry["BKG_JPsi"] = "B_0\\ BAKGROUND";

  Leg_entry["SGN_Rho"] = "SIGNAL\\ \\rho (770)";
  Leg_entry["BKG_Rho"] = "BACKGROUND\\ \\rho (770)";

  Leg_entry["SGN_K0s"] = "SIGNAL \\ K_0^s";
  Leg_entry["BKG_K0s"] = "BACKGROUND";

  Leg_entry["SGN_Psi2S"] = "SIGNAL #Psi(2S)"; 
  Leg_entry["BKG_Psi2S"] = "J\\Psi - BKG";

  Leg_entry["SGN_B0"] = "B_0\\ SIGNAL";
  Leg_entry["BKG_B0"] = "DATA (PRE-CUT)";
  Leg_entry["BKG_B0_pstcut"] = "DATA (POST-CUT)";

  return Leg_entry[category];
}


TString pngName(const TString& histo1_name, const TString& category2){
  TString png_name = "./FitResults/" + histo1_name ;
  if (category2 != "") png_name += "_" + category2;
  
  return png_name + ".png";
}

TString pdfName(const TString& histo1_name, const TString& category2){
  TString pdf_name = "./FitResults/" + histo1_name ;
  if (category2 != "") pdf_name += "_" + category2;
  
  return pdf_name + ".pdf";
}

void histo_SetUp(TH1* histo, const TString& category, const TString& x_name,bool fill = true , bool norm = true){
  //AXIS LABEL                                                                                                                                                                  
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.035);
  histo->GetXaxis()->SetLabelSize(0.03);
  TString y_name = Form("%s %4.3f %s", "dN/", histo->GetXaxis()->GetBinWidth(1),"GeV");
  histo->GetYaxis()->SetTitle(y_name);
  histo->GetYaxis()->SetTitleSize(0.035);
  histo->GetYaxis()->SetLabelSize(0.03);

  if (category == "??"){
    histo->GetXaxis()->SetTitleOffset(1.);
    histo->GetYaxis()->SetTitleOffset(1.);
  }

  //WIDTH & COLOR
  gStyle->SetLineWidth(3);
  histo->SetLineWidth(4);
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



int draw_many_histo(std::vector<TString> histos, std::vector<TString> categories, const TString & x_name, const TString MyPngName ,bool stack = true){
  TFile* input_file = open_file();
  //FILL THE STACK
  THStack* Stk = new THStack("hStack",";"+x_name+";counts");;
  
  //LEGEND
  UInt_t NH = histos.size();
  auto legend = new TLegend(0.62,1 - 0.08*NH,.89,.89);
  
  for (UInt_t i = 0; i < NH; i++){
    TH1* h = (TH1*)input_file->Get(histos[i]);
    histo_SetUp(h, categories[i], "", true, false);
    
    Stk->Add(h);
    legend->AddEntry(h, CategoryLegend(categories[i]) ,"f");
  }
  
  Stk->SetMaximum(1.1*Stk->GetMaximum());
 
  //DRAW
  TString png_name;
  if(MyPngName == "") png_name = pngName("STK_"+histos[1], "");
  else png_name = pngName(MyPngName,"");

  TCanvas* c1 = new TCanvas("c1","canvas", 1248, 1024);

  TString DRAWopt = "HIST";
  if(!stack) DRAWopt = "nostack HIST";
  Stk->Draw(DRAWopt);
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
  TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);

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


int Fit_sidebandsB0(TString model){

	
	// LOAD TREE 
	TFile* input_file = open_file("DAT");
	TTree* TreeBKG = (TTree*)input_file->Get("B0_Psi2Ssignal");
	
	// B0 mass RANGE
	RooRealVar B0m("M_B0", "B0 mass sidebands", 4.9, 5.6);

	B0m.setRange("SB1", 4.9, 5.17944); // SIDEBAND 1
	B0m.setRange("SB2", 5.38, 5.6); // SIDEBAND 2
	B0m.setRange("FULLregion", 4.9, 5.6); // FULLregion RANGE 
	B0m.setRange("SGNregion", 5.17944, 5.37937); // SGNregion RANGE 

	TString CutOpt = "((M_B0 < 5.17944) && (M_B0 > 4.9)) || ((M_B0 > 5.38) && (M_B0 < 5.6))";
	RooDataSet B0_mass ("B0_mass","B0_mass", B0m, Import(*TreeBKG), Cut(CutOpt)); 
	// IV polynomial (a0 + a1 * x + a2 *x^2 + a3* x^3 + a4 * x^4)
	// Parameters with starting value, lower bound, upper bound:
	RooRealVar a0("a0", "",  10,  0.,  100000);
	RooRealVar a1("a1", "",  10,  -100000,  0.);
	RooRealVar a2("a2", "",  10,  0.,  100000);
	RooRealVar a3("a3", "",  10,  -100000,  0.);
	RooRealVar a4("a4", "",  10, 0., 100000);

	RooPolynomial poly2("poly2", "poly2", B0m, RooArgList(a0, a1, a2));
	RooPolynomial poly4("poly4", "poly4", B0m, RooArgList(a0, a1, a2, a3, a4));
	RooProdPdf BKGpoly("BKGpoly", "BKG polynomial",poly4);

	// Constant 
	RooRealVar tau("tau", "", -8.0, -20.0 , -5.);
	RooExponential expo("expo", "", B0m, tau);
	RooProdPdf BKGdesc("BKGdesc", "BKGdesc", expo); 
	RooPolynomial Const("Const", "", B0m);
	RooProdPdf BKGbase("BKGbase", "BKGbase", Const); 
	RooRealVar b0("b0", "",  0.5,  -1.,  1.);
	RooRealVar b1("b1", "",  0.2,  -1.,  1.);
	RooRealVar b2("b2", "",  0.5,  -1.,  1.);
	RooRealVar b3("b3", "",  0.2,  -1.,  1.);
	RooChebychev Cheb("Cheb", "", B0m, RooArgSet(b0, b1, b2, b3) );
	RooPolynomial Line("Line", "", B0m, RooArgList(b0, b1));
	RooProdPdf BKGline("BKGbase", "BKGbase", Line); 

	// Create my PDF
	RooRealVar A("A", "", 15.0, 10.0 , 30.);
	RooRealVar Flex("Flex", "", 4., 3. , 5.3);
	RooRealVar C("C", "", 0.5, 0 , 1);
	
	RooGenericPdf Fermi("Fermi", "", "1./(1. + exp(( @0 - @1)*@2)) + @3", RooArgList(B0m, Flex, A, C));

	// New PDF
	RooRealVar Slope("Slope", "", -10.0, -30.0 , -5.);
	RooRealVar Min("Min", "", 4.5, 4. , 5.);
	RooRealVar EXP("EXP", "", 4.);
	
	RooGenericPdf myPDF("myPDF", "", "pow((@0 - @1), @4) * exp(( @0 - @1)*@2) + @3", RooArgList(B0m, Min,  Slope, C, EXP));

	// Buld PDF
	RooRealVar bkg_yield("nBKG", "", 4000, 3000, 100000);
	RooRealVar f("f", "", 0., .5);
	//RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Cheb, bkg_yield);
	//RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Fermi, bkg_yield);
	RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", myPDF, bkg_yield);
	//RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", RooArgList(expo, Const), f);

	//RooFitResult *ResBKGmodelSB2 = BKGmodel.fitTo(B0_mass, Range("SB2"), Save());
	//RooFitResult *ResBKGmodelSB1 = BKGmodel.fitTo(B0_mass, Range("SB1"), Save());
	RooFitResult *ResBKGmodel= BKGmodel.fitTo(B0_mass, Range("SB1,SB2"), Save());
	//ResBKGmodelSB2->Print();
	//ResBKGmodelSB1->Print();
	ResBKGmodel->Print("v");

	RooAbsReal* IntSreg = BKGmodel.createIntegral(B0m, NormSet(B0m), Range("SGNregion"));
	double Ib = IntSreg->getVal();
	//int Bb = h_BKG_B0->GetEntries();
	//double Bs = Bb * (Ib)/(1. - Ib); 
	//RooArgSet nset(B0m);

	// LEGEND
	auto legend = new TLegend(0.50, 0.80,.89,.89);
	legend->SetTextSize(0.025);
	legend->SetBorderSize(0);
	

	//DRAW
	TCanvas* c1 = new TCanvas("c1","canvas", 1480, 1024);
	RooPlot* mesframe = B0m.frame(Bins(70));

	B0_mass.plotOn(mesframe);
//	BKGmodel.plotOn(mesframe, Components(Const), LineColor(kRed) ,LineStyle(kDashed), Range("FULLregion"), RooFit::NormRange("SB1,SB2"));
	BKGmodel.plotOn(mesframe, Range("FULLregion"), RooFit::NormRange("SB1,SB2"));
	//BKGmodel.plotOn(mesframe, LineColor(kRed) ,LineStyle(kDashed), Range("SB2"), RooFit::NormRange("SB1,SB2"));
	mesframe->Draw("SAME");
	cout << " ---> Chi^2 = " << mesframe->chiSquare() << endl;
	//cout << " NORMALIZATION  = " << BKGmodel.createIntegral(B0m, NormSet(B0m))->getVal() << endl;


	//legend->AddEntry(h_BKG_B0   , CategoryLegend("BKG_B0") , "f");
	//legend->AddEntry(h_BKG_B0cut, CategoryLegend("BKG_B0_pstcut"), "f");
	legend->AddEntry("BKGmodel", "FIT : "+ model  + "\t \\Chi^2 = " + Form("%4.3f",mesframe->chiSquare() ), "l");
	legend->Draw();

	gStyle->SetOptStat(0);


	TString png_name = pngName("FitB0sidebands"+model, "");
	c1->SaveAs(png_name);
	
	input_file->Close();

	return 0;

}



int Fit_mass_X(TString func = "GAUS+CB"){

	// load HISTO
	TString histo_name = "SGN_Psi2S_M";
	TFile* input_file = open_file();
	TH1F* h = (TH1F*)input_file->Get(histo_name);
	if ( !h ){
		std::cout<< "null pointer for histogram named " << histo_name << std::endl;
		exit(-1);}


	// optimal FUNC
	TF1 * fitFunc(0); 
	TF1 * fitSubFunc_G1(0); 
	TF1 * fitSubFunc_CB(0); 

	double Mlow = 3.63 , Mhigh = 3.75;
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
		fitter.Config().ParSettings(0).SetLimits(100, 400);
		fitter.Config().ParSettings(1).SetLimits(3.680,3.690);
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.05);
		fitter.Config().ParSettings(3).SetLimits(0, 100);
		fitter.Config().ParSettings(4).SetLimits(3.680,3.690);
		fitter.Config().ParSettings(5).SetLimits(0.01, 1.);
	}
	if (func == "CB"){
		double par0[] = {130., 5.27, 0.02, 1, 3};
		fitter.Config().ParSettings(0).SetLimits(0, 500);
		fitter.Config().ParSettings(1).SetLimits(3.680, 3.690);
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.1);
		fitter.Config().ParSettings(3).SetLimits(0., 5.);
		fitter.Config().ParSettings(4).SetLimits(0., 100.);
	}
	if (func == "GAUS+CB"){
		double par0[] = {180., 5.27, 0.02, 150., 5.27, 0.02, 1, 3};
		// ... gauss
		fitter.Config().ParSettings(0).SetLimits(20, 300);   // ampl
		fitter.Config().ParSettings(1).SetLimits(3.68, 3.69); // mean 
		//fitter.Config().ParSettings(1).Fix();
		fitter.Config().ParSettings(2).SetLimits(0.0, 0.05); // sigma
		// ... crystal ball
		fitter.Config().ParSettings(3).SetLimits(0, 100);   // ampl
		fitter.Config().ParSettings(4).SetLimits(3.680,3.690); // mean
		//fitter.Config().ParSettings(4).Fix();
		fitter.Config().ParSettings(5).SetLimits(0.0, 0.1);  // sigma
		fitter.Config().ParSettings(6).SetLimits(0., 5.);  // N
		fitter.Config().ParSettings(7).SetLimits(0., 1.);  // alpha		
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
	fitFunc->SetLineWidth(3);
	if (func == "GAUS+CB"){
		fitSubFunc_G1 = new TF1("fitSubFunc_G1","gaus",Mlow, Mhigh); // [0] norm [1] mean [2] sigma
		fitSubFunc_G1->SetParameters(*FitParams, *(FitParams+1), *(FitParams+2));
		fitSubFunc_G1->SetLineColor(kGreen); fitSubFunc_G1->SetLineWidth(3);
		fitSubFunc_CB = new TF1("fitSubFunc_CB", "gaus", Mlow, Mhigh);
		fitSubFunc_CB->SetParameters(*(FitParams+3), *(FitParams+4), *(FitParams+5), *(FitParams+6), *(FitParams+7));
		fitSubFunc_CB->SetLineColor(kOrange); fitSubFunc_CB->SetLineWidth(3);
	}

	// SIDEBANDS
	float meanXm = (*(FitParams+1) + *(FitParams+4))*0.5;	
	float wG1 = fitSubFunc_G1->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float wCB = fitSubFunc_CB->Integral(Mlow,Mhigh)/fitFunc->Integral(Mlow,Mhigh);
	float SIGMA = 2 * sqrt( (*(FitParams+2) * *(FitParams+2) * wG1*wG1  )  + (*(FitParams+5) * *(FitParams+5) * wCB*wCB));
	float SIGMA_ERR = 4./SIGMA * sqrt( (*(FitParams+2) * wG1 )*(*(FitParams+2) * wG1 ) * (wG1  * *(FitParErr+2))*(wG1  * *(FitParErr+2))  + 
												  (*(FitParams+5) * wCB)*(*(FitParams+5) * wCB) * (wCB * *(FitParErr+5))*(wCB * *(FitParErr+5))); 
	std::cout << "\nMEAN " << meanXm << "\tSIGMA " << SIGMA << " +/- "  << SIGMA_ERR << std::endl;
	float nSIG = 4;
	float MXnearLeft = meanXm - nSIG * SIGMA; 
	float MXnearRight = meanXm + nSIG * SIGMA; 
	std::cout << " Psi(2S) Signal Region : [ " << MXnearLeft << ", " << MXnearRight << "]" << std::endl;

	TLine* NearR = new TLine(MXnearRight ,0., MXnearRight ,20);
	TLine* NearL = new TLine(MXnearLeft  ,0., MXnearLeft  ,20);
	NearR->SetLineColor(kBlack); NearL->SetLineColor(kBlack); 
	NearR->SetLineWidth(3); NearL->SetLineWidth(3) ; 

	double SGNefficiency= h->Integral(h->FindBin(MXnearLeft), h->FindBin(MXnearRight));
	SGNefficiency = SGNefficiency/h->Integral();
	std::cout << "SGN EFFICIENCY: " << SGNefficiency << std::endl;

	//LEGEND
	auto legend = new TLegend(0.18,0.70,.45,.89);
	legend->SetBorderSize(0);
	gStyle->SetLegendTextSize(0.03);

	// DRAW results 
	TString png_name = pngName(histo_name + "_Fit" + func, "");
	TString png_name_log = pngName(histo_name + "_FitLog_" + func, "");
	TString pdf_name = pdfName(histo_name + "_Fit" + func, "");
	TString pdf_name_log = pdfName(histo_name + "_FitLog_" + func, "");
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);

	histo_SetUp(h, "SGN_Psi2S", "M(J/#Psi #pi^{+}#pi^{-})", false, false);
	h->SetTitle("");

	h->Draw(); 
	fitSubFunc_G1->Draw("same"); fitSubFunc_CB->Draw("same");
	fitFunc->Draw("same");
	NearR->Draw("same"); NearL->Draw("same");
	legend->AddEntry(h, CategoryLegend("SGN_Psi2S") ,"f");
	legend->AddEntry("fitFunc", "FIT " ,"l");
	legend->AddEntry("fitSubFunc_G1", "GAUSIAN" ,"l");
	legend->AddEntry("fitSubFunc_CB", "CRYSTAL BALL" ,"l");
	legend->Draw();

	gStyle->SetOptStat(0);

	c1->SaveAs(png_name);
	c1->SaveAs(pdf_name);
	c1->SetLogy();
	c1->SaveAs(png_name_log);
	c1->SaveAs(pdf_name_log);
	input_file->Close();
	
	//Write on file
	std::string outFileParams = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toPsi2SK0s/results/SGNfit/Psi2Sparams.txt";	
	std::ofstream outFile;
	outFile.open(outFileParams);
	if(!outFile.is_open()) std::cout << "ERROR cannot open file" << outFileParams << std::endl;
	outFile << "# Psi(2S) fit parameters ..." << std::endl;	
	outFile << func << std::endl;
	outFile << "MPsi2SnearLeft \t" << MXnearLeft << "\n";
	outFile << "MPsi2SnearRight \t" << MXnearRight << "\n";
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

	// SIGNAL REGION
	float nSIGs = 3.;
	float MB_sgnL = meanBm - nSIGs *SIGMA;
	float MB_sgnR = meanBm + nSIGs *SIGMA;
	std::cout << "SGN INTERVAL +/-" << nSIGs << " SIGMA [" << MB_sgnL << " ," << MB_sgnR << "]" <<std::endl;

	//LEGEND
	auto legend = new TLegend(0.11,0.75,.41,.89);
	legend->SetBorderSize(0);

	// DRAW results 
	//TString png_name = pngName(histo_name + "_Fit" + func, "");
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	//TString png_name_log = pngName(histo_name + "_FitLOG" + func, "");

	histo_SetUp(h, "SGN_B0", "",  false, false);
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

	//c1->SaveAs(png_name);
	c1->SetLogy();
	//c1->SaveAs(png_name_log);

	input_file->Close();

	std::string outFileParams = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/B0params.txt";	
	std::ofstream outFile;
	outFile.open(outFileParams);
	if(!outFile.is_open()) std::cout << "ERROR cannot open file" << outFileParams << std::endl;
	outFile << "# B0 fit parameters ..." << std::endl;	
	outFile << func << std::endl;
	outFile << "MB_nearLeft"  << "\t" << MB_nearLeft  << "\n";
	outFile << "MB_farLeft"   << "\t" << MB_farLeft   << "\n";
	outFile << "MB_nearRight" << "\t" << MB_nearRight << "\n";
	outFile << "MB_farRight"  << "\t" << MB_farRight  << "\n";
	outFile << "SGN CONTAMINATION : " << SGNcontamination << "\n"; 
	outFile << "# SIGNAL REGION\n";
	outFile << "MB_sgnL"      << "\t" << MB_sgnL      << "\n"; 
	outFile << "MB_sgnR"      << "\t" << MB_sgnR      << "\n"; 
	outFile << "# PAR-NAME \t # VALUE \t # ERROR"<< "\n";
	for (unsigned int p = 0; p < result.NFreeParameters(); p++){
		//result.GetParameterName(p)
		outFile << fitter.Config().ParSettings(p).Name() << "\t" << *(FitParams + p) << "\t" << *(FitParErr + p) << "\n"; 
	}	
	outFile << "Chi2/Ndf" << "\t" << result.Chi2()/result.Ndf() << "\n";
	outFile.close();

	return 0;
}//Fit_mass_B0()
