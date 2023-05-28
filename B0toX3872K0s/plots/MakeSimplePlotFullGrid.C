#include "TGraphSmooth.h"
#include "TMultiGraph.h"
#include"TPaveText.h"
#include "TPaletteAxis.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TChain.h"
#include "TH1F.h"
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TLatex.h"
#include <iostream>
#include <fstream>

using namespace std;
  

//double lumi = 4860.;

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"
#include <iostream>

void FPCanvasStyle(TCanvas* pad, std::string left_text="", std::string right_text="", int iPosX=0, TString extraText="",
                   bool outOfFrame=false, bool is2DCOLZ=false)
{
  //
  // Global variables
  //

  TString leftText     = left_text;
  float leftTextFont   = 62;  // default is helvetic-bold

  float extraTextFont = 52;  // default is helvetica-italics

  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float rightTextSize     = 0.6;
  float rightTextOffset   = 0.2;
  float leftTextSize      = 0.6;
  float leftTextOffset    = 0.1;  // only used in outOfFrame version


  float relPosX    = 0.045;
  float relPosY    = 0.035;
  float relExtraDX = 1.1;
  float relExtraDY = 1.2;



  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;

  bool drawLogo      = false;

  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.17);
  pad->SetTopMargin(0.08);
  pad->SetRightMargin(0.05);
  if(is2DCOLZ)
    {
      TGaxis::SetMaxDigits(4);
      pad->SetTopMargin(0.07);
      pad->SetRightMargin(0.17);
      pad->SetLeftMargin(0.15);
      for(auto obj : *pad->GetListOfPrimitives())
        {
	  auto obj_name = TString(obj->ClassName());
	  if(obj_name.Contains("2"))
            {
                TPaletteAxis* palette =
		  (TPaletteAxis*)((TH2*)gDirectory->Get(obj->GetName()))->GetListOfFunctions()->FindObject("palette");
                palette->SetX1NDC(0.835);
                palette->SetX2NDC(0.875);
                palette->SetY1NDC(0.13);
                palette->SetY2NDC(0.93);
                gPad->Modified();
                gPad->Update();
            }
        }
    }
    
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();
   TString rightText = right_text;
   
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*leftTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(rightTextSize*t);    
  latex.DrawLatex(1-r,1-t+rightTextOffset*t,rightText);

  if(iPosX==0)
    leftText += "#scale[0.76]{#bf{ "+extraText+"}}";
  if( outOfFrame )
    {
      latex.SetTextFont(leftTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(leftTextSize*t);    
      latex.DrawLatex(l,1-t+rightTextOffset*t,leftText);
    }
  
  pad->cd();

  float posX_=0;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
	{
	  posX_ =   l + 0.045*(1-l-r)*W/H;
	  posY_ = 1-t - 0.045*(1-t-b);
	  float xl_0 = posX_;
	  float yl_0 = posY_ - 0.15;
 	  float xl_1 = posX_ + 0.15*H/W;
	  float yl_1 = posY_;
	  TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
	  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
	  pad_logo->Draw();
	  pad_logo->cd();
	  CMS_logo->Draw("X");
	  pad_logo->Modified();
	  pad->cd();
	}
      else
	{
	  latex.SetTextFont(leftTextFont);
	  latex.SetTextSize(leftTextSize*t);
	  latex.SetTextAlign(align_);
	  latex.DrawLatex(posX_, posY_, leftText);
	  if( extraText != "" ) 
	    {
	      latex.SetTextFont(extraTextFont);
	      latex.SetTextAlign(align_);
	      latex.SetTextSize(extraTextSize*t);
	      latex.DrawLatex(posX_, posY_- relExtraDY*leftTextSize*t, extraText);
	    }
	}
    }
  else if(extraText != "" && iPosX != 0)
    {
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, extraText);      
    }

  pad->Update();
    
  return;
}




double makeInterpolation(TGraph* g,std::string ctau){

  double*  xvals=g->GetX();
  double*  yvals=g->GetY();

  double    deltaylow=99999999;
  double   deltayup=99999999;
  double yobs=1.;

  double yup=999.;
  double xup=999.;
  double ylow=999.;
  double xlow=999.;

  for(int i=0;i<12; i++){
    double deltay=yvals[i]-yobs;
    std::cout<<i<<" xval: "<<xvals[i]<<" yval: "<<yvals[i]<<" deltay: "<<deltay<<std::endl;
    if(deltay>=0 && deltay<deltayup){
      deltayup=deltay;
      yup=yvals[i];
      xup=xvals[i];
      }else if(deltay<0 && abs(deltay)<abs(deltaylow)){
      deltaylow=deltay;
      xlow=xvals[i];
      ylow=yvals[i];
      
    }
  }
  std::cout<<yup<<" "<<ylow<<" "<<xup<<" "<<xlow<<std::endl;
  
  double m=(yup-ylow)/(xup-xlow);
  
  double    xobs=(1-ylow+m*xlow)/m;
  std::cout<<" m: "<< m << "ul: "<<xobs<< std::endl;                                                                                                                                                                 
  TCanvas* c = new TCanvas("c","c",1);
  c->SetLogy();
  
  g->Draw("ape");
  c->SaveAs(TString::Format("~/www/DisplacedPhotons/DisplacedPhotons-TDR/BTL-Studies-2022/Interpolation_FIT_ctau%scm.png",ctau.c_str()));
  return xobs;

}



TGraph* makeROCcurve(TH1F* sigHist, TH1F* bkgHist){

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

  //TGraph *gout;
  //TGraphSmooth *gs = new TGraphSmooth("normal");
  //  grout = gs->SmoothKern(g,"normal",2.0);
  //  gout = gs->SmoothSuper(g,"",0,0);

  return g;
}

void fitGausBkg(TH1F* h, std::string res){
  TCanvas *c_ratio = new TCanvas("c_ratio","",700,800);
  c_ratio->cd();
  double par[6];
  TF1 *G1 = new TF1 ("G1","gaus",0,2.0);
  TF1 *G2 = new TF1 ("G2","gaus",0,2.0);

  //G1->FixParameter(1,0.0);
  //  h->Fit(G1,"R");
  // h->Fit(G2,"R");
    
  //  G1->GetParameters(&par[0]);
  //  G2->GetParameters(&par[3]);

  TF1 *f= new TF1 ("f","gaus(0)+gaus(3)",0,2.0);
  f->SetParameters(par);
  h->Fit(f,"R");
  f->SetLineColor(2);

  // 6 parameters of the fit: (p0 to p5)
  double Ncons  = f->GetParameter(0);
  double Gmean= f->GetParameter(1);
  double Gwdth = f->GetParameter(2);
  h->Draw("hist");
  f->Draw("same");
  //  legend2->Draw("same");
  FPCanvasStyle(c_ratio, "CMS", "14 TeV", 0,"#it{Phase-2 Simulation}", 1,0);
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/1obeta_DY_Fit_%s.png",res.c_str()));
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/1obeta_DY_Fit_%s.pdf",res.c_str()));

}
void makePlot1obeta(std::string smass, double* s_eff30, double* b_eff30, double* s_eff90, double* b_eff90, double _1obeta_cut30, double _1obeta_cut90, double* rms_DY_30ps, double* rms_DY_90ps){

  TFile* file_S;
  TFile* file_B;
  file_S = new TFile(TString::Format("ntuples_%s_V0.root",smass.c_str()), "READ");
  file_B = new TFile("ntuples_DY_V2.root", "READ");
  TTree* tree_S = (TTree*) file_S->Get("demo/tree");
  TTree* tree_B = (TTree*) file_B->Get("demo/tree");

  TH1F* h1_1ob_S=new TH1F("h1_1ob_S", "", 300,0.9,3.);
  tree_S->Draw("genpar_1obeta>>h1_1ob_S", "genpar_pt>1&&(abs(genpar_id)==1000015)");
  TH1F* h1_1ob_B=new TH1F("h1_1ob_B", "", 300,0.9,3.);
  tree_B->Draw("genpar_1obeta>>h1_1ob_B", "genpar_pt>1&&(abs(genpar_id)==11||abs(genpar_id)==13||abs(genpar_id)==15)");

  TH1F* h1_1ob30_S=new TH1F("h1_1ob30_S", "", 300,0.9,3.);
  tree_S->Draw("genpar_1obeta_smear30>>h1_1ob30_S", "genpar_pt>1&&(abs(genpar_id)==1000015)");
  TH1F* h1_1ob30_B=new TH1F("h1_1ob30_B", "", 300,0.9,3.);
  tree_B->Draw("genpar_1obeta_smear30>>h1_1ob30_B", "genpar_pt>1&&(abs(genpar_id)==11||abs(genpar_id)==13||abs(genpar_id)==15)");

  TH1F* h1_1ob90_S=new TH1F("h1_1ob90_S", "", 300,0.9,3.);
  tree_S->Draw("genpar_1obeta_smear90>>h1_1ob90_S", "genpar_pt>1&&(abs(genpar_id)==1000015)");
  TH1F* h1_1ob90_B=new TH1F("h1_1ob90_B", "", 300,0.9,3.);
  tree_B->Draw("genpar_1obeta_smear90>>h1_1ob90_B", "genpar_pt>1&&(abs(genpar_id)==11||abs(genpar_id)==13||abs(genpar_id)==15)");



  //------fix legend-----//
  TLegend* legend1;
  legend1 = new TLegend(0.35, 0.6, 0.89, 0.9);
  legend1->SetFillColor(kWhite);
  legend1->SetTextFont(42);
  legend1->SetTextSize(1.5*legend1->GetTextSize());
  legend1->SetHeader(TString::Format("HSCP - #tilde{#tau} (%s GeV) ",smass.c_str()));
  //  legend1->AddEntry(h1_1ob_S,"HSCP - GEN", "LF");
  //  legend1->AddEntry(h1_1ob_B,"DY - GEN", "LF");
  legend1->AddEntry(h1_1ob30_S,"HSCP - MTD(30ps)", "LF");
  legend1->AddEntry(h1_1ob30_B,"DY - MTD(30ps)", "LF");
  legend1->AddEntry(h1_1ob90_S,"HSCP - BTL(90ps) + ETL(30ps)", "LF");
  legend1->AddEntry(h1_1ob90_B,"DY - BTL(90ps) + ETL(30ps)", "LF");

  TCanvas *c_ratio = new TCanvas("c_ratio","",700,800);
  c_ratio->cd();
  //  c_ratio->SetLogy();

  TH1F* h1_test = new TH1F("h1_test", "", 300,0.9,3.);
  h1_test->GetYaxis()->SetTitle("Events [A.U.]");
  h1_test->GetXaxis()->SetTitle("1/#beta");
  //  h1_test->GetXaxis()->SetRangeUser(0.9,1.7);
  //  h1_1ob_S->SetFillStyle(kDashed);
  h1_1ob_S->SetLineStyle(kDashed);

  h1_1ob30_S->SetLineColor(kRed);
  h1_1ob30_B->SetLineColor(kRed);
  h1_1ob30_S->SetLineStyle(kDashed);

  h1_1ob90_S->SetLineColor(kBlue);
  h1_1ob90_B->SetLineColor(kBlue);
  h1_1ob90_S->SetLineStyle(kDashed);

  h1_test->Draw("hist");
  //  h1_1ob_B->DrawNormalized("histsame");
  // h1_1ob_S->DrawNormalized("histsame");
  h1_1ob30_S->Scale(1./h1_1ob30_S->Integral());
  h1_1ob30_B->Scale(1./h1_1ob30_B->Integral());
  h1_1ob30_S->Draw("histsame");
  h1_1ob30_B->Draw("histsame");

  h1_1ob90_S->Scale(1./h1_1ob90_S->Integral());
  h1_1ob90_B->Scale(1./h1_1ob90_B->Integral());
  h1_1ob90_S->Draw("histsame");
  h1_1ob90_B->Draw("histsame");

  legend1->Draw("same");
  FPCanvasStyle(c_ratio, "CMS", "14 TeV", 0,"#it{Phase-2 Simulation}", 1,0);
  c_ratio->SaveAs("~/www/HSCP-TDR/1obeta.png");
  c_ratio->SaveAs("~/www/HSCP-TDR/1obeta.pdf");
  

  c_ratio->SetLogy();
  h1_test->GetYaxis()->SetRangeUser(0.0001, 10);
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/1obeta_%s_LOG.png",smass.c_str()));
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/1obeta_%s_LOG.pdf",smass.c_str()));
  
  //  fitGausBkg(h1_1ob30_B, "30");
  //  fitGausBkg(h1_1ob90_B, "90");
  

  TGraph* g_GEN = makeROCcurve(h1_1ob_S,h1_1ob_B);
  TGraph* g_30 = makeROCcurve(h1_1ob30_S,h1_1ob30_B);
  TGraph* g_90 = makeROCcurve(h1_1ob90_S,h1_1ob90_B);
  TLegend* legend2;
  legend2 = new TLegend(0.35, 0.45, 0.89, 0.65);
  legend2->SetFillColor(kWhite);
  legend2->SetTextFont(42);
  legend2->SetTextSize(1.5*legend2->GetTextSize());
  legend2->SetHeader(TString::Format("HSCP - #tilde{#tau} (%s GeV) ",smass.c_str()));
  //  legend2->AddEntry(g_GEN,"ROC - GEN", "LF");
  legend2->AddEntry(g_30,"ROC - MTD(30ps)", "LF");
  legend2->AddEntry(g_90,"ROC - BTL(90ps) + ETL(30ps)", "LF");

  g_30->SetLineColor(kRed);
  g_90->SetLineColor(kBlue);
  TH1F* h1_test2 = new TH1F("h1_test2", "", 1000,0.0001,1.);
  h1_test2->GetYaxis()->SetTitle("Background Efficiency");
  h1_test2->GetXaxis()->SetTitle("Signal Efficiency");
  h1_test2->Draw("HIST");
  //  g_90->GetXaxis()->SetRangeUser(0.0000001,5);
  //  g_90->GetXaxis()->SetTitle("Background Efficiency");
  //  g_90->GetYaxis()->SetTitle("Signal Efficiency");
  g_90->Draw("LSAME");
  g_30->Draw("LSAME");
  //  g_GEN->Draw("LSAME");
  legend2->Draw("same");
  FPCanvasStyle(c_ratio, "CMS", "14 TeV", 0,"#it{Phase-2 Simulation}", 1,0);
  c_ratio->SetLogy(0);
  //  c_ratio->SetLogx(0);
  c_ratio->SetLogx();
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/ROC_%s.png",smass.c_str()));
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/ROC_%s.pdf",smass.c_str()));



  //compute efficiencies
  *s_eff30 = h1_1ob30_S->Integral(h1_1ob30_S->FindBin(_1obeta_cut30), h1_1ob30_S->GetNbinsX());
  *b_eff30 = h1_1ob30_B->Integral(h1_1ob30_B->FindBin(_1obeta_cut30), h1_1ob30_B->GetNbinsX());
  *s_eff90 = h1_1ob90_S->Integral(h1_1ob90_S->FindBin(_1obeta_cut90), h1_1ob90_S->GetNbinsX());
  *b_eff90 = h1_1ob90_B->Integral(h1_1ob90_B->FindBin(_1obeta_cut90), h1_1ob90_B->GetNbinsX());

  *rms_DY_30ps = h1_1ob30_B->GetRMS(1);
  *rms_DY_90ps = h1_1ob90_B->GetRMS(1);
}



void makePlotMass(std::string smass, double min, double max){

  TFile* file_S;
  file_S = new TFile(TString::Format("ntuples_%s_V0.root",smass.c_str()), "READ");
  TTree* tree_S = (TTree*) file_S->Get("demo/tree");

  TH1F* h1_GEN=new TH1F("h1_GEN", "", 750,100,1600);
  TH1F* h1_smear30=new TH1F("h1_smear30", "", 750,100,1600);
  TH1F* h1_smear50=new TH1F("h1_smear50", "", 750,100,1600);
  TH1F* h1_smear70=new TH1F("h1_smear70", "", 750,100,1600);
  TH1F* h1_smear90=new TH1F("h1_smear90", "", 750,100,1600);
  tree_S->Draw("genpar_mass>>h1_GEN", "genpar_pt>1&&(abs(genpar_id)==1000015)");
  tree_S->Draw("genpar_mass_smear30>>h1_smear30", "genpar_pt>1&&(abs(genpar_id)==1000015)");
  tree_S->Draw("genpar_mass_smear50>>h1_smear50", "genpar_pt>1&&(abs(genpar_id)==1000015)");
  tree_S->Draw("genpar_mass_smear70>>h1_smear70", "genpar_pt>1&&(abs(genpar_id)==1000015)");
  tree_S->Draw("genpar_mass_smear90>>h1_smear90", "genpar_pt>1&&(abs(genpar_id)==1000015)");


  //------fix legend-----//
  TLegend* legend1;
  legend1 = new TLegend(0.56, 0.55, 0.89, 0.9);
  legend1->SetFillColor(kWhite);
  legend1->SetTextFont(42);
  legend1->SetTextSize(1.5*legend1->GetTextSize());
  legend1->SetHeader("HSCP - #tilde{#tau} (432 GeV) ");
  //  legend1->AddEntry(h1_GEN,"HSCP - GEN", "LF");
  legend1->AddEntry(h1_smear30,"HSCP - MTD(30ps)", "LF");
  legend1->AddEntry(h1_smear50,"HSCP - BTL(50ps) + ETL(30ps)", "LF");
  legend1->AddEntry(h1_smear70,"HSCP - BTL(70ps) + ETL(30ps)", "LF");
  legend1->AddEntry(h1_smear90,"HSCP - BTL(90ps) + ETL(30ps)", "LF");


  TCanvas *c_ratio = new TCanvas("c_ratio","",700,800);
  c_ratio->cd();
  //  c_ratio->SetLogy();

  TH1F* h1_test = new TH1F("h1_test", "", 750,100,1600);
  h1_smear90->GetYaxis()->SetTitle("Events [A.U.]");
  h1_smear90->GetXaxis()->SetTitle("mass[GeV]");
  h1_smear90->GetXaxis()->SetRangeUser(min,max);
  //  h1_smear90->GetYaxis()->SetRangeUser(0.00001,1);
  //  h1_1ob_S->SetFillStyle(kDashed);
  h1_GEN->SetLineStyle(kDashed);

  h1_smear30->SetLineColor(kRed);
  //  h1_smear30->SetLineStyle(kDashed);
  h1_smear50->SetLineColor(kOrange-3);
  //  h1_smear50->SetLineStyle(kDashed);
  h1_smear70->SetLineColor(kGreen+2);
  //  h1_smear70->SetLineStyle(kDashed);
  h1_smear90->SetLineColor(kBlue);
  //  h1_smear90->SetLineStyle(kDashed);

  h1_smear90->DrawNormalized("hist");
  //  h1_GEN->DrawNormalized("histsame");

  h1_smear30->DrawNormalized("histsame");
  h1_smear50->DrawNormalized("histsame");
  h1_smear70->DrawNormalized("histsame");
  h1_smear90->DrawNormalized("histsame");


  legend1->Draw("same");
  FPCanvasStyle(c_ratio, "CMS", "14 TeV", 0,"#it{Phase-2 Simulation}", 1,0);
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/mass_%s.png",smass.c_str()));
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/mass_%s.pdf",smass.c_str()));
  

  c_ratio->SetLogy();
  h1_test->GetYaxis()->SetRangeUser(0.0001, 10);
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/mass_%s_LOG.png",smass.c_str()));
  c_ratio->SaveAs(TString::Format("~/www/HSCP-TDR/mass_%s_LOG.pdf",smass.c_str()));
  
}

void  writeDatacard(std::string smass,double sig,std::string res,std::string lumi){
  
  ofstream datacard;
  datacard.open(TString::Format("datacards/datacard_%sGeV_Res%sps_Lumi%sfbinv.txt",smass.c_str(),res.c_str(),lumi.c_str()));
  datacard << "# Simple counting experiment\n";
  datacard << "# Simplified version of the HSCP analysis \n";
  datacard << "imax 1\n";
  datacard << "jmax 1\n";
  datacard << "kmax 1\n";
  datacard << "# we have just one channel, in which we observe 0 events\n";
  datacard << "bin b1\n";
  datacard << "observation 0\n";
  datacard << "bin                 b1 b1\n";
  datacard << "process         signal  background \n";
  datacard << "process          0      1  \n";
  datacard <<   TString::Format("rate           %f   0.00000000001 \n",sig);
  datacard << "  lumi    lnN      1.26  1.26   \n";
  datacard.close();
}

void makeAllPlots1oBeta(){

  double s_eff30[10];
  double b_eff30[10];
  double s_eff90[10];
  double b_eff90[10];
  double rms30;
  double rms90;
  double cut30=1.05;
  double cut90=1.10;
  std::vector<std::string>smasses;
  smasses.push_back("200");
  smasses.push_back("308");
  smasses.push_back("432");
  smasses.push_back("557");
  smasses.push_back("651");
  smasses.push_back("745");
  smasses.push_back("871");
  smasses.push_back("1029");
  smasses.push_back("1218");
  smasses.push_back("1409");
  double masses[10]={200,308,432,557,651,745,871,1029,1218,1409};
  for(int i=0;i<10;i++) makePlot1obeta(smasses[i],&s_eff30[i],&b_eff30[i],&s_eff90[i],&b_eff90[i],cut30,cut90, &rms30, &rms90);
  std::cout<<"--------------> "<<rms30<< " "<<rms90<<std::endl;

  //  for(int i=0;i<10;i++)  std::cout<<s_eff30[i]<<" "<<b_eff30[i]<<" "<<s_eff90[i]<<" "<<b_eff90[i]<<std::endl;
  TCanvas *c_ratio = new TCanvas("c_ratio","",700,800);
  c_ratio->cd();
  TGraph* g_30 = new TGraph(10, masses,s_eff30);
  TGraph* g_90 = new TGraph(10, masses,s_eff90);
  TLegend* legend2;
  legend2 = new TLegend(0.2, 0.2, 0.89, 0.45);
  legend2->SetFillColor(kWhite);
  legend2->SetTextFont(42);
  legend2->SetTextSize(1.5*legend2->GetTextSize());
  legend2->AddEntry(g_30,"HSCP - MTD(30ps) -  1/#beta > 1.05", "L");                                                                                                      
  legend2->AddEntry(g_90,"HSCP - BTL(90ps) + ETL(30ps) -  1/#beta > 1.10", "L");                                                                                                         
  g_30->SetLineColor(kRed);
  g_90->SetLineColor(kBlue);
  g_30->GetXaxis()->SetTitle("Mass[GeV]");
  g_30->GetYaxis()->SetTitle("Signal Efficiency");
  g_30->GetYaxis()->SetRangeUser(0.,1.1);
  g_30->Draw("AL");
  g_90->Draw("LSAME");

  legend2->Draw("same");
  FPCanvasStyle(c_ratio, "CMS", "14 TeV", 0,"#it{Phase-2 Simulation}", 1,0);
  c_ratio->SaveAs("~/www/HSCP-TDR/SignalEff_vs_mass.png");
  c_ratio->SaveAs("~/www/HSCP-TDR/SignalEff_vs_mass.pdf");


  //writing datacards

  double lumi=3000;
  double xsec[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.}; //xsec in fb*0.001


  double signorm_30[10];
  double signorm_90[10];
  for(int i=0;i<10;i++){
    std::cout<<s_eff30[i]<<" "<<s_eff90[i]<<std::endl;
    signorm_30[i]=0.001*xsec[i]*lumi*s_eff30[i];
    signorm_90[i]=0.001*xsec[i]*lumi*s_eff90[i];
    writeDatacard(smasses[i],signorm_30[i],"30","3000");
    writeDatacard(smasses[i],signorm_90[i],"90","3000");
  }


}

void makeAllPlotsMass(){

  makePlotMass("200", 100,450);
  makePlotMass("308", 150,550);
  makePlotMass("432",300,700);
  makePlotMass("557", 400, 800);
  makePlotMass("651", 500, 900);
  makePlotMass("745", 600,1000);
  makePlotMass("871", 600, 1100);
  makePlotMass("1029", 850,1300);
  makePlotMass("1218", 1000, 1500);
  makePlotMass("1409", 1200, 1600);
}




double getExpLimAndUnc(std::string smass, double* err,std::string res,std::string lumi){

  TFile* file; 
  file= TFile::Open(TString::Format("higgsCombine%sGeV_Res%sps_Lumi%sfbinv.AsymptoticLimits.mH120.root", smass.c_str(),res.c_str(),lumi.c_str()), "READ");
  
  TTree *t = (TTree *) file->Get("limit");
  Double_t  limit, limitErr = 0;
  Float_t quantileExpected;
  double lim_exp=0;
  if(t!= NULL){
    t->SetBranchAddress("limit", &limit);
    t->SetBranchAddress("quantileExpected", &quantileExpected);

    for (int i = 0, n = t->GetEntries(); i < n; ++i) {
      t->GetEntry(i);
      if(i==2)     lim_exp=limit;
      if(i==1)    err[1]=limit;
    }
  }

  std::cout<<smass<<" "<<lim_exp<<std::endl;
  return lim_exp;

}



TGraph* plotLimit(std::string res,std::string lumi){
  //get expected limits                                                                                                                                    
  double exp_lim_200[2];
  exp_lim_200[0]=getExpLimAndUnc("200",exp_lim_200,res,lumi);
  double exp_lim_308[2];
  exp_lim_308[0]=getExpLimAndUnc("308",exp_lim_308,res,lumi);
  double exp_lim_432[2];
  exp_lim_432[0]=getExpLimAndUnc("432",exp_lim_432,res,lumi);
  double exp_lim_557[2];
  exp_lim_557[0]=getExpLimAndUnc("557",exp_lim_557,res,lumi);
  double exp_lim_651[2];
  exp_lim_651[0]=getExpLimAndUnc("651",exp_lim_651,res,lumi);
  double exp_lim_745[2];
  exp_lim_745[0]=getExpLimAndUnc("745",exp_lim_745,res,lumi);
  double exp_lim_871[2];
  exp_lim_871[0]=getExpLimAndUnc("871",exp_lim_871,res,lumi);
  double exp_lim_1029[2];
  exp_lim_1029[0]=getExpLimAndUnc("1029",exp_lim_1029,res,lumi);
  double exp_lim_1218[2];
  exp_lim_1218[0]=getExpLimAndUnc("1218",exp_lim_1218,res,lumi);
  double exp_lim_1409[2];
  exp_lim_1409[0]=getExpLimAndUnc("1409",exp_lim_1409,res,lumi);

  double masses[10]={200,308,432,557,651,745,871,1029,1218,1409};
  double exp_lim[10]={exp_lim_200[0],exp_lim_308[0],exp_lim_432[0],exp_lim_557[0],exp_lim_651[0],exp_lim_745[0],exp_lim_871[0],exp_lim_1029[0],exp_lim_1218[0],exp_lim_1409[0]};
  double xsec_injected=0.001;//fb

  for(int i=0;i<10;i++)exp_lim[i]=exp_lim[i]*0.001;

  TGraph* g = new TGraph(10,masses,exp_lim);
  TCanvas* c = new TCanvas("c","c",1);
  c->SetTitle(TString::Format("UL vs Mass for (%sps)",res.c_str()));
  c->cd();
  c->SetLogy();
  g->SetLineColor(kRed);
  g->SetMarkerColor(kRed);
  g->GetYaxis()->SetRangeUser(0.00001,0.01);

  g->GetYaxis()->SetTitle("95% C.L. UL on cross section [fb]");
  g->GetXaxis()->SetTitle("Mass [GeV]");
  g->Draw("APL");
  c->SaveAs(TString::Format("~/www/HSCP-TDR/Limits_Res%sps_Lumi%s.png",res.c_str(),lumi.c_str()));

  return g;
}



void plotAllLimits(){

  TGraph* g_30 = plotLimit("30", "3000");
  TGraph* g_90 = plotLimit("90", "3000");

  TCanvas *c_ratio = new TCanvas("c_ratio","",700,800);
  c_ratio->cd();
  TLegend* legend2;
  legend2 = new TLegend(0.2, 0.2, 0.89, 0.4);  legend2->SetFillColor(kWhite);
  legend2->SetTextFont(42);
  legend2->SetHeader("Integrate Luminosity: 3000 fb^{-1}");
  legend2->AddEntry(g_30,"HSCP - MTD(30ps) -  1/#beta > 1.05", "L");                                                                                        
  legend2->AddEntry(g_90,"HSCP - BTL(90ps) + ETL(30ps) -  1/#beta > 1.10", "L"); 



  g_30->SetLineColor(kRed);
  g_90->SetLineColor(kBlue);
  g_90->GetYaxis()->SetRangeUser(0., 0.005);
  g_90->GetXaxis()->SetTitle("Mass[GeV]");
  g_90->GetYaxis()->SetTitle("95% C.L. UL on cross section [fb]");

  g_90->Draw("AL");
  g_30->Draw("LSAME");

  legend2->Draw("same");
  FPCanvasStyle(c_ratio, "CMS", "14 TeV", 0,"#it{Phase-2 Simulation}", 1,0);
  c_ratio->SaveAs("~/www/HSCP-TDR/Limits_vs_mass.png");
  c_ratio->SaveAs("~/www/HSCP-TDR/Limits_vs_mass.pdf");
  c_ratio->SetLogy();
  g_90->GetYaxis()->SetRangeUser(0.0001, 0.01);
  g_90->GetXaxis()->SetTitle("Mass[GeV]");
  g_90->GetYaxis()->SetTitle("95% C.L. UL on cross section [fb]");

  g_90->Draw("AL");
  g_30->Draw("LSAME");
  legend2->Draw("same");

  FPCanvasStyle(c_ratio, "CMS", "14 TeV", 0,"#it{Phase-2 Simulation}", 1,0);
  c_ratio->SaveAs("~/www/HSCP-TDR/Limits_vs_mass_LOG.png");
  c_ratio->SaveAs("~/www/HSCP-TDR/Limits_vs_mass_LOG.pdf");
}
