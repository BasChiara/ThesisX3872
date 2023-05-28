#include "RootPlots.h"

#include <iostream>
#include <cmath>

namespace myRootLib{

void histoSetUp(TH1* histo, const TString& x_label, const TString& y_label, Color_t color, bool fill, bool norm){
	//AXIS LABEL
	// ... x axis
	histo->GetXaxis()->SetTitle(x_label);
	histo->GetXaxis()->SetTitleSize(0.04);
	histo->GetXaxis()->SetLabelSize(0.04);
	// ... y axis
	histo->GetYaxis()->SetTitle(y_label);
	histo->GetYaxis()->SetTitleSize(0.04);
	histo->GetYaxis()->SetLabelSize(0.04);

	//WIDTH & COLOR                                                                                                                                                               
	histo->SetLineWidth(3);
	histo->SetLineColor(color);
	if(fill)histo->SetFillColorAlpha(color, 0.3);
	//NORMALIZATION
	if(norm) histo->Scale(1./histo->Integral());
}



// h1 has to be the DATA
// h2 has to be the MONTE CARLO
TCanvas* RatioPlot(TH1* h1, TH1* h2){
	
	TCanvas* c = new TCanvas("c", "", 1024, 1248);
	// create upper TPad
	TPad* up_pad = new TPad("up_pad", "", 0., 0.30, 1.,1.); // xlow, ylow, xup, yup (mother pad reference system)
	up_pad->SetBottomMargin(0); // Upper and lower plot are joined
	up_pad->Draw();
	up_pad->cd();
	h1->SetStats(0);          // No statistics on upper plot
   	h2->Draw("HISTE");               
   	h1->Draw("PE0 same");   
	
	// Avoid the first label (0) to be clipped.
	TAxis *Yaxis = h2->GetYaxis();
	Yaxis->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
	Yaxis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	Yaxis->SetLabelSize(30);
	Yaxis->SetTitleOffset(1.);

	//create lower TPad
	c->cd();
	TPad* ratio_pad = new TPad("ratio_pad", "", 0., 0., 1.,0.30);
	ratio_pad->SetTopMargin(0);
   	ratio_pad->SetBottomMargin(0.4);
   	ratio_pad->SetGridy(); // vertical grid
   	ratio_pad->Draw();
   	ratio_pad->cd();

	TH1F* h_ratio = (TH1F*)h1->Clone("h_ratio");
	h_ratio->SetLineColor(kBlack);
   	h_ratio->SetMinimum(0.);
   	h_ratio->SetMaximum(3.);
   	h_ratio->Sumw2();
	h_ratio->SetStats(0);
	h_ratio->Divide(h2);
	// ratio plot style ...
	h_ratio->SetMarkerStyle(20);
	Yaxis = h_ratio->GetYaxis();
	Yaxis->SetTitle("Data/MC ratio");
	Yaxis->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
	Yaxis->ChangeLabel(-1, -1, -1, -1, -1, -1, " ");
   	Yaxis->SetTitleFont(43);
	Yaxis->SetTitleSize(30);
   	Yaxis->SetTitleOffset(1.75);
	Yaxis->SetLabelSize(0.1);
	TAxis *Xaxis = h_ratio->GetXaxis();
	Xaxis->SetTitle(h2->GetXaxis()->GetTitle());
   ////Xaxis->SetLabelFont(43);
   	Xaxis->SetLabelSize(0.1);
   	Xaxis->SetTitleFont(43);
	Xaxis->SetTitleSize(40);
	Xaxis->SetTitleOffset(3.5);
	c->Update();
	h_ratio->Draw("PE");
	up_pad->cd();
	return c;

}//RatioPlot()




};
