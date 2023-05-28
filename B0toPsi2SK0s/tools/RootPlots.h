//RootPlots.h
//
//   my ROOT library

#include "TH1F.h"
#include "TH2F.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"

namespace myRootLib{

void histoSetUp(TH1* histo, const TString& x_label, const TString& y_label, Color_t color, bool fill = true , bool norm = true);
TCanvas* RatioPlot(TH1* h1, TH1* h2);

};
