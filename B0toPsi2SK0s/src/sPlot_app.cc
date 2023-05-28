#include "../tools/RootPlots.h"

#include "TFile.h"
#include <iostream>

int main (int argc, char** argv){

    TFile* inFile = new TFile("outRoot/sPlot_MCpreUL.root");
    if(!inFile->IsOpen()){
        std::cout << " [ERROR] cannot open " << inFile->GetName() << std::endl;
        return -1;
    }

    TH1F* h_MC = (TH1F*)inFile->Get("BDTout_sgnMC");
    TH1F* h_Dat = (TH1F*)inFile->Get("BDTout_sgn__BDTout");
    h_Dat->SetTitle("");

    TCanvas* c = myRootLib::RatioPlot(h_Dat, h_MC);

    c->SaveAs("/eos/user/c/cbasile/www/B0toX3872K0s/Psi2S/MVAcuts/BEST_CUT/RatioBDTout_sPlot.png");
    inFile->Close();
    return 0;
}