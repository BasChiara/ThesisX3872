#include "PrepAnalysis.C"
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

int RunAnalysis(const TString& which_analysis) {

  //================ Loading files

  TString path = "/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKs_2022Mar30/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220330_155143/0000/xNANO_mc_2022Mar30_";

  TString tree_path = "";
  int file_num[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 11, 14, 16, 18, 19, 20, 21, 23, 24, 25, 26, 29, 30, 34, 37, 38, 42, 43, 44, 47, 48, 49, 50};

  int Nfiles = sizeof(file_num)/ sizeof(*file_num);
  std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
  
  //"BParkNANO_mc_10215.root"
  //================ Creating chain                                                               

  TChain* chain =new TChain("Events");

  for (int f = 0; f < Nfiles; f++){
 
    tree_path = path + std::to_string(file_num[f]) + ".root/Events";
    chain->Add(tree_path);
    //"BParkNANO_mc_10215.root/Events"
  }
  
  cout<<" Number of events: " <<chain->GetEntries()<<std::endl;

  //================ Run analysis                                                                                                                         
  PrepAnalysis tree( chain );
  tree.Loop(which_analysis);

  return 0;
}// RunAnalysis()
