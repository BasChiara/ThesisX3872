#include "CheckGenLev.C"
#include <iostream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

//USARE SOLO PER FARE ANALISI GEN LEVEL!!
//--> path di Chiara
//  TString path ="/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKs_2022Apr22/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220422_134814/0000/xNANO_mc_2022Apr22_"; 
//int j = 0, FilesToSkip[] = {22, 32, 33, 41, 0};
//
//--> path di Livia
//  TString path ="/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKsLivia_2022Apr26/BPH_BToXKs_XToJPsiRho_JPsiToMuMu/crab_BdToX3872Ks/220426_134635/0000/xNANO_mc_2022Apr26_"xNANO_mc_2022Apr26_";
//


int RunAnalysis() {

  //================ Loading files

   const string dir_name = "/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKsLivia_2022Apr26/BPH_BToXKs_XToJPsiRho_JPsiToMuMu/crab_BdToX3872Ks/220426_134635/0000/";
   const char *directory = "/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKsLivia_2022Apr26/BPH_BToXKs_XToJPsiRho_JPsiToMuMu/crab_BdToX3872Ks/220426_134635/0000";

   if (directory == NULL) {
      std::cout << "ERROR IN ACCESSING THE DIRECTORY" << std::endl;       
      return 1;
   }
   
   DIR* dir_pointer = opendir(directory);//point to the directory
   struct dirent* file = NULL; // entry in the directory
   struct stat file_stats;

   TString tree_path, tree_name = "/Events";
   int Nfiles = 0;
   //================ Creating chain                                                               

   TChain* chain =new TChain("Events");

   while((file = readdir (dir_pointer))){
      if(file == NULL){
	 std::cout << "ERROR null pointer to file" << std::endl;
	 return 1;
      }

      if (strcmp(file->d_name, "xNANO_mc_") < 0) continue; // skip "." and ".."
      tree_path = dir_name + file->d_name;
      stat(tree_path, &file_stats);
      if (file_stats.st_size < 2280000) continue;
      std::cout << file_stats.st_size << "\t" << file->d_name << std::endl;
      Nfiles ++;

      tree_path = dir_name + file->d_name + tree_name; 
      chain->Add(tree_path);
   }

   std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
   cout<<" Number of events: " <<chain->GetEntries()<<std::endl;

   //================ Run analysis                                                                                                                         
   CheckGenLev tree( chain );
   tree.Loop();


   closedir (dir_pointer);

   return 0;
}// RunAnalysis()
