#include "MassAnalysis.C"
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

int RunAnalysis( ) {


	//================ Loading files


	//std::string DirPath = "/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/B0ToXKs_2022Apr29/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220429_084035/0000/";

	std::string DirPath = "/eos/cms/store/group/phys_egamma/soffi/XKs/";
	int Nfiles = 0; 
	TString tree_path = "";
	const TString treeName = "/Events";

	//================ Creating chain                                                               

	TChain* chain =new TChain("Events");


	struct dirent* file = NULL; // entry in the directory
	struct stat file_stats;
	const char* directory;
	DIR* dir_pointer = NULL;

	directory = DirPath.c_str();
	dir_pointer = opendir(directory);//point to the directory

	while((file = readdir (dir_pointer))){
		if(file == NULL){
			std::cout << "ERROR null pointer to file" << std::endl;
			exit(-1);
		}

		if (strcmp(file->d_name, "xNANO_") < 0) continue; // skip "." and ".." and "log"

		Nfiles ++;
		//std::cout << file->d_name << std::endl;
		tree_path = DirPath + "/" + file->d_name + treeName; 
		chain->Add(tree_path);
	}

	std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
	cout<<" Number of events: " <<chain->GetEntries()<<std::endl;

  //================ Run analysis                                                                                                                         
  MassAnalysis tree( chain );
  tree.Loop();

  return 0;
}// RunAnalysis()
