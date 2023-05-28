#include <iostream>
#include <string>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <dirent.h>
#include <sys/stat.h>
using namespace std;


TChain* LoadTree(){

	//================ Loading files
	const string dir_name = "/pnfs/roma1.infn.it/data/cms/store/user/crovelli/Charmonium/crab_data_Run2017C/220520_202615/0000/"; 
	const char* directory = "/pnfs/roma1.infn.it/data/cms/store/user/crovelli/Charmonium/crab_data_Run2017C/220520_202615/0000";

	if (directory == NULL) {
		std::cout << "ERROR IN ACCESSING THE DIRECTORY" << std::endl;       
		return NULL;	
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
			return NULL;
		}

		if (strcmp(file->d_name, "xNANO_mc_") < 0) continue; // skip "." and ".." and "log"
		
		Nfiles ++;
		std::cout << file->d_name << std::endl;
		tree_path = dir_name + file->d_name + tree_name; 
	}

	std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;

	return chain;

}//LoadTree()


int TestNtupla(){

	TChain* chain = LoadTree();

	return 0;
}// TestNtupla()

int RunAnalysis(const TString& which_analysis) {

  //================ Loading files

  TString path = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017D/Run2017D_000";

  TString tree_path = "";
  int file_num[] = {0, 1};

  int Nfiles = sizeof(file_num)/ sizeof(*file_num);
  std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
  
  //================ Creating chain                                                               

  TChain* chain =new TChain("Events");

  for (int f = 0; f < Nfiles; f++){
 		
	  tree_path = path + std::to_string(file_num[f]) + ".root" + "/Events";
	  //tree_path = path + std::to_string(file_num[f]) + ".root/Events";
	  chain->Add(tree_path);
  }
  
  cout<<" Number of events: " <<chain->GetEntries()<<std::endl;

  //================ Run analysis                                                                                                                         

  return 0;
}// RunAnalysis()

int RunPrePostFitAnalysis() {

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	//PrePostFit tree( chain );

	//tree.Loop();

	return 0;
}// RunPrePostFitAnalysis()


int RunTriggerSel() {

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	//TriggerSelection tree( chain );

	//tree.Loop();

	return 0;
}// RunTriggerSel()


int RunSgnBkgAnalysis(){ 

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	//SGNvsBKGvariables tree( chain );

	//tree.Loop();

	return 0;
}//RunSgnBkgAnalysis() 
