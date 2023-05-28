#include "./src/B0toX3872K0s.C"
#include "./src/PreSelDATA2017.C"
#include "./src/PrePostFit.C"
#include "./src/TriggerSelection.C"
#include "./src/HLTapply.C"
#include "./src/SGNvsBKGvariables.C"

#include <iostream>
#include <string>
#include <TStyle.h>
#include <TCanvas.h>
#include <dirent.h>
#include <sys/stat.h>
using namespace std;


std::string MCdata2017 = "/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/B0ToXKs_2022Apr29/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220429_084035/";
std::string Data2017_B = "/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/Xdata2017_2022May17/Charmonium/crab_data_Run2017B/220517_110819/";
std::string Data2017_C = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017C/Run2017C_";
std::string Data2017_D = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017D/Run2017D_";

std::string Data2017_B_0001 = "/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/Xdata2017_2022May17/Charmonium/crab_data_Run2017B/220517_110819/0001";
std::string Data2017_D_0000 = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017D/Run2017D_0000.root";
std::string Data2017_D_0001 = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017D/Run2017D_0001.root";

TChain* LoadTree(TString dataset = "MC"){

	bool LxP = false, T2 = false;	
	
	string LxPlusDirName = "";
	string LxPlusDirPath = "";
	string T2DirName = "";
	int Ndirectory = 1;

	if (dataset ==    "MC"  ){
		LxPlusDirPath =  MCdata2017; 
		LxP = true;
	}
	if (dataset == "data17B" or dataset == "data17"){ 
		LxPlusDirPath = Data2017_B; 
		LxP = true; Ndirectory = 2;
	}
	if (dataset == "data17D" or dataset == "data17"){
		T2DirName = Data2017_D; 
		T2 = true; Ndirectory = 2;
	}
	

	int Nfiles = 0;
	TString tree_path, tree_name = "/Events";
	TChain* chain =new TChain("Events");

	//================ LOADING FILES FROM LXPLUS
	if (LxP){

		struct dirent* file = NULL; // entry in the directory
		struct stat file_stats;
		const char* directory;
		DIR* dir_pointer = NULL;

		if (dataset == "MC") std::cout << " ===== READING MC data ===== " << std::endl;;
		if (dataset == "data17B" or dataset == "data17") std::cout << " ===== READING RUN2 2017(B) data ===== " << std::endl;
		
		for (int d = 0; d < Ndirectory; d++ ){
			LxPlusDirName = LxPlusDirPath + Form("%.4d", d);
			//std::cout << LxPlusDirName << std::endl;
			directory = LxPlusDirName.c_str();
			dir_pointer = opendir(directory);//point to the directory

			while((file = readdir (dir_pointer))){
				if(file == NULL){
					std::cout << "ERROR null pointer to file" << std::endl;
					return NULL;
				}

				if (strcmp(file->d_name, "xNANO_") < 0) continue; // skip "." and ".." and "log"

				Nfiles ++;
				//std::cout << file->d_name << std::endl;
				tree_path = LxPlusDirName + "/" + file->d_name + tree_name; 
				chain->Add(tree_path);
			}
		}
	}
	cout<<" Number of events: " <<chain->GetEntries()<<std::endl;
	//================ LOADING FILES FROM T2 
	if (T2){
		if (dataset == "data17D" or dataset == "data17")	std::cout << " ===== READING RUN2 2017(D) data ===== " << std::endl;
		for (int d = 0; d < Ndirectory; d++ ){
			tree_path = T2DirName + Form("%.4d", d) + ".root" + tree_name; 
			chain->Add(tree_path);	
			//std::cout << tree_path  << std::endl;
		}
	}

	std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
	cout<<" Number of events: " <<chain->GetEntries()<<std::endl;

	return chain;

}//LoadTree()


int JustCheck(TString dataset = "MC"){

	
	TChain* chain = LoadTree(dataset);
	return 0;

}


int RunPrePostFitAnalysis() {

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	PrePostFit tree( chain );

	tree.Loop();

	return 0;
}// RunPrePostFitAnalysis()


int RunTriggerSel() {

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	TriggerSelection tree( chain );

	tree.Loop();

	return 0;
}// RunTriggerSel()

int RunHLTapply(TString dataset ="MC") {

	TChain* chain = LoadTree(dataset);

	//================ Run analysis                                                                                                                         
	HLTapply tree( chain , dataset);

	tree.Loop();

	return 0;
}// RunHLTapply()

int RunSgnBkgAnalysis(){ 

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	SGNvsBKGvariables tree( chain );

	tree.Loop();

	return 0;
}//RunSgnBkgAnalysis() 
