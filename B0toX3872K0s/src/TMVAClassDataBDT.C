#include <iostream>
#include <map>
#include <string>
 
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
 
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"


TString outFile_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/"; 
TString outTree_name = "TMVAoutData";
TString inFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/TMVAinputs.root";
TString SGNtreeName = "inputSIGNAL";
TString inFileBKG_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BKG_data17.root";
TString BKGtreeName = "B0sidebands";
TString SEED = "100"; // use the default SEED

int TMVAClassDataBDT( TString myMethodList = "" )
{

// This loads the library
	TMVA::Tools::Instance();
	// Default MVA methods to be trained + tested
	std::map<std::string,int> Use;

	// Boosted Decision Trees
	Use["BDT"]             = 1; // uses Adaptive Boost
	Use["BDTbest"]         = 1; // uses 5 best performing
	Use["BDTseed"]         = 1; // sto impazzendo con lo SplitSeed
	Use["BDTfinal"]        = 1; // uses Gradient Boost
	// Neural Networks (all are feed-forward Multilayer Perceptrons)
	Use["MLP"]             = 0; // Recommended ANN
	Use["DNN_GPU"]         = 1; // CUDA-accelerated DNN training.
	Use["DNN_CPU"]         = 1; // Multi-core accelerated DNN.
	Use["DNNbest"]         = 1; // Multi-core accelerated DNN.
	// ---------------------------------------------------------------//
	std::cout << std::endl;
	std::cout << "==> Start TMVAClassification" << std::endl;
	//                                                                                                                                                                           
	//--------------------------------------------------------------------------------------------------//

	// Select methods 
	if (myMethodList != "") {
		for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

		std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
		for (UInt_t i=0; i<mlist.size(); i++) {
			std::string regMethod(mlist[i]);

			if (Use.find(regMethod) == Use.end()) {
				std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
				for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
				std::cout << std::endl;
				return 1;
			}
			Use[regMethod] = 1;
		}
	}
	// Here the preparation phase begins

	// Read training and test data

	TFile *inputSGN(0);
	TFile *inputBKG(0);
	TString fname_S = inFileSGN_path;//"../results/TMVAinputs.root";
	if (!gSystem->AccessPathName( fname_S )) {
		inputSGN = TFile::Open( fname_S ); // check if file in local directory exists
	}
	if (!inputSGN) {
		std::cout << "ERROR: could not open SIGNAL data file" << std::endl;
		exit(1);
	}
	TString fname_B = inFileBKG_path;//"../results/TMVAinputs.root";
	if (!gSystem->AccessPathName( fname_B )) {
		inputBKG = TFile::Open( fname_B ); // check if file in local directory exists
	}
	if (!inputBKG) {
		std::cout << "ERROR: could not open BACKGROUND data file" << std::endl;
		exit(1);
	}
	
	std::cout << "\n--- TMVAClassification:        Using input file: (SIGNAL) " << inputSGN->GetName() << std::endl;
	std::cout << "\n--- TMVAClassification:        Using input file: (BACKGROUND) " << inputBKG->GetName() << std::endl;

	// Register the training and test trees

	TTree *signalTree     = (TTree*)inputSGN->Get(SGNtreeName);//TreeS
	TTree *background     = (TTree*)inputBKG->Get(BKGtreeName);//TreeB
 
   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
	outFile_path.Append(outTree_name);
	if (Use["BDT"] )outFile_path.Append("BDTtuning");
	if (Use["BDTseed"] ){outFile_path.Append("_S"); outFile_path.Append(SEED);}
	if (Use["BDTbest"] ) outFile_path.Append("BDTbest"); 
	if (Use["DNN_CPU"] or Use["DNN_GPU"] or Use["MLP"]) outFile_path.Append("NNtuning");
	if (Use["DNNbest"]) outFile_path.Append("DNNbest");
	outFile_path.Append(".root");
   TString outfileName(outFile_path); 
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   
	// Create the factory object. Later you can choose the methods
	// whose performance you'd like to investigate. The factory is
	// the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   //
   // All TMVA output can be suppressed by removing the "!" (not) in front of the "Silent" argument in the option string

   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
   "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  	
	TString outDataset("dataset");
	if (Use["BDT"] or Use["BDTbest"]) outDataset.Append("BDT"); 
	if (Use["DNN_CPU"] or Use["DNNbest"]) outDataset.Append("DNN"); 
   TMVA::DataLoader *dataloader=new TMVA::DataLoader(outDataset);
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
   
   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]


   dataloader->AddVariable( "pTM_B0", 'F' );
   dataloader->AddVariable( "SVprob", 'F' );
   dataloader->AddVariable( "LxySign_B0", 'F' );
   dataloader->AddVariable( "CosAlpha_B0", 'F' );
   dataloader->AddVariable( "DR_Pi1B0", 'F' );
   dataloader->AddVariable( "pT_Pi1", 'F' );
   dataloader->AddVariable( "pT_Rho", 'F' );
   dataloader->AddVariable( "D0_Rho", 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   
   dataloader->AddSpectator( "M_Rho",  "\\ M(\\rho)", "units", 'F' );
   //dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

	//global event weights per tree (see below for setting event-wise weights)
	Double_t signalWeight     = 1.0;
	Double_t backgroundWeight = 1.0;

	// You can add an arbitrary number of signal or background trees
	dataloader->AddSignalTree    ( signalTree,     signalWeight );
	dataloader->AddBackgroundTree( background, backgroundWeight );

	// To give different trees for training and testing, do as follows:
	//
	//     dataloader->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
	//     dataloader->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

	// Use the following code instead of the above two or four lines to add signal and background training and test events "by hand"
	//


	// Apply additional cuts on the signal and background samples (can be different)
	TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

	// Tell the dataloader how to use the training and testing events
	//
	// If no numbers of events are given, half of the events in the tree are used
	// for training, and the other half for testing:
	//
	//    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
	//
	// To also specify the number of testing events, use:
	//
	//    dataloader->PrepareTrainingAndTestTree( mycut,
	//         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
	
	TString PrepareOpt("nTrain_Signal=478:nTrain_Background=2568:SplitMode=Random:SplitSeed=");
	PrepareOpt.Append(SEED);
	PrepareOpt.Append(":NormMode=None:!V");
	std::cout << " ---- TEST TRAINING SPLITTING OPTIONS ----" << std::endl;
	std::cout << " ----  " << PrepareOpt << std::endl;
	dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,PrepareOpt); 
			//"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );

	// ### Book MVA methods
	//
	// Please lookup the various method configuration options in the corresponding cxx files, eg:
	// src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
	// it is possible to preset ranges in the option string in which the cut optimisation should be done:
	// "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

	// Boosted Decision Trees
	if (Use["BDTG"]) // Gradient Boost
		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
				"!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );


	
	if (Use["BDTbest"]){  // Adaptive Boost

	
		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT200_S25_D4_b03_nC35",
		 "!H:!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=35");

		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT200_S25_D5_b03_nC30",
		 "!H:!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30");	

		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT250_S25_D4_b03_nC35",
		 "!H:!V:NTrees=250:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=35");	

		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT250_S25_D5_b03_nC30",
		 "!H:!V:NTrees=250:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30");	

		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT300_S25_D4_b03_nC35",
		 "!H:!V:NTrees=300:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=35");	

		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT300_S25_D3_b03_nC30",
		 "!H:!V:NTrees=300:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30");	

		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT200_S25_D3_b03_nC30",
		 "!H:!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30");	

	}

	TString BDTname = "BDT";
	TString BDTopt("!H:!V:NTrees=");
	int NT = 3;
	TString NTrees[] = {"200", "250", "300", "350"}; 
	int NS = 2;
	TString NodeSize[] = {"2.5%", "5%"};
	int Nd = 3;
	TString MaxDepth[] = {"3","4","5"};
	int Nb = 3;
	TString Beta[] = {"0.3", "0.5", "0.7"};
	int NC = 3;
	TString NCuts[] = {"20", "30", "35"};

	if (Use["BDT"]){  // Adaptive Boost
		for (int t = 0; t <NT; t++ ){
			for (int s = 0; s <NS; s++ ){
				for (int d = 0; d <Nd; d++ ){
					for (int b = 0; b <Nb; b++ ){
						for (int c = 0; c <NC; c++ ){

							BDTname.Append("_nT");			BDTname.Append(NTrees[t]);
							BDTname.Append("_S");			BDTname.Append(NodeSize[s]);
							BDTname.Append("_D");			BDTname.Append(MaxDepth[d]);
							BDTname.Append("_b");         BDTname.Append(Beta[b]);
							BDTname.Append("_nC");			BDTname.Append(NCuts[c]);

							BDTopt.Append(NTrees[t]); 
							BDTopt.Append(":MinNodeSize=");	BDTopt.Append(NodeSize[s]);
							BDTopt.Append(":MaxDepth=");		BDTopt.Append( MaxDepth[d]);
							BDTopt.Append(":BoostType=AdaBoost:AdaBoostBeta=");	BDTopt.Append(Beta[b]);
							BDTopt.Append(":UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=");	BDTopt.Append(NCuts[c]);

							std::cout << BDTopt << std::endl;					

							factory->BookMethod( dataloader, TMVA::Types::kBDT, BDTname, BDTopt);	

							BDTname = "BDT"; 
							BDTopt = "!H:!V:NTrees=";
						}	
					}
				}
			}
		}
	}
	if (Use["BDTseed"]){
		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT100_b08_nC30",
				"!H:!V:NTrees=100:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.8:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30");

		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT50_b05_nC30",
				"!H:!V:NTrees=50:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30");

	}
	if (Use["BDTfinal"]) // Bagging
		factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT_nT100_b08_nC30",
		 "!H:!V:NTrees=100:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.8:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=30");	


	if (Use["DNN_CPU"] or Use["DNN_GPU"]) {

		// General Options.
		TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
				"WeightInitialization=XAVIERUNIFORM");

		// General layout.
		int NLyt = 4;
		TString layoutString[] = {"Layout=TANH|64,TANH|64,TANH|32,LINEAR", "Layout=TANH|64,TANH|128,TANH|128,TANH|32,LINEAR","Layout=TANH|64,TANH|32,TANH|32,TANH|16,LINEAR","Layout=RELU|32,RELU|64,RELU|64,RELU|64,RELU|32,LINEAR", "Layout=TANH|128,RELU|64,RELU|64,RELU|32,RELU|8,LINEAR"};
		// Parameters grid
		int Nrates = 2;
		TString LRates[] = {"1e-3", "1e-4", "1e-2"};
		int NmE = 1;
		TString MaxEpoch[] = {"1000", "150", "500"};
		int Np = 3;
		TString momentum[] = {"0.9", "0.8", "0.7"};
		int NBs = 2;
		TString BatchSize[] = {"30", "50", "100"};
		//Define Training strategy. One could define multiple strategy string separated by the "|" delimiter
		//
		TString trainingStrategyString = "TrainingStrategy=LearningRate=";
		
		//TString trainingStrategyString = "TrainingStrategy=LearningRate=1e-1,MaxEpoch=100,DropConfig=0.0+0.1+0.1+0.1+0.1+0.1,ConvergenceSteps=20|LearningRate=1e-2,MaxEpoch=1000,DropConfig=0.0,ConvergenceSteps=100|LearningRate=1e-3,MaxEpoch=1200,DropConfig=0.0,ConvergenceSteps=100";


		TString NN_name = "DNN";

		// Cuda implementation.
		if (Use["DNN_GPU"]) {
			TString gpuOptions = dnnOptions +":"+ layoutString[1]+ ":" + trainingStrategyString + ":Architecture=CPU";
			factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_GPU", gpuOptions);
		}
		// Multi-core CPU implementation.
		if (Use["DNN_CPU"]) {
			for(int l = 0; l < NLyt; l++){
				for(int r = 0; r < Nrates; r++){
					for(int e = 0; e < NmE; e++){
						for(int p = 0; p < Np; p++){
							for(int b = 0; b < NBs; b++){

								NN_name.Append("_CPU_L");	NN_name.Append(std::to_string(l));// layout 
								NN_name.Append("_LR"); 		NN_name.Append(LRates[r]);// learning rate 
								NN_name.Append("_E"); 		NN_name.Append(MaxEpoch[e]);// MAX epoch 
								NN_name.Append("_P"); 		NN_name.Append(momentum[p]);// MAX epoch 
								NN_name.Append("_BS"); 		NN_name.Append(BatchSize[b]);// MAX epoch 

								dnnOptions.Append (":"); dnnOptions.Append (layoutString[l]);
								dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString); dnnOptions.Append (LRates[r]);
								dnnOptions.Append (",MaxEpoch=");  dnnOptions.Append(MaxEpoch[e]);
								dnnOptions.Append (",Momentum=");  dnnOptions.Append(momentum[p]);
								dnnOptions.Append (",BatchSize=");  dnnOptions.Append(BatchSize[b]);
								dnnOptions.Append ("ConvergenceSteps=80,TestRepetitions=1,Optimizer=ADAM,WeightDecay=1e-4,Regularization=None,DropConfig=0.0");

								TString cpuOptions = dnnOptions + ":Architecture=CPU";
								factory->BookMethod(dataloader, TMVA::Types::kDL, NN_name , cpuOptions);



								dnnOptions = "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM";
								NN_name = "DNN";
								trainingStrategyString = "TrainingStrategy=LearningRate=";
							}
						}	
					}
				}
			}
		}
	}


	if(Use["DNNbest"]){

	
		TString dnnInitOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM:");
		TString dnnFinOptions ("ConvergenceSteps=80,TestRepetitions=1,Optimizer=ADAM,WeightDecay=1e-4,Regularization=None,DropConfig=0.0");
		TString layoutString[] = {"Layout=TANH|64,TANH|64,TANH|32,LINEAR", "Layout=TANH|64,TANH|128,TANH|128,TANH|32,LINEAR","Layout=TANH|64,TANH|32,TANH|32,TANH|16,LINEAR","Layout=RELU|32,RELU|64,RELU|64,RELU|64,RELU|32,LINEAR", "Layout=TANH|128,RELU|64,RELU|64,RELU|32,RELU|8,LINEAR"};
		TString cpuOptions("");		

		cpuOptions = dnnInitOptions + layoutString[2] + ":TrainingStrategy=LearningRate=1e-4,MaxEpoch=1000,Momentum=0.9,BatchSize=30" + dnnFinOptions; 
		factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU_L2_LR1e4_E1000_P09_BS30", cpuOptions);
		
		cpuOptions = dnnInitOptions + layoutString[2] + ":TrainingStrategy=LearningRate=1e-4,MaxEpoch=1000,Momentum=0.7,BatchSize=30" + dnnFinOptions; 
		factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU_L2_LR1e4_E1000_P07_BS30", cpuOptions);

		cpuOptions = dnnInitOptions + layoutString[0] + ":TrainingStrategy=LearningRate=1e-4,MaxEpoch=1000,Momentum=0.8,BatchSize=30" + dnnFinOptions; 
		factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU_L0_LR1e4_E1000_P08_BS30", cpuOptions);

		cpuOptions = dnnInitOptions + layoutString[0] + ":TrainingStrategy=LearningRate=1e-4,MaxEpoch=1000,Momentum=0.7,BatchSize=30" + dnnFinOptions; 
		factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU_L0_LR1e4_E1000_P07_BS30", cpuOptions);

		cpuOptions = dnnInitOptions + layoutString[1] + ":TrainingStrategy=LearningRate=1e-4,MaxEpoch=1000,Momentum=0.8,BatchSize=30" + dnnFinOptions; 
		factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU_L1_LR1e4_E1000_P08_BS30", cpuOptions);


	}


	// Now you can tell the factory to train, test, and evaluate the MVAs
	//
	// Train MVAs using the set of training events
	factory->TrainAllMethods();

	// Evaluate all MVAs using the set of test events
	factory->TestAllMethods();

	// Evaluate and compare performance of all configured MVAs
	factory->EvaluateAllMethods();

	// --------------------------------------------------------------

	// Save the output
	outputFile->Close();

	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;

	delete factory;
	delete dataloader;
	// Launch the GUI for the root macros
	if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

	return 0;

}//TMVAClassification()


int main( int argc, char** argv )
{
	// Select methods (don't look at this code - not of interest)
	TString methodList;
	for (int i=1; i<argc; i++) {
		TString regMethod(argv[i]);
		if(regMethod=="-b" || regMethod=="--batch") continue;
		if (!methodList.IsNull()) methodList += TString(",");
		methodList += regMethod;
	}
	return TMVAClassDataBDT(methodList);
}
