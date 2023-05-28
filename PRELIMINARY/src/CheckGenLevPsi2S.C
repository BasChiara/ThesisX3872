#define CheckGenLevPsi2S_cxx
#include "../include/CheckGenLevPsi2S.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void CheckGenLevPsi2S::Loop()
{


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   const ULong64_t nbreak = 1000000;


   //BRANCHES
   fChain->SetBranchStatus("*",0); 

   fChain->SetBranchStatus("HLT_DoubleMu4_JpsiTrk_Displaced",1);
   fChain->SetBranchStatus("HLT_Dimuon25_Jpsi",1);
   fChain->SetBranchStatus("HLT_DoubleMu4_JpsiTrkTrk_Displaced",1);
   fChain->SetBranchStatus("HLT_Dimuon18_PsiPrime",1);
   fChain->SetBranchStatus("HLT_DoubleMu4_PsiPrimeTrk_Displaced",1);
   fChain->SetBranchStatus("nGenPart",1);
   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_genPartIdxMother",1);
   fChain->SetBranchStatus("GenPart_eta",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   fChain->SetBranchStatus("GenPart_pt",1);
   fChain->SetBranchStatus("GenPart_mass",1);


   
  //HISTOGRAMS

  //----Generator----//
  
  Int_t Nbins = 40;
  TH1F h_Gpt_mu("gen_pt_mu","", Nbins, 0.,20);
  TH1F h_GptL_mu("gen_Lpt_mu","", Nbins, 0.,20);
  TH1F h_GptSL_mu("gen_SLpt_mu","", Nbins, 0.,20);
  TH1F h_Gpt_pi("gen_pt_pi","", Nbins, 0.,6);
  TH1F h_Gpt_JPsi("gen_pt_JPsi","", Nbins, 5.,35);
  TH1F h_Gpt_PiPi("gen_pt_PiPi","", Nbins, 0.,10.);
  TH1F h_Gpt_Psi("gen_pt_Psi","", Nbins, 5., 50.);
  TH1F h_Gpt_Ks("gen_pt_Ks","", Nbins, 0.,20);
  TH1F h_Gpt_B("gen_pt_B","", Nbins, 0.,50);

  TH1F h_Geta_mu("gen_eta_mu","", Nbins, -3.5,3.5);
  TH1F h_Geta_pi("gen_eta_pi","", Nbins, -3.5,3.5);
  TH1F h_Geta_JPsi("gen_eta_JPsi","", Nbins, -3.5,3.5);
  TH1F h_Geta_PiPi("gen_eta_PiPi","", Nbins, -3.5,3.5);
  TH1F h_Geta_Psi("gen_eta_Psi","", Nbins, -3.5,3.5);
  TH1F h_Geta_Ks("gen_eta_Ks","", Nbins, -3.5,3.5);
  TH1F h_Geta_B("gen_eta_B","", Nbins, -3.5,3.5);

  TH1F h_Gm_mu("gen_m_mu","", Nbins, .05, .15);
  TH1F h_Gm_pi("gen_m_pi","", Nbins, .1, .2);
  TH1F h_Gm_JPsi("gen_m_JPsi","", Nbins, 2.9,3.2);
  TH1F h_Gm_PiPi("gen_m_PiPi","", Nbins, 0.2 ,1.);
  TH1F h_Gm_Psi("gen_m_Psi","", Nbins, 3.6, 3.8);
  TH1F h_Gm_Ks("gen_m_Ks","", Nbins, 0.45, 0.55);
  TH1F h_Gm_B("gen_m_B","", Nbins, 5., 5.5 );

    //----Multiplicity----//
  TH1F h_nB0("nB0", "", 20, 0, 20);
  TH1F h_nK0s("nK0s", "", 10, 0, 10);
  TH1F h_nJPsi("nJPsi", "", 3, 0, 3);
  TH1F h_nPiPi("nPiPi", "", 40, 500, 10000);
  TH1F h_nMuon("nMuon", "", 15, 0, 15);
  TH1F h_nTracks("nTracks", "", 50, 100, 800);
   

   const ULong64_t percentToPrint = (ULong64_t)(nentries/10.);
   float perc;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //if (!HLT_DoubleMu4_PsiPrimeTrk_Displaced) continue;

      if ((jentry+1) % percentToPrint ==0) std::cout << "--> " << Form("%.0f",(float)(jentry+1)/nentries* 100.) << " \%"<< std::endl;
      // ------- GENERATOR -------//
      GenPartFillP4();
      GenPart_FillKinHist( &GenP4_Mum, &h_Gpt_mu, &h_Geta_mu, &h_Gm_mu);
      GenPart_FillKinHist( &GenP4_Mup, &h_Gpt_mu, &h_Geta_mu, &h_Gm_mu);
      GenPart_FillKinHist( &GenP4_JPsi ,&h_Gpt_JPsi, &h_Geta_JPsi, &h_Gm_JPsi);
      GenPart_FillKinHist( &GenP4_Pim, &h_Gpt_pi, &h_Geta_pi, &h_Gm_pi);
      GenPart_FillKinHist( &GenP4_Pip, &h_Gpt_pi, &h_Geta_pi, &h_Gm_pi);
      GenPart_FillKinHist( &GenP4_PiPi, &h_Gpt_PiPi, &h_Geta_PiPi, &h_Gm_PiPi);
      GenPart_FillKinHist( &GenP4_Psi2S, &h_Gpt_Psi ,&h_Geta_Psi, &h_Gm_Psi);
      GenPart_FillKinHist( &GenP4_K0s, &h_Gpt_Ks, &h_Geta_Ks, &h_Gm_Ks);
      GenPart_FillKinHist( &GenP4_B0,&h_Gpt_B, &h_Geta_B, &h_Gm_B);
      if(GenP4_Mum.Pt() > GenP4_Mup.Pt()){
			h_GptL_mu.Fill(GenP4_Mum.Pt());
			h_GptSL_mu.Fill(GenP4_Mup.Pt());
		}else{
			h_GptL_mu.Fill(GenP4_Mup.Pt());
			h_GptSL_mu.Fill(GenP4_Mum.Pt());
		}
    }
    TFile gen_outfile("./outRoot/CheckGenLev_" + tags_ + ".root","RECREATE"); // create the output file
    std::cout << " ---> [OUT] output saved in file " << gen_outfile.GetName() << std::endl;
    h_Gpt_mu.Write();
    h_GptL_mu.Write();
    h_GptSL_mu.Write();
    h_Gpt_pi.Write();
    h_Gpt_JPsi.Write();
    h_Gpt_PiPi.Write();
    h_Gpt_Psi.Write();
    h_Gpt_Ks.Write();
    h_Gpt_B.Write();

    h_Geta_mu.Write();
    h_Geta_pi.Write();
    h_Geta_JPsi.Write();
    h_Geta_PiPi.Write();
    h_Geta_Psi.Write();
    h_Geta_Ks.Write();
    h_Geta_B.Write();

    h_Gm_mu.Write();
    h_Gm_pi.Write();
    h_Gm_JPsi.Write();
    h_Gm_PiPi.Write();
    h_Gm_Psi.Write();
    h_Gm_Ks.Write();
    h_Gm_B.Write();
      
    gen_outfile.Close();

  
}//Loop()


void CheckGenLevPsi2S::GenPart_FillKinHist(ROOT::Math::PtEtaPhiMVector* GenVec, TH1* h_pt, TH1* h_eta, TH1* h_mass){ 

	h_pt->Fill(GenVec->Pt());
	h_eta->Fill(GenVec->Eta());
	h_mass->Fill(GenVec->M());

}//GenPart_FillKinHist()

void CheckGenLevPsi2S::GenPartFillP4(){

	Int_t MumIdx = -1, MupIdx = -1, PimIdx = -1 , PipIdx = -1, JPsiIdx = -1, Psi2SIdx = -1, K0sIdx = -1, B0Idx = -1;
   for (UInt_t g = 0; g < nGenPart; g++){
		// MUONS
		if( (abs(GenPart_pdgId[g]) ==   isMum) &&
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isJPsi) && 
			 (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isPsi2S) &&
			 (abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]]) == isB0) )
			{
			JPsiIdx = GenPart_genPartIdxMother[g];
			Psi2SIdx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]];
			B0Idx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]];
			if (GenPart_pdgId[g] ==   isMum) MumIdx = g;
			if (GenPart_pdgId[g] ==  -isMum) MupIdx = g;
			}
		// PIONS
      if( (abs(GenPart_pdgId[g]) ==   isPip) && 
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isPsi2S) &&
			 (abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]) == isB0) )
			{
			if (GenPart_pdgId[g] == -isPip) PimIdx = g;
			if (GenPart_pdgId[g] ==  isPip) PipIdx = g;
			}
      if( (GenPart_pdgId[g] ==  isK0s) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == abs(isB0))) K0sIdx = g;


   }//on generated ptl

	if(0){
		std::cout << " IDX Mu-  / Mu+  " << MumIdx  << "\t" << MupIdx << std::endl; 
		std::cout << " IDX Pi-  / Pi+  " << PimIdx  << "\t" << PipIdx << std::endl; 
		std::cout << " IDX JPsi " << JPsiIdx << std::endl; 
		std::cout << " IDX Psi2S/ K0s  " << Psi2SIdx << "\t" << K0sIdx << std::endl; 
		std::cout << " IDX B0  " << B0Idx << std::endl; 
	}
	
	if(MumIdx > 0){   
      GenP4_Mum.SetPt(GenPart_pt[MumIdx]);   GenP4_Mum.SetEta(GenPart_eta[MumIdx]);   
		GenP4_Mum.SetPhi(GenPart_phi[MumIdx]); GenP4_Mum.SetM(mMuon);
	} else std::cout << "No generated mu- is found" << std::endl;  

	if(MupIdx > 0){   
      GenP4_Mup.SetPt(GenPart_pt[MupIdx]);   GenP4_Mup.SetEta(GenPart_eta[MupIdx]);
		GenP4_Mup.SetPhi(GenPart_phi[MupIdx]); GenP4_Mup.SetM(mMuon); 
	} else std::cout << "No generated mu+ is found" << std::endl;  

	if(JPsiIdx > 0){   
      GenP4_JPsi.SetPt(GenPart_pt[JPsiIdx]);    GenP4_JPsi.SetEta(GenPart_eta[JPsiIdx]);
		GenP4_JPsi.SetPhi(GenPart_phi[JPsiIdx]);  GenP4_JPsi.SetM(GenPart_mass[JPsiIdx]); 
	} else std::cout << "No generated JPsi is found" << std::endl;  

	if(PimIdx > 0){   
      GenP4_Pim.SetPt(GenPart_pt[PimIdx]);   GenP4_Pim.SetEta(GenPart_eta[PimIdx]);
		GenP4_Pim.SetPhi(GenPart_phi[PimIdx]); GenP4_Pim.SetM(mPion);
	} else std::cout << "No generated pi- is found" << std::endl;  

	if(PipIdx > 0){   
      GenP4_Pip.SetPt(GenPart_pt[PipIdx]);   GenP4_Pip.SetEta(GenPart_eta[PipIdx]);
		GenP4_Pip.SetPhi(GenPart_phi[PipIdx]); GenP4_Pip.SetM(mPion); 
	} else std::cout << "No generated pi+ is found" << std::endl;  

	if((PipIdx > 0) && (PimIdx > 0)){   
      GenP4_PiPi = GenP4_Pip + GenP4_Pim;
	} else std::cout << "No generated PiPi is found" << std::endl;  
	
	if(Psi2SIdx > 0){   
      GenP4_Psi2S.SetPt(GenPart_pt[Psi2SIdx]);   GenP4_Psi2S.SetEta(GenPart_eta[Psi2SIdx]);
		GenP4_Psi2S.SetPhi(GenPart_phi[Psi2SIdx]); GenP4_Psi2S.SetM(GenPart_mass[Psi2SIdx]); 
	} else std::cout << "No generated X(3872) is found" << std::endl;  

	if(K0sIdx > 0){   
      GenP4_K0s.SetPt(GenPart_pt[K0sIdx]);   GenP4_K0s.SetEta(GenPart_eta[K0sIdx]);
		GenP4_K0s.SetPhi(GenPart_phi[K0sIdx]); GenP4_K0s.SetM(GenPart_mass[K0sIdx]); 
	} else std::cout << "No generated K0s is found" << std::endl;  

	if(B0Idx >= 0){   
      GenP4_B0.SetPt(GenPart_pt[B0Idx]);   GenP4_B0.SetEta(GenPart_eta[B0Idx]);
		GenP4_B0.SetPhi(GenPart_phi[B0Idx]); GenP4_B0.SetM(GenPart_mass[B0Idx]); 
	} else std::cout << "No generated B0 is found" << std::endl;  

	
}//GenPartFillP4()