#define PrepAnalysisPsi2S_cxx
#include "../include/PrepAnalysisPsi2S.h"
#include <TH2F.h>
#include <TH1I.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>

#include <TMath.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <Math/GenVector/VectorUtil.h>
#include <Math/GenVector/PtEtaPhiM4D.h>

//#define GRAPHYCS_

void PrepAnalysisPsi2S::Loop(const TString& which_analysis)
{
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  const Long64_t nbreak = nentries + 200;//100000;
  const ULong64_t entry_CP = 2000;
  Int_t code = analysis_code(which_analysis);
  
  //BRANCHES
  fChain->SetBranchStatus("*",0); 

 
  fChain->SetBranchStatus("nGenPart",1);
  fChain->SetBranchStatus("GenPart_pdgId",1);
  fChain->SetBranchStatus("GenPart_genPartIdxMother",1);
  fChain->SetBranchStatus("GenPart_eta",1);
  fChain->SetBranchStatus("GenPart_phi",1);
  fChain->SetBranchStatus("GenPart_pt",1);
  fChain->SetBranchStatus("GenPart_mass",1);

  if(code == 0){
    fChain->SetBranchStatus("nB0",1);
    fChain->SetBranchStatus("nK0s",1);
    fChain->SetBranchStatus("nJPsiToMuMu",1);
    fChain->SetBranchStatus("npipi",1);    
    fChain->SetBranchStatus("nMuon",1);
    fChain->SetBranchStatus("nProbeTracks",1);
  }

  if(code == 1){
    
    fChain->SetBranchStatus("nMuon",1);
    fChain->SetBranchStatus("Muon_eta",1);
    fChain->SetBranchStatus("Muon_phi",1);
    fChain->SetBranchStatus("Muon_charge",1);
    fChain->SetBranchStatus("Muon_pt",1);
    fChain->SetBranchStatus("Muon_softId",1);
    fChain->SetBranchStatus("Muon_looseId",1);
    fChain->SetBranchStatus("Muon_isGlobal",1);
    fChain->SetBranchStatus("Muon_genPartIdx",1);
    fChain->SetBranchStatus("Muon_genPartFlav",1);
  } //JPsi branches
  
  
  if(code == 2){

    fChain->SetBranchStatus("nProbeTracks",1);
    fChain->SetBranchStatus("ProbeTracks_eta",1);
    fChain->SetBranchStatus("ProbeTracks_phi",1);
    fChain->SetBranchStatus("ProbeTracks_pt",1);
    fChain->SetBranchStatus("ProbeTracks_charge",1);
    fChain->SetBranchStatus("ProbeTracks_genPartIdx",1);
    fChain->SetBranchStatus("ProbeTracks_genPartFlav",1);
    fChain->SetBranchStatus("ProbeTracks_isMatchedToMuon",1); // is Matched to any Muon!
    fChain->SetBranchStatus("ProbeTracks_isMatchedToLooseMuon",1);
    fChain->SetBranchStatus("ProbeTracks_isMatchedToSoftMuon",1);
  }
  
  if(code == 3){
    fChain->SetBranchStatus("nK0s",1);
    fChain->SetBranchStatus("K0s_fitted_phi",1);
    fChain->SetBranchStatus("K0s_fitted_eta",1);
    fChain->SetBranchStatus("K0s_fitted_pt",1);
    fChain->SetBranchStatus("K0s_prefit_mass",1);
    fChain->SetBranchStatus("K0s_fitted_mass_womc",1);
    fChain->SetBranchStatus("K0s_fitted_mass",1);  
  }

 
  //VARIABLES

  const Int_t isMu = 13, isPi = 211, isPsi2S = 9920443, isKs = 310, isB0 = 511 ;
  double Delta_Rm, Delta_Rp, DRp_min, DRm_min;
  int DRm_min_idx = 0, DRp_min_idx = 0, DRm_min2_idx = 0, DRp_min2_idx = 0;
  bool isOppositeCharge = false;

  //---- JPsi ---- // 
  const Int_t isJPsi = 443;
  double JPsi_mass;
  ULong64_t nGenMuon = 0;

  const double m_mu = 0.105658;
  ROOT::Math::RhoEtaPhiVector Mum_vect3D(1., 0.,0.), Mup_vect3D(1., 0.,0.);
  ROOT::Math::PtEtaPhiMVector gen_Mum_vect4D(0.,0.,0., m_mu), gen_Mup_vect4D(0., 0., 0., m_mu);
  double score_pt_mum, score_pt_mup;
  int MuMulty =0, TotMu = 0, TotJPsi = 0; 
  double TotSoftId =0., TotGlobal = 0.;
  bool isRecoMum = false , isRecoMup = false ;
   
  //---- Rho ---- // 
  const Int_t isRho = 113;
  double Rho_mass;
  ULong64_t TotRhoToPiPi = 0, TotTracks = 0;
  const double m_son = 0.1395704; //pion mass
  ROOT::Math::PtEtaPhiMVector gen_Sonm_vect4D(0.,0.,0., m_son), gen_Sonp_vect4D(0.,0.,0., m_son);
  ROOT::Math::PtEtaPhiMVector Sonm_vect4D(0.,0.,0., m_son), Sonp_vect4D(0.,0.,0., m_son); 
  double score_pt_pim, score_pt_pip;
  bool isRecoPim = false , isRecoPip = false ;
  int PiMulty = 0;
  int isSgnBkgNone = -2; //Sgn = 1 Bkg = 0 None = -1
  

  //---- K0s ----//
  const int isK0s = 310;
  const double m_K0 = 0.497648;
  ROOT::Math::PtEtaPhiMVector gen_K0_vect4D(0.,0.,0., m_K0), K0_vect4D(0.,0.,0., m_K0);
  double DRKmin = 100., score_pt_K0;
  int TotK0s = 0, DRKmin_idx = 0;  
  bool isRecoK0 = false;
  

  //HISTOGRAMS

  //----Generator----//
  
  Int_t Nbins = 40;
  TH1F h_Gpt_mu("gen_pt_mu","", Nbins, 0.,20);
  TH1F h_GptL_mu("gen_Lpt_mu","", Nbins, 0.,20);
  TH1F h_GptSL_mu("gen_SLpt_mu","", Nbins, 0.,20);
  TH1F h_Gpt_pi("gen_pt_pi","", Nbins, 0.,6);
  TH1F h_Gpt_JPsi("gen_pt_JPsi","", Nbins, 5.,35);
  TH1F h_Gpt_Rho("gen_pt_Rho","", Nbins, 0.,10.);
  TH1F h_Gpt_Psi2S("gen_pt_Psi2S","", Nbins, 5., 50.);
  TH1F h_Gpt_Ks("gen_pt_Ks","", Nbins, 0.,20);
  TH1F h_Gpt_B("gen_pt_B","", Nbins, 0.,50);

  TH1F h_Geta_mu("gen_eta_mu","", Nbins, -3.5,3.5);
  TH1F h_Geta_pi("gen_eta_pi","", Nbins, -3.5,3.5);
  TH1F h_Geta_JPsi("gen_eta_JPsi","", Nbins, -3.5,3.5);
  TH1F h_Geta_Rho("gen_eta_Rho","", Nbins, -3.5,3.5);
  TH1F h_Geta_Psi2S("gen_eta_Psi2S","", Nbins, -3.5,3.5);
  TH1F h_Geta_Ks("gen_eta_Ks","", Nbins, -3.5,3.5);
  TH1F h_Geta_B("gen_eta_B","", Nbins, -3.5,3.5);

  TH1F h_Gm_mu("gen_m_mu","", Nbins, .05, .15);
  TH1F h_Gm_pi("gen_m_pi","", Nbins, .1, .2);
  TH1F h_Gm_JPsi("gen_m_JPsi","", Nbins, 2.9,3.2);
  TH1F h_Gm_Rho("gen_m_Rho","", Nbins, 0.2 ,1.);
  TH1F h_Gm_Psi2S("gen_m_Psi2S","", Nbins, 3.8, 4.);
  TH1F h_Gm_Ks("gen_m_Ks","", Nbins, 0.45, 0.55);
  TH1F h_Gm_B("gen_m_B","", Nbins, 5., 5.5 );


  //----Multiplicity----//
  TH1F h_nB0("nB0", "", 20, 0, 20);
  TH1F h_nK0s("nK0s", "", 10, 0, 10);
  TH1F h_nJPsi("nJPsi", "", 3, 0, 3);
  TH1F h_nPiPi("nPiPi", "", 40, 500, 10000);
  TH1F h_nMuon("nMuon", "", 15, 0, 15);
  TH1F h_nTracks("nTracks", "", 50, 100, 800);

  //---- JPsi ---- //
  int nbins = 50;
  
  TH2F  h_DRmin_vs_DpT_mu("DRmin_vs_DpT_mu", "", nbins, 0.0, 0.05 , nbins, 0., 1.);

  // ... dal generatore ...//
  TH1F  h_pt_mu_fromJPsi("pT_mu_fromJPsi" ,"", nbins, 0., 20.);
  TH1F  h_eta_mu_fromJPsi("Eta_mu_fromJPsi" ,"", nbins, -3.1, 3.1);
  // ... other muons ...//
  TH1F  h_pt_mu_Disc("pT_mu_Disc" ,"", nbins, 0., 20.);
  TH1F  h_eta_mu_Disc("Eta_mu_Disc" ,"", nbins, -3.1, 3.1);
  // ... those having a gen track associated ...//
  TH1F  h_pt_mu_MatchGen("pT_mu_MatchGen" ,"", nbins, 0., 20.);
  TH1F  h_eta_mu_MatchGen("Eta_mu_MatchGen" ,"", nbins, -3.1, 3.1);
  TH1F  h_N_mu_MatchGen("N_mu_MatchGen", "", 3,0, 3);
  // ... mia selezione ...//
  TH1F  h_pt_mu_Rec("pT_mu_Rec" ,"", nbins, 0., 20.);
  TH1F  h_eta_mu_Rec("Eta_mu_Rec" ,"", nbins, -3.1, 3.1);  
  TH1F  h_MuMulty("MuMulty","",10, 0, 10); 
  // --> Quality check
  TH1F  h_rec_muCharge("rec_muCharge","", 3, -0.25, 1.25);
  TH1F  h_rec_muTrk("rec_muTrk", "", 3, -0.25, 1.25);
  TH1F  h_rec_mumu("rec_mumu", "", 3, -0.25, 1.25);
  TH1F  h_rec_Nmu("rec_Nmu", "", 3, 0, 3);
  
  TH1F  h_MuSoftId_Reco("MuSoftId_Reco", "", 2, -0.5, 1.5);
  TH1F  h_MuLooseId_Reco("MuLooseId_Reco", "", 2, -0.5, 1.5);
  TH1F  h_MuIsGlobal_Reco("MuIsGlobal_Reco", "", 2, -0.5, 1.5);
  TH1F  h_MuSoftId_Disc("MuSoftId_Disc", "", 2, -0.5, 1.5);
  TH1F  h_MuLooseId_Disc("MuLooseId_Disc", "", 2, -0.5, 1.5);
  TH1F  h_MuIsGlobal_Disc("MuIsGlobal_Disc", "", 2, -0.5, 1.5);

  // ... leading & subleading mu ...//
  TH1F  h_lead_mu_DRmin("pT_leading_mu", "", nbins, 0., 20);
  TH1F  h_sublead_mu_DRmin("pT_subleading_mu", "", nbins, 0., 20);
  // ... inv. mass ...//
  TH1F  h_mumu_mass_Gen("Mmumu_Gen","", nbins, 2.9, 3.3);
  TH1F  h_mumu_mass_DRmin("Mmumu_DRmin","", nbins, 2.9, 3.3);
  
  //---- Rho ---- // 
  nbins = 50;
  
  TH2F  h_DRmin_vs_DpT("DRmin_vs_DpT", "", 50, 0.0, 0.1 , 50, 0., 1.5);

  // ... dal generatore ...//
  TH1F  h_pt_pi_fromRho("pT_pi_fromRho" ,"", nbins, 0., 5.);
  TH1F  h_eta_pi_fromRho("Eta_pi_fromRho" ,"", nbins, -3.1, 3.1);
  // ... tutte le tracce ...//
  TH1F  h_pt_trk("pT_trk" ,"", nbins, 0., 5.);
  TH1F  h_eta_trk("Eta_trk" ,"", nbins, -3.1, 3.1);
  // ... mia selezione ...//
  TH1F  h_pt_pi_Rec("pT_pi_Rec" ,"", nbins, 0., 5.);
  TH1F  h_eta_pi_Rec("Eta_pi_Rec" ,"", nbins, -3.1, 3.1);
  TH1F  h_PiMulty("PiMulty","",80, 0, 800); 
  
  TH1F  h_rec_piCharge("rec_piCharge","", 3, -0.25, 1.25);
  TH1F  h_rec_piTrk("rec_piTrk", "", 3, -0.25, 1.25);
  TH1F  h_rec_pipi("rec_pipi", "", 3, -0.25, 1.25);
  TH1F  h_rec_Npi("rec_Npi", "", 3, 0, 3);
  //---> Quality Check
  TH1F  h_TrkIsSoftMu_Reco("TrkIsSoftMu_Reco", "", 2, -0.5, 1.5);
  TH1F  h_TrkIsLooseMu_Reco("TrkIsLooseMu_Reco", "", 2, -0.5, 1.5);
  TH1F  h_TrkIsMu_Reco("TrkIsMu_Reco", "", 2, -0.5, 1.5);
  TH1F  h_TrkIsSoftMu_Disc("TrkIsSoftMu_Disc", "", 2, -0.5, 1.5);
  TH1F  h_TrkIsLooseMu_Disc("TrkIsLooseMu_Disc", "", 2, -0.5, 1.5);
  TH1F  h_TrkIsMu_Disc("TrkIsMu_Disc", "", 2, -0.5, 1.5);

  // ... inv mass ...//
  double low = 0.4 , high = 1.;
  TH1F  h_pipi_mass_Gen("Mpipi_Gen","", nbins,low , high);
  TH1F  h_pipi_mass_DRmin("Mpipi_DRmin","", nbins, low, high);
  TH1F  h_pipi_mass_Sel2("Mpipi_Sel2","", nbins, low, high);


  //---- K0 short ---- //
  nbins = 50;

  TH2F  h_DRmin_vs_DpT_K0("DRmin_vs_DpT_K0", "", 50, 0.0, 0.1 , 50, 0., 1.5);

  // ... dal generatore ...//
  TH1F  h_pt_K_fromB("pT_K_fromB" ,"", nbins, 0., 20.);
  TH1F  h_eta_K_fromB("Eta_K_fromB" ,"", nbins, -3.1, 3.1);
  // ... tutte le tracce ...//
  TH1F  h_pt_K_Disc("pT_Kdisc" ,"", nbins, 0., 20.);  
  TH1F  h_eta_K_Disc("Eta_Kdisc" ,"", nbins, -3.1, 3.1);
  // ... mia selezione ...//
  TH1F  h_pt_K_Rec("pT_K_Rec" ,"", nbins, 0., 20.);
  TH1F  h_eta_K_Rec("Eta_K_Rec" ,"", nbins, -3.1, 3.1);
 
  TH1F  h_rec_K("rec_K", "", 3, -0.25, 1.25);
  TH1F  h_K0sMulty("K0sMulty","",10, 0, 10); 

  // ... inv mass ...//
  low = 0.45; high = 0.55;
  TH1F  h_K0s_mass_Gen("MK0s_Gen","", nbins, low, high);
  TH1F  h_K0s_mass_prefit("MK0s_DRmin_prefit","", nbins, low, high);
  TH1F  h_K0s_mass_postfit_womc("MK0s_DRmin_fit_womc","", nbins, low, high);
  TH1F  h_K0s_mass_fitted("MK0s_DRmin_fit","", nbins, low, high);

  //..... Loop .....//

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (jentry == nbreak) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    if((jentry +1 ) % entry_CP == 0) std::cout<< "--> Processing event number " << jentry + 1 << std::endl;
    
    if (code == 0){
      // ------- GENERATOR -------//
		GenPartFillP4();
      GenPart_FillKinHist( &GenP4_Mum, &h_Gpt_mu, &h_Geta_mu, &h_Gm_mu);
      GenPart_FillKinHist( &GenP4_Mup, &h_Gpt_mu, &h_Geta_mu, &h_Gm_mu);
      GenPart_FillKinHist( &GenP4_JPsi ,&h_Gpt_JPsi, &h_Geta_JPsi, &h_Gm_JPsi);
      GenPart_FillKinHist( &GenP4_Pim, &h_Gpt_pi, &h_Geta_pi, &h_Gm_pi);
      GenPart_FillKinHist( &GenP4_Pip, &h_Gpt_pi, &h_Geta_pi, &h_Gm_pi);
      GenPart_FillKinHist( &GenP4_Rho, &h_Gpt_Rho, &h_Geta_Rho, &h_Gm_Rho);
      GenPart_FillKinHist( &GenP4_Psi2S, &h_Gpt_Psi2S ,&h_Geta_Psi2S, &h_Gm_Psi2S);
      GenPart_FillKinHist( &GenP4_K0s, &h_Gpt_Ks, &h_Geta_Ks, &h_Gm_Ks);
      GenPart_FillKinHist( &GenP4_B0,&h_Gpt_B, &h_Geta_B, &h_Gm_B);
		if(GenP4_Mum.Pt() > GenP4_Mup.Pt()){
			h_GptL_mu.Fill(GenP4_Mum.Pt());
			h_GptSL_mu.Fill(GenP4_Mup.Pt());
		}else{
			h_GptL_mu.Fill(GenP4_Mup.Pt());
			h_GptSL_mu.Fill(GenP4_Mum.Pt());
		}
		//PrintDecayChain(); 
      TFile gen_outfile("./plots/analysisPsi2S_Gen_UL17.root","RECREATE"); // create the output file
      h_Gpt_mu.Write();
		h_GptL_mu.Write();
		h_GptSL_mu.Write();
      h_Geta_mu.Write();
      h_Gm_mu.Write();
      h_Gpt_pi.Write();
      h_Geta_pi.Write();
      h_Gm_pi.Write();
      h_Gpt_JPsi.Write();
      h_Geta_JPsi.Write();
      h_Gm_JPsi.Write();
      h_Gpt_Rho.Write();
      h_Geta_Rho.Write();
      h_Gm_Rho.Write();
      h_Gpt_Psi2S.Write();
      h_Geta_Psi2S.Write();
      h_Gm_Psi2S.Write();
      h_Gpt_Ks.Write();
      h_Geta_Ks.Write();
      h_Gm_Ks.Write();
      h_Gpt_B.Write();
      h_Geta_B.Write();
      h_Gm_B.Write();
      
      gen_outfile.Close();
      
      //----- MULTIPLICITY----//
      
      h_nB0.Fill(nB0);
      h_nK0s.Fill(nK0s);
      h_nJPsi.Fill(nJPsiToMuMu);
      h_nPiPi.Fill(npipi);
      h_nMuon.Fill(nMuon);
      h_nTracks.Fill(nProbeTracks);
      
      TFile mul_outfile("./plots/analysisPsi2S_Multiplicity_UL17.root","RECREATE"); // create the output file
      h_nB0.Write();
      h_nK0s.Write();
      h_nJPsi.Write();
      h_nPiPi.Write();
      h_nMuon.Write();
      h_nTracks.Write();
      
      mul_outfile.Close();
    
    } //----- GENERATOR & MULTIPLICITY -------//

	 if(code == 1){

		 //--- GENERATORE

		 gen_Daugh_ptetaphiM(isJPsi, &gen_Mum_vect4D, &gen_Mup_vect4D); //costruisco un quadrivettore dei muoni generati
		 h_pt_mu_fromJPsi.Fill(gen_Mum_vect4D.Pt());
		 h_eta_mu_fromJPsi.Fill(gen_Mum_vect4D.Eta());
		 h_pt_mu_fromJPsi.Fill(gen_Mup_vect4D.Pt());
		 h_eta_mu_fromJPsi.Fill(gen_Mup_vect4D.Eta());
		 h_mumu_mass_Gen.Fill((gen_Mum_vect4D + gen_Mup_vect4D).M());

		 DRp_min = 100.;
		 DRm_min = 100.;

		 MuMulty = 0;
		 TotMu = 0;
		 isRecoMum = false;
		 isRecoMup = false;

		 for(ULong64_t m = 0; m < nMuon; m++){
			 if (Muon_softId[m] == false) continue;
			 MuMulty++;	
			 //Delta R-
			 find_minDeltaR(m, gen_Mum_vect4D, &DRm_min, &DRm_min_idx, code);
			 //Delta R+
			 find_minDeltaR(m, gen_Mup_vect4D, &DRp_min, &DRp_min_idx, code);

		 }//on reconstructed Mu 
		 h_MuMulty.Fill(MuMulty);
		 // Delta pT @ DRmin
		 score_pt_mum = score_PT(DRm_min_idx, &gen_Mum_vect4D, code);
		 h_DRmin_vs_DpT_mu.Fill(DRm_min,score_pt_mum);
		 score_pt_mup = score_PT(DRp_min_idx, &gen_Mup_vect4D, code);
		 h_DRmin_vs_DpT_mu.Fill(DRp_min,score_pt_mup);

		 //check on mu-
		 isRecoMum = isReconstructedTrack(DRm_min_idx, DRm_min, score_pt_mum, code);
		 if (isRecoMum){
			 h_pt_mu_Rec.Fill(Muon_pt[DRm_min_idx]);
			 h_eta_mu_Rec.Fill(Muon_eta[DRm_min_idx]);
			 h_rec_muTrk.Fill(1.);

			 if(Muon_charge[DRm_min_idx] < 0 ){ h_rec_muCharge.Fill(1.);
			 }else{ h_rec_muCharge.Fill(0.);}

			 TotMu += 1;
			 if (Muon_softId[DRm_min_idx]) TotSoftId += 1.;
			 if (Muon_isGlobal[DRm_min_idx]) TotGlobal +=1.;
		 }else{ h_rec_muTrk.Fill(0.);}

		 // check sul mu+
		 isRecoMup = isReconstructedTrack(DRp_min_idx, DRp_min, score_pt_mup, code);
		 if (isRecoMup){
			 h_pt_mu_Rec.Fill(Muon_pt[DRp_min_idx]);
			 h_eta_mu_Rec.Fill(Muon_eta[DRp_min_idx]);
			 h_rec_muTrk.Fill(1.);
			 TotMu += 1;
			 if (Muon_softId[DRp_min_idx]) TotSoftId += 1.;
			 if (Muon_isGlobal[DRp_min_idx]) TotGlobal +=1.;

			 if(Muon_charge[DRp_min_idx] > 0 ){ h_rec_muCharge.Fill(1.);
			 }else{ h_rec_muCharge.Fill(0.);}

		 }else{ h_rec_muTrk.Fill(0.);}

		 h_rec_Nmu.Fill(TotMu);

		 //NON-signal muons
		 OtherMuons(DRm_min_idx, DRm_min, DRp_min_idx, DRp_min, isRecoMum, isRecoMup, &h_pt_mu_Disc, &h_eta_mu_Disc, &h_pt_mu_MatchGen, &h_eta_mu_MatchGen, &h_N_mu_MatchGen);

		 //check sulla coppia mu+mu-
		 isOppositeCharge = (Muon_charge[DRm_min_idx]* Muon_charge[DRp_min_idx]) < 0;
		 if ( isRecoMum && isRecoMup && isOppositeCharge){
			 h_mumu_mass_DRmin.Fill(MyInvMass(DRm_min_idx,DRp_min_idx, code));
			 h_rec_mumu.Fill(1.);
			 TotJPsi += 1;
		 }else{h_rec_mumu.Fill(0.);}

		 MuTrkQualityCheck(DRm_min_idx, DRp_min_idx, isRecoMum, isRecoMup, &h_MuSoftId_Reco, &h_MuSoftId_Disc, &h_MuLooseId_Reco, &h_MuLooseId_Disc, &h_MuIsGlobal_Reco, &h_MuIsGlobal_Disc);

		 //LEADING AND SUBLEADING BEST MUONS 
		 if (Muon_pt[DRm_min_idx] > Muon_pt[DRp_min_idx]){
			 h_lead_mu_DRmin.Fill(Muon_pt[DRm_min_idx]);
			 h_sublead_mu_DRmin.Fill(Muon_pt[DRp_min_idx]);
		 }else {
			 h_lead_mu_DRmin.Fill(Muon_pt[DRp_min_idx]);
			 h_sublead_mu_DRmin.Fill(Muon_pt[DRm_min_idx]);
		 }

	 }// --- JPsi ANALYSIS ---


	 if (code == 2) {


		 //...GENERATORE...
		 gen_Daugh_ptetaphiM(isRho, &gen_Sonm_vect4D, &gen_Sonp_vect4D); //costruisco un quadrivettore dei pioni della rho generati
		 h_pt_pi_fromRho.Fill(gen_Sonm_vect4D.Pt());
		 h_eta_pi_fromRho.Fill(gen_Sonm_vect4D.Eta());
		 h_pt_pi_fromRho.Fill(gen_Sonp_vect4D.Pt());
		 h_eta_pi_fromRho.Fill(gen_Sonp_vect4D.Eta());
		 h_pipi_mass_Gen.Fill((gen_Sonm_vect4D + gen_Sonp_vect4D).M());

		 DRp_min = 100.;
		 DRm_min = 100.;

		 PiMulty = 0;
		 TotTracks = 0;

		 for(ULong64_t m = 0; m < nProbeTracks; m++){
			 if (ProbeTracks_isMatchedToMuon[m]) continue;
			 //salvo tutte le tracce
			 PiMulty++;
			 h_pt_trk.Fill(ProbeTracks_pt[m]);
			 h_eta_trk.Fill(ProbeTracks_eta[m]);

			 //Delta R-
			 find_minDeltaR(m, gen_Sonm_vect4D, &DRm_min, &DRm_min_idx, code);
			 //Delta R+
			 find_minDeltaR(m, gen_Sonp_vect4D, &DRp_min, &DRp_min_idx, code);
		 } // on tracks
		 h_PiMulty.Fill(PiMulty);

		 // Delta pT @ DR minimo
		 score_pt_pim = score_PT(DRm_min_idx, &gen_Sonm_vect4D, code);
		 h_DRmin_vs_DpT.Fill(DRm_min,score_pt_pim);
		 score_pt_pip = score_PT(DRp_min_idx, &gen_Sonp_vect4D, code);
		 h_DRmin_vs_DpT.Fill(DRp_min,score_pt_pip);

		 //check sulla traccia -
		 isRecoPim = isReconstructedTrack(DRm_min_idx, DRm_min, score_pt_pim, code);
		 if (isRecoPim && (gen_Sonm_vect4D.Pt() > 0.5)){
			 h_pt_pi_Rec.Fill(ProbeTracks_pt[DRm_min_idx]);
			 h_eta_pi_Rec.Fill(ProbeTracks_eta[DRm_min_idx]);
			 h_rec_piTrk.Fill(1.);

			 TotTracks += 1;
			 if(ProbeTracks_charge[DRm_min_idx] < 0 ){ h_rec_piCharge.Fill(1.);
			 }else{ h_rec_piCharge.Fill(0.);}
		 }else if(gen_Sonm_vect4D.Pt() > 0.5) { h_rec_piTrk.Fill(0.);}
		 //check sulla traccia +
		 isRecoPip = isReconstructedTrack(DRp_min_idx, DRp_min, score_pt_pip, code);
			 if (isRecoPip && (gen_Sonp_vect4D.Pt() > 0.5)){
				 h_pt_pi_Rec.Fill(ProbeTracks_pt[DRp_min_idx]);
				 h_eta_pi_Rec.Fill(ProbeTracks_eta[DRp_min_idx]);
				 h_rec_piTrk.Fill(1.);
				 TotTracks += 1;

				 if(ProbeTracks_charge[DRp_min_idx] > 0 ){ h_rec_piCharge.Fill(1.);
				 }else{ h_rec_piCharge.Fill(0.);}
			 }else if (gen_Sonp_vect4D.Pt() > 0.5){ h_rec_piTrk.Fill(0.);}
		 h_rec_Npi.Fill(TotTracks);
		 PiTrkQualityCheck( DRm_min_idx, DRp_min_idx, isRecoPim, isRecoPip, &h_TrkIsSoftMu_Reco, &h_TrkIsSoftMu_Disc, &h_TrkIsLooseMu_Reco, &h_TrkIsLooseMu_Disc,&h_TrkIsMu_Reco, &h_TrkIsMu_Disc);

		 //check sulla coppia pi+pi-
		 if ( isRecoPim && isRecoPip){
			 h_pipi_mass_DRmin.Fill(MyInvMass(DRm_min_idx,DRp_min_idx, code));
			 h_rec_pipi.Fill(1.);
			 TotRhoToPiPi += 1;
		 }else{ h_rec_pipi.Fill(0.);}

		 isRecoPim = (DRm_min < 0.1) && (score_pt_pim < 0.5);
		 isRecoPip = (DRp_min < 0.1) && (score_pt_pip < 0.5);
		 if(isRecoPim && isRecoPip) h_pipi_mass_Sel2.Fill(MyInvMass(DRm_min_idx,DRp_min_idx, code));

	 }// --- RHO ANALYSIS ---




	 if(code == 3){
		 //...GENERATORE...
		 gen_K0s_ptetaphi(isB0, &gen_K0_vect4D); //costruisco un quadrivettore dei K0short dal B0 generati
		 h_pt_K_fromB.Fill(gen_K0_vect4D.Pt());
		 h_eta_K_fromB.Fill(gen_K0_vect4D.Eta());
		 h_K0s_mass_Gen.Fill(gen_K0_vect4D.M());
		 h_K0sMulty.Fill(nK0s);
		 DRKmin = 100.;

		 for(UInt_t k = 0; k < nK0s; k++){
			 //... salvo tutti i K0 ...//
			 h_pt_K_Disc.Fill(K0s_fitted_pt[k]);
			 h_eta_K_Disc.Fill(K0s_fitted_eta[k]);

			 //DeltaR
			 find_minDeltaR(k, gen_K0_vect4D, &DRKmin, &DRKmin_idx, code);

		 }//on all K0s

		 // Delta pT @ DR minimo
		 score_pt_K0 = score_PT(DRKmin_idx, &gen_K0_vect4D, code);
		 h_DRmin_vs_DpT_K0.Fill(DRKmin,score_pt_K0);
		 //check on K
		 isRecoK0 = isReconstructedTrack(DRKmin_idx, DRKmin, score_pt_K0, code);
		 if (isRecoK0){
			 h_pt_K_Rec.Fill(K0s_fitted_pt[DRKmin_idx]);
			 h_eta_K_Rec.Fill(K0s_fitted_eta[DRKmin_idx]);
			 h_rec_K.Fill(1.);

			 h_K0s_mass_prefit.Fill(K0s_prefit_mass[DRKmin_idx]);
			 h_K0s_mass_postfit_womc.Fill(K0s_fitted_mass_womc[DRKmin_idx]);
			 h_K0s_mass_fitted.Fill(K0s_fitted_mass[DRKmin_idx]);

			 TotK0s++;
		 }else h_rec_K.Fill(0.);

	 }// --- K0s ANALYSIS --- //










  }//on events
  
  if(code == 1){
    
    TotSoftId /= h_pt_mu_Rec.GetEntries();
    TotGlobal /= h_pt_mu_Rec.GetEntries();

    std::cout << " Total number of reconstructed JPsi " << TotJPsi << std::endl;
    std::cout << " Total number of muon track with true SoftId " << TotSoftId << "\t with TotGlobal "<< TotGlobal << std::endl;

    //h_pt_mu_Disc.Add(&h_pt_mu_Rec, -1);
    //h_eta_mu_Disc.Add(&h_eta_mu_Rec, -1);
    
    TFile outfile("./plots/analysis_MuonsREDO.root","RECREATE"); // create the output file
    h_DRmin_vs_DpT_mu.Write();
    h_MuMulty.Write(); 
 
    h_pt_mu_fromJPsi.Write();
    h_eta_mu_fromJPsi.Write();
    h_pt_mu_Disc.Write();
    h_eta_mu_Disc.Write();
    h_pt_mu_MatchGen.Write();
    h_eta_mu_MatchGen.Write();
    h_N_mu_MatchGen.Write();
    h_pt_mu_Rec.Write();
    h_eta_mu_Rec.Write();

    h_rec_muCharge.Write();
    h_rec_muTrk.Write();
    h_rec_mumu.Write();
    h_rec_Nmu.Write();

    h_MuSoftId_Reco.Write();
    h_MuLooseId_Reco.Write();
    h_MuIsGlobal_Reco.Write();
    h_MuSoftId_Disc.Write();
    h_MuLooseId_Disc.Write();
    h_MuIsGlobal_Disc.Write();

    h_mumu_mass_Gen.Write();
    h_mumu_mass_DRmin.Write();

    h_lead_mu_DRmin.Write();
    h_sublead_mu_DRmin.Write();
    outfile.Close();
  }

  if(code == 2){
    
    std::cout << " Total number of reconstructed Rho to pipi " << TotRhoToPiPi << std::endl;
    
    h_pt_trk.Add(&h_pt_pi_Rec, -1);
    h_eta_trk.Add(&h_eta_pi_Rec, -1);
    
    TFile outfile("./plots/analysis_TracksREDO.root","RECREATE"); // create the output file
    h_DRmin_vs_DpT.Write();

    h_pt_trk.Write();
    h_pt_pi_fromRho.Write();
    h_pt_pi_Rec.Write();
    h_eta_trk.Write();
    h_eta_pi_fromRho.Write();
    h_eta_pi_Rec.Write();
	 h_PiMulty.Write();
        
    h_rec_piCharge.Write();
    h_rec_piTrk.Write();
    h_rec_pipi.Write();
    h_rec_Npi.Write();

    h_TrkIsSoftMu_Reco.Write();
    h_TrkIsLooseMu_Reco.Write();
    h_TrkIsMu_Reco.Write();
    h_TrkIsSoftMu_Disc.Write();
    h_TrkIsLooseMu_Disc.Write();
    h_TrkIsMu_Disc.Write();

    h_pipi_mass_Gen.Write();
    h_pipi_mass_DRmin.Write();
    h_pipi_mass_Sel2.Write();

    outfile.Close();
  }


if(code == 3){
  h_pt_K_Disc.Add(&h_pt_K_Rec, -1);
  h_eta_K_Disc.Add(&h_eta_K_Rec, -1);
  
  std::cout << " Total number of reconstructed K0 " << TotK0s << std::endl;

  TFile outfile("./plots/analysis_K0REDO.root","RECREATE"); // create the output file
  h_DRmin_vs_DpT_K0.Write();

  h_pt_K_fromB.Write();
  h_pt_K_Disc.Write();
  h_pt_K_Rec.Write();
  h_eta_K_fromB.Write();
  h_eta_K_Disc.Write();
  h_eta_K_Rec.Write();

  h_rec_K.Write();
  h_K0sMulty.Write();

  h_K0s_mass_Gen.Write();
  h_K0s_mass_prefit.Write();
  h_K0s_mass_postfit_womc.Write();
  h_K0s_mass_fitted.Write();

}




}// Loop()



void PrepAnalysisPsi2S::gen_K0s_ptetaphi(const int &Mother_pdgId, ROOT::Math::PtEtaPhiMVector* V){
  
  const int isK0s = 310;
  for(ULong64_t i = 0; i < nGenPart ; i++){
    if(abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) != Mother_pdgId) continue;
    
    V->SetPt(GenPart_pt[i]);    
    V->SetEta(GenPart_eta[i]);	
    V->SetPhi(GenPart_phi[i]);     

      
    }//loop on gen ptls
}// gen_K0s_ptetaphi 


void PrepAnalysisPsi2S::GenPart_FillKinHist(ROOT::Math::PtEtaPhiMVector* GenVec, TH1* h_pt, TH1* h_eta, TH1* h_mass){ 

	h_pt->Fill(GenVec->Pt());
	h_eta->Fill(GenVec->Eta());
	h_mass->Fill(GenVec->M());

}//GenPart_FillKinHist()

void PrepAnalysisPsi2S::GenPartFillP4(){

	Int_t MumIdx = -1, MupIdx = -1, PimIdx = -1 , PipIdx = -1, JPsiIdx = -1, Psi2SIdx = -1, K0sIdx = -1, RhoIdx = -1, B0Idx = -1;
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
			 (abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]) == isB0) )//&&
			 //(abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]]) == isB0) )
			{
			RhoIdx = GenPart_genPartIdxMother[g];
			if (GenPart_pdgId[g] == -isPip) PimIdx = g;
			if (GenPart_pdgId[g] ==  isPip) PipIdx = g;
			}
      if( (GenPart_pdgId[g] ==  isK0s) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == abs(isB0))) K0sIdx = g;


   }//on generated ptl

	if(0){
		std::cout << " IDX Mu-  / Mu+  " << MumIdx  << "\t" << MupIdx << std::endl; 
		std::cout << " IDX Pi-  / Pi+  " << PimIdx  << "\t" << PipIdx << std::endl; 
		std::cout << " IDX JPsi / Rho  " << JPsiIdx << "\t" << RhoIdx << std::endl; 
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

	if(RhoIdx > 0){   
      GenP4_Rho.SetPt(GenPart_pt[RhoIdx]);   GenP4_Rho.SetEta(GenPart_eta[RhoIdx]);
		GenP4_Rho.SetPhi(GenPart_phi[RhoIdx]); GenP4_Rho.SetM(GenPart_mass[RhoIdx]); 
	} else std::cout << "No generated rho is found" << std::endl;  
	
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


void PrepAnalysisPsi2S::gen_Daugh_ptetaphiM( const Int_t& Mother_pdgId,ROOT::Math::PtEtaPhiMVector* V_neg,ROOT::Math::PtEtaPhiMVector* V_pos ){

	int Daugh_neg_pdgId =0;
	if (Mother_pdgId == 443) Daugh_neg_pdgId = 13; //JPsitoMuMu
	if (Mother_pdgId == 113) Daugh_neg_pdgId = -211; //RhoToPiPi
	//  std::cout << "Number of generated ptls "<< nGenPart << std::endl;
	for(ULong64_t i = 0; i < nGenPart ; i++){
		if(GenPart_pdgId[GenPart_genPartIdxMother[i]] == Mother_pdgId){

			if(GenPart_pdgId[i] == Daugh_neg_pdgId){ //check ptl-                                                                                                                               
				//    std::cout << i << " particle is ptl-" << std::endl;                                                                                                              
				V_neg->SetEta(GenPart_eta[i]);
				V_neg->SetPhi(GenPart_phi[i]);
				V_neg->SetPt(GenPart_pt[i]);
				//	V_neg->SetPt(GenPart_mass[i]);

			}else if (GenPart_pdgId[i] == -Daugh_neg_pdgId) {  //check ptl+

				//std::cout << i << " particle is ptl+" << std::endl;                                                                                                             
				V_pos->SetEta(GenPart_eta[i]);
				V_pos->SetPhi(GenPart_phi[i]);
				V_pos->SetPt(GenPart_pt[i]);
				//V_pos->SetPt(GenPart_mass[i]);
			}
		}

	}//loop on gen ptls    
}// gen_Daugh_ptetaphiM()



void PrepAnalysisPsi2S::find_minDeltaR(const int& Trk_idx, const ROOT::Math::PtEtaPhiMVector& gen4V, double* DR_min, int* DR_min_idx, const int analysis_code){

  ROOT::Math::PtEtaPhiMVector rec4V(0.,0.,0.,0.);
  if (analysis_code == 1){
    rec4V.SetPt(Muon_pt[Trk_idx]);
    rec4V.SetEta(Muon_eta[Trk_idx]);
    rec4V.SetPhi(Muon_phi[Trk_idx]);
  
  }
  if (analysis_code == 2){
    rec4V.SetPt(ProbeTracks_pt[Trk_idx]);
    rec4V.SetEta(ProbeTracks_eta[Trk_idx]);
    rec4V.SetPhi(ProbeTracks_phi[Trk_idx]);
  
  }
  if (analysis_code == 3){
    rec4V.SetPt(K0s_fitted_pt[Trk_idx]);
    rec4V.SetEta(K0s_fitted_eta[Trk_idx]);
    rec4V.SetPhi(K0s_fitted_phi[Trk_idx]);

  } 


  //DeltaR
  double Delta_R = ROOT::Math::VectorUtil::DeltaR (gen4V, rec4V);
  //...cerco il minimo                                                                                                                                
  if (Delta_R < *DR_min ){
    *DR_min = Delta_R;
    *DR_min_idx = Trk_idx;
  }

}//find_minDeltaR()


double PrepAnalysisPsi2S::score_PT(const int track_idx, const ROOT::Math::PtEtaPhiMVector* gen4V, const int analysis_code){

  double rec_pt = 0; // dipende se mu o pi
  if (analysis_code == 1) rec_pt = Muon_pt[track_idx];
  if (analysis_code == 2) rec_pt = ProbeTracks_pt[track_idx];
  if (analysis_code == 3) rec_pt = K0s_fitted_pt[track_idx];
  
  double score = TMath::Abs(gen4V->Pt() - rec_pt)/gen4V->Pt();

  return score;
}

bool PrepAnalysisPsi2S::isReconstructedTrack(const int& DR_min_idx ,const double& DR_min, const double& DpT_DRmin, const int& analysis_code){
  float DR_threshold = 0.03;
  float DpT_threshold = 0.5;
  if (analysis_code == 1) DpT_threshold = 100.; // for the JPsi analysis no constrain on DeltaPT
  bool track_score = false;

  if ((DR_min < DR_threshold)&& (DpT_DRmin < DpT_threshold)) track_score = true;

  return track_score;

}//isReconstructedTrack

void PrepAnalysisPsi2S::OtherMuons(const UInt_t& DR_neg_idx, const double& DRmin_neg, const UInt_t& DR_pos_idx,const double& DRmin_pos, const bool& isReco_neg, const bool& isReco_pos, TH1* h_pt_mu, TH1* h_eta_mu, TH1* h_pt_mu_MatchedGen, TH1* h_eta_mu_MatchedGen, TH1* h_N_mu_MatchedGen){
  int NmatchedGen = 0;
  if (!isReco_neg && DRmin_neg > 0.05){
    h_pt_mu->Fill(Muon_pt[DR_neg_idx]);
    h_eta_mu->Fill(Muon_eta[DR_neg_idx]);
  }
  if (!isReco_pos && DRmin_pos > 0.05){
    h_pt_mu->Fill(Muon_pt[DR_neg_idx]);
    h_eta_mu->Fill(Muon_eta[DR_neg_idx]);
  }
  for(UInt_t m = 0; m < nMuon; m++){
    if ( (m == DR_neg_idx) || (m == DR_pos_idx)) continue;
    if (Muon_softId[m] == false) continue;
    h_pt_mu->Fill(Muon_pt[m]);
    h_eta_mu->Fill(Muon_eta[m]);

    if (Muon_genPartIdx[m] != -1){
      //std::cout << Muon_genPartFlav[m] << std::endl;
      h_pt_mu_MatchedGen->Fill(Muon_pt[m]);
      h_eta_mu_MatchedGen->Fill(Muon_eta[m]);
      NmatchedGen++;
    }
  }
  h_N_mu_MatchedGen->Fill(NmatchedGen);
}//flagMuons()

void PrepAnalysisPsi2S::MuTrkQualityCheck(const UInt_t& DR_neg_idx, const UInt_t& DR_pos_idx, const bool& isReco_neg, const bool& isReco_pos, TH1* h_SoftId_Reco, TH1* h_SoftId_Disc, TH1* h_LooseId_Reco, TH1* h_LooseId_Disc,TH1* h_isGlobal_Reco, TH1* h_isGlobal_Disc){

  for(UInt_t m = 0; m < nMuon; m++ ){
    if( Muon_softId[m] == false ) continue;
    if( ((m == DR_neg_idx) && isReco_neg) || ((m == DR_pos_idx) && isReco_pos) ){
      h_SoftId_Reco->Fill((float)Muon_softId[m]);
      h_LooseId_Reco->Fill((float)Muon_looseId[m]);
      h_isGlobal_Reco->Fill((float)Muon_isGlobal[m]);
    }else{
      h_SoftId_Disc->Fill((float)Muon_softId[m]);
      h_LooseId_Disc->Fill((float)Muon_looseId[m]);
      h_isGlobal_Disc->Fill((float)Muon_isGlobal[m]);
    }
  }

}//MuTrackQualityCheck()



void PrepAnalysisPsi2S::PiTrkQualityCheck(const UInt_t& DR_neg_idx, const UInt_t& DR_pos_idx, const bool& isReco_neg, const bool& isReco_pos, TH1* h_isSoftMu_Reco, TH1* h_isSoftMu_Disc, TH1* h_isLooseMu_Reco, TH1* h_isLooseMu_Disc,TH1* h_isMu_Reco, TH1* h_isMu_Disc){

  for(UInt_t p = 0; p < nProbeTracks; p++ ){

    if( ((p == DR_neg_idx) && isReco_neg) || ((p == DR_pos_idx) && isReco_pos) ){
      h_isSoftMu_Reco->Fill((float)ProbeTracks_isMatchedToSoftMuon[p]);
      h_isLooseMu_Reco->Fill((float)ProbeTracks_isMatchedToLooseMuon[p]);
      h_isMu_Reco->Fill((float)ProbeTracks_isMatchedToMuon[p]);
    }else{
      h_isSoftMu_Disc->Fill((float)ProbeTracks_isMatchedToSoftMuon[p]);
      h_isLooseMu_Disc->Fill((float)ProbeTracks_isMatchedToLooseMuon[p]);
      h_isMu_Disc->Fill((float)ProbeTracks_isMatchedToMuon[p]);
    }
  }

}//PiTrackQualityCheck()


Int_t PrepAnalysisPsi2S::analysis_code(const TString& which_analysis){

  Int_t code = -1;
  if (which_analysis == "Intro") {
    code = 0;
    std::cout << "\n... Starting Generator & Multipicity analysis... \n" << std::endl;
  }
  if (which_analysis == "JPsi") {
    code = 1;
    std::cout << "\n... Starting JPsi analysis... \n" << std::endl;
  }
  if (which_analysis == "Rho"){
    code = 2;
    std::cout << "\n... Starting Rho analysis... \n" << std::endl;
  }
  if (which_analysis == "K0s"){
    code = 3;
    std::cout << "\n... Starting K-short analysis... \n" << std::endl;
  }
  return code;

}


double PrepAnalysisPsi2S::MyInvMass(const int idx1, const int idx2, const int& analysis_code){

  double mass = 0;
  ROOT::Math::PtEtaPhiMVector particle_1(0.,0.,0.,0.), particle_2(0.,0.,0.,0.), Ptot(0.,0.,0.,0.);
  
  if (analysis_code == 1){
    const double mmu = 0.105658; 
    particle_1.SetPt(Muon_pt[idx1]);
    particle_1.SetEta(Muon_eta[idx1]);
    particle_1.SetPhi(Muon_phi[idx1]);
    particle_1.SetM( mmu );
    particle_2.SetPt(Muon_pt[idx2]);
    particle_2.SetEta(Muon_eta[idx2]);
    particle_2.SetPhi(Muon_phi[idx2]);
    particle_2.SetM( mmu );
  }
  
  if (analysis_code == 2){
    const double mpi = 0.1395704; 
    particle_1.SetPt(ProbeTracks_pt[idx1]);
    particle_1.SetEta(ProbeTracks_eta[idx1]);
    particle_1.SetPhi(ProbeTracks_phi[idx1]);
    particle_1.SetM( mpi );
    particle_2.SetPt(ProbeTracks_pt[idx2]);
    particle_2.SetEta(ProbeTracks_eta[idx2]);
    particle_2.SetPhi(ProbeTracks_phi[idx2]);
    particle_2.SetM( mpi );
  }

  Ptot = particle_1 + particle_2;
  mass = Ptot.M();
  return mass;
}


void PrepAnalysisPsi2S::PrintDecayChain(){

	int GenPtl_PDGid =-1; 	
	std::cout << "#GEN-IDX \t  GEN PTL ID \t MOTHER IDX " << std::endl;
	for (ULong64_t g = 0; g < nGenPart; g++){
		GenPtl_PDGid = GenPart_pdgId[g];
	  
		std::cout<< "  "  << g << "  \t   " <<  GenPart_pdgId[g] << "  \t   " << GenPart_genPartIdxMother[g]<< std::endl;
// << "\t <--" << GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] << " <-- " << GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]] << std::endl;


		}
}
