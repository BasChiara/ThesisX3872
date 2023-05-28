#define MassAnalysis_cxx
#include "MassAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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

void MassAnalysis::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  const ULong64_t nbreak = nentries + 100;
  const ULong64_t entry_CP = 200;
  
  //BRANCHES                                                                                                                                                                                                        
  fChain->SetBranchStatus("*",0);
  // MC truth
  fChain->SetBranchStatus("nGenPart",1);
  fChain->SetBranchStatus("GenPart_pdgId",1);
  fChain->SetBranchStatus("GenPart_genPartIdxMother",1);
  fChain->SetBranchStatus("GenPart_eta",1);
  fChain->SetBranchStatus("GenPart_phi",1);
  fChain->SetBranchStatus("GenPart_pt",1);
  fChain->SetBranchStatus("GenPart_mass",1);
  // Muons
  fChain->SetBranchStatus("nMuon",1);
  fChain->SetBranchStatus("Muon_pt",1);
  fChain->SetBranchStatus("Muon_eta",1);
  fChain->SetBranchStatus("Muon_phi",1);
  fChain->SetBranchStatus("Muon_charge",1);
  fChain->SetBranchStatus("Muon_softId",1);
  //Tracks
  fChain->SetBranchStatus("nProbeTracks",1);
  fChain->SetBranchStatus("ProbeTracks_eta",1);
  fChain->SetBranchStatus("ProbeTracks_phi",1);
  fChain->SetBranchStatus("ProbeTracks_pt",1);
  fChain->SetBranchStatus("ProbeTracks_charge",1);
  fChain->SetBranchStatus("ProbeTracks_isMatchedToMuon",1);
  //Kaons
  fChain->SetBranchStatus("nK0s",1);
  fChain->SetBranchStatus("K0s_fitted_eta",1);
  fChain->SetBranchStatus("K0s_fitted_phi",1);
  fChain->SetBranchStatus("K0s_fitted_pt",1);
  
  //VARIABLES
  const int isMuon= 13, isPion = 211, isJPsi = 443, isRho = 113, isK0s = 310, isB0 = 511;
  const unsigned int JPsi_code = 1, Rho_code = 2, K0s_code = 3; 
  const double m_mu = 0.105658, m_pi = 0.1395704, m_K0s = 0.4976;
  const double mJPsiPDG = 3.096900;

  double Delta_Rm, Delta_Rp, DRp_min, DRm_min, Delta_Rk, DRk_min;
  int RecoMum_idx = 0, RecoMup_idx = 0, RecoPip_idx = 0, RecoPim_idx = 0, RecoK0s_idx = 0;
  ROOT::Math::PtEtaPhiMVector gen_Mum_vect4D(0.,0.,0., m_mu), gen_Mup_vect4D(0., 0., 0., m_mu);
  ROOT::Math::PtEtaPhiMVector Mum_vect4D(0.,0.,0., m_mu), Mup_vect4D(0., 0., 0., m_mu);
  ROOT::Math::PtEtaPhiMVector gen_Pim_vect4D(0.,0.,0., m_pi), gen_Pip_vect4D(0., 0., 0., m_pi);
  ROOT::Math::PtEtaPhiMVector Pim_vect4D(0.,0.,0., m_pi), Pip_vect4D(0., 0., 0., m_pi);
  ROOT::Math::PtEtaPhiMVector gen_K0s_vect4D(0.,0.,0., m_K0s);
  ROOT::Math::PtEtaPhiMVector K0s_vect4D(0.,0.,0., m_K0s);
  bool isRecoMum = false , isRecoMup = false ;
  bool isRecoPim = false , isRecoPip = false ;
  bool isRecoJPsi = false , isRecoRho = false;
  bool isRecoK0s = false;
  int Mum_idx[15], Mup_idx[15], totMum = 0, totMup = 0;
  int Pim_idx[1000], Pip_idx[1000], totPim = 0, totPip = 0;
  int nJPsi = 0, nJPsi_ev = 0, nRho = 0, nPip_CutOn = 0, nPim_CutOn = 0, nX3872 = 0, nK0s = 0, nB0 = 0;
  const double limMuMuPi_low = 3.1, limMuMuPi_high = 3.9;
  const double limXPt_low = 10. , limXPt_high = 14.;
  double Mcheck, PtCheck;
  const bool PTcut = false, Mcut = true;

  //HISTOGRAM
  int nbins = 50;
  float Mlow = 2.7 , Mhigh = 3.5;
  TH1F h_JPsi_M("JPsi_M", "", nbins, Mlow, Mhigh);
  TH1F h_MuMu_M("MuMu_M", "", nbins, Mlow, Mhigh);
  
  Mlow = .4 ; Mhigh = 1.;
  TH1F h_Rho_M("Rho_M", "", nbins, Mlow, Mhigh);
  TH1F h_PiPi_M("PiPi_M", "", nbins, Mlow, Mhigh);
  
  Mlow = 3.7 ; Mhigh = 4.05;
  TH1F h_X_M("X_M", "", nbins, Mlow, Mhigh);
  TH1F h_X_M_GEN("X_M_GEN", "", nbins, Mlow, Mhigh);
  TH1F h_JPsiPiPi_M("JPsiPiPi_M", "", nbins, Mlow, Mhigh);
  TH1F h_MuMuPiPi_M("MuMuPiPi_M", "", nbins, Mlow, Mhigh);
  TH1F h_JPsiPi_M("JPsiPi_M", "", nbins, 3.0, 4.);
  TH1F h_NPi_CutOn("NPi_CutOn", "", nbins, 0, 10000);
  
  TH1F h_X_pT("X_pT", "", nbins, 5., 25.);
  TH1F h_JPsiPiPi_pT("JPsiPiPi_pT", "", nbins, 5., 25.);
  //.... cut on |M(mu mu pi) - M(mu mu) - M(JPsi)| ....//
  TH1F h_Mdiff_JPsiPi_SGN("Mdiff_JPsiPi_SGN", "", nbins, 3.2, 4.2);
  TH1F h_Mdiff_JPsiPi_BKG("Mdiff_JPsiPi_BKG", "", nbins, 3.2, 4.2);
  //.... cut on |M(mu mu pi pi) - M(mu mu) - M(JPsi)| ....//    
  TH1F h_Mdiff_JPsiPiPi_SGN("Mdiff_JPsiPiPi_SGN", "", nbins, 3.6, 4.2);
  TH1F h_Mdiff_JPsiPiPi_BKG("Mdiff_JPsiPiPi_BKG", "", nbins, 3.6, 4.2);

	// B meson
	Mlow = 5.0, Mhigh = 5.6;
	TH1F h_B0_M_GEN("B0_M_GEN", "", nbins, Mlow, Mhigh);
	TH1F h_B0_M("B0_M", "", nbins, Mlow, Mhigh);
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
	  // if (Cut(ientry) < 0) continue;
	  if(jentry == nbreak) break;
	  if((jentry +1 ) % entry_CP == 0) std::cout<< "--> Processing event number " << jentry + 1 << std::endl;

	  //MUON Reco
	  gen_Daugh_ptetaphiM(isJPsi, &gen_Mum_vect4D, &gen_Mup_vect4D); //costruisco un quadrivettore dei muoni generati

	  totMum = 0;
	  totMup = 0;

	  //Find DRmin
	  RecoMum_idx = find_minDeltaR(gen_Mum_vect4D, JPsi_code);
	  isRecoMum = isReconstructedTrack(RecoMum_idx, gen_Mum_vect4D, JPsi_code);
	  RecoMup_idx = find_minDeltaR(gen_Mup_vect4D, JPsi_code);
	  isRecoMup = isReconstructedTrack(RecoMup_idx, gen_Mup_vect4D, JPsi_code);
	  if (isRecoMum && isRecoMup) nJPsi++;
	  //... separo mu- e mu+
	  for(UInt_t m =0; m<nMuon; m++ ){
		  if( !Muon_softId[m] )continue;
		  if(Muon_charge[m] < 0.){ 
			  Mum_idx[totMum] = m;
			  totMum +=1;
		  }else{
			  Mup_idx[totMup]= m;
			  totMup += 1;
		  }
	  }//on Muons


	  //PION Reco
	  gen_Daugh_ptetaphiM(isRho, &gen_Pim_vect4D, &gen_Pip_vect4D); //costruisco un quadrivettore dei pioni generati
	  h_X_M_GEN.Fill( (gen_Mum_vect4D+gen_Mup_vect4D + gen_Pim_vect4D+gen_Pip_vect4D).M() );

	  totPim = 0;
	  totPip = 0;

	  //Find DRmin
	  RecoPim_idx = find_minDeltaR(gen_Pim_vect4D, Rho_code);
	  isRecoPim = isReconstructedTrack(RecoPim_idx, gen_Pim_vect4D, Rho_code);
	  RecoPip_idx = find_minDeltaR(gen_Pip_vect4D, Rho_code);
	  isRecoPip = isReconstructedTrack(RecoPip_idx, gen_Pip_vect4D, Rho_code);
	  if ( isRecoPim && isRecoPip ) nRho++;
	  //... separo pi+ e pi-
	  for (UInt_t p = 0; p < nProbeTracks; p++){
		  if( ProbeTracks_isMatchedToMuon[p]  )continue;	
		  if(ProbeTracks_charge[p] < 0){	
			  Pim_idx[totPim] = p;
			  totPim++;
		  }else{
			  Pip_idx[totPip]= p;
			  totPip++;
		  }
	  }//on Tracks

	  //KAON Reco
	  gen_K0s_ptetaphiM(&gen_K0s_vect4D);
	  //std::cout << "GEN  -> Pt Eta Phi  " << gen_K0s_vect4D.Pt() << " " << gen_K0s_vect4D.Eta() << std::endl;
	  RecoK0s_idx = find_minDeltaR(gen_K0s_vect4D, K0s_code);
	  isRecoK0s = isReconstructedTrack(RecoK0s_idx, gen_K0s_vect4D, K0s_code);
	  if (isRecoK0s) nK0s++;
		

	  for(int mum = 0; mum < totMum;  mum++){ //mu-
		  FillPtEtaPhiM(Mum_idx[mum], &Mum_vect4D, isMuon);
		  for(int mup = 0;mup < totMup; mup++){//mu+
			  FillPtEtaPhiM(Mup_idx[mup], &Mup_vect4D, isMuon);

			  //mu-mu+
			  isRecoJPsi = ((Mum_idx[mum] == RecoMum_idx) && isRecoMum) && ((Mup_idx[mup] == RecoMup_idx) && isRecoMup);
			  if( isRecoJPsi ) h_JPsi_M.Fill((Mum_vect4D+Mup_vect4D).M()); //-->JPsi
			  else h_MuMu_M.Fill((Mum_vect4D+Mup_vect4D).M());//M(MuMu)

			  for(int pim = 0; pim < totPim; pim++){//pi-
				  FillPtEtaPhiM(Pim_idx[pim], &Pim_vect4D, isPion);

				  Mcheck = ( Mum_vect4D+Mup_vect4D + Pim_vect4D).M(); // 3.1 < M(MuMu Pi-) < 3.75
				  if(Mcut) if ((Mcheck < limMuMuPi_low) || (Mcheck > limMuMuPi_high)) continue;
				  if (isRecoJPsi) nPim_CutOn++;

				  for(int pip = 0; pip < totPip; pip++){//pi+
					  FillPtEtaPhiM(Pip_idx[pip], &Pip_vect4D, isPion);

					  Mcheck = ( Mum_vect4D+Mup_vect4D + Pip_vect4D).M();// 3.1 < M(MuMu Pi+) < 3.75
					  if(Mcut) if ((Mcheck < limMuMuPi_low) || (Mcheck > limMuMuPi_high)) continue;
					  if (isRecoJPsi) nPip_CutOn++;

					  isRecoRho = ((Pim_idx[pim] == RecoPim_idx) && isRecoPim ) && ((Pip_idx[pip] == RecoPip_idx) && isRecoPip );
					  PtCheck = (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).Pt();
					  if (PTcut )if(PtCheck < limXPt_low || PtCheck > limXPt_high) continue;

					  if(isRecoJPsi){
						  if (isRecoRho) h_Rho_M.Fill( (Pim_vect4D+Pip_vect4D).M()  );
						  else h_PiPi_M.Fill( (Pim_vect4D+Pip_vect4D).M()  );
					  }

					  if( isRecoRho && isRecoJPsi ){ //M(JPsi-Rho)---> X(3872)
						  nX3872++;
						  h_X_M.Fill( (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).M() );
						  h_X_pT.Fill( (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).Pt());

						  h_Mdiff_JPsiPiPi_SGN.Fill( TMath::Abs( (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).M() - (Mum_vect4D+Mup_vect4D).M() + mJPsiPDG ) );
						  h_Mdiff_JPsiPi_SGN.Fill( TMath::Abs( (Mum_vect4D+Mup_vect4D + Pim_vect4D).M() - (Mum_vect4D+Mup_vect4D).M() + mJPsiPDG ) );
						  h_Mdiff_JPsiPi_SGN.Fill( TMath::Abs( (Mum_vect4D+Mup_vect4D + Pip_vect4D).M() - (Mum_vect4D+Mup_vect4D).M() + mJPsiPDG ) );

						  h_JPsiPi_M.Fill((Mum_vect4D+Mup_vect4D + Pim_vect4D).M());
						  h_JPsiPi_M.Fill((Mum_vect4D+Mup_vect4D + Pip_vect4D).M());

					  }
					  if (isRecoJPsi && !isRecoRho){ 

						  h_JPsiPiPi_M.Fill((Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).M()); //M(JPsi-PiPi)
						  h_Mdiff_JPsiPiPi_BKG.Fill( TMath::Abs( (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).M() - (Mum_vect4D+Mup_vect4D).M() + mJPsiPDG ) );
						  h_Mdiff_JPsiPi_BKG.Fill( TMath::Abs( (Mum_vect4D+Mup_vect4D + Pim_vect4D).M() - (Mum_vect4D+Mup_vect4D).M() + mJPsiPDG ) );
						  h_Mdiff_JPsiPi_BKG.Fill( TMath::Abs( (Mum_vect4D+Mup_vect4D + Pip_vect4D).M() - (Mum_vect4D+Mup_vect4D).M() + mJPsiPDG ) );
						  h_JPsiPiPi_pT.Fill( (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).Pt() );

					  }
					  if (!isRecoJPsi || !isRecoRho) h_MuMuPiPi_M.Fill( (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D).M() ); //M(MuMu-PiPi)


				  }//on pi+
			  }//on pi-
			  nPip_CutOn/= totPim;
			  h_NPi_CutOn.Fill(nPip_CutOn + nPim_CutOn);

		  }//on mu+
	  }//on mu-
		
	  h_B0_M_GEN.Fill((gen_Mum_vect4D+gen_Mup_vect4D + gen_Pim_vect4D+gen_Pip_vect4D + gen_K0s_vect4D).M());
	  if ( isRecoMum && isRecoMup  &&  isRecoPim && isRecoPip  &&  isRecoK0s){
		  nB0++;
		  FillPtEtaPhiM(RecoMum_idx, &Mum_vect4D, isMuon);
		  FillPtEtaPhiM(RecoMup_idx, &Mup_vect4D, isMuon);
		  FillPtEtaPhiM(RecoPim_idx, &Pim_vect4D, isPion);
		  FillPtEtaPhiM(RecoPip_idx, &Pip_vect4D, isPion);
		  FillPtEtaPhiM(RecoK0s_idx, &K0s_vect4D, isK0s);
		  h_B0_M.Fill( (Mum_vect4D+Mup_vect4D + Pim_vect4D+Pip_vect4D + K0s_vect4D).M() );
	  }
	
  }//Entries

  std::cout << "Number of reco JPsi    : " << nJPsi << std::endl;
  std::cout << "Number of reco Rho     : "  << nRho << std::endl;
  std::cout << "Number of reco X(3872) : "  << nX3872 << std::endl;
  std::cout << "Number of reco K0s     : "  << nK0s << std::endl;
  std::cout << "Number of reco B0      : "  << nB0 << std::endl;
  TString OutFilePath = "./plots/invMass";
  if(Mcut) OutFilePath.Append("CutON.root");
  else     OutFilePath.Append("CutOFF.root");
  TFile M_outfile(OutFilePath, "RECREATE"); // create the output file
  h_JPsi_M.Write();
  h_MuMu_M.Write();
  h_X_M.Write();
  h_X_M_GEN.Write();
  h_JPsiPiPi_M.Write();
  h_MuMuPiPi_M.Write();

  h_JPsiPi_M.Write();

  h_Rho_M.Write();
  h_PiPi_M.Write();
  h_NPi_CutOn.Write();

  h_Mdiff_JPsiPiPi_SGN.Write();
  h_Mdiff_JPsiPiPi_BKG.Write();
  h_Mdiff_JPsiPi_SGN.Write();
  h_Mdiff_JPsiPi_BKG.Write();

  h_X_pT.Write();
  h_JPsiPiPi_pT.Write();

	h_B0_M_GEN.Write();
	h_B0_M.Write();

  M_outfile.Close();
  std::cout << " Output --> " << OutFilePath << std::endl;

}//Loop()


void MassAnalysis::gen_Daugh_ptetaphiM( const Int_t& Mother_pdgId,ROOT::Math::PtEtaPhiMVector* V_neg,ROOT::Math::PtEtaPhiMVector* V_pos ){

  int Daugh_neg_pdgId =0;
  int Nanny_pdgId =9920443;
  if (Mother_pdgId == 443) Daugh_neg_pdgId = 13; //JPsitoMuMu
  if (Mother_pdgId == 113) Daugh_neg_pdgId = -211; //RhoToPiPi 
  for(ULong64_t i = 0; i < nGenPart ; i++){
    if( (GenPart_pdgId[GenPart_genPartIdxMother[i]] == Mother_pdgId) && (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[i]]] == Nanny_pdgId) ){

      if(GenPart_pdgId[i] == Daugh_neg_pdgId){ //check ptl-
        V_neg->SetEta(GenPart_eta[i]);
        V_neg->SetPhi(GenPart_phi[i]);
        V_neg->SetPt(GenPart_pt[i]);
      }else if (GenPart_pdgId[i] == -Daugh_neg_pdgId) {  //check ptl+
        V_pos->SetEta(GenPart_eta[i]);
        V_pos->SetPhi(GenPart_phi[i]);
        V_pos->SetPt(GenPart_pt[i]);
      }
    }

  }//loop on gen ptls                                                                                                                                                                                               
}// gen_Daugh_ptetaphiM

void MassAnalysis::gen_K0s_ptetaphiM(ROOT::Math::PtEtaPhiMVector* V_gen){
	int Daugh_neg_pdgId = 310;
	int Mother_pdgId = 511;
	
	for(ULong64_t i = 0; i < nGenPart ; i++){
		if( abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == Mother_pdgId){ 

        V_gen->SetEta(GenPart_eta[i]);
        V_gen->SetPhi(GenPart_phi[i]);
        V_gen->SetPt(GenPart_pt[i]);
		}
	}

}//gen_K0s_ptetaphiM()

void MassAnalysis::FillPtEtaPhiM(const int& Ptl_Idx, ROOT::Math::PtEtaPhiMVector* rec4V, const int PDGid){

  if (PDGid == 13){
    rec4V->SetPt(Muon_pt[Ptl_Idx]);
    rec4V->SetEta(Muon_eta[Ptl_Idx]);
    rec4V->SetPhi(Muon_phi[Ptl_Idx]);
  }//muons
  if (PDGid == 211){
    rec4V->SetPt(ProbeTracks_pt[Ptl_Idx]);
    rec4V->SetEta(ProbeTracks_eta[Ptl_Idx]);
    rec4V->SetPhi(ProbeTracks_phi[Ptl_Idx]);
  }//pions
  if (PDGid == 310){
    rec4V->SetPt(K0s_fitted_pt[Ptl_Idx]);
    rec4V->SetEta(K0s_fitted_eta[Ptl_Idx]);
    rec4V->SetPhi(K0s_fitted_phi[Ptl_Idx]);
  }//pions

}//FillPtEtaPhiM()

int MassAnalysis::find_minDeltaR(const ROOT::Math::PtEtaPhiMVector& gen4V, const int analysis_code){

  ROOT::Math::PtEtaPhiMVector rec4V(0.,0.,0.,0.);
  double Delta_R, DRmin = 1000.;
  int DR_min_idx = -1;

  if (analysis_code == 1){
    for (ULong64_t m = 0; m < nMuon; m++){
      if( !Muon_softId[m] )continue;
      FillPtEtaPhiM(m, &rec4V, 13);
      Delta_R = ROOT::Math::VectorUtil::DeltaR (gen4V, rec4V);
      if (Delta_R < DRmin ){
	DRmin = Delta_R;
	DR_min_idx = m;
      }
    }
  }//Muons
  if (analysis_code == 2){
    for(ULong64_t p = 0; p < nProbeTracks; p++){
      if( ProbeTracks_isMatchedToMuon[p]  )continue;
      FillPtEtaPhiM(p, &rec4V, 211); 
      Delta_R = ROOT::Math::VectorUtil::DeltaR (gen4V, rec4V);
      if (Delta_R < DRmin ){
        DRmin = Delta_R;
        DR_min_idx = p;
      }
    }

  }//Pions   
  if (analysis_code == 3){
    for(ULong64_t k = 0; k < nK0s; k++){
      FillPtEtaPhiM(k, &rec4V, 310); 
		//std::cout << "RECO PT " << rec4V.Pt() << std::endl;
      Delta_R = ROOT::Math::VectorUtil::DeltaR (gen4V, rec4V);
      if (Delta_R < DRmin ){
        DRmin = Delta_R;
        DR_min_idx = k;
      }
    }

  }//Kaons
  
  return DR_min_idx;
  
}//find_minDeltaR() 


double MassAnalysis::score_PT(const int track_idx, const ROOT::Math::PtEtaPhiMVector& gen4V, const int analysis_code){

  double rec_pt = 0; // dipende se mu o pi                                                                                                                                                                          
  if (analysis_code == 1) rec_pt = Muon_pt[track_idx];
  if (analysis_code == 2) rec_pt = ProbeTracks_pt[track_idx];
  if (analysis_code == 3) rec_pt = K0s_fitted_pt[track_idx];

  double score = TMath::Abs(gen4V.Pt() - rec_pt)/gen4V.Pt();

  return score;
}//score_PT()

bool MassAnalysis::isReconstructedTrack(const int& DR_min_idx , const ROOT::Math::PtEtaPhiMVector& gen4V, const int& analysis_code){
  float DR_threshold = 0.03;
  float DpT_threshold = 0.5;
  bool track_score = false;
  int ptlPDGid = -1;

  if (analysis_code == 1) ptlPDGid = 13;
  if (analysis_code == 2) ptlPDGid = 211;
  if (analysis_code == 3) ptlPDGid = 310;

  ROOT::Math::PtEtaPhiMVector rec4V(0.,0.,0.,0.);
  FillPtEtaPhiM(DR_min_idx, &rec4V, ptlPDGid);

  double DR_min = ROOT::Math::VectorUtil::DeltaR (gen4V, rec4V);
  double DpT_DRmin = score_PT(DR_min_idx, gen4V, analysis_code);
  if ((DR_min < DR_threshold)&& (DpT_DRmin < DpT_threshold)) track_score = true;
  return track_score;

}//isReconstructedTrack()

