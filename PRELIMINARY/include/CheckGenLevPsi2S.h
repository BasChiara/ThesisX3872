//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 28 18:06:53 2022 by ROOT version 6.24/06
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKsLivia_2022Apr26/BPH_BToXKs_XToJPsiRho_JPsiToMuMu/crab_BdToX3872Ks/220426_134635/0000/xNANO_mc_2022Apr26_380.root
//////////////////////////////////////////////////////////

#ifndef CheckGenLevPsi2S_h
#define CheckGenLevPsi2S_h
#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TMath.h>
#include "Math/Vector4D.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
// Header file for the classes stored in the TTree if any.

class CheckGenLevPsi2S {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nB0;
   Float_t         B0_K0s_matchTrack1_D0sign[1];   //[nB0]
   Float_t         B0_K0s_matchTrack1_eta[1];   //[nB0]
   Float_t         B0_K0s_matchTrack1_phi[1];   //[nB0]
   Float_t         B0_K0s_matchTrack1_pt[1];   //[nB0]
   Float_t         B0_K0s_matchTrack2_D0sign[1];   //[nB0]
   Float_t         B0_K0s_matchTrack2_eta[1];   //[nB0]
   Float_t         B0_K0s_matchTrack2_phi[1];   //[nB0]
   Float_t         B0_K0s_matchTrack2_pt[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_eta[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_mass[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_phi[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_pt[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_svprob[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxX[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxXE[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxY[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxYE[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxZ[1];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxZE[1];   //[nB0]
   Float_t         B0_K0s_nmcFitted_mass[1];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi1eta[1];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi1phi[1];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi1pt[1];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi2eta[1];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi2phi[1];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi2pt[1];   //[nB0]
   Float_t         B0_MuMu_DCA[1];   //[nB0]
   Float_t         B0_MuMu_LxySign[1];   //[nB0]
   Float_t         B0_MuMu_cosAlpha[1];   //[nB0]
   Float_t         B0_MuMu_fitted_eta[1];   //[nB0]
   Float_t         B0_MuMu_fitted_mass[1];   //[nB0]
   Float_t         B0_MuMu_fitted_phi[1];   //[nB0]
   Float_t         B0_MuMu_fitted_pt[1];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxX[1];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxXE[1];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxY[1];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxYE[1];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxZ[1];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxZE[1];   //[nB0]
   Float_t         B0_MuMu_mu1_dr[1];   //[nB0]
   Float_t         B0_MuMu_mu2_dr[1];   //[nB0]
   Float_t         B0_MuMu_prefit_mu1_eta[1];   //[nB0]
   Float_t         B0_MuMu_prefit_mu1_phi[1];   //[nB0]
   Float_t         B0_MuMu_prefit_mu1_pt[1];   //[nB0]
   Float_t         B0_MuMu_prefit_mu2_eta[1];   //[nB0]
   Float_t         B0_MuMu_prefit_mu2_phi[1];   //[nB0]
   Float_t         B0_MuMu_prefit_mu2_pt[1];   //[nB0]
   Float_t         B0_MuMu_sv_prob[1];   //[nB0]
   Float_t         B0_PVEx[1];   //[nB0]
   Float_t         B0_PVEy[1];   //[nB0]
   Float_t         B0_PVEz[1];   //[nB0]
   Float_t         B0_PVx[1];   //[nB0]
   Float_t         B0_PVy[1];   //[nB0]
   Float_t         B0_PVz[1];   //[nB0]
   Float_t         B0_PiPi_pi1_d0sig[1];   //[nB0]
   Float_t         B0_PiPi_pi2_d0sig[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_eta[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_phi[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_pt[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_vx[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_vy[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_vz[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_eta[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_phi[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_pt[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_vx[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_vy[1];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_vz[1];   //[nB0]
   Float_t         B0_cosAlpha_PV[1];   //[nB0]
   Float_t         B0_decayVtxX[1];   //[nB0]
   Float_t         B0_decayVtxXE[1];   //[nB0]
   Float_t         B0_decayVtxY[1];   //[nB0]
   Float_t         B0_decayVtxYE[1];   //[nB0]
   Float_t         B0_decayVtxZ[1];   //[nB0]
   Float_t         B0_decayVtxZE[1];   //[nB0]
   Float_t         B0_finalFit_JPsi_mass[1];   //[nB0]
   Float_t         B0_finalFit_Rho_mass[1];   //[nB0]
   Float_t         B0_finalFit_X_mass[1];   //[nB0]
   Float_t         B0_finalFit_eta[1];   //[nB0]
   Float_t         B0_finalFit_k0s_eta[1];   //[nB0]
   Float_t         B0_finalFit_k0s_phi[1];   //[nB0]
   Float_t         B0_finalFit_k0s_pt[1];   //[nB0]
   Float_t         B0_finalFit_mass[1];   //[nB0]
   Float_t         B0_finalFit_mu1_eta[1];   //[nB0]
   Float_t         B0_finalFit_mu1_phi[1];   //[nB0]
   Float_t         B0_finalFit_mu1_pt[1];   //[nB0]
   Float_t         B0_finalFit_mu2_eta[1];   //[nB0]
   Float_t         B0_finalFit_mu2_phi[1];   //[nB0]
   Float_t         B0_finalFit_mu2_pt[1];   //[nB0]
   Float_t         B0_finalFit_phi[1];   //[nB0]
   Float_t         B0_finalFit_pi1_eta[1];   //[nB0]
   Float_t         B0_finalFit_pi1_phi[1];   //[nB0]
   Float_t         B0_finalFit_pi1_pt[1];   //[nB0]
   Float_t         B0_finalFit_pi2_eta[1];   //[nB0]
   Float_t         B0_finalFit_pi2_phi[1];   //[nB0]
   Float_t         B0_finalFit_pi2_pt[1];   //[nB0]
   Float_t         B0_finalFit_pt[1];   //[nB0]
   Float_t         B0_fitted_mass_womc[1];   //[nB0]
   Float_t         B0_lxySign_PV[1];   //[nB0]
   Float_t         B0_svchi2[1];   //[nB0]
   Float_t         B0_svprob[1];   //[nB0]
   Int_t           B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrkTrk_Displaced[1];   //[nB0]
   Int_t           B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[1];   //[nB0]
   Int_t           B0_K0s_matchTrack1_fired_DoubleMu4_PsiPrimeTrk_Displaced[1];   //[nB0]
   Int_t           B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrkTrk_Displaced[1];   //[nB0]
   Int_t           B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[1];   //[nB0]
   Int_t           B0_K0s_matchTrack2_fired_DoubleMu4_PsiPrimeTrk_Displaced[1];   //[nB0]
   Int_t           B0_MuMu_mu1_fired_Dimuon18_PsiPrime[1];   //[nB0]
   Int_t           B0_MuMu_mu1_fired_Dimuon25_Jpsi[1];   //[nB0]
   Int_t           B0_MuMu_mu1_fired_DoubleMu4_JpsiTrkTrk_Displaced[1];   //[nB0]
   Int_t           B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[1];   //[nB0]
   Int_t           B0_MuMu_mu1_fired_DoubleMu4_PsiPrimeTrk_Displaced[1];   //[nB0]
   Int_t           B0_MuMu_mu2_fired_Dimuon18_PsiPrime[1];   //[nB0]
   Int_t           B0_MuMu_mu2_fired_Dimuon25_Jpsi[1];   //[nB0]
   Int_t           B0_MuMu_mu2_fired_DoubleMu4_JpsiTrkTrk_Displaced[1];   //[nB0]
   Int_t           B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[1];   //[nB0]
   Int_t           B0_MuMu_mu2_fired_DoubleMu4_PsiPrimeTrk_Displaced[1];   //[nB0]
   Int_t           B0_PiPi_p1_fired_DoubleMu4_JpsiTrkTrk_Displaced[1];   //[nB0]
   Int_t           B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[1];   //[nB0]
   Int_t           B0_PiPi_p1_fired_DoubleMu4_PsiPrimeTrk_Displaced[1];   //[nB0]
   Int_t           B0_PiPi_p2_fired_DoubleMu4_JpsiTrkTrk_Displaced[1];   //[nB0]
   Int_t           B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[1];   //[nB0]
   Int_t           B0_PiPi_p2_fired_DoubleMu4_PsiPrimeTrk_Displaced[1];   //[nB0]
   Int_t           B0_dimuon_idx[1];   //[nB0]
   Int_t           B0_dipion_idx[1];   //[nB0]
   Int_t           B0_k0short_idx[1];   //[nB0]
   Int_t           B0_mu1_idx[1];   //[nB0]
   Int_t           B0_mu2_idx[1];   //[nB0]
   Int_t           B0_pi1_idx[1];   //[nB0]
   Int_t           B0_pi2_idx[1];   //[nB0]
   Int_t           B0_pv_idx[1];   //[nB0]
   UInt_t          nJPsiToMuMu;
   UInt_t          nK0s;
   UInt_t          npipi;
   UInt_t          nGenPart;
   Float_t         GenPart_eta[500];   //[nGenPart]
   Float_t         GenPart_mass[500];   //[nGenPart]
   Float_t         GenPart_phi[500];   //[nGenPart]
   Float_t         GenPart_pt[500];   //[nGenPart]
   Float_t         GenPart_vx[500];   //[nGenPart]
   Float_t         GenPart_vy[500];   //[nGenPart]
   Float_t         GenPart_vz[500];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[500];   //[nGenPart]
   Int_t           GenPart_pdgId[500];   //[nGenPart]
   Int_t           GenPart_status[500];   //[nGenPart]
   Int_t           GenPart_statusFlags[500];   //[nGenPart]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         genWeight;
   UInt_t          nPSWeight;
   Float_t         PSWeight[4];   //[nPSWeight]
   UInt_t          nMuon;
   Float_t         Muon_eta[4];   //[nMuon]
   Float_t         Muon_mass[4];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[4];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[4];   //[nMuon]
   Float_t         Muon_phi[4];   //[nMuon]
   Float_t         Muon_pt[4];   //[nMuon]
   Float_t         Muon_ptErr[4];   //[nMuon]
   Float_t         Muon_vx[4];   //[nMuon]
   Float_t         Muon_vy[4];   //[nMuon]
   Float_t         Muon_vz[4];   //[nMuon]
   Int_t           Muon_charge[4];   //[nMuon]
   Int_t           Muon_pdgId[4];   //[nMuon]
   Bool_t          Muon_isGlobal[4];   //[nMuon]
   Bool_t          Muon_looseId[4];   //[nMuon]
   Bool_t          Muon_softId[4];   //[nMuon]
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nProbeTracks;
   Float_t         ProbeTracks_d0sig[400];   //[nProbeTracks]
   Float_t         ProbeTracks_eta[400];   //[nProbeTracks]
   Float_t         ProbeTracks_mass[400];   //[nProbeTracks]
   Float_t         ProbeTracks_phi[400];   //[nProbeTracks]
   Float_t         ProbeTracks_pt[400];   //[nProbeTracks]
   Float_t         ProbeTracks_vx[400];   //[nProbeTracks]
   Float_t         ProbeTracks_vy[400];   //[nProbeTracks]
   Float_t         ProbeTracks_vz[400];   //[nProbeTracks]
   Int_t           ProbeTracks_charge[400];   //[nProbeTracks]
   Int_t           ProbeTracks_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced[400];   //[nProbeTracks]
   Int_t           ProbeTracks_fired_HLT_DoubleMu4_JpsiTrk_Displaced[400];   //[nProbeTracks]
   Int_t           ProbeTracks_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced[400];   //[nProbeTracks]
   Int_t           ProbeTracks_nValidHits[400];   //[nProbeTracks]
   Int_t           ProbeTracks_pdgId[400];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToLooseMuon[400];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToMuon[400];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToSoftMuon[400];   //[nProbeTracks]
   Int_t           HLT_Dimuon25_Jpsi;
   Int_t           HLT_Dimuon20_Jpsi_Barrel_Seagulls;
   Int_t           HLT_DoubleMu4_JpsiTrk_Displaced;
   Int_t           HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   Int_t           HLT_DoubleMu4_3_Jpsi_Displaced;
   Int_t           HLT_DoubleMu4_3_Jpsi;
   Int_t           HLT_DoubleMu4_Jpsi_Displaced;
   Int_t           HLT_Dimuon18_PsiPrime;
   Int_t           HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
   Int_t           HLT_DoubleMu4_PsiPrimeTrk_Displaced;
   Int_t           HLT_Dimuon0_Jpsi3p5_Muon2;
   Int_t           HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;
   Int_t           HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;
   Int_t           HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
   Int_t           HLT_DoubleMu4_3_Bs;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[1];   //[nTrigObj]
   Float_t         TrigObj_eta[1];   //[nTrigObj]
   Float_t         TrigObj_phi[1];   //[nTrigObj]
   Int_t           TrigObj_id[1];   //[nTrigObj]
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[1];   //[nSV]
   Float_t         SV_dlenSig[1];   //[nSV]
   Float_t         SV_dxy[1];   //[nSV]
   Float_t         SV_dxySig[1];   //[nSV]
   Float_t         SV_pAngle[1];   //[nSV]
   Int_t           SV_charge[1];   //[nSV]
   Int_t           Muon_genPartIdx[4];   //[nMuon]
   Int_t           Muon_genPartFlav[4];   //[nMuon]
   Float_t         SV_chi2[1];   //[nSV]
   Float_t         SV_eta[1];   //[nSV]
   Float_t         SV_mass[1];   //[nSV]
   Float_t         SV_ndof[1];   //[nSV]
   Float_t         SV_phi[1];   //[nSV]
   Float_t         SV_pt[1];   //[nSV]
   Float_t         SV_x[1];   //[nSV]
   Float_t         SV_y[1];   //[nSV]
   Float_t         SV_z[1];   //[nSV]
   UChar_t         SV_ntracks[1];   //[nSV]
   Int_t           ProbeTracks_genPartIdx[400];   //[nProbeTracks]
   Int_t           ProbeTracks_genPartFlav[400];   //[nProbeTracks]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nB0;   //!
   TBranch        *b_B0_K0s_matchTrack1_D0sign;   //!
   TBranch        *b_B0_K0s_matchTrack1_eta;   //!
   TBranch        *b_B0_K0s_matchTrack1_phi;   //!
   TBranch        *b_B0_K0s_matchTrack1_pt;   //!
   TBranch        *b_B0_K0s_matchTrack2_D0sign;   //!
   TBranch        *b_B0_K0s_matchTrack2_eta;   //!
   TBranch        *b_B0_K0s_matchTrack2_phi;   //!
   TBranch        *b_B0_K0s_matchTrack2_pt;   //!
   TBranch        *b_B0_K0s_mcFitted_eta;   //!
   TBranch        *b_B0_K0s_mcFitted_mass;   //!
   TBranch        *b_B0_K0s_mcFitted_phi;   //!
   TBranch        *b_B0_K0s_mcFitted_pt;   //!
   TBranch        *b_B0_K0s_mcFitted_svprob;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxX;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxXE;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxY;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxYE;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxZ;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxZE;   //!
   TBranch        *b_B0_K0s_nmcFitted_mass;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi1eta;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi1phi;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi1pt;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi2eta;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi2phi;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi2pt;   //!
   TBranch        *b_B0_MuMu_DCA;   //!
   TBranch        *b_B0_MuMu_LxySign;   //!
   TBranch        *b_B0_MuMu_cosAlpha;   //!
   TBranch        *b_B0_MuMu_fitted_eta;   //!
   TBranch        *b_B0_MuMu_fitted_mass;   //!
   TBranch        *b_B0_MuMu_fitted_phi;   //!
   TBranch        *b_B0_MuMu_fitted_pt;   //!
   TBranch        *b_B0_MuMu_fitted_vtxX;   //!
   TBranch        *b_B0_MuMu_fitted_vtxXE;   //!
   TBranch        *b_B0_MuMu_fitted_vtxY;   //!
   TBranch        *b_B0_MuMu_fitted_vtxYE;   //!
   TBranch        *b_B0_MuMu_fitted_vtxZ;   //!
   TBranch        *b_B0_MuMu_fitted_vtxZE;   //!
   TBranch        *b_B0_MuMu_mu1_dr;   //!
   TBranch        *b_B0_MuMu_mu2_dr;   //!
   TBranch        *b_B0_MuMu_prefit_mu1_eta;   //!
   TBranch        *b_B0_MuMu_prefit_mu1_phi;   //!
   TBranch        *b_B0_MuMu_prefit_mu1_pt;   //!
   TBranch        *b_B0_MuMu_prefit_mu2_eta;   //!
   TBranch        *b_B0_MuMu_prefit_mu2_phi;   //!
   TBranch        *b_B0_MuMu_prefit_mu2_pt;   //!
   TBranch        *b_B0_MuMu_sv_prob;   //!
   TBranch        *b_B0_PVEx;   //!
   TBranch        *b_B0_PVEy;   //!
   TBranch        *b_B0_PVEz;   //!
   TBranch        *b_B0_PVx;   //!
   TBranch        *b_B0_PVy;   //!
   TBranch        *b_B0_PVz;   //!
   TBranch        *b_B0_PiPi_pi1_d0sig;   //!
   TBranch        *b_B0_PiPi_pi2_d0sig;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_eta;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_phi;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_pt;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_vx;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_vy;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_vz;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_eta;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_phi;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_pt;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_vx;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_vy;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_vz;   //!
   TBranch        *b_B0_cosAlpha_PV;   //!
   TBranch        *b_B0_decayVtxX;   //!
   TBranch        *b_B0_decayVtxXE;   //!
   TBranch        *b_B0_decayVtxY;   //!
   TBranch        *b_B0_decayVtxYE;   //!
   TBranch        *b_B0_decayVtxZ;   //!
   TBranch        *b_B0_decayVtxZE;   //!
   TBranch        *b_B0_finalFit_JPsi_mass;   //!
   TBranch        *b_B0_finalFit_Rho_mass;   //!
   TBranch        *b_B0_finalFit_X_mass;   //!
   TBranch        *b_B0_finalFit_eta;   //!
   TBranch        *b_B0_finalFit_k0s_eta;   //!
   TBranch        *b_B0_finalFit_k0s_phi;   //!
   TBranch        *b_B0_finalFit_k0s_pt;   //!
   TBranch        *b_B0_finalFit_mass;   //!
   TBranch        *b_B0_finalFit_mu1_eta;   //!
   TBranch        *b_B0_finalFit_mu1_phi;   //!
   TBranch        *b_B0_finalFit_mu1_pt;   //!
   TBranch        *b_B0_finalFit_mu2_eta;   //!
   TBranch        *b_B0_finalFit_mu2_phi;   //!
   TBranch        *b_B0_finalFit_mu2_pt;   //!
   TBranch        *b_B0_finalFit_phi;   //!
   TBranch        *b_B0_finalFit_pi1_eta;   //!
   TBranch        *b_B0_finalFit_pi1_phi;   //!
   TBranch        *b_B0_finalFit_pi1_pt;   //!
   TBranch        *b_B0_finalFit_pi2_eta;   //!
   TBranch        *b_B0_finalFit_pi2_phi;   //!
   TBranch        *b_B0_finalFit_pi2_pt;   //!
   TBranch        *b_B0_finalFit_pt;   //!
   TBranch        *b_B0_fitted_mass_womc;   //!
   TBranch        *b_B0_lxySign_PV;   //!
   TBranch        *b_B0_svchi2;   //!
   TBranch        *b_B0_svprob;   //!
   TBranch        *b_B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_K0s_matchTrack1_fired_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_K0s_matchTrack2_fired_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu1_fired_Dimuon18_PsiPrime;   //!
   TBranch        *b_B0_MuMu_mu1_fired_Dimuon25_Jpsi;   //!
   TBranch        *b_B0_MuMu_mu1_fired_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu1_fired_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu2_fired_Dimuon18_PsiPrime;   //!
   TBranch        *b_B0_MuMu_mu2_fired_Dimuon25_Jpsi;   //!
   TBranch        *b_B0_MuMu_mu2_fired_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu2_fired_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p1_fired_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p1_fired_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p2_fired_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p2_fired_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_B0_dimuon_idx;   //!
   TBranch        *b_B0_dipion_idx;   //!
   TBranch        *b_B0_k0short_idx;   //!
   TBranch        *b_B0_mu1_idx;   //!
   TBranch        *b_B0_mu2_idx;   //!
   TBranch        *b_B0_pi1_idx;   //!
   TBranch        *b_B0_pi2_idx;   //!
   TBranch        *b_B0_pv_idx;   //!
   TBranch        *b_nJPsiToMuMu;   //!
   TBranch        *b_nK0s;   //!
   TBranch        *b_npipi;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_vx;   //!
   TBranch        *b_GenPart_vy;   //!
   TBranch        *b_GenPart_vz;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nProbeTracks;   //!
   TBranch        *b_ProbeTracks_d0sig;   //!
   TBranch        *b_ProbeTracks_eta;   //!
   TBranch        *b_ProbeTracks_mass;   //!
   TBranch        *b_ProbeTracks_phi;   //!
   TBranch        *b_ProbeTracks_pt;   //!
   TBranch        *b_ProbeTracks_vx;   //!
   TBranch        *b_ProbeTracks_vy;   //!
   TBranch        *b_ProbeTracks_vz;   //!
   TBranch        *b_ProbeTracks_charge;   //!
   TBranch        *b_ProbeTracks_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_ProbeTracks_fired_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_ProbeTracks_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_ProbeTracks_nValidHits;   //!
   TBranch        *b_ProbeTracks_pdgId;   //!
   TBranch        *b_ProbeTracks_isMatchedToLooseMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToSoftMuon;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi;   //!
   TBranch        *b_HLT_Dimuon20_Jpsi_Barrel_Seagulls;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi;   //!
   TBranch        *b_HLT_DoubleMu4_Jpsi_Displaced;   //!
   TBranch        *b_HLT_Dimuon18_PsiPrime;   //!
   TBranch        *b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls;   //!
   TBranch        *b_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;   //!
   TBranch        *b_HLT_DoubleMu4_3_Bs;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_SV_charge;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_SV_ntracks;   //!
   TBranch        *b_ProbeTracks_genPartIdx;   //!
   TBranch        *b_ProbeTracks_genPartFlav;   //!

   CheckGenLevPsi2S(TTree *tree=0, const TString& tag = "MC");
   virtual ~CheckGenLevPsi2S();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void     GenPart_FillKinHist(ROOT::Math::PtEtaPhiMVector* GenVec, TH1* h_pt, TH1* h_eta, TH1* h_mass);
   virtual void     GenPartFillP4();
   private:

   // I/O
   TString tags_;

   // PDG particle id & mass
	int isMum = 13, isPip = 211, isJPsi = 443, isK0s = 310, isPsi2S = 100443, isB0 = 511;
	const float mMuon = 0.105658, mPion = 0.1395704, mK0s = 0.497648; 
	// generator ptls 4-vectors
	ROOT::Math::PtEtaPhiMVector GenP4_Mum, GenP4_Mup; // muons
	ROOT::Math::PtEtaPhiMVector GenP4_JPsi;           //JPsi
	ROOT::Math::PtEtaPhiMVector GenP4_Pim, GenP4_Pip; // pions
	ROOT::Math::PtEtaPhiMVector GenP4_PiPi;           // PiPi 
	ROOT::Math::PtEtaPhiMVector GenP4_Psi2S;          // Psi(2S)
	ROOT::Math::PtEtaPhiMVector GenP4_K0s;            // K0s
	ROOT::Math::PtEtaPhiMVector GenP4_B0;             // B0


};

#endif

#ifdef CheckGenLevPsi2S_cxx
CheckGenLevPsi2S::CheckGenLevPsi2S(TTree *tree, const TString& tag) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKsLivia_2022Apr26/BPH_BToXKs_XToJPsiRho_JPsiToMuMu/crab_BdToX3872Ks/220426_134635/0000/xNANO_mc_2022Apr26_380.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKsLivia_2022Apr26/BPH_BToXKs_XToJPsiRho_JPsiToMuMu/crab_BdToX3872Ks/220426_134635/0000/xNANO_mc_2022Apr26_380.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);

   tags_ = tag;
}

CheckGenLevPsi2S::~CheckGenLevPsi2S()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CheckGenLevPsi2S::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CheckGenLevPsi2S::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CheckGenLevPsi2S::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nB0", &nB0, &b_nB0);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_D0sign", B0_K0s_matchTrack1_D0sign, &b_B0_K0s_matchTrack1_D0sign);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_eta", B0_K0s_matchTrack1_eta, &b_B0_K0s_matchTrack1_eta);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_phi", B0_K0s_matchTrack1_phi, &b_B0_K0s_matchTrack1_phi);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_pt", B0_K0s_matchTrack1_pt, &b_B0_K0s_matchTrack1_pt);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_D0sign", B0_K0s_matchTrack2_D0sign, &b_B0_K0s_matchTrack2_D0sign);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_eta", B0_K0s_matchTrack2_eta, &b_B0_K0s_matchTrack2_eta);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_phi", B0_K0s_matchTrack2_phi, &b_B0_K0s_matchTrack2_phi);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_pt", B0_K0s_matchTrack2_pt, &b_B0_K0s_matchTrack2_pt);
   fChain->SetBranchAddress("B0_K0s_mcFitted_eta", B0_K0s_mcFitted_eta, &b_B0_K0s_mcFitted_eta);
   fChain->SetBranchAddress("B0_K0s_mcFitted_mass", B0_K0s_mcFitted_mass, &b_B0_K0s_mcFitted_mass);
   fChain->SetBranchAddress("B0_K0s_mcFitted_phi", B0_K0s_mcFitted_phi, &b_B0_K0s_mcFitted_phi);
   fChain->SetBranchAddress("B0_K0s_mcFitted_pt", B0_K0s_mcFitted_pt, &b_B0_K0s_mcFitted_pt);
   fChain->SetBranchAddress("B0_K0s_mcFitted_svprob", B0_K0s_mcFitted_svprob, &b_B0_K0s_mcFitted_svprob);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxX", B0_K0s_mcFitted_vtxX, &b_B0_K0s_mcFitted_vtxX);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxXE", B0_K0s_mcFitted_vtxXE, &b_B0_K0s_mcFitted_vtxXE);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxY", B0_K0s_mcFitted_vtxY, &b_B0_K0s_mcFitted_vtxY);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxYE", B0_K0s_mcFitted_vtxYE, &b_B0_K0s_mcFitted_vtxYE);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxZ", B0_K0s_mcFitted_vtxZ, &b_B0_K0s_mcFitted_vtxZ);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxZE", B0_K0s_mcFitted_vtxZE, &b_B0_K0s_mcFitted_vtxZE);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_mass", B0_K0s_nmcFitted_mass, &b_B0_K0s_nmcFitted_mass);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi1eta", B0_K0s_nmcFitted_pi1eta, &b_B0_K0s_nmcFitted_pi1eta);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi1phi", B0_K0s_nmcFitted_pi1phi, &b_B0_K0s_nmcFitted_pi1phi);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi1pt", B0_K0s_nmcFitted_pi1pt, &b_B0_K0s_nmcFitted_pi1pt);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi2eta", B0_K0s_nmcFitted_pi2eta, &b_B0_K0s_nmcFitted_pi2eta);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi2phi", B0_K0s_nmcFitted_pi2phi, &b_B0_K0s_nmcFitted_pi2phi);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi2pt", B0_K0s_nmcFitted_pi2pt, &b_B0_K0s_nmcFitted_pi2pt);
   fChain->SetBranchAddress("B0_MuMu_DCA", B0_MuMu_DCA, &b_B0_MuMu_DCA);
   fChain->SetBranchAddress("B0_MuMu_LxySign", B0_MuMu_LxySign, &b_B0_MuMu_LxySign);
   fChain->SetBranchAddress("B0_MuMu_cosAlpha", B0_MuMu_cosAlpha, &b_B0_MuMu_cosAlpha);
   fChain->SetBranchAddress("B0_MuMu_fitted_eta", B0_MuMu_fitted_eta, &b_B0_MuMu_fitted_eta);
   fChain->SetBranchAddress("B0_MuMu_fitted_mass", B0_MuMu_fitted_mass, &b_B0_MuMu_fitted_mass);
   fChain->SetBranchAddress("B0_MuMu_fitted_phi", B0_MuMu_fitted_phi, &b_B0_MuMu_fitted_phi);
   fChain->SetBranchAddress("B0_MuMu_fitted_pt", B0_MuMu_fitted_pt, &b_B0_MuMu_fitted_pt);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxX", B0_MuMu_fitted_vtxX, &b_B0_MuMu_fitted_vtxX);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxXE", B0_MuMu_fitted_vtxXE, &b_B0_MuMu_fitted_vtxXE);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxY", B0_MuMu_fitted_vtxY, &b_B0_MuMu_fitted_vtxY);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxYE", B0_MuMu_fitted_vtxYE, &b_B0_MuMu_fitted_vtxYE);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxZ", B0_MuMu_fitted_vtxZ, &b_B0_MuMu_fitted_vtxZ);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxZE", B0_MuMu_fitted_vtxZE, &b_B0_MuMu_fitted_vtxZE);
   fChain->SetBranchAddress("B0_MuMu_mu1_dr", B0_MuMu_mu1_dr, &b_B0_MuMu_mu1_dr);
   fChain->SetBranchAddress("B0_MuMu_mu2_dr", B0_MuMu_mu2_dr, &b_B0_MuMu_mu2_dr);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu1_eta", B0_MuMu_prefit_mu1_eta, &b_B0_MuMu_prefit_mu1_eta);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu1_phi", B0_MuMu_prefit_mu1_phi, &b_B0_MuMu_prefit_mu1_phi);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu1_pt", B0_MuMu_prefit_mu1_pt, &b_B0_MuMu_prefit_mu1_pt);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu2_eta", B0_MuMu_prefit_mu2_eta, &b_B0_MuMu_prefit_mu2_eta);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu2_phi", B0_MuMu_prefit_mu2_phi, &b_B0_MuMu_prefit_mu2_phi);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu2_pt", B0_MuMu_prefit_mu2_pt, &b_B0_MuMu_prefit_mu2_pt);
   fChain->SetBranchAddress("B0_MuMu_sv_prob", B0_MuMu_sv_prob, &b_B0_MuMu_sv_prob);
   fChain->SetBranchAddress("B0_PVEx", B0_PVEx, &b_B0_PVEx);
   fChain->SetBranchAddress("B0_PVEy", B0_PVEy, &b_B0_PVEy);
   fChain->SetBranchAddress("B0_PVEz", B0_PVEz, &b_B0_PVEz);
   fChain->SetBranchAddress("B0_PVx", B0_PVx, &b_B0_PVx);
   fChain->SetBranchAddress("B0_PVy", B0_PVy, &b_B0_PVy);
   fChain->SetBranchAddress("B0_PVz", B0_PVz, &b_B0_PVz);
   fChain->SetBranchAddress("B0_PiPi_pi1_d0sig", B0_PiPi_pi1_d0sig, &b_B0_PiPi_pi1_d0sig);
   fChain->SetBranchAddress("B0_PiPi_pi2_d0sig", B0_PiPi_pi2_d0sig, &b_B0_PiPi_pi2_d0sig);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_eta", B0_PiPi_prefit_pi1_eta, &b_B0_PiPi_prefit_pi1_eta);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_phi", B0_PiPi_prefit_pi1_phi, &b_B0_PiPi_prefit_pi1_phi);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_pt", B0_PiPi_prefit_pi1_pt, &b_B0_PiPi_prefit_pi1_pt);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_vx", B0_PiPi_prefit_pi1_vx, &b_B0_PiPi_prefit_pi1_vx);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_vy", B0_PiPi_prefit_pi1_vy, &b_B0_PiPi_prefit_pi1_vy);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_vz", B0_PiPi_prefit_pi1_vz, &b_B0_PiPi_prefit_pi1_vz);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_eta", B0_PiPi_prefit_pi2_eta, &b_B0_PiPi_prefit_pi2_eta);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_phi", B0_PiPi_prefit_pi2_phi, &b_B0_PiPi_prefit_pi2_phi);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_pt", B0_PiPi_prefit_pi2_pt, &b_B0_PiPi_prefit_pi2_pt);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_vx", B0_PiPi_prefit_pi2_vx, &b_B0_PiPi_prefit_pi2_vx);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_vy", B0_PiPi_prefit_pi2_vy, &b_B0_PiPi_prefit_pi2_vy);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_vz", B0_PiPi_prefit_pi2_vz, &b_B0_PiPi_prefit_pi2_vz);
   fChain->SetBranchAddress("B0_cosAlpha_PV", B0_cosAlpha_PV, &b_B0_cosAlpha_PV);
   fChain->SetBranchAddress("B0_decayVtxX", B0_decayVtxX, &b_B0_decayVtxX);
   fChain->SetBranchAddress("B0_decayVtxXE", B0_decayVtxXE, &b_B0_decayVtxXE);
   fChain->SetBranchAddress("B0_decayVtxY", B0_decayVtxY, &b_B0_decayVtxY);
   fChain->SetBranchAddress("B0_decayVtxYE", B0_decayVtxYE, &b_B0_decayVtxYE);
   fChain->SetBranchAddress("B0_decayVtxZ", B0_decayVtxZ, &b_B0_decayVtxZ);
   fChain->SetBranchAddress("B0_decayVtxZE", B0_decayVtxZE, &b_B0_decayVtxZE);
   fChain->SetBranchAddress("B0_finalFit_JPsi_mass", B0_finalFit_JPsi_mass, &b_B0_finalFit_JPsi_mass);
   fChain->SetBranchAddress("B0_finalFit_Rho_mass", B0_finalFit_Rho_mass, &b_B0_finalFit_Rho_mass);
   fChain->SetBranchAddress("B0_finalFit_X_mass", B0_finalFit_X_mass, &b_B0_finalFit_X_mass);
   fChain->SetBranchAddress("B0_finalFit_eta", B0_finalFit_eta, &b_B0_finalFit_eta);
   fChain->SetBranchAddress("B0_finalFit_k0s_eta", B0_finalFit_k0s_eta, &b_B0_finalFit_k0s_eta);
   fChain->SetBranchAddress("B0_finalFit_k0s_phi", B0_finalFit_k0s_phi, &b_B0_finalFit_k0s_phi);
   fChain->SetBranchAddress("B0_finalFit_k0s_pt", B0_finalFit_k0s_pt, &b_B0_finalFit_k0s_pt);
   fChain->SetBranchAddress("B0_finalFit_mass", B0_finalFit_mass, &b_B0_finalFit_mass);
   fChain->SetBranchAddress("B0_finalFit_mu1_eta", B0_finalFit_mu1_eta, &b_B0_finalFit_mu1_eta);
   fChain->SetBranchAddress("B0_finalFit_mu1_phi", B0_finalFit_mu1_phi, &b_B0_finalFit_mu1_phi);
   fChain->SetBranchAddress("B0_finalFit_mu1_pt", B0_finalFit_mu1_pt, &b_B0_finalFit_mu1_pt);
   fChain->SetBranchAddress("B0_finalFit_mu2_eta", B0_finalFit_mu2_eta, &b_B0_finalFit_mu2_eta);
   fChain->SetBranchAddress("B0_finalFit_mu2_phi", B0_finalFit_mu2_phi, &b_B0_finalFit_mu2_phi);
   fChain->SetBranchAddress("B0_finalFit_mu2_pt", B0_finalFit_mu2_pt, &b_B0_finalFit_mu2_pt);
   fChain->SetBranchAddress("B0_finalFit_phi", B0_finalFit_phi, &b_B0_finalFit_phi);
   fChain->SetBranchAddress("B0_finalFit_pi1_eta", B0_finalFit_pi1_eta, &b_B0_finalFit_pi1_eta);
   fChain->SetBranchAddress("B0_finalFit_pi1_phi", B0_finalFit_pi1_phi, &b_B0_finalFit_pi1_phi);
   fChain->SetBranchAddress("B0_finalFit_pi1_pt", B0_finalFit_pi1_pt, &b_B0_finalFit_pi1_pt);
   fChain->SetBranchAddress("B0_finalFit_pi2_eta", B0_finalFit_pi2_eta, &b_B0_finalFit_pi2_eta);
   fChain->SetBranchAddress("B0_finalFit_pi2_phi", B0_finalFit_pi2_phi, &b_B0_finalFit_pi2_phi);
   fChain->SetBranchAddress("B0_finalFit_pi2_pt", B0_finalFit_pi2_pt, &b_B0_finalFit_pi2_pt);
   fChain->SetBranchAddress("B0_finalFit_pt", B0_finalFit_pt, &b_B0_finalFit_pt);
   fChain->SetBranchAddress("B0_fitted_mass_womc", B0_fitted_mass_womc, &b_B0_fitted_mass_womc);
   fChain->SetBranchAddress("B0_lxySign_PV", B0_lxySign_PV, &b_B0_lxySign_PV);
   fChain->SetBranchAddress("B0_svchi2", B0_svchi2, &b_B0_svchi2);
   fChain->SetBranchAddress("B0_svprob", B0_svprob, &b_B0_svprob);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrkTrk_Displaced", B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrkTrk_Displaced, &b_B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced", B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_fired_DoubleMu4_PsiPrimeTrk_Displaced", B0_K0s_matchTrack1_fired_DoubleMu4_PsiPrimeTrk_Displaced, &b_B0_K0s_matchTrack1_fired_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrkTrk_Displaced", B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrkTrk_Displaced, &b_B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced", B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_fired_DoubleMu4_PsiPrimeTrk_Displaced", B0_K0s_matchTrack2_fired_DoubleMu4_PsiPrimeTrk_Displaced, &b_B0_K0s_matchTrack2_fired_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu1_fired_Dimuon18_PsiPrime", B0_MuMu_mu1_fired_Dimuon18_PsiPrime, &b_B0_MuMu_mu1_fired_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("B0_MuMu_mu1_fired_Dimuon25_Jpsi", B0_MuMu_mu1_fired_Dimuon25_Jpsi, &b_B0_MuMu_mu1_fired_Dimuon25_Jpsi);
   fChain->SetBranchAddress("B0_MuMu_mu1_fired_DoubleMu4_JpsiTrkTrk_Displaced", B0_MuMu_mu1_fired_DoubleMu4_JpsiTrkTrk_Displaced, &b_B0_MuMu_mu1_fired_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced", B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu1_fired_DoubleMu4_PsiPrimeTrk_Displaced", B0_MuMu_mu1_fired_DoubleMu4_PsiPrimeTrk_Displaced, &b_B0_MuMu_mu1_fired_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu2_fired_Dimuon18_PsiPrime", B0_MuMu_mu2_fired_Dimuon18_PsiPrime, &b_B0_MuMu_mu2_fired_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("B0_MuMu_mu2_fired_Dimuon25_Jpsi", B0_MuMu_mu2_fired_Dimuon25_Jpsi, &b_B0_MuMu_mu2_fired_Dimuon25_Jpsi);
   fChain->SetBranchAddress("B0_MuMu_mu2_fired_DoubleMu4_JpsiTrkTrk_Displaced", B0_MuMu_mu2_fired_DoubleMu4_JpsiTrkTrk_Displaced, &b_B0_MuMu_mu2_fired_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced", B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu2_fired_DoubleMu4_PsiPrimeTrk_Displaced", B0_MuMu_mu2_fired_DoubleMu4_PsiPrimeTrk_Displaced, &b_B0_MuMu_mu2_fired_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p1_fired_DoubleMu4_JpsiTrkTrk_Displaced", B0_PiPi_p1_fired_DoubleMu4_JpsiTrkTrk_Displaced, &b_B0_PiPi_p1_fired_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced", B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p1_fired_DoubleMu4_PsiPrimeTrk_Displaced", B0_PiPi_p1_fired_DoubleMu4_PsiPrimeTrk_Displaced, &b_B0_PiPi_p1_fired_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p2_fired_DoubleMu4_JpsiTrkTrk_Displaced", B0_PiPi_p2_fired_DoubleMu4_JpsiTrkTrk_Displaced, &b_B0_PiPi_p2_fired_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced", B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p2_fired_DoubleMu4_PsiPrimeTrk_Displaced", B0_PiPi_p2_fired_DoubleMu4_PsiPrimeTrk_Displaced, &b_B0_PiPi_p2_fired_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("B0_dimuon_idx", B0_dimuon_idx, &b_B0_dimuon_idx);
   fChain->SetBranchAddress("B0_dipion_idx", B0_dipion_idx, &b_B0_dipion_idx);
   fChain->SetBranchAddress("B0_k0short_idx", B0_k0short_idx, &b_B0_k0short_idx);
   fChain->SetBranchAddress("B0_mu1_idx", B0_mu1_idx, &b_B0_mu1_idx);
   fChain->SetBranchAddress("B0_mu2_idx", B0_mu2_idx, &b_B0_mu2_idx);
   fChain->SetBranchAddress("B0_pi1_idx", B0_pi1_idx, &b_B0_pi1_idx);
   fChain->SetBranchAddress("B0_pi2_idx", B0_pi2_idx, &b_B0_pi2_idx);
   fChain->SetBranchAddress("B0_pv_idx", B0_pv_idx, &b_B0_pv_idx);
   fChain->SetBranchAddress("nJPsiToMuMu", &nJPsiToMuMu, &b_nJPsiToMuMu);
   fChain->SetBranchAddress("nK0s", &nK0s, &b_nK0s);
   fChain->SetBranchAddress("npipi", &npipi, &b_npipi);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_vx", GenPart_vx, &b_GenPart_vx);
   fChain->SetBranchAddress("GenPart_vy", GenPart_vy, &b_GenPart_vy);
   fChain->SetBranchAddress("GenPart_vz", GenPart_vz, &b_GenPart_vz);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_vx", Muon_vx, &b_Muon_vx);
   fChain->SetBranchAddress("Muon_vy", Muon_vy, &b_Muon_vy);
   fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   fChain->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nProbeTracks", &nProbeTracks, &b_nProbeTracks);
   fChain->SetBranchAddress("ProbeTracks_d0sig", ProbeTracks_d0sig, &b_ProbeTracks_d0sig);
   fChain->SetBranchAddress("ProbeTracks_eta", ProbeTracks_eta, &b_ProbeTracks_eta);
   fChain->SetBranchAddress("ProbeTracks_mass", ProbeTracks_mass, &b_ProbeTracks_mass);
   fChain->SetBranchAddress("ProbeTracks_phi", ProbeTracks_phi, &b_ProbeTracks_phi);
   fChain->SetBranchAddress("ProbeTracks_pt", ProbeTracks_pt, &b_ProbeTracks_pt);
   fChain->SetBranchAddress("ProbeTracks_vx", ProbeTracks_vx, &b_ProbeTracks_vx);
   fChain->SetBranchAddress("ProbeTracks_vy", ProbeTracks_vy, &b_ProbeTracks_vy);
   fChain->SetBranchAddress("ProbeTracks_vz", ProbeTracks_vz, &b_ProbeTracks_vz);
   fChain->SetBranchAddress("ProbeTracks_charge", ProbeTracks_charge, &b_ProbeTracks_charge);
   fChain->SetBranchAddress("ProbeTracks_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced", ProbeTracks_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_ProbeTracks_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("ProbeTracks_fired_HLT_DoubleMu4_JpsiTrk_Displaced", ProbeTracks_fired_HLT_DoubleMu4_JpsiTrk_Displaced, &b_ProbeTracks_fired_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("ProbeTracks_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced", ProbeTracks_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_ProbeTracks_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("ProbeTracks_nValidHits", ProbeTracks_nValidHits, &b_ProbeTracks_nValidHits);
   fChain->SetBranchAddress("ProbeTracks_pdgId", ProbeTracks_pdgId, &b_ProbeTracks_pdgId);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToLooseMuon", ProbeTracks_isMatchedToLooseMuon, &b_ProbeTracks_isMatchedToLooseMuon);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToMuon", ProbeTracks_isMatchedToMuon, &b_ProbeTracks_isMatchedToMuon);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToSoftMuon", ProbeTracks_isMatchedToSoftMuon, &b_ProbeTracks_isMatchedToSoftMuon);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon20_Jpsi_Barrel_Seagulls", &HLT_Dimuon20_Jpsi_Barrel_Seagulls, &b_HLT_Dimuon20_Jpsi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced, &b_HLT_DoubleMu4_3_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi", &HLT_DoubleMu4_3_Jpsi, &b_HLT_DoubleMu4_3_Jpsi);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced, &b_HLT_DoubleMu4_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime, &b_HLT_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon10_PsiPrime_Barrel_Seagulls", &HLT_Dimuon10_PsiPrime_Barrel_Seagulls, &b_HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2, &b_HLT_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", &HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi, &b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05, &b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("SV_charge", SV_charge, &b_SV_charge);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("SV_ntracks", SV_ntracks, &b_SV_ntracks);
   fChain->SetBranchAddress("ProbeTracks_genPartIdx", ProbeTracks_genPartIdx, &b_ProbeTracks_genPartIdx);
   fChain->SetBranchAddress("ProbeTracks_genPartFlav", ProbeTracks_genPartFlav, &b_ProbeTracks_genPartFlav);
   Notify();
}

Bool_t CheckGenLevPsi2S::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CheckGenLevPsi2S::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CheckGenLevPsi2S::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CheckGenLevPsi2S_cxx
