//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 31 10:21:20 2022 by ROOT version 6.24/06
// from TTree Events/Events
// found on file: /eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKs_2022Mar30/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220330_155143/0000/xNANO_mc_2022Mar30_1.root
//////////////////////////////////////////////////////////

#ifndef PrepAnalysisPsi2S_h
#define PrepAnalysisPsi2S_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TString.h>
#include <TH1.h>
#include <TMath.h>
#include "Math/Vector4D.h"
#include <Math/GenVector/VectorUtil.h>
#include <Math/GenVector/PtEtaPhiM4D.h>

// Header file for the classes stored in the TTree if any.

class PrepAnalysisPsi2S {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nB0;
   Float_t         B0_cosAlpha_PV[20];   //[nB0]
   Float_t         B0_dxySignMu1[20];   //[nB0]
   Float_t         B0_dxySignMu2[20];   //[nB0]
   Float_t         B0_dxySignTr1[20];   //[nB0]
   Float_t         B0_dxySignTr2[20];   //[nB0]
   Float_t         B0_dzSignMu1[20];   //[nB0]
   Float_t         B0_dzSignMu2[20];   //[nB0]
   Float_t         B0_dzSignTr1[20];   //[nB0]
   Float_t         B0_dzSignTr2[20];   //[nB0]
   Float_t         B0_fit_eta[20];   //[nB0]
   Float_t         B0_fit_mass[20];   //[nB0]
   Float_t         B0_fit_phi[20];   //[nB0]
   Float_t         B0_fit_pt[20];   //[nB0]
   Float_t         B0_fitted_JPsi_mass[20];   //[nB0]
   Float_t         B0_fitted_Rho_mass[20];   //[nB0]
   Float_t         B0_fitted_X_mass[20];   //[nB0]
   Float_t         B0_fitted_mass_womc[20];   //[nB0]
   Float_t         B0_lxySign_PV[20];   //[nB0]
   Float_t         B0_svchi2[20];   //[nB0]
   Float_t         B0_svprob[20];   //[nB0]
   Int_t           B0_dilepton_idx[20];   //[nB0]
   Int_t           B0_dipion_idx[20];   //[nB0]
   Int_t           B0_k0short_idx[20];   //[nB0]
   Int_t           B0_mu1_idx[20];   //[nB0]
   Int_t           B0_mu2_idx[20];   //[nB0]
   Int_t           B0_pv_idx[20];   //[nB0]
   Int_t           B0_trk1Rho_idx[20];   //[nB0]
   Int_t           B0_trk2Rho_idx[20];   //[nB0]
   UInt_t          nJPsiToMuMu;
   Float_t         JPsiToMuMu_eta[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_fitted_mass[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_fitted_massErr[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_mass[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_phi[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_pt[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_svchi2[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_svndof[2];   //[nJPsiToMuMu]
   Float_t         JPsiToMuMu_svprob[2];   //[nJPsiToMuMu]
   Int_t           JPsiToMuMu_charge[2];   //[nJPsiToMuMu]
   Int_t           JPsiToMuMu_l1_idx[2];   //[nJPsiToMuMu]
   Int_t           JPsiToMuMu_l2_idx[2];   //[nJPsiToMuMu]
   Int_t           JPsiToMuMu_pdgId[2];   //[nJPsiToMuMu]
   UInt_t          nK0s;
   Float_t         K0s_fitted_eta[8];   //[nK0s]
   Float_t         K0s_fitted_mass[8];   //[nK0s]
   Float_t         K0s_fitted_mass_womc[8];   //[nK0s]
   Float_t         K0s_fitted_phi[8];   //[nK0s]
   Float_t         K0s_fitted_pt[8];   //[nK0s]
   Float_t         K0s_prefit_mass[8];   //[nK0s]
   Float_t         K0s_svprob[8];   //[nK0s]
   Float_t         K0s_trk1_eta[8];   //[nK0s]
   Float_t         K0s_trk1_phi[8];   //[nK0s]
   Float_t         K0s_trk1_pt[8];   //[nK0s]
   Float_t         K0s_trk2_eta[8];   //[nK0s]
   Float_t         K0s_trk2_phi[8];   //[nK0s]
   Float_t         K0s_trk2_pt[8];   //[nK0s]
   Int_t           K0s_trk1_charge[8];   //[nK0s]
   Int_t           K0s_trk2_charge[8];   //[nK0s]
   UInt_t          npipi;
   Float_t         pipi_eta[20000];   //[npipi]
   Float_t         pipi_mass[20000];   //[npipi]
   Float_t         pipi_phi[20000];   //[npipi]
   Float_t         pipi_pt[20000];   //[npipi]
   Int_t           pipi_charge[20000];   //[npipi]
   Int_t           pipi_pdgId[20000];   //[npipi]
   Int_t           pipi_trk1_idx[20000];   //[npipi]
   Int_t           pipi_trk2_idx[20000];   //[npipi]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[150];   //[nGenPart]
   Float_t         GenPart_mass[150];   //[nGenPart]
   Float_t         GenPart_phi[150];   //[nGenPart]
   Float_t         GenPart_pt[150];   //[nGenPart]
   Float_t         GenPart_vx[150];   //[nGenPart]
   Float_t         GenPart_vy[150];   //[nGenPart]
   Float_t         GenPart_vz[150];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[150];   //[nGenPart]
   Int_t           GenPart_pdgId[150];   //[nGenPart]
   Int_t           GenPart_status[150];   //[nGenPart]
   Int_t           GenPart_statusFlags[150];   //[nGenPart]
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
   UInt_t          nMuon;
   Float_t         Muon_eta[15];   //[nMuon]
   Float_t         Muon_mass[15];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[15];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[15];   //[nMuon]
   Float_t         Muon_phi[15];   //[nMuon]
   Float_t         Muon_pt[15];   //[nMuon]
   Float_t         Muon_ptErr[15];   //[nMuon]
   Float_t         Muon_vx[15];   //[nMuon]
   Float_t         Muon_vy[15];   //[nMuon]
   Float_t         Muon_vz[15];   //[nMuon]
   Int_t           Muon_charge[15];   //[nMuon]
   Int_t           Muon_fired_HLT_Dimuon0_Jpsi3p5_Muon2[15];   //[nMuon]
   Int_t           Muon_fired_HLT_Dimuon10_PsiPrime_Barrel_Seagulls[15];   //[nMuon]
   Int_t           Muon_fired_HLT_Dimuon18_PsiPrime[15];   //[nMuon]
   Int_t           Muon_fired_HLT_Dimuon20_Jpsi_Barrel_Seagulls[15];   //[nMuon]
   Int_t           Muon_fired_HLT_Dimuon25_Jpsi[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu4_3_Jpsi[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu4_3_Jpsi_Displaced[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu4_JpsiTrk_Displaced[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu4_Jpsi_Displaced[15];   //[nMuon]
   Int_t           Muon_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced[15];   //[nMuon]
   Int_t           Muon_highQualityTrack[15];   //[nMuon]
   Int_t           Muon_isTriggering[15];   //[nMuon]
   Int_t           Muon_pdgId[15];   //[nMuon]
   Bool_t          Muon_isGlobal[15];   //[nMuon]
   Bool_t          Muon_looseId[15];   //[nMuon]
   Bool_t          Muon_softId[15];   //[nMuon]
   Float_t         Pileup_nTrueInt;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nProbeTracks;
   Float_t         ProbeTracks_DCASig[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_dzTrg[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_eta[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_mass[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_phi[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_pt[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_vx[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_vy[1000];   //[nProbeTracks]
   Float_t         ProbeTracks_vz[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_charge[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_isLostTrk[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_isPacked[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_nValidHits[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_pdgId[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToLooseMuon[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToMuon[1000];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToSoftMuon[1000];   //[nProbeTracks]
   UChar_t         HLT_Dimuon25_Jpsi;
   UChar_t         HLT_Dimuon20_Jpsi_Barrel_Seagulls;
   UChar_t         HLT_DoubleMu4_JpsiTrk_Displaced;
   UChar_t         HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   UChar_t         HLT_DoubleMu4_3_Jpsi_Displaced;
   UChar_t         HLT_DoubleMu4_3_Jpsi;
   UChar_t         HLT_DoubleMu4_Jpsi_Displaced;
   UChar_t         HLT_Dimuon18_PsiPrime;
   UChar_t         HLT_Dimuon10_PsiPrime_Barrel_Seagulls;
   UChar_t         HLT_DoubleMu4_PsiPrimeTrk_Displaced;
   UChar_t         HLT_Dimuon0_Jpsi3p5_Muon2;
   UChar_t         HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;
   UChar_t         HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;
   UChar_t         HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
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
   Float_t         SV_dlen[3];   //[nSV]
   Float_t         SV_dlenSig[3];   //[nSV]
   Float_t         SV_pAngle[3];   //[nSV]
   Int_t           Muon_genPartIdx[15];   //[nMuon]
   Int_t           Muon_genPartFlav[15];   //[nMuon]
   Float_t         SV_chi2[3];   //[nSV]
   Float_t         SV_eta[3];   //[nSV]
   Float_t         SV_mass[3];   //[nSV]
   Float_t         SV_ndof[3];   //[nSV]
   Float_t         SV_phi[3];   //[nSV]
   Float_t         SV_pt[3];   //[nSV]
   Float_t         SV_x[3];   //[nSV]
   Float_t         SV_y[3];   //[nSV]
   Float_t         SV_z[3];   //[nSV]
   Int_t           ProbeTracks_genPartIdx[1000];   //[nProbeTracks]
   Int_t           ProbeTracks_genPartFlav[1000];   //[nProbeTracks]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nB0;   //!
   TBranch        *b_B0_cosAlpha_PV;   //!
   TBranch        *b_B0_dxySignMu1;   //!
   TBranch        *b_B0_dxySignMu2;   //!
   TBranch        *b_B0_dxySignTr1;   //!
   TBranch        *b_B0_dxySignTr2;   //!
   TBranch        *b_B0_dzSignMu1;   //!
   TBranch        *b_B0_dzSignMu2;   //!
   TBranch        *b_B0_dzSignTr1;   //!
   TBranch        *b_B0_dzSignTr2;   //!
   TBranch        *b_B0_fit_eta;   //!
   TBranch        *b_B0_fit_mass;   //!
   TBranch        *b_B0_fit_phi;   //!
   TBranch        *b_B0_fit_pt;   //!
   TBranch        *b_B0_fitted_JPsi_mass;   //!
   TBranch        *b_B0_fitted_Rho_mass;   //!
   TBranch        *b_B0_fitted_X_mass;   //!
   TBranch        *b_B0_fitted_mass_womc;   //!
   TBranch        *b_B0_lxySign_PV;   //!
   TBranch        *b_B0_svchi2;   //!
   TBranch        *b_B0_svprob;   //!
   TBranch        *b_B0_dilepton_idx;   //!
   TBranch        *b_B0_dipion_idx;   //!
   TBranch        *b_B0_k0short_idx;   //!
   TBranch        *b_B0_mu1_idx;   //!
   TBranch        *b_B0_mu2_idx;   //!
   TBranch        *b_B0_pv_idx;   //!
   TBranch        *b_B0_trk1Rho_idx;   //!
   TBranch        *b_B0_trk2Rho_idx;   //!
   TBranch        *b_nJPsiToMuMu;   //!
   TBranch        *b_JPsiToMuMu_eta;   //!
   TBranch        *b_JPsiToMuMu_fitted_mass;   //!
   TBranch        *b_JPsiToMuMu_fitted_massErr;   //!
   TBranch        *b_JPsiToMuMu_mass;   //!
   TBranch        *b_JPsiToMuMu_phi;   //!
   TBranch        *b_JPsiToMuMu_pt;   //!
   TBranch        *b_JPsiToMuMu_svchi2;   //!
   TBranch        *b_JPsiToMuMu_svndof;   //!
   TBranch        *b_JPsiToMuMu_svprob;   //!
   TBranch        *b_JPsiToMuMu_charge;   //!
   TBranch        *b_JPsiToMuMu_l1_idx;   //!
   TBranch        *b_JPsiToMuMu_l2_idx;   //!
   TBranch        *b_JPsiToMuMu_pdgId;   //!
   TBranch        *b_nK0s;   //!
   TBranch        *b_K0s_fitted_eta;   //!
   TBranch        *b_K0s_fitted_mass;   //!
   TBranch        *b_K0s_fitted_mass_womc;   //!
   TBranch        *b_K0s_fitted_phi;   //!
   TBranch        *b_K0s_fitted_pt;   //!
   TBranch        *b_K0s_prefit_mass;   //!
   TBranch        *b_K0s_svprob;   //!
   TBranch        *b_K0s_trk1_eta;   //!
   TBranch        *b_K0s_trk1_phi;   //!
   TBranch        *b_K0s_trk1_pt;   //!
   TBranch        *b_K0s_trk2_eta;   //!
   TBranch        *b_K0s_trk2_phi;   //!
   TBranch        *b_K0s_trk2_pt;   //!
   TBranch        *b_K0s_trk1_charge;   //!
   TBranch        *b_K0s_trk2_charge;   //!
   TBranch        *b_npipi;   //!
   TBranch        *b_pipi_eta;   //!
   TBranch        *b_pipi_mass;   //!
   TBranch        *b_pipi_phi;   //!
   TBranch        *b_pipi_pt;   //!
   TBranch        *b_pipi_charge;   //!
   TBranch        *b_pipi_pdgId;   //!
   TBranch        *b_pipi_trk1_idx;   //!
   TBranch        *b_pipi_trk2_idx;   //!
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
   TBranch        *b_Muon_fired_HLT_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_Muon_fired_HLT_Dimuon10_PsiPrime_Barrel_Seagulls;   //!
   TBranch        *b_Muon_fired_HLT_Dimuon18_PsiPrime;   //!
   TBranch        *b_Muon_fired_HLT_Dimuon20_Jpsi_Barrel_Seagulls;   //!
   TBranch        *b_Muon_fired_HLT_Dimuon25_Jpsi;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu4_3_Jpsi;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu4_3_Jpsi_Displaced;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu4_Jpsi_Displaced;   //!
   TBranch        *b_Muon_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced;   //!
   TBranch        *b_Muon_highQualityTrack;   //!
   TBranch        *b_Muon_isTriggering;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nProbeTracks;   //!
   TBranch        *b_ProbeTracks_DCASig;   //!
   TBranch        *b_ProbeTracks_dzTrg;   //!
   TBranch        *b_ProbeTracks_eta;   //!
   TBranch        *b_ProbeTracks_mass;   //!
   TBranch        *b_ProbeTracks_phi;   //!
   TBranch        *b_ProbeTracks_pt;   //!
   TBranch        *b_ProbeTracks_vx;   //!
   TBranch        *b_ProbeTracks_vy;   //!
   TBranch        *b_ProbeTracks_vz;   //!
   TBranch        *b_ProbeTracks_charge;   //!
   TBranch        *b_ProbeTracks_isLostTrk;   //!
   TBranch        *b_ProbeTracks_isPacked;   //!
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
   TBranch        *b_SV_pAngle;   //!
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
   TBranch        *b_ProbeTracks_genPartIdx;   //!
   TBranch        *b_ProbeTracks_genPartFlav;   //!

   PrepAnalysisPsi2S(TTree *tree=0);
   virtual ~PrepAnalysisPsi2S();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const TString& which_analysis);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //...... my functions ......//
   virtual Int_t    analysis_code(const TString& which_analysis);
   
   virtual void     gen_K0s_ptetaphi(const int &Mother_pdgId, ROOT::Math::PtEtaPhiMVector* V); 
   virtual void     gen_Daugh_ptetaphiM(const Int_t& Mother_pdgId ,ROOT::Math::PtEtaPhiMVector* V_neg,ROOT::Math::PtEtaPhiMVector* V_pos );
	virtual void     GenPartFillP4();
   virtual void     GenPart_FillKinHist(ROOT::Math::PtEtaPhiMVector* GenVec, TH1* h_pt, TH1* h_eta, TH1* h_mass);
   

   virtual void     find_minDeltaR(const int& Trk_idx, const ROOT::Math::PtEtaPhiMVector& gen4V, double* DR_min, int* DR_min_idx, const int analysis_code);
   virtual double   score_PT(const int track_idx, const ROOT::Math::PtEtaPhiMVector* gen4V, const int analysis_code = 1);
   virtual bool     isReconstructedTrack(const int& DR_min_idx, const double& DR_min, const double& DpT_DRmin, const int& analysis_code);
   virtual void     OtherMuons(const UInt_t& DR_neg_idx,const double& DRmin_neg, const UInt_t& DR_pos_idx,const double& DRmin_pos, const bool& isReco_neg, const bool& isReco_pos, TH1* h_pt_mu, TH1* h_eta_mu, TH1* h_pt_mu_MatchedGen, TH1* h_eta_mu_MatchedGen, TH1* h_N_mu_MatchedGen);
   virtual void     MuTrkQualityCheck(const UInt_t& DR_neg_idx, const UInt_t& DR_pos_idx, const bool& isReco_neg, const bool& isReco_pos, TH1* h_SoftId_Reco, TH1* h_SoftId_Disc, TH1* h_LooseId_Reco, TH1* h_LooseId_Disc,TH1* h_isGlobal_Reco, TH1* h_isGlobal_Disc);
   virtual void     PiTrkQualityCheck(const UInt_t& DR_neg_idx, const UInt_t& DR_pos_idx, const bool& isReco_neg, const bool& isReco_pos, TH1* h_isSoftMu_Reco, TH1* h_isSoftMu_Disc, TH1* h_isLooseMu_Reco, TH1* h_isLooseMu_Disc,TH1* h_isMu_Reco, TH1* h_isMu_Disc);
   
   virtual double   MyInvMass(const int idx1, const int idx2, const int& analysis_code);
	
	virtual void PrintDecayChain();

	private:
	
	// PDG particle id & mass
	int isMum = 13, isPip = 211, isJPsi = 443, isRho = 113, isK0s = 310, isPsi2S= 100443, isB0 = 511;
	const float mMuon = 0.105658, mPion = 0.1395704, mK0s = 0.497648; 
	// generator ptls 4-vectors
	ROOT::Math::PtEtaPhiMVector GenP4_Mum, GenP4_Mup; // muons
	ROOT::Math::PtEtaPhiMVector GenP4_JPsi;           //JPsi
	ROOT::Math::PtEtaPhiMVector GenP4_Pim, GenP4_Pip; // pions
	ROOT::Math::PtEtaPhiMVector GenP4_Rho;            // Rho 
	ROOT::Math::PtEtaPhiMVector GenP4_Psi2S;          // Psi2S
	ROOT::Math::PtEtaPhiMVector GenP4_K0s;            // K0s
	ROOT::Math::PtEtaPhiMVector GenP4_B0;             // K0s



};

#endif

#ifdef PrepAnalysisPsi2S_cxx
PrepAnalysisPsi2S::PrepAnalysisPsi2S(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKs_2022Mar30/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220330_155143/0000/xNANO_mc_2022Mar30_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/group/phys_egamma/crovelli/LowPtEle/B0ToXKs_2022Mar30/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220330_155143/0000/xNANO_mc_2022Mar30_1.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

PrepAnalysisPsi2S::~PrepAnalysisPsi2S()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PrepAnalysisPsi2S::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PrepAnalysisPsi2S::LoadTree(Long64_t entry)
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

void PrepAnalysisPsi2S::Init(TTree *tree)
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
   fChain->SetBranchAddress("B0_cosAlpha_PV", B0_cosAlpha_PV, &b_B0_cosAlpha_PV);
   fChain->SetBranchAddress("B0_dxySignMu1", B0_dxySignMu1, &b_B0_dxySignMu1);
   fChain->SetBranchAddress("B0_dxySignMu2", B0_dxySignMu2, &b_B0_dxySignMu2);
   fChain->SetBranchAddress("B0_dxySignTr1", B0_dxySignTr1, &b_B0_dxySignTr1);
   fChain->SetBranchAddress("B0_dxySignTr2", B0_dxySignTr2, &b_B0_dxySignTr2);
   fChain->SetBranchAddress("B0_dzSignMu1", B0_dzSignMu1, &b_B0_dzSignMu1);
   fChain->SetBranchAddress("B0_dzSignMu2", B0_dzSignMu2, &b_B0_dzSignMu2);
   fChain->SetBranchAddress("B0_dzSignTr1", B0_dzSignTr1, &b_B0_dzSignTr1);
   fChain->SetBranchAddress("B0_dzSignTr2", B0_dzSignTr2, &b_B0_dzSignTr2);
   fChain->SetBranchAddress("B0_fit_eta", B0_fit_eta, &b_B0_fit_eta);
   fChain->SetBranchAddress("B0_fit_mass", B0_fit_mass, &b_B0_fit_mass);
   fChain->SetBranchAddress("B0_fit_phi", B0_fit_phi, &b_B0_fit_phi);
   fChain->SetBranchAddress("B0_fit_pt", B0_fit_pt, &b_B0_fit_pt);
   fChain->SetBranchAddress("B0_fitted_JPsi_mass", B0_fitted_JPsi_mass, &b_B0_fitted_JPsi_mass);
   fChain->SetBranchAddress("B0_fitted_Rho_mass", B0_fitted_Rho_mass, &b_B0_fitted_Rho_mass);
   fChain->SetBranchAddress("B0_fitted_X_mass", B0_fitted_X_mass, &b_B0_fitted_X_mass);
   fChain->SetBranchAddress("B0_fitted_mass_womc", B0_fitted_mass_womc, &b_B0_fitted_mass_womc);
   fChain->SetBranchAddress("B0_lxySign_PV", B0_lxySign_PV, &b_B0_lxySign_PV);
   fChain->SetBranchAddress("B0_svchi2", B0_svchi2, &b_B0_svchi2);
   fChain->SetBranchAddress("B0_svprob", B0_svprob, &b_B0_svprob);
   fChain->SetBranchAddress("B0_dilepton_idx", B0_dilepton_idx, &b_B0_dilepton_idx);
   fChain->SetBranchAddress("B0_dipion_idx", B0_dipion_idx, &b_B0_dipion_idx);
   fChain->SetBranchAddress("B0_k0short_idx", B0_k0short_idx, &b_B0_k0short_idx);
   fChain->SetBranchAddress("B0_mu1_idx", B0_mu1_idx, &b_B0_mu1_idx);
   fChain->SetBranchAddress("B0_mu2_idx", B0_mu2_idx, &b_B0_mu2_idx);
   fChain->SetBranchAddress("B0_pv_idx", B0_pv_idx, &b_B0_pv_idx);
   fChain->SetBranchAddress("B0_trk1Rho_idx", B0_trk1Rho_idx, &b_B0_trk1Rho_idx);
   fChain->SetBranchAddress("B0_trk2Rho_idx", B0_trk2Rho_idx, &b_B0_trk2Rho_idx);
   fChain->SetBranchAddress("nJPsiToMuMu", &nJPsiToMuMu, &b_nJPsiToMuMu);
   fChain->SetBranchAddress("JPsiToMuMu_eta", JPsiToMuMu_eta, &b_JPsiToMuMu_eta);
   fChain->SetBranchAddress("JPsiToMuMu_fitted_mass", JPsiToMuMu_fitted_mass, &b_JPsiToMuMu_fitted_mass);
   fChain->SetBranchAddress("JPsiToMuMu_fitted_massErr", JPsiToMuMu_fitted_massErr, &b_JPsiToMuMu_fitted_massErr);
   fChain->SetBranchAddress("JPsiToMuMu_mass", JPsiToMuMu_mass, &b_JPsiToMuMu_mass);
   fChain->SetBranchAddress("JPsiToMuMu_phi", JPsiToMuMu_phi, &b_JPsiToMuMu_phi);
   fChain->SetBranchAddress("JPsiToMuMu_pt", JPsiToMuMu_pt, &b_JPsiToMuMu_pt);
   fChain->SetBranchAddress("JPsiToMuMu_svchi2", JPsiToMuMu_svchi2, &b_JPsiToMuMu_svchi2);
   fChain->SetBranchAddress("JPsiToMuMu_svndof", JPsiToMuMu_svndof, &b_JPsiToMuMu_svndof);
   fChain->SetBranchAddress("JPsiToMuMu_svprob", JPsiToMuMu_svprob, &b_JPsiToMuMu_svprob);
   fChain->SetBranchAddress("JPsiToMuMu_charge", JPsiToMuMu_charge, &b_JPsiToMuMu_charge);
   fChain->SetBranchAddress("JPsiToMuMu_l1_idx", JPsiToMuMu_l1_idx, &b_JPsiToMuMu_l1_idx);
   fChain->SetBranchAddress("JPsiToMuMu_l2_idx", JPsiToMuMu_l2_idx, &b_JPsiToMuMu_l2_idx);
   fChain->SetBranchAddress("JPsiToMuMu_pdgId", JPsiToMuMu_pdgId, &b_JPsiToMuMu_pdgId);
   fChain->SetBranchAddress("nK0s", &nK0s, &b_nK0s);
   fChain->SetBranchAddress("K0s_fitted_eta", K0s_fitted_eta, &b_K0s_fitted_eta);
   fChain->SetBranchAddress("K0s_fitted_mass", K0s_fitted_mass, &b_K0s_fitted_mass);
   fChain->SetBranchAddress("K0s_fitted_mass_womc", K0s_fitted_mass_womc, &b_K0s_fitted_mass_womc);
   fChain->SetBranchAddress("K0s_fitted_phi", K0s_fitted_phi, &b_K0s_fitted_phi);
   fChain->SetBranchAddress("K0s_fitted_pt", K0s_fitted_pt, &b_K0s_fitted_pt);
   fChain->SetBranchAddress("K0s_prefit_mass", K0s_prefit_mass, &b_K0s_prefit_mass);
   fChain->SetBranchAddress("K0s_svprob", K0s_svprob, &b_K0s_svprob);
   fChain->SetBranchAddress("K0s_trk1_eta", K0s_trk1_eta, &b_K0s_trk1_eta);
   fChain->SetBranchAddress("K0s_trk1_phi", K0s_trk1_phi, &b_K0s_trk1_phi);
   fChain->SetBranchAddress("K0s_trk1_pt", K0s_trk1_pt, &b_K0s_trk1_pt);
   fChain->SetBranchAddress("K0s_trk2_eta", K0s_trk2_eta, &b_K0s_trk2_eta);
   fChain->SetBranchAddress("K0s_trk2_phi", K0s_trk2_phi, &b_K0s_trk2_phi);
   fChain->SetBranchAddress("K0s_trk2_pt", K0s_trk2_pt, &b_K0s_trk2_pt);
   fChain->SetBranchAddress("K0s_trk1_charge", K0s_trk1_charge, &b_K0s_trk1_charge);
   fChain->SetBranchAddress("K0s_trk2_charge", K0s_trk2_charge, &b_K0s_trk2_charge);
   fChain->SetBranchAddress("npipi", &npipi, &b_npipi);
   fChain->SetBranchAddress("pipi_eta", pipi_eta, &b_pipi_eta);
   fChain->SetBranchAddress("pipi_mass", pipi_mass, &b_pipi_mass);
   fChain->SetBranchAddress("pipi_phi", pipi_phi, &b_pipi_phi);
   fChain->SetBranchAddress("pipi_pt", pipi_pt, &b_pipi_pt);
   fChain->SetBranchAddress("pipi_charge", pipi_charge, &b_pipi_charge);
   fChain->SetBranchAddress("pipi_pdgId", pipi_pdgId, &b_pipi_pdgId);
   fChain->SetBranchAddress("pipi_trk1_idx", pipi_trk1_idx, &b_pipi_trk1_idx);
   fChain->SetBranchAddress("pipi_trk2_idx", pipi_trk2_idx, &b_pipi_trk2_idx);
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
   fChain->SetBranchAddress("Muon_fired_HLT_Dimuon0_Jpsi3p5_Muon2", Muon_fired_HLT_Dimuon0_Jpsi3p5_Muon2, &b_Muon_fired_HLT_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("Muon_fired_HLT_Dimuon10_PsiPrime_Barrel_Seagulls", Muon_fired_HLT_Dimuon10_PsiPrime_Barrel_Seagulls, &b_Muon_fired_HLT_Dimuon10_PsiPrime_Barrel_Seagulls);
   fChain->SetBranchAddress("Muon_fired_HLT_Dimuon18_PsiPrime", Muon_fired_HLT_Dimuon18_PsiPrime, &b_Muon_fired_HLT_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("Muon_fired_HLT_Dimuon20_Jpsi_Barrel_Seagulls", Muon_fired_HLT_Dimuon20_Jpsi_Barrel_Seagulls, &b_Muon_fired_HLT_Dimuon20_Jpsi_Barrel_Seagulls);
   fChain->SetBranchAddress("Muon_fired_HLT_Dimuon25_Jpsi", Muon_fired_HLT_Dimuon25_Jpsi, &b_Muon_fired_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi, &b_Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi, &b_Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05, &b_Muon_fired_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu4_3_Jpsi", Muon_fired_HLT_DoubleMu4_3_Jpsi, &b_Muon_fired_HLT_DoubleMu4_3_Jpsi);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu4_3_Jpsi_Displaced", Muon_fired_HLT_DoubleMu4_3_Jpsi_Displaced, &b_Muon_fired_HLT_DoubleMu4_3_Jpsi_Displaced);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced", Muon_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_Muon_fired_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu4_JpsiTrk_Displaced", Muon_fired_HLT_DoubleMu4_JpsiTrk_Displaced, &b_Muon_fired_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu4_Jpsi_Displaced", Muon_fired_HLT_DoubleMu4_Jpsi_Displaced, &b_Muon_fired_HLT_DoubleMu4_Jpsi_Displaced);
   fChain->SetBranchAddress("Muon_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced", Muon_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_Muon_fired_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("Muon_highQualityTrack", Muon_highQualityTrack, &b_Muon_highQualityTrack);
   fChain->SetBranchAddress("Muon_isTriggering", Muon_isTriggering, &b_Muon_isTriggering);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nProbeTracks", &nProbeTracks, &b_nProbeTracks);
   fChain->SetBranchAddress("ProbeTracks_DCASig", ProbeTracks_DCASig, &b_ProbeTracks_DCASig);
   fChain->SetBranchAddress("ProbeTracks_dzTrg", ProbeTracks_dzTrg, &b_ProbeTracks_dzTrg);
   fChain->SetBranchAddress("ProbeTracks_eta", ProbeTracks_eta, &b_ProbeTracks_eta);
   fChain->SetBranchAddress("ProbeTracks_mass", ProbeTracks_mass, &b_ProbeTracks_mass);
   fChain->SetBranchAddress("ProbeTracks_phi", ProbeTracks_phi, &b_ProbeTracks_phi);
   fChain->SetBranchAddress("ProbeTracks_pt", ProbeTracks_pt, &b_ProbeTracks_pt);
   fChain->SetBranchAddress("ProbeTracks_vx", ProbeTracks_vx, &b_ProbeTracks_vx);
   fChain->SetBranchAddress("ProbeTracks_vy", ProbeTracks_vy, &b_ProbeTracks_vy);
   fChain->SetBranchAddress("ProbeTracks_vz", ProbeTracks_vz, &b_ProbeTracks_vz);
   fChain->SetBranchAddress("ProbeTracks_charge", ProbeTracks_charge, &b_ProbeTracks_charge);
   fChain->SetBranchAddress("ProbeTracks_isLostTrk", ProbeTracks_isLostTrk, &b_ProbeTracks_isLostTrk);
   fChain->SetBranchAddress("ProbeTracks_isPacked", ProbeTracks_isPacked, &b_ProbeTracks_isPacked);
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
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", &TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", &TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", &TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_id", &TrigObj_id, &b_TrigObj_id);
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
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
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
   fChain->SetBranchAddress("ProbeTracks_genPartIdx", ProbeTracks_genPartIdx, &b_ProbeTracks_genPartIdx);
   fChain->SetBranchAddress("ProbeTracks_genPartFlav", ProbeTracks_genPartFlav, &b_ProbeTracks_genPartFlav);
   Notify();
}

Bool_t PrepAnalysisPsi2S::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PrepAnalysisPsi2S::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PrepAnalysisPsi2S::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PrepAnalysisPsi2S_cxx
