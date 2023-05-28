#define CheckGenLev_cxx
#include "CheckGenLev.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void CheckGenLev::Loop()
{


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   const ULong64_t nbreak = 100000;
   const ULong64_t entry_CP = 2000;


   //BRANCHES
   fChain->SetBranchStatus("*",0); 


   fChain->SetBranchStatus("nGenPart",1);
   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_genPartIdxMother",1);
   fChain->SetBranchStatus("GenPart_eta",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   fChain->SetBranchStatus("GenPart_pt",1);
   fChain->SetBranchStatus("GenPart_mass",1);

   // VARIABLES
   const Int_t isMu = 13, isJPsi = 443, isRho = 113, isPi = 211, isX3872 = 9920443, isKs = 310, isB0 = 511 ;


  //HISTOGRAMS

  //----Generator----//
  
  Int_t Nbins = 40;
  TH1F h_Gpt_mu("gen_pt_mu","", Nbins, 0.,20);
  TH1F h_Gpt_pi("gen_pt_pi","", Nbins, 0.,6);
  TH1F h_Gpt_JPsi("gen_pt_JPsi","", Nbins, 5.,25);
  TH1F h_Gpt_Rho("gen_pt_Rho","", Nbins, 0.,8.);
  TH1F h_Gpt_X("gen_pt_X","", Nbins, 5.,35.);
  TH1F h_Gpt_Ks("gen_pt_Ks","", Nbins, 0.,20);
  TH1F h_Gpt_B("gen_pt_B","", Nbins, 0.,30);

  TH1F h_Geta_mu("gen_eta_mu","", Nbins, -3.5,3.5);
  TH1F h_Geta_pi("gen_eta_pi","", Nbins, -3.5,3.5);
  TH1F h_Geta_JPsi("gen_eta_JPsi","", Nbins, -3.5,3.5);
  TH1F h_Geta_Rho("gen_eta_Rho","", Nbins, -3.5,3.5);
  TH1F h_Geta_X("gen_eta_X","", Nbins, -3.5,3.5);
  TH1F h_Geta_Ks("gen_eta_Ks","", Nbins, -3.5,3.5);
  TH1F h_Geta_B("gen_eta_B","", Nbins, -3.5,3.5);


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if ( (jentry + 1) % entry_CP == 0) std::cout << "=> EV " << jentry + 1<< std::endl;
      // ------- GENERATOR -------//
      GenPart_FillKinHist( isMu, isJPsi, &h_Gpt_mu, &h_Geta_mu);
      GenPart_FillKinHist( isPi, isRho, &h_Gpt_pi, &h_Geta_pi);
      GenPart_FillKinHist( isJPsi, isX3872, &h_Gpt_JPsi, &h_Geta_JPsi);
      GenPart_FillKinHist( isRho, isX3872, &h_Gpt_Rho, &h_Geta_Rho);
      GenPart_FillKinHist( isX3872 , isB0, &h_Gpt_X ,&h_Geta_X);
      GenPart_FillKinHist( isKs, isB0, &h_Gpt_Ks, &h_Geta_Ks);
      GenPart_FillKinHist( isB0, 1 ,&h_Gpt_B, &h_Geta_B);
      
      TFile gen_outfile("./plots/CheckGenLev_Livia.root","RECREATE"); // create the output file
      h_Gpt_mu.Write();
      h_Geta_mu.Write();
      h_Gpt_pi.Write();
      h_Geta_pi.Write();
      h_Gpt_JPsi.Write();
      h_Geta_JPsi.Write();
      h_Gpt_Rho.Write();
      h_Geta_Rho.Write();
      h_Gpt_X.Write();
      h_Geta_X.Write();
      h_Gpt_Ks.Write();
      h_Geta_Ks.Write();
      h_Gpt_B.Write();
      h_Geta_B.Write();
      
      gen_outfile.Close();

      // if (Cut(ientry) < 0) continue;
   }
}//Loop()


void CheckGenLev::GenPart_FillKinHist(const int& GenPtl_PDGid, const int & GenPtl_MotherPDGid, TH1* h_pt, TH1* h_eta){ 
  
  for (ULong64_t g = 0; g < nGenPart; g++){
    if ( TMath::Abs(GenPart_pdgId[g]) != GenPtl_PDGid) continue; //conto anche anti ptl
    if (GenPtl_PDGid != 511){ //se guardo i B la mother ptl non mi interessa
      if ( TMath::Abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) != GenPtl_MotherPDGid) continue;
      h_pt->Fill(GenPart_pt[g]);
      h_eta->Fill(GenPart_eta[g]);
    }else{
      h_pt->Fill(GenPart_pt[g]);
      h_eta->Fill(GenPart_eta[g]);
    }
  }
}//GenPart_FillKinHist()

