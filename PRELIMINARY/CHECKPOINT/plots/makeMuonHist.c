{
  gStyle->SetLegendTextSize(0.025);
  
  //---- Matching Histo
  draw_matching("rec_muCharge", "Muon charge reconstruction efficiency");
  draw_matching("rec_muTrk", "Muon track reconstruction efficiency");
  draw_matching("rec_mumu", "Muon pairs recostruction efficiency");

  draw_NReco("rec_Nmu");
  draw_NReco("N_mu_MatchGen", "OtherP");
  
  //---- Track Quality check
  draw_QC_histo("MuIsGlobal_Reco", "MuIsGlobal_Disc", "isGlobal");
  draw_QC_histo("MuSoftId_Reco", "MuSoftId_Disc", "Soft ID");
  draw_QC_histo("MuLooseId_Reco", "MuLooseId_Disc", "Loose ID");

  //---- pT 
  draw_two_histograms("pT_mu_Disc", "disc", "pT_mu_fromJPsi", "gen", "p_T\\ [GeV]");
  draw_two_histograms("pT_mu_fromJPsi", "gen", "pT_mu_Rec", "rec", "p_T\\ [GeV]");
  draw_two_histograms("pT_mu_Disc", "disc", "pT_mu_Rec", "rec", "p_T\\ [GeV]");
  draw_many_histo({"pT_mu_Disc", "pT_mu_Rec", "pT_mu_MatchGen"}, {"disc", "rec", "OtherP"}, "\\ p_T [GeV]");

  //--- Eta
  draw_two_histograms("Eta_mu_Disc", "disc", "Eta_mu_fromJPsi", "gen", "\\eta");
  draw_two_histograms("Eta_mu_fromJPsi", "gen", "Eta_mu_Rec", "rec", "\\eta");
  draw_two_histograms("Eta_mu_Disc", "disc", "Eta_mu_Rec", "rec", "\\eta");
  draw_many_histo({"Eta_mu_Disc", "Eta_mu_Rec", "Eta_mu_MatchGen"}, {"disc", "rec", "OtherP"}, "\\eta");

  //--- INVARIANT MASS
  draw_two_histograms("Mmumu_Gen", "gen", "Mmumu_DRmin", "rec", "m(\\mu^+ \\mu^-)[GeV]");

  //--- (2D)DeltaR MINIMO & DPT + PROIEZIONI
  draw_2Dhisto("DRmin_vs_DpT_mu");

  //--- Leading & Subleading muon
  draw_lead_pt("pT_leading_mu", "pT_subleading_mu");
    
}
