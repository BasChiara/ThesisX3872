{
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //---- Matching Histo
  draw_matching("rec_piCharge", "Track charge reconstruction efficiency");
  draw_matching("rec_piTrk", "Pion track reconstruction efficiency");
  draw_matching("rec_pipi", "Pion pairs recostruction efficiency");
  draw_single_histogram("PiMulty", "n(#pi)", "rec");
  
  draw_NReco("rec_Npi");
  
  //---- Quality Check
  draw_QC_histo("TrkIsSoftMu_Reco","TrkIsSoftMu_Disc","Is Matched to Soft Muon");
  draw_QC_histo("TrkIsLooseMu_Reco","TrkIsLooseMu_Disc","Is Matched to Loose Muon");
  draw_QC_histo("TrkIsMu_Reco","TrkIsMu_Disc","Is Matched to Muon");

  //---- pT 
  draw_two_histograms("pT_trk", "disc", "pT_pi_fromRho", "gen", "p_T\\ [GeV]", 0);
  draw_two_histograms("pT_pi_fromRho", "gen", "pT_pi_Rec", "rec", "p_T\\ [GeV]", 0);
  draw_two_histograms("pT_trk", "disc", "pT_pi_Rec", "rec", "p_T\\ [GeV]", 0);

  //--- Eta
  draw_two_histograms("Eta_trk", "disc", "Eta_pi_fromRho", "gen", "\\eta", 0);
  draw_two_histograms("Eta_pi_fromRho", "gen", "Eta_pi_Rec", "rec", "\\eta", 1);
  draw_two_histograms("Eta_trk", "disc", "Eta_pi_Rec", "rec", "\\eta", 0);

  //--- INVARIANT MASS
  draw_two_histograms("Mpipi_Gen", "gen", "Mpipi_DRmin", "rec", "m(\\pi^+ \\pi^-)[GeV]", 1);
  draw_two_histograms("Mpipi_DRmin", "rec", "Mpipi_Sel2", "sel2", "m(\\pi^+ \\pi^-)[GeV]", 0);

  //--- (2D)DeltaR MINIMO & DPT + PROIEZIONI
  draw_2Dhisto("DRmin_vs_DpT");
}
