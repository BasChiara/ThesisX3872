{
  gStyle->SetLegendTextSize(0.03);

  //---- Matching Histo
  draw_matching("rec_K", "Pion track reconstruction efficiency");
  draw_single_histogram("K0sMulty", "rec", "K^{0}_{s} candidates");
   
  //---- pT 
  draw_two_histograms("pT_Kdisc", "disc", "pT_K_fromB", "gen", "p_{T} [GeV]", 0);
  draw_two_histograms("pT_K_fromB", "gen", "pT_K_Rec", "rec", "p_{T}  [GeV]", 0);
  draw_two_histograms("pT_Kdisc", "disc", "pT_K_Rec", "rec", "p_{T}  [GeV]", 0);

  //--- Eta
  draw_two_histograms("Eta_Kdisc", "disc", "Eta_K_fromB", "gen", "\\eta", 0);
  draw_two_histograms("Eta_K_fromB", "gen", "Eta_K_Rec", "rec", "\\eta", 1);
  draw_two_histograms("Eta_Kdisc", "disc", "Eta_K_Rec", "rec", "\\eta", 0);

  //--- INVARIANT MASS
  draw_single_histogram("MK0s_DRmin_fit_womc", "rec", "m(K^{0}_{s})[GeV]");
  draw_two_histograms("MK0s_Gen", "gen", "MK0s_DRmin_fit_womc", "rec", "m(K^{0}_{s})[GeV]", 1);
  draw_two_histograms("MK0s_DRmin_prefit", "prefit", "MK0s_DRmin_fit_womc", "womc", "m(K^{0}_{s})[GeV]", 1, false);
  draw_many_histo({"MK0s_DRmin_fit", "MK0s_DRmin_prefit", "MK0s_DRmin_fit_womc"}, {"rec", "prefit", "womc"},"M(K^{0}_{s})" );
  //--- (2D)DeltaR MINIMO & DPT + PROIEZIONI
  draw_2Dhisto("DRmin_vs_DpT_K0");
}
