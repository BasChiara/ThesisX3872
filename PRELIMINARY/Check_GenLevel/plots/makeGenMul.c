{
  gStyle->SetLegendTextSize(0.04);
  // HISTO TRANSVERSE MOMENTUM
  draw_pT_histo("gen_pt_mu", "mu");
  draw_pT_histo("gen_pt_pi", "pi");
  draw_pT_histo("gen_pt_JPsi", "JPsi");
  draw_pT_histo("gen_pt_Rho", "Rho");
  draw_pT_histo("gen_pt_X", "X");
  draw_pT_histo("gen_pt_Ks", "Ks");
  draw_pT_histo("gen_pt_B", "B0");

  //HISTO PSEUDORAPIDITY
  draw_Eta_histo("gen_eta_mu", "mu");
  draw_Eta_histo("gen_eta_pi", "pi");
  draw_Eta_histo("gen_eta_JPsi", "JPsi");
  draw_Eta_histo("gen_eta_Rho", "Rho");
  draw_Eta_histo("gen_eta_X", "X");
  draw_Eta_histo("gen_eta_Ks", "Ks");
  draw_Eta_histo("gen_eta_B", "B0");

}
