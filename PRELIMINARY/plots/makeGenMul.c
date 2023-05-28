{
	gStyle->SetLineWidth(3);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	// HISTO TRANSVERSE MOMENTUM
  draw_pT_histo("gen_pt_mu", "mu");
  draw_pT_histo("gen_pt_pi", "pi");
  draw_pT_histo("gen_pt_JPsi", "JPsi");
  draw_pT_histo("gen_pt_Rho", "Rho");
  draw_pT_histo("gen_pt_X", "X");
  draw_pT_histo("gen_pt_Ks", "Ks");
  draw_pT_histo("gen_pt_B", "B0");
  draw_two_histograms("gen_SLpt_mu", "muSL", "gen_Lpt_mu", "muL", "p_{T}(#mu) [GeV]");

  //HISTO PSEUDORAPIDITY
  draw_Eta_histo("gen_eta_mu", "mu");
  draw_Eta_histo("gen_eta_pi", "pi");
  draw_Eta_histo("gen_eta_JPsi", "JPsi");
  draw_Eta_histo("gen_eta_Rho", "Rho");
  draw_Eta_histo("gen_eta_X", "X");
  draw_Eta_histo("gen_eta_Ks", "Ks");
  draw_Eta_histo("gen_eta_B", "B0");

  //HISTO MULTIPLICITY
  draw_Mul("nB0", "B0");
  draw_Mul("nK0s", "Ks");
  draw_Mul("nJPsi", "JPsi");
  draw_Mul("nPiPi", "PiPi");
  draw_Mul("nMuon", "mu");
  draw_Mul("nTracks", "pi");

/*  //HISTO MASS
  draw_mass_histo("gen_m_mu", "mu");
  draw_mass_histo("gen_m_pi", "pi");
  draw_mass_histo("gen_m_JPsi", "JPsi");
  draw_mass_histo("gen_m_Rho", "Rho");
  draw_mass_histo("gen_m_X", "X");
  draw_mass_histo("gen_m_Ks", "Ks");
  draw_mass_histo("gen_m_B", "B0");
*/
}
