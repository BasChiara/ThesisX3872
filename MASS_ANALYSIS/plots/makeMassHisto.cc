{
	
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //JPsi
  draw2StackedHistos("MuMu_M", "JPsiBKG", "JPsi_M", "JPsiSGN", "M(#mu^{+} #mu^{-}) [GeV]");
  
  //Pions
  //draw_two_histograms("PiPi_M", "PiPiBKG", "Rho_M", "RhoSGN", "M(\\pi^+ \\pi^-) [GeV]");
  draw2StackedHistos("PiPi_M", "PiPiBKG", "Rho_M", "RhoSGN", "M(#pi^{+} #pi^{-}) [GeV]");
  draw_single_histogram("NPi_CutOn", "number of pion per event", "Pions");
  
  //X
  draw_single_histogram("X_M", "M(J/#psi#rho) [GeV]", "XSGN");
  //draw_single_histogram("MuMuPiPi_M", "M(\\mu^+ \\mu^-\\pi^+ \\pi^-) [GeV]", "XBBKG");
  draw_single_histogram("JPsiPi_M", "M(J/#psi #pi) [GeV]", "JPsiPi");
  //draw_single_histogram("JPsiPiPi_M", "M(J\\Psi\\pi^+ \\pi^-) [GeV]", "XBKG");
  //draw_two_histograms("X_M", "XSGN", "MuMuPiPi_M", "XBBKG", "M(#mu^{+} #mu^{-}#pi^{+} #pi^{-}) [GeV]");
  draw_two_histograms("X_M_GEN", "XGEN", "X_M", "XSGN", "M(#mu^{+} #mu^{-}#pi^{+} #pi^{-}) [GeV]", true, false);
  draw2StackedHistos("MuMuPiPi_M", "XBBKG", "X_M", "XSGN", "M(#mu^{+} #mu^{-}#pi^{+} #pi^{-}) [GeV]");
  draw2StackedHistos("JPsiPiPi_M", "XBKG", "X_M", "XSGN", "M(J/#psi#pi^{+} #pi^{-}) [GeV]");
  
  draw_two_histograms("Mdiff_JPsiPiPi_BKG", "XBKG", "Mdiff_JPsiPiPi_SGN", "XSGN", "|M(#mu^{+} #mu^{-}#pi^{+} #pi^{-}) - M(#mu^{+} #mu^{-}) + M_{J/#psi}^{PDG}| GeV", true);
  draw_two_histograms("Mdiff_JPsiPi_BKG", "XBKG", "Mdiff_JPsiPi_SGN", "XSGN", "|M(#mu^{+} #mu^{-}#pi^{+/-}) - M(#mu^{+} #mu^{-}) + M_{J/#psi}^{PDG}| GeV", true);

	// B0
  draw_two_histograms("B0_M_GEN", "BGEN", "B0_M", "BSGN", "M(#mu^{+} #mu^{-}#pi^{+} #pi^{-} K_{s}^{0}) [GeV]", true, false);
	

}
