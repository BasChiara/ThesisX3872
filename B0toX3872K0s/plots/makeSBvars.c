{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
// B0

draw2StackedHistos("BKG_B0_M" , "BKG_B0", "SGN_B0_M", "SGN_B0", "M(J\\Psi \\pi^+ \\pi^- K_0^s)[GeV]");
draw_many_histo({"BKGrk_B0_M", "BKGk_B0_M", "BKGr_B0_M", "SGN_B0_M"}, {"BKG_RK_B0", "BKG_K0s_B0", "BKG_Rho_B0", "SGN_B0"}, "M(J/#psi #pi^{+} #pi^{-} K_{0}^{s})[GeV]", "BKGdetail_B0_mass");

draw_many_2Dhisto({"BKGk_B0vsX_M", "BKG_B0vsX_M", "SGN_B0vsX_M"}, {"BKG_K0s_B0", "BKG_RK_B0", "SGN_B0"}, "M(J\\Psi \\pi^+ \\pi^-) [GeV]","M(J\\Psi \\pi^+ \\pi^- K_0^s)[GeV]", "");
draw_many_2Dhisto({"BKGr_RhovsK0s_M", "BKGk_RhovsK0s_M","BKGrk_RhovsK0s_M", "SGN_RhovsK0s_M"}, {"BKG_Rho_B0", "BKG_K0s_B0","BKG_RK_B0", "SGN_B0"}, "\\ M(K_0^s)[GeV]", "M(\\pi^+ \\pi^-) [GeV]","RhovsK0s_M");


// Rho
draw2StackedHistos("BKG_Rho_M" , "BKG_Rho", "SGN_Rho_M", "SGN_Rho", "M(\\pi^+ \\pi^-)[GeV]");
draw_two_histograms("BKGb_Rho_M" , "BKG_B0_Rho", "SGN_Rho_M", "SGN_Rho", "M(\\pi^+ \\pi^-)[GeV]");

// pT
draw_two_histograms("SGN_B0_pT", "SGN_B0", "BKG_B0_pT", "BKG_B0", "\\ p_T(B_0)/M(B_0)", false);

//draw_single_histogram("TOTfit_SGN_JPsi_pT", "TOT_SGN_JPsi", "\\ p_T(J\\Psi)/p_T(B_0)", false);

draw_two_histograms("SGN_Pi1_pT", "SGN_Pi", "BKG_Pi1_pT", "BKG_Pi", "\\ p_T(\\pi_1)/p_T(B_0)", false);
draw_two_histograms("SGN_Pi1_pT", "SGN_B0", "BKGb_Pi1_pT", "BKG_B0", "\\ p_T(\\pi_1)/p_T(B_0)", false);

draw_two_histograms("SGN_Pi2_pT", "SGN_Pi", "BKG_Pi2_pT", "BKG_Pi", "\\ p_T(\\pi_2)/p_T(B_0)", false);

draw_two_histograms("SGN_Rho_pT", "SGN_Rho", "BKG_Rho_pT", "BKG_Rho", "\\ p_T(\\rho)/p_T(B_0)", false);
draw_two_histograms("SGN_Rho_pT", "SGN_B0", "BKGb_Rho_pT", "BKG_B0", "\\ p_T(\\rho)/p_T(B_0)", false);

draw_two_histograms("SGN_X3872_pT", "SGN_X3872", "BKG_X3872_pT", "BKG_X3872", "\\ p_T(X(3872))/p_T(B_0)", false);

draw_two_histograms("SGN_K0s_pT", "SGN_K0s", "BKG_K0s_pT", "BKG_K0s", "\\ p_T(K_0^s)/p_T(B_0)", false);
draw_two_histograms("SGN_K0s_pT", "SGN_K0s", "BKGb_K0s_pT", "BKG_B0_K0s", "\\ p_T(K_0^s)/p_T(B_0)", false);


// DR(pi,B0)
draw_two_histograms("SGN_DR_Pi1B0_Rho", "SGN_Pi", "BKG_DR_Pi1B0_Rho", "BKG_Pi", "\\Delta R(\\pi_1, B_0)", false );
draw_two_histograms("SGN_DR_Pi1B0_Rho", "SGN_B0", "BKGb_DR_Pi1B0_Rho", "BKG_B0", "\\Delta R(\\pi_1, B_0)", false );
draw_two_histograms("SGN_DR_Pi2B0_Rho", "SGN_Pi", "BKG_DR_Pi2B0_Rho", "BKG_Pi", "\\Delta R(\\pi_2, B_0)", false );
draw_two_histograms("SGN_DR_Pi1B0_K0s", "SGN_K0trk", "BKG_DR_Pi1B0_K0s", "BKG_K0trk", "\\Delta R(\\pi_1, B_0)", false );
draw_two_histograms("SGN_DR_Pi2B0_K0s", "SGN_K0trk", "BKG_DR_Pi2B0_K0s", "BKG_K0trk", "\\Delta R(\\pi_2, B_0)", false );

// Rho vertex
draw_two_histograms("SGN_Rho_D0", "SGN_Rho", "BKG_Rho_D0", "BKG_Rho", "D0 sig.",false);
draw_two_histograms("SGN_Rho_D0", "SGN_B0", "BKGb_Rho_D0", "BKG_B0", "D0 sig.",false);
draw_two_histograms("SGN_Rho_D0max", "SGN_Pi", "BKG_Rho_D0max", "BKG_Pi", "D0 MAX [cm]",false);
draw_two_histograms("SGN_Rho_D0max", "SGN_Pi", "BKGb_Rho_D0max", "BKG_B0", "D0 MAX [cm]",false);
draw_two_histograms("SGN_Rho_D0min", "SGN_Pi", "BKG_Rho_D0min", "BKG_Pi", "D0 min [cm]",false);
draw_two_histograms("SGN_Rho_D0min", "SGN_Pi", "BKGb_Rho_D0min", "BKG_B0", "D0 min [cm]",false);

// K0s vertex
draw_two_histograms("SGN_K0s_SVp", "SGN_K0s", "BKG_K0s_SVp", "BKG_K0s", "Prob(SV)", false);
draw_two_histograms("SGN_K0s_D0", "SGN_K0s", "BKG_K0s_D0", "BKG_K0s", "D0 sig.",false);
draw_two_histograms("SGN_K0s_D0", "SGN_K0s", "BKGb_K0s_D0", "BKG_B0_K0s", "D0 sig.",false);
draw_many_2Dhisto({"SGN_K0s_D0vsPt", "BKG_K0s_D0vsPt"}, {"SGN_K0s", "BKG_K0s"}, "\\ p_T(K_0^s)/p_T(B_0)", "D0 sig.", "");

draw_two_histograms("SGN_LeadTrk_pT", "SGN_Pi", "BKG_LeadTrk_pT", "BKG_Pi", "\\ p_T(\\pi \\ leading)/p_T(B_0)", false);

// B0 vertex
draw_two_histograms("SGN_B0_LxySign", "SGN_B0", "BKG_B0_LxySign", "BKG_B0", "\\frac{L_{xy}}{\\sigma_{xy}}", false);
draw_two_histograms("SGN_B0_SVchi2", "SGN_B0", "BKG_B0_SVchi2", "BKG_B0", "\\Chi^2(SV)", false);
draw_two_histograms("SGN_B0_SVp", "SGN_B0", "BKG_B0_SVp", "BKG_B0", "Prob(SV)", false);
draw_two_histograms("SGN_B0_cosA", "SGN_B0", "BKG_B0_cosA", "BKG_B0", "cos(\\alpha)", false, true);

draw_two_2Dhisto("SGN_B0_rVSz_decayV" , "SGN_B0", "BKG_B0_rVSz_decayV", "BKG_B0", "z_{decayV} [cm]", "r_{decayV} [cm]");

//DELTA phi B0 rest frame
draw_two_histograms("BKG_Dphi_B0RF", "BKG_Dphi", "SGN_Dphi_B0RF", "SGN_Dphi", "\\Delta\\phi(X,K_0^s)", false, true);

}
