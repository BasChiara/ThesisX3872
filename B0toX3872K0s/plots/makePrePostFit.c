{
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	draw_many_histo({"MCmatch_TOTfit_JPsi_M","MCmatch_VTXfit_JPsi_M", "MCmatch_prefit_JPsi_M" }, {"TOT", "VTX", "PRE"}, "M(\\mu^+ \\mu^-)", true, false);

	draw_two_histograms("MCmatch_prefit_JPsi_M", "PRE", "MCmatch_VTXfit_JPsi_M", "VTX", "M(\\mu^+\\mu^-)", false, false);

	draw_two_histograms("MCmatch_prefit_Rho_M", "PRE", "MCmatch_TOTfit_Rho_M", "TOT", "M(\\pi^+\\pi^-)", false, false);
	draw_single_histogram("MCmatch_gen_Rho_M", "M(\\pi^+\\pi^-)", "PRE");
	draw_single_histogram("genPi_pt", "p_T(\\pi)", "PRE");
	draw_single_histogram("genPi_eta", "\\eta (\\pi)", "PRE");

	draw_two_histograms("MCmatch_TOTfit_X3872_M", "TOTX", "MCmatch_prefit_X3872_M", "PREX", "M(#mu^{+} #mu^{-} #pi^{+} #pi^{-})", true, true);
	draw_two_histograms("MCmatch_TOTfit_B0_M", "TOTB", "MCmatch_prefit_B0_M", "PREB", "M(#mu^{+} #mu^{-}#pi^{+} #pi^{-} K_{s}^{0})", true, true);

	draw_two_histograms("MCmatch_prefit_K0s_M", "PRE", "MCmatch_VTXfit_K0s_M", "VTX", "\\ M(K^0_s)", false, false);

	draw_many_histo({"MCmatch_TOTfit_K0s_M","MCmatch_VTXfit_K0s_M", "MCmatch_prefit_K0s_M" }, {"TOT", "VTX", "PRE"}, "\\ M(K^0_s)", true, false);

	draw_many_histo({"MCmatch_TOTfit_B0_M","MCmatch_WOMCfit_B0_M", "MCmatch_prefit_B0_M" }, {"TOT", "VTX", "PRE"}, "M(\\mu^+\\mu^- \\pi^+\\pi^- K^0_s)", true, false);






}

