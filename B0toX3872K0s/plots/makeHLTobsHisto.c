{

draw_multiple_histo({"FiredMCmatch_Mu_pT", "FiredMCmatch_Mu_eta", "FiredMCmatch_Mu_dr"}, {"pT", "eta", "DCA"}, {"\\ p_T \\ [GeV]", "\\eta", "DCA\\ p_\\mu -BS [cm]"}, "MuTrigger_L0");
draw_multiple_histo({"FiredMCmatch_MuMu_pT", "FiredMCmatch_MuMu_M", "FiredMCmatch_MuMu_DCA"}, {"pT", "M", "DCA"}, {"\\ p_T \\ [GeV]", "M(\\mu^+ \\mu^-)", "DCA \\ \\mu^+ \\mu^-"}, "MuTrigger_L1");
draw_multiple_histo({"FiredMCmatch_MuMu_LxySign", "FiredMCmatch_MuMu_cosAlpha", "FiredMCmatch_MuMu_SVp"}, {"sign", "cosA", "SVp"}, {"\\frac{L_{xy}}{\\sigma}","cos(\\alpha)","SV probability"} ,"MuTrigger_L2");

draw_Nfired("FiredMCmatch_K0s_vs_RhoTrk_N", "Fired tracks mother");

draw_multiple_histo({"FiredMCmatch_Trk_pT", "FiredMCmatch_Trk_eta", "FiredMCmatch_Trk_d0Sign"},{"pT", "eta", "DCA"}  , {"\\ p_T(trk) \\ [GeV]", "\\eta (trk)", "DCA trk-BS "}  ,"TrkTrigger_L1");

draw_single_histogram("FiredMCmatch_Trk_pT", "p_T(trk) \\ [GeV]", "pT");
draw_single_histogram("FiredMCmatch_Trk_d0Sign", "DCA trk-BS ", "DCA");

}
