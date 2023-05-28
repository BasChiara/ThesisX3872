{
//--> Delta R vs Delta pT
draw_2Dhisto("DRminVSDpT_Mu", "mu");
draw_2Dhisto("DRminVSDpT_Pi", "pi");
draw_2Dhisto("DRminVSDpT_K0s", "K0s");
draw_2Dhisto("DRminVSDpT_aiut_Pi", "pi");
draw_2Dhisto("DRminVSDpT_aiut_K0s", "K0s");

//--> JPsi
draw_Ncandidates("Ncandidate_JPsi", "JPsi", "");
draw2StackedHistos("PREfit_BKG_JPsi_M" , "PRE_BKG_JPsi", "PREfit_SGN_JPsi_M", "PRE_SGN_JPsi", "M(\\mu^+ \\mu^-)[GeV]");
draw2StackedHistos("VTXfit_BKG_JPsi_M" , "VTX_BKG_JPsi", "VTXfit_SGN_JPsi_M", "VTX_SGN_JPsi", "M(\\mu^+ \\mu^-)[GeV]");

//-->Rho
draw2StackedHistos("PREfit_BKG_Rho_M" , "PRE_BKG_Rho", "PREfit_SGN_Rho_M", "PRE_SGN_Rho", "M(\\pi^+ \\pi^-)[GeV]");

//--> KVTXfit_SGN_K0s_M0s 
draw2StackedHistos("VTXfit_BKG_K0s_M" , "VTX_BKG_K0s", "VTXfit_SGN_K0s_M", "VTX_SGN_K0s", "\\ M(K_0^s)[GeV]");

//--> X(3872)
draw_Ncandidates("Ncandidate_X3872", "X3872", "");
draw2StackedHistos("PREfit_BKG_X3872_M" , "PRE_BKG_X3872", "PREfit_SGN_X3872_M", "PRE_SGN_X3872", "M(\\mu^+ \\mu^- \\pi^+ \\pi^-)[GeV]");
draw2StackedHistos("TOTfit_BKG_X3872_M" , "TOT_BKG_X3872", "TOTfit_SGN_X3872_M", "TOT_SGN_X3872", "M(\\mu^+ \\mu^- \\pi^+ \\pi^-)[GeV]");

//--> B0
draw_Ncandidates("Ncandidate_B0", "B0", "");
draw2StackedHistos("TOTfit_BKG_B0_M" , "TOT_BKG_B0", "TOTfit_SGN_B0_M", "TOT_SGN_B0", "M(\\mu^+ \\mu^- \\pi^+ \\pi^- K_0^s)[GeV]");
draw2StackedHistos("VTXfit_BKG_B0_M" , "VTX_BKG_B0", "VTXfit_SGN_B0_M", "VTX_SGN_B0", "M(\\mu^+ \\mu^- \\pi^+ \\pi^- K_0^s)[GeV]");
draw2StackedHistos("PREfit_BKG_B0_M" , "PRE_BKG_B0", "PREfit_SGN_B0_M", "PRE_SGN_B0", "M(\\mu^+ \\mu^- \\pi^+ \\pi^- K_0^s)[GeV]");

//--> DEBUG
draw_ptlMCmissed("MCmatch_WhichPTLFail", "MC", "Particles with no Gen-Level matching");
draw_ptlMCmissed("B0cand_WhichPtlMCmiss", "B0", "B0 candidates daughters which are non-signal");
}
