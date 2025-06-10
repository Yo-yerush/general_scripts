source("P:/yonatan/scripts/methylome/venn_diagram_for_DMRs.R")

venn.di("mto1_mto3_wt2",
        name1 = "mto1_vs_wt_2",
        name2 = "mto3_vs_wt_2")

venn.di("mto1_mto3_wt2_35s",
        name1 = "mto1_vs_wt_2",
        name2 = "mto3_vs_wt_2",
        name3 = "p35s_vs_ev_35s")

venn.di("mto1_mto3_wt2_35s_sse",
        name1 = "mto1_vs_wt_2",
        name2 = "mto3_vs_wt_2",
        name3 = "p35s_vs_ev_35s",
        name4 = "sse_high_vs_sse_ev",
        name5 = "sse_low_vs_sse_ev")

venn.di("sse_ev_35s",
        name1 = "sse_high_vs_ev_35s",
        name2 = "sse_low_vs_ev_35s",
        name3 = "p35s_vs_ev_35s")

venn.di("sse_35s",
        name1 = "sse_high_vs_sse_ev",
        name2 = "sse_low_vs_sse_ev",
        name3 = "p35s_vs_ev_35s")

venn.di("ev_35s_mto_wt_2",
        name1 = "mto1_vs_wt_2",
        name2 = "mto3_vs_wt_2",
        name3 = "ev_35s_vs_wt_2")

venn.di("ev_35s_mto1_wt_2",
        name1 = "mto1_vs_wt_2",
        name2 = "ev_35s_vs_wt_2")

venn.di("p35s_ev_35s_mto_wt_2",
        name1 = "mto1_vs_wt_2",
        name2 = "mto3_vs_wt_2",
        name3 = "p35s_vs_ev_35s",
        name4 = "ev_35s_vs_wt_2")

venn.di("mto1_sse_low",
        name1 = "mto1_vs_wt_2",
        name2 = "sse_low_vs_sse_ev")

venn.di("mto1_mto3_sse_low",
        name1 = "mto1_vs_wt_2",
        name2 = "mto3_vs_wt_2",
        name3 = "sse_low_vs_sse_ev")

#venn.di("mto1_mto3_35s",
#        name1 = "mto1_vs_wt",
#        name2 = "mto3_vs_wt",
#        name3 = "p35s_vs_ev_35s")

#venn.di("mto1_mto3",
#        name1 = "mto1_vs_wt",
#        name2 = "mto3_vs_wt")

#venn.di("mto1_mto3_35s_sse",
#        name1 = "mto1_vs_wt",
#        name2 = "mto3_vs_wt",
#        name3 = "p35s_vs_ev_35s",
#        name4 = "sse_high_vs_sse_ev",
#        name5 = "sse_low_vs_sse_ev")

#venn.di("mto1_mto3_35s_sse_2",
#        name1 = "mto1_vs_wt",
#        name2 = "mto3_vs_wt",
#        name3 = "p35s_vs_ev_35s",
#        name4 = "sse_high_vs_ev_35s",
#        name5 = "sse_low_vs_ev_35s")
