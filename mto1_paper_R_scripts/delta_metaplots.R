library(dplyr)
library(tidyr)
library(ggplot2)

main_dir <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/mto1_vs_wt/metaPlots/TEs/metaPlot_tables"

delta_by_cntx <- function(cntx) {
    wt_cntx <- rbind(
        read.csv(paste0(main_dir, "/wt.", cntx, ".up.stream.csv")),
        read.csv(paste0(main_dir, "/wt.", cntx, ".gene.body.csv")),
        read.csv(paste0(main_dir, "/wt.", cntx, ".down.stream.csv"))
    )

    mto1_cntx <- rbind(
        read.csv(paste0(main_dir, "/mto1.", cntx, ".up.stream.csv")),
        read.csv(paste0(main_dir, "/mto1.", cntx, ".gene.body.csv")),
        read.csv(paste0(main_dir, "/mto1.", cntx, ".down.stream.csv"))
    )

    mto1_cntx$Proportion - wt_cntx$Proportion
}

cntx_bind <- cbind(
    data.frame(pos = 1:60),
    CG = delta_by_cntx("CG"),
    CHG = delta_by_cntx("CHG"),
    CHH = delta_by_cntx("CHH")
) %>%
    pivot_longer(cols = c(CG, CHG, CHH), names_to = "context", values_to = "delta")

min_value <- -0.1 # min(v.cntx.stream$delta)
max_value <- 0.1 # max(v.cntx.stream$delta)
# q1_value <- min_value + ((max_value - min_value) / 3)
# q2_value <- min_value + ((max_value - min_value) / 3) + ((max_value - min_value) / 3)
# # middle_value = round(mean(c(min_value, max_value)), 2)
# q_value <- min_value + ((max_value - min_value) / 4)


#legend_labels <- c("wt", "\nmto1")
legend_labels <- c("CG", "\nCHG", "\n\nCHH")

plot_colors <- c("#3d53b4","#3b8f3e", "#bb4949")

main_title <- "TEs"
breaks_and_labels <- list(breaks = c(1, 20, 40, 60), labels = c("    -2kb", "TSS", "TTS", "+2kb    "))

plot_out <- ggplot(data = cntx_bind, aes(x = pos, y = delta, color = context, group = context)) +
    geom_vline(xintercept = c(20, 40), colour = "gray", linetype = "solid", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "gray30", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.5) +
    theme_classic() +
    labs(
        title = main_title,
        x = "",
        y = "Δ methylation"
    ) +
    theme(
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9)
    ) +
    scale_x_continuous(breaks = breaks_and_labels$breaks, labels = breaks_and_labels$labels, expand = expansion(add = c(0, 0))) +
    scale_y_continuous(breaks = c(0, 0.05), limits = c(-0.005, 0.055)) +
    scale_color_manual(values = plot_colors) +
      annotate("text",
               x = 3,
               y = 0.054,
               label = legend_labels,
               hjust = 0, vjust = 0.75, size = 2.75, 
               color = plot_colors, fontface = "bold")


svg(paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/mto1_paper/TEs_delta_metaPlot.svg"), width = 1.88, height = 1.94, family = "serif")
print(plot_out)
dev.off()
