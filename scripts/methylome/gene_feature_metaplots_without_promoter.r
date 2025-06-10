library(ggplot2)
library(dplyr)

setwd("C:/Users/yonatany/OneDrive - Migal/Desktop/")

var1 <- "wt"
var2 <- "mto1"

region_names <- c("fiveUTR", "CDS", "introns", "threeUTR")

# Plotting the data
# bind data frames for metaPlot
binSize <- 10
pos.end <- binSize * length(region_names)

var1_metaPlot = list()
var2_metaPlot <- list()

for (region_name in region_names) {
    for (cntx in c("CG", "CHG", "CHH")) {
        var1_metaPlot[[region_name]][[cntx]] <- read.csv(paste0("mto1_vs_wt/metaPlots/Gene_features/features_metaPlot_tables/", var1, ".", cntx, ".", region_name, ".features.csv"))
        var2_metaPlot[[region_name]][[cntx]] <- read.csv(paste0("mto1_vs_wt/metaPlots/Gene_features/features_metaPlot_tables/", var2, ".", cntx, ".", region_name, ".features.csv"))
    }
}

for (cntx in c("CG", "CHG", "CHH")) {
    # bind data frames by feature
    var1_proportions <- c()
    var2_proportions <- c()
    for (region_name in region_names) {
        var1_proportions <- c(var1_proportions, var1_metaPlot[[region_name]][[cntx]]$Proportion)
        var2_proportions <- c(var2_proportions, var2_metaPlot[[region_name]][[cntx]]$Proportion)
    }
    v.cntx <- rbind(
        data.frame(pos = 1:pos.end, Proportion = var1_proportions, V = "V1"),
        data.frame(pos = 1:pos.end, Proportion = var2_proportions, V = "V2")
    )

    # plot configuration
    min_value <- min(v.cntx$Proportion)
    max_value <- max(v.cntx$Proportion)
    q1_value <- min_value + ((max_value - min_value) / 3)
    q2_value <- min_value + ((max_value - min_value) / 3) + ((max_value - min_value) / 3)
    # middle_value = round(mean(c(min_value, max_value)), 2)
    q_value <- min_value + ((max_value - min_value) / 4)

    legend_labels <- c(paste0(" ", var1), paste0(" ", var2))

    region_names_plot <- c("5'UTR", "CDS", "Intron", "3'UTR")
    main_title <- "Protein Coding Genes"
    br.max <- binSize * length(region_names)
    breaks_x <- seq(0, br.max, by = binSize)
    breaks_vline <- breaks_x[-c(1, length(breaks_x))]
    breaks_labels <- breaks_x[1:length(region_names)] + binSize / 2
    # breaks_and_labels <- list(breaks = seq(0,100,by=binSize), labels = region_names)

    plot_out <- ggplot(data = v.cntx, aes(x = pos, y = Proportion, color = V, group = V)) +
        geom_vline(xintercept = breaks_vline, colour = "gray", linetype = "solid", size = 0.6) +
        geom_line(size = 0.65) + # , aes(linetype = V)) +
        scale_color_manual(values = c("V1" = "gray50", "V2" = "#bf6828")) +
        theme_classic() +
        labs(title = main_title, x = "", y = paste0(cntx, " Methylation")) +
        theme(
            legend.position = "none",
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8, face = "bold")
        ) +
        scale_x_continuous(breaks = breaks_labels, labels = region_names_plot, minor_breaks = breaks_x, expand = expansion(add = c(0, 0))) +
        scale_y_continuous(
            breaks = c(min_value, q1_value, q2_value, max_value),
            labels = c(
                round(min_value, 3), round(q1_value, 3),
                round(q2_value, 3), round(max_value, 3)
            )
            # ) +
            # annotate("text",
            #    x = max(breaks_vline),
            #    y = max_value * c(0.95, 0.85),
            #    #y = c(max_value * 0.95, max_value * 0.85),
            #    label = legend_labels,
            #    hjust = 0, vjust = 0.75, size = 3.25,
            #    color = c("gray40", "#bf6828"), fontface = "bold"
        )


    svg(paste0("Genes_features_", cntx, "_metaPlot_", var2, "_vs_", var1, ".svg"), width = 2.5, height = 2, family = "serif")
    # par(mar = c(2,2,1,2))
    print(plot_out)
    dev.off()
}

message(paste("Processed average metaPlot for", length(Genes), "Protein Coding Genes\n"))
