library(scales)
library(ggplot2)
library(dplyr)



ec_squere_plots <- function(treatment=NULL, by_folder="by_DEseq2") {
    res_dir <- paste0("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/EC_kegg_maps/", by_folder, "/with_EC")


    output_dir <- paste0(res_dir, "/square_plots/", treatment)
    dir.create(paste0(res_dir, "/square_plots/"), showWarnings = F)
    dir.create(output_dir, showWarnings = F)


    # col_fun index
    col_fun_up <- scales::col_numeric(palette = c("white", "red"), domain = c(0, 2))
    col_fun_down <- scales::col_numeric(palette = c("blue", "white"), domain = c(-2, 0))

    map_ids <- c(
        "03020",
        "04122",
        "04712",
        "00270",
        "00300",
        "00310",
        "00340",
        "00920",
        "00020",
        "00500",
        "00520",
        "00260",
        "00290",
        "00330",
        "00250",
        "00620"
    )

    for (n.map in map_ids) {
        pathway_name <- paste0("ath", n.map)

        file_full_name <- list.files(
            path = paste0(res_dir, "/", treatment),
            pattern = paste0(".*_", pathway_name, "_.*csv"),
            full.names = T
        )

        # check if the file exists
        if (length(file_full_name) == 1) {
            dir.create(paste0(output_dir, "/", pathway_name), showWarnings = F)

            path_res <- read.csv(file_full_name) %>% filter(pValue < 0.05)

            # change max and min limits to +-2
            path_res$log2FC[path_res$log2FC > 2] <- 2
            path_res$log2FC[path_res$log2FC < -2] <- -2

            # add colors
            path_res$col[path_res$log2FC > 0] <- col_fun_up(path_res$log2FC[path_res$log2FC > 0])
            path_res$col[path_res$log2FC < 0] <- col_fun_down(path_res$log2FC[path_res$log2FC < 0])


            # group by EC and plot it
            unique_ecs <- unique(path_res$ec)
            for (ec.i in unique_ecs) {
                ec_res <- path_res %>% filter(ec == ec.i)

                log2fc_df <- data.frame(
                    xx = 1:nrow(ec_res),
                    yy = ec_res$log2FC,
                    zz = ec_res$col
                )

                p <- ggplot(log2fc_df, aes(x = xx, y = 1)) +
                    geom_tile(aes(fill = I(zz)), color = "black", size = 1) +
                    theme_void() +
                    theme(
                        legend.position = "none",
                        plot.margin = margin(-1, -3.5, -1, -3.5),
                        panel.spacing = unit(0, "lines")
                    )

                svg(paste0(output_dir, "/", pathway_name, "/", ec.i, ".svg"),
                    width = 1.305,
                    height = 0.5
                )
                print(p)
                dev.off()
            }
        }
    }
}

ec_squere_plots("mto1")
ec_squere_plots("mto3")
ec_squere_plots("dCGS")
ec_squere_plots("SSE_high")
ec_squere_plots("SSE_low")
ec_squere_plots("high_vs_low")


######################### old color_fun shit
# {
#    # add color gradient
#    pathway_up <- path_res %>% filter(log2FC > 0)
#    pathway_down <- path_res %>% filter(log2FC < 0)
#
#    # make NULL col_fun
#    col_fun_up <- function(x) NULL
#    col_fun_down <- function(x) NULL
#
#    # make col_fun (when its not enpty...)
#    if (nrow(pathway_up) != 0) {
#        col_fun_up <- scales::col_numeric(palette = c("white", "red"), domain = range(pathway_up$log2FC))
#    }
#
#    if (nrow(pathway_down) != 0) {
#        col_fun_down <- scales::col_numeric(palette = c("blue", "white"), domain = range(pathway_down$log2FC))
#    }
#
#    # adjust color if lower than first quantile
#    q_up <- min(pathway_up$log2FC) + ((max(pathway_up$log2FC) - min(pathway_up$log2FC)) * 0.25)
#    q_down <- max(pathway_up$log2FC) + ((min(pathway_up$log2FC) - max(pathway_up$log2FC)) * 0.25)
#    q_down <- pathway_down %>%
#        .$log2FC %>%
#        min() * 0.25
#
#    path_res$col[path_res$log2FC > 0] <- col_fun_up(pathway_up$log2FC)
#    path_res$col[path_res$log2FC < 0] <- col_fun_down(pathway_down$log2FC)
#
#    path_res$col[path_res$log2FC <= q_up & path_res$log2FC > 0] <- col_fun_up(q_up)
#    path_res$col[path_res$log2FC >= q_down & path_res$log2FC < 0] <- col_fun_down(q_down)
#    #############
# }
