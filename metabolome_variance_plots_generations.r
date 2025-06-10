variance_plots_metabolome <- function(genotype, control, results_path, save_plots = FALSE, save_svg = TRUE, save_pdf = FALSE) {
    library(tidyverse)
    library(ggplot2)
    library(ggbreak)
    library(car) # For statistical testing of variances

    # results_path = "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/תוצאות ממיכל דפני/תוצאות_260224/variance_plots"
    # genotype = "SSE"
    # control = "EV"
    # file_name = "sse_met_generations"

    ### override ggplot2 default colors
    scale_fill_discrete <- function(...) {
        scale_fill_manual(..., values = c("#567c54", "gray80", "gray50", "gray30"))
    }

    # Set file path
    file_path <- paste0(results_path, "/", file_name, "_GCMS.csv")

    # Read the data (using fileEncoding to handle Hebrew characters)
    data <- read.csv(file_path, check.names = FALSE, fileEncoding = "UTF-8")
    data$X = gsub(" $", "", data$X)

    # remove total AA levels
    data <- data[-grep("TOTAL|total", data$X), ]

    ###### ###### ###### ###### ###### ######

    # Convert to long format for easier analysis
    data_long <- data %>%
        pivot_longer(cols = -X, names_to = "sample", values_to = "value") # %>% mutate(value = log(value))

    # Identify sample types ('genotype'-1, 'genotype'-2, 'genotype'-3, 'control')
    data_long <- data_long %>%
        mutate(group_type = case_when(
            grepl(paste0(genotype, ".*1"), sample, ignore.case = TRUE) ~ paste0(genotype, "_T1"),
            grepl(paste0(genotype, ".*2"), sample, ignore.case = TRUE) ~ paste0(genotype, "_T2"),
            grepl(paste0(genotype, ".*3"), sample, ignore.case = TRUE) ~ paste0(genotype, "_T3"),
            grepl(paste0(genotype, ".*4"), sample, ignore.case = TRUE) ~ paste0(genotype, "_T4"),
            grepl(paste0(genotype, ".*5"), sample, ignore.case = TRUE) ~ paste0(genotype, "_T5"),
            grepl(control, sample, ignore.case = TRUE) ~ control,
            TRUE ~ genotype
        ))

#    # ANOVA test
#    AOV_fun <- function(data_long, compound) {
#        subset_data <- data_long %>% filter(X == compound)
#        stat <- aov(value ~ group_type, data = subset_data)
#        # Perform Tukey's HSD test
#        tukey <- TukeyHSD(stat, "group_type", conf.level = 0.95)
#        Tukey.levels <- tukey[["group_type"]][, 4]
#        Tukey.labels <- data.frame(multcompView::multcompLetters(Tukey.levels)["Letters"])
#
#        Tukey.labels$group_type <- rownames(Tukey.labels)
#        Tukey.labels$X <- compound
#        row.names(Tukey.labels) <- NULL
#
#        return(select(Tukey.labels, X, group_type, Letters))
#    }
#
#    aov_letters_results <- purrr::map_df(
#        unique(data_long$X), function(compound) {
#            AOV_fun(data_long, compound)
#        }
#    ) %>%
#        select(X, group_type, Letters)
#
#    # Print results
#    print("T-test results for each compound:")
#    print(t_test_results)

    ###### ###### ###### ###### ###### ######

    # Create a grouped version (all 'genotype' vs 'control')
    variance_data_grouped <- data_long %>%
        # mutate(group_type = ifelse(grepl(genotype, sample_type), genotype, sample_type)) %>%
        group_by(X, group_type) %>%
        summarise(
            variance = var(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        # add 't.test' results for the plot
        # left_join(t_test_results, by = "X") %>%
        mutate(y_pos = mean + sd + 200)

    ## p-value and stars for each compound
    # single_label_data <- variance_data_grouped %>%
    #    group_by(X) %>%
    #    summarise(
    #        # The same p_value across EV/dCGS if your t-test is genotype vs. control
    #        p_value = first(p_value),
    #        stars = first(stars),
    #        # Put the label above the highest bar (mean+sd) among the groups:
    #        y_pos = max(mean + sd, na.rm = TRUE) + 200, # <— small offset, tweak as needed
    #        .groups = "drop"
    #    )

    ###### ###### ###### ###### ###### ######

    # Plot variance comparison
    p1 <- ggplot(variance_data_grouped, aes(x = X, y = log(variance), fill = group_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.55),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(),
            panel.grid.major.x = element_blank()
        ) +
        labs(
            title = bquote("Variance of AA Level: " * italic(.(genotype)) * " vs. " * .(control)),
            x = "",
            y = "log(Variance)",
            fill = ""
        )

    # if (genotype == "mto1") {
    #  p1 <- p1 +
    #    scale_y_break(c(6e5, 3.4e6)) +
    #    scale_y_break(c(3.5e6, 3.9e6))
    # } else if (genotype == "mto3") {
    #  p1 <- p1 +
    #    scale_y_break(c(2e6, 1.175e7)) # +
    #  # scale_y_break(c(3.5e6, 3.9e6))
    # } else if (genotype == "dCGS") {
    #  p1 <- p1 +
    #    scale_y_break(c(1.75e6, 3e6)) # +
    #  # scale_y_break(c(3.5e6, 3.9e6))
    # }

    # Plot mean comparison
    p2 <- ggplot(variance_data_grouped, aes(x = X, y = mean, fill = group_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
            position = position_dodge(width = 0.9),
            width = 0.25
        ) +
        #geom_text(
        #    data = single_label_data,
        #    aes(x = X, y = y_pos, label = stars),
        #    inherit.aes = FALSE,
        #    position = position_dodge(width = 0.9)
        #) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.55),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(),
            panel.grid.major.x = element_blank(),
            # legend.position = c(1, 1), # topright legend
            # legend.justification = c(1, 1)
        ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
        labs(
            title = bquote("AA Levels: " * italic(.(genotype)) * " vs. " * .(control)),
            x = "",
            y = "AA levels (nmol/gr DW)",
            fill = ""
        )

    # Plot SD comparison
    max_p3 <- max(variance_data_grouped$sd, na.rm = TRUE) * 1.05
    p3 <- ggplot(variance_data_grouped, aes(x = X, y = sd, fill = group_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.55),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(),
            panel.grid.major.x = element_blank(),
            # legend.position = c(1, 1), # topright legend
            # legend.justification = c(1, 1)
        ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max_p3)) +
        labs(
            title = bquote("Standard Deviation of AA Levels: " * italic(.(genotype)) * " vs. " * .(control)),
            x = "",
            y = "Standard Deviation",
            fill = ""
        )

    ###### ###### ###### ###### ###### ######

    # Statistical comparison - test if variances differ between 'genotype' and 'control'
    # Function to test variance equality
    test_variance_equality <- function(compound) {
        subset_data <- data_long %>%
            filter(X == compound) %>%
            mutate(group_type = ifelse(grepl(genotype, group_type), genotype, group_type))

        # Skip if insufficient data
        if (n_distinct(subset_data$group_type) < 2) {
            return(data.frame(compound = compound, p_value = NA, significant = FALSE))
        }

        test_result <- tryCatch(
            {
                levene_test <- leveneTest(value ~ group_type, data = subset_data)
                data.frame(
                    compound = compound,
                    p_value = levene_test$`Pr(>F)`[1],
                    significant = levene_test$`Pr(>F)`[1] < 0.05
                )
            },
            error = function(e) {
                data.frame(compound = compound, p_value = NA, significant = FALSE)
            }
        )

        return(test_result)
    }

    # Apply test to all compounds
    all_compounds <- unique(data$X)
    variance_test_results <- map_df(all_compounds, test_variance_equality)

    # Create plot showing 'genotype' vs 'control' variance ratio
    variance_ratio <- variance_data_grouped %>%
        select(X, group_type, variance) %>%
        pivot_wider(names_from = group_type, values_from = variance) %>%
        mutate(
            ratio = !!sym(genotype) / !!sym(control),
            log2_ratio = log2(ratio),
            significant = X %in% (variance_test_results %>% filter(significant == TRUE) %>% pull(compound))
        )

    # Plot ratio
    p4 <- ggplot(variance_ratio, aes(x = reorder(X, log2_ratio), y = log2_ratio, fill = significant)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red3")) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.55),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(),
            panel.grid.major.x = element_blank() # ,
            # legend.position = c(0, 1), # topleft legend
            # legend.justification = c(0, 1)
        ) +
        labs(
            title = bquote("Log2 Ratio of Variance (" * italic(.(genotype)) * "/" * .(control) * ")"),
            subtitle = "Red bars indicate statistically significant difference (Levene's test)",
            x = "",
            y = "Log2(Variance Ratio)",
            fill = "Significant"
        )

    ###### ###### ###### ###### ###### ######

    # Save plots
    if (save_plots) {
        if (save_svg) {
            svg(paste0(results_path, "/", genotype, "_amino_acid_variance.svg"), width = 6, height = 4, family = "serif")
            print(p1)
            dev.off()

            svg(paste0(results_path, "/", genotype, "_amino_acid_mean.svg"), width = 6, height = 4, family = "serif")
            print(p2)
            dev.off()

            svg(paste0(results_path, "/", genotype, "_amino_acid_sd.svg"), width = 6, height = 4, family = "serif")
            print(p3)
            dev.off()

            svg(paste0(results_path, "/", genotype, "_amino_acid_variance_ratio.svg"), width = 6, height = 4, family = "serif")
            print(p4)
            dev.off()
        }

        if (save_pdf) {
            pdf(paste0(results_path, "/", genotype, "_amino_acid_plots.pdf"), width = 6, height = 4, family = "serif")
            print(p1)
            print(p2)
            print(p3)
            print(p4)
            dev.off()
        }
    }

    # Print significant compounds
    significant_compounds <- variance_test_results %>% filter(significant == TRUE)
    print(paste0("Compounds with significantly different variance between ", genotype, " and ", control, ":"))
    print(significant_compounds)

    return(list(
        variance = p1,
        mean = p2,
        sd = p3,
        variance_ratio = p4
    ))
}
