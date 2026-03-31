metabolome_boxplot <- function(comp.name,
                               df,
                               exp = "mutant",
                               ctrl = "wt",
                               title_x = NULL,
                               text_x = NA, # 'NULL' to remove it
                               log_norm = TRUE,
                               title_y = ifelse(log_norm, "log(nmol/gr)", "nmol/gr"), # 'NULL' to remove it
                               group_lines = TRUE,
                               p_as_star = TRUE,
                               x_cex = 1,
                               box_col = "gray60",
                               jitter_col = "#496b40",
                               jitter_size = 1,
                               width = 2.01,
                               height = 2.53,
                               path_for_save_file = NULL,
                               image_formate = "pdf") {
  library(dplyr)
  library(ggplot2)
  library(ggsignif)
  library(tidyr)
  library(stringr)
  library(multcompView)
  library(rlang)
  library(ggpubr)
  library(rstatix)

  ##############

  names(df)[1] <- "X"
  data <- data.frame(colnames(df), unlist(df[df$X == comp.name, ], use.names = F))[-1, ]

  names(data) <- data.frame(colnames(df), unlist(df[df$X == comp.name, ], use.names = F))[1, ]
  names(data) <- c("X", "Y")
  data[, 2] <- as.numeric(data[, 2])

  line_wt <- grep(ctrl, data$X)
  data$X[line_wt] <- gsub("\\.[0-9]+|\\_[0-9]+", "", data$X[line_wt])

  if (group_lines) {
    line_treatment <- grep(exp, data$X)
    data$X[line_treatment] <- gsub("\\.[0-9]+|_[0-9]+", "", data$X[line_treatment])

    data <- data[c(line_wt, line_treatment), ]
    n_exp_groups <- 1
  } else {
    # if there is technical triplicates
    # (finished with 'exp.1.1' instead of 'exp.1')
    trip_char_vec <- paste0("^", exp, "\\.[0-9]+\\.[0-9]+|^", exp, "_[0-9]+\\.[0-9]+")
    triplicates_exp <- grepl(trip_char_vec, data$X)
    if (sum(triplicates_exp) > 1) {
      data$X[triplicates_exp] <- gsub("\\.[0-9]+$", "", data$X[triplicates_exp])
    }
    line_exp <- grep(paste0("^", exp, "\\.|", "^", exp, "_"), data$X)
    # exp_labels <- seq_along(line_exp)
    # data$X[line_exp] <- paste0(exp, "-", exp_labels)
    data$X[line_exp] <- gsub("\\.|_", "-", data$X[line_exp])

    n_exp_groups <- length(unique(data$X[line_exp]))
    data <- data[c(line_wt, line_exp), ]
  }

  title_main <- comp.name
  level.order <- unique(data$X)

  ############### normelize y-values
  if (log_norm) {
    data$Y <- log1p(data$Y)
  }

  ############### pValue vector
  stat.test <- data %>%
    t_test(Y ~ X, ref.group = ctrl, var.equal = T) %>% ################## problem
    add_significance("p") %>%
    add_xy_position(x = "X")
  stat.test$xmin <- 1
  stat.test$xmax <- stat.test$xmax + 1

  p_labels <- ifelse(p_as_star, "p.signif", "p")
  p_labels_size <- ifelse(p_as_star, 4, 2)


  y_max <- max(stat.test$y.position, na.rm = T) * 1.1
  y_min <- min(data$Y, na.rm = T)
  if (log_norm == T & y_max > 8.5) {
    y_min <- y_min * 0.9
  } else {
    y_min <- y_min * 0.75
  }

  #########################

  if (is.null(text_x)) {
     axis_text_x <- element_blank()
     ticks_len_x <- unit(0, "cm")
  } else {
     axis_text_x <- element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "bold", size = x_cex)
     ticks_len_x <- unit(0.1, "cm")
  }

  axis_title_x <- if (is.null(title_x)) {
     element_blank()
  } else {
     element_text(size = 12, face = "bold")
  }

  axis_title_y <- if (is.null(title_y)) {
     element_blank()
  } else {
     element_text(size = 12, face = "bold")
  }

  #########################

  if (nchar(comp.name) >= 40) {
    comp.name <- paste(substring(comp.name, 1, 35), "....", substring(comp.name, nchar(comp.name) - 1, nchar(comp.name)), sep = "")
    title_main <- comp.name
  }

  #########################
  comp.plot <- data %>%
    ggplot(aes(x = factor(x = X, level = level.order), y = Y)) +
    geom_boxplot(fill = c("gray100", rep(box_col, n_exp_groups)), colour = "black") +
    theme_classic() +
    theme( # panel.spacing = unit(2, "lines"),
      text = element_text(family = "serif"),
      axis.text.x = axis_text_x,
      axis.text.y = element_text(face = "bold", size = 8),
      axis.title.y = axis_title_y,
      axis.title.x = axis_title_x,
      axis.line = element_blank(),
      axis.ticks = element_line("black", size = 0.75),
      axis.ticks.length.x = ticks_len_x,
      axis.ticks.length.y = unit(ifelse(group_lines, -0.125, -0.1), "cm")
      # plot.title = element_text(face="bold.italic")
      # strip.text = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      # panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    labs(title = title_main, x = title_x, y = title_y) +

    # geom_jitter(color="#383838", size=0.625, alpha=0.9, width = 0.18)  +
    geom_jitter(color = jitter_col, alpha = 0.9, size = jitter_size) +
    # geom_text(aes(label = ifelse(pValue < 0.05, ifelse(pValue < 0.01, ifelse(pValue < 0.001, '***', '**'), '*'), '')),
    #          vjust = -0.5, size = 4) +
    geom_rect(aes(xmin = 0.5, xmax = length(level.order) + 0.5, ymin = y_min, ymax = y_max), # ymax = max(data[,2])*1.05),
      fill = "transparent", color = "black", size = 0.75
    ) +
    stat_pvalue_manual(stat.test, label = p_labels, label.size = p_labels_size, hide.ns = T, family = "serif", face = "italic") + # , tip.length = 0.01)
    # annotate("text", label= p_value_vector, size = 8)
    scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max))


  if (is.null(path_for_save_file)) {
    comp.plot
  } else {
    # remove "/" if its in the end of "path_for_save_file" row
    if (substr(path_for_save_file, nchar(path_for_save_file) - 1 + 1, nchar(path_for_save_file)) == "/") {
      path_for_save_file <- substr(path_for_save_file, 1, nchar(path_for_save_file) - 1)
    }

    dir.create(paste0(path_for_save_file, "/GCMS_", exp, "_", Sys.Date()), showWarnings = F)
    ggsave(paste0(path_for_save_file, "/", plot_dir_name, "/GCMS_", exp, "_", Sys.Date(), "/", comp.name, ".", image_formate),
      plot = comp.plot, width = width, height = height
    )
  }
}
