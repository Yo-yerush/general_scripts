GCMS_met23_plot_save <- function(comp.name,
                                 exp = "mto1",
                                 ctrl = "wt",
                                 title_x = "",
                                 title_y = "log(nmol/gr)",
                                 file_for_metabo,
                                 group_lines = TRUE,
                                 p_as_star = TRUE,
                                 width = 2.01,
                                 height = 2.53,
                                 x_cex = 1,
                                 box_col = "gray60",
                                 jitter_col = "#496b40",
                                 log_norm = TRUE,
                                 path_for_save_file = "C:/Users/yonye/Migal/Rachel Amir Team - General/yonatan/methionine/GCMS_23/") {
  library(dplyr)
  library(ggplot2)
  library(ggsignif)
  library(tidyr)
  library(stringr)
  library(multcompView)
  library(rlang)
  library(ggpubr)
  library(rstatix)

  # remove "/" if its in the end of "path_for_save_file" row
  if (substr(path_for_save_file, nchar(path_for_save_file) - 1 + 1, nchar(path_for_save_file)) == "/") {
    path_for_save_file <- substr(path_for_save_file, 1, nchar(path_for_save_file) - 1)
  }

  # if there is a pthway name in the file from lc_pipline
  plot_dir_name <- "plots"

  table <- file_for_metabo
  # names(table) = gsub("\\.[0-9+]", "", names(table))
  names(table)[1] <- "X"

  ###
  ###

  comp.name <- comp.name

  ###
  ###

  y.axis.hight <- 1.1
  y.label.hight <- 1

  bar.width <- 0.5
  error.width <- 0.1

  ##############



  data <- data.frame(colnames(table), unlist(table[table$X == comp.name, ], use.names = F))[-1, ]

  names(data) <- data.frame(colnames(table), unlist(table[table$X == comp.name, ], use.names = F))[1, ]
  names(data) <- c("X", "Y")
  data[, 2] <- as.numeric(data[, 2])

  line_wt <- grep(ctrl, data$X)
  data$X[line_wt] <- gsub("\\.[0-9]+", "", data$X[line_wt])

  if (group_lines) {
    line_treatment <- grep(exp, data$X)
    data$X[line_treatment] <- gsub("\\.[0-9]+", "", data$X[line_treatment])

    data <- data[c(line_wt, line_treatment), ]
    n_exp_groups <- 1
  } else {
    # if there is technical triplicates
    # (finished with 'exp.1.1' instead of 'exp.1')
    triplicates_exp <- grepl(paste0("^", exp, "\\.[0-9]+\\.[0-9]+"), data$X)
    if (sum(triplicates_exp) > 1) {
       data$X[triplicates_exp] <- gsub("\\.[0-9]+$", "", data$X[triplicates_exp])
    }
    line_exp <- grep(paste0("^", exp, "\\."), data$X)
    # exp_labels <- seq_along(line_exp)
    # data$X[line_exp] <- paste0(exp, "-", exp_labels)
    data$X[line_exp] <- gsub("\\.", "-", data$X[line_exp])

    n_exp_groups <- length(unique(data$X[line_exp]))
    data <- data[c(line_wt, line_exp), ]
  }

  ### normelize y-values
  if (log_norm) {
    data$Y <- log1p(data$Y)
  }


  title_main <- comp.name
  #  X = names(data)[1]
  #  Y = names(data)[2]
  level.order <- unique(data$X)


  ############### pValue vector
  stat.test <- data %>%
    t_test(Y ~ X, ref.group = ctrl, var.equal = T) %>% ################## problem
    add_significance("p") %>%
    add_xy_position(x = "X")
  stat.test$xmin <- 1
  stat.test$xmax <- stat.test$xmax + 1

  p_labels <- ifelse(p_as_star, "p.signif", "p")
  p_labels_size <- ifelse(p_as_star, 4, 2)


  y_max = max(stat.test$y.position, na.rm = T) * 1.06
  y_min <- min(data$Y, na.rm = T)
  if (log_norm == T & y_max > 8.5) {
    y_min <- y_min * 0.9
  } else {
    y_min <- y_min * 0.8
  }

  # p.value_fun <- function(line, ctrl = data$methionine[line_wt]) {
  #  t.test_p.value = t.test(line,ctrl, var.equal = T)$p.value
  #
  #  if (t.test_p.value <= 0.05 & t.test_p.value > 0.01) {
  #    p_value_lable = "*"
  #  } else if (t.test_p.value <= 0.01 & t.test_p.value > 0.001) {
  #    p_value_lable = "**"
  #  } else if (t.test_p.value <= 0.001) {
  #    p_value_lable = "***"
  #  } else if (t.test_p.value > 0.05) {
  #    p_value_lable = "ns"
  #  }
  #  return(p_value_lable)
  # }

  # pValue vector (arranged)
  # p_value_vector = c("",
  #                   p.value_fun(data$methionine[line_1]),
  #                   p.value_fun(data$methionine[line_2]),
  #                   p.value_fun(data$methionine[line_3]))




  #########################


  if (nchar(comp.name) >= 40) {
    comp.name <- paste(substring(comp.name, 1, 35), "....", substring(comp.name, nchar(comp.name) - 1, nchar(comp.name)), sep = "")
    title_main <- comp.name
  }

  # c("firebrick4","grey70")
  #########################
  comp.plot <- data %>%
    ggplot(aes(x = factor(x = X, level = level.order), y = Y)) +
    geom_boxplot(fill = c("gray100", rep(box_col, n_exp_groups)), colour = "black") +
    theme_classic() +
    theme( # panel.spacing = unit(2, "lines"),
      text = element_text(family = "serif"),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, face = "bold", size = x_cex),
      axis.text.y = element_text(face = "bold", size = 8),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.line = element_blank(),
      axis.ticks = element_line("black", size = 0.75),
      axis.ticks.length.x = unit(0.1, "cm"),
      axis.ticks.length.y = unit(ifelse(group_lines, -0.125, -0.1), "cm")
      # plot.title = element_text(face="bold.italic")
      # strip.text = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      # panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    labs(title = title_main, x = title_x, y = title_y) +

    # geom_jitter(color="#383838", size=0.625, alpha=0.9, width = 0.18)  +
    geom_jitter(color = jitter_col, alpha = 0.9, size = 1) +
    # geom_text(aes(label = ifelse(pValue < 0.05, ifelse(pValue < 0.01, ifelse(pValue < 0.001, '***', '**'), '*'), '')),
    #          vjust = -0.5, size = 4) +
    geom_rect(aes(xmin = 0.5, xmax = length(level.order) + 0.5, ymin = y_min, ymax = y_max), # ymax = max(data[,2])*1.05),
      fill = "transparent", color = "black", size = 0.75
    ) +
    stat_pvalue_manual(stat.test, label = p_labels, label.size = p_labels_size, hide.ns = T, family = "serif", face = "italic") + # , tip.length = 0.01)
    # annotate("text", label= p_value_vector, size = 8)
    scale_y_continuous(expand = c(0,0), limits = c(y_min, y_max))

  dir.create(paste0(path_for_save_file, "/", plot_dir_name), showWarnings = F)
  dir.create(paste0(path_for_save_file, "/", plot_dir_name, "/GCMS_", exp, "_", Sys.Date()), showWarnings = F)

  ggsave(paste0(path_for_save_file, "/", plot_dir_name, "/GCMS_", exp, "_", Sys.Date(), "/", comp.name, ".png"),
    plot = comp.plot, width = width, height = height
  )
}
