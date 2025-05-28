yo_theme_base = function (base_size = 12, base_family = "", base_line_size = base_size/22, 
    base_rect_size = base_size/22) {
    half_line <- base_size/2
    theme_bw(
        base_size = base_size, base_family = base_family,
        base_line_size = base_line_size, base_rect_size = base_rect_size
    ) %+replace%
        theme(
            axis.text = element_text(colour = "black", size = rel(0.8)),
            axis.ticks = element_line(colour = "black", linewidth = rel(0.5)),
            panel.border = element_rect(
                fill = NA, colour = "black",
                linewidth = rel(1)
            ),
            panel.grid = element_blank(),
            #panel.grid.major = element_line(linewidth = rel(0.1)),
            #panel.grid.minor = element_line(linewidth = rel(0.05)),
            #strip.background = element_rect(fill = "black"),
            strip.text = element_text(
                colour = "white", size = rel(0.8),
                margin = margin(
                    0.8 * half_line, 0.8 * half_line,
                    0.8 * half_line, 0.8 * half_line
                )
            ),
            complete = TRUE
        )
    }

#yo_theme_base <- function(base_size = 16, base_family = "") {
    #library(ggthemes)
#
    # yo_theme_foundation <- function(base_size = 12, base_family = "") {
    #    thm <- theme_gray(base_size = base_size, base_family = base_family)
    #    thm$rect$colour <- "red"
    #    #thm$rect$linetype = 0
    #    for (i in names(thm)) {
    #        if ("colour" %in% names(thm[[i]])) {
    #            thm[[i]]["colour"] <- list(NULL)
    #        }
    #        if ("fill" %in% names(thm[[i]])) {
    #            thm[[i]]["fill"] <- list(NULL)
    #        }
    #    }
    #    thm + theme(
    #        panel.border = element_rect(fill = NA), legend.background = element_rect(colour = "blue"),
    #        line = element_line(colour = "black"),
    #        rect = element_rect(
    #            fill = "white",
    #            colour = "black"
    #        ),
    #        text = element_text(colour = "black")
    #    )
    # }
    #
    # yo_theme_foundation() +
    #theme(
    #    line = element_line(
    #        colour = "black",
    #        lineend = "round", linetype = "solid"
    #    ), rect = element_rect(
    #        fill = "white",
    #        colour = "black", linetype = "solid"
    #    ), text = element_text(
    #        colour = "black",
    #        face = "plain", family = base_family, size = base_size,
    #        vjust = 0.5, hjust = 0.5, lineheight = 1
    #    ), panel.grid = element_blank(),
    #    strip.background = element_rect(colour = NA), legend.key = element_rect(colour = NA),
    #    title = element_text(size = rel(1)), plot.title = element_text(
    #        size = rel(1.2),
    #        face = "bold"
    #    ), strip.text = element_text(), axis.ticks.length = unit(
    #        0.5,
    #        "lines"
    #    )
    #)

