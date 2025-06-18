  ####
# This multiplot function script taken from the R Cookbook website:
# http://www.cookbook-r.com/
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
# modifued by Yo - added heiths and widths parameters
 ####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# an example for 3 rows and 6 columns plots with different heights and widths:
# layout = matrix(1:18, nrow = 3, ncol = 6, byrow = TRUE)
# heights = c(2, 1, 1)  # First row twice as tall
# widths = c(1, 1, 2, 1, 1, 1)  # Third column twice as wide

multiplot <- function (..., plotlist = NULL, file, cols = 1, layout = NULL, heights = NULL, widths = NULL) 
{
    library(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
            ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    # Set default heights if not provided
    if (is.null(heights)) {
        heights <- rep(1, nrow(layout))  # Equal heights by default
    }
    
    # Set default widths if not provided
    if (is.null(widths)) {
        widths <- rep(1, ncol(layout))  # Equal widths by default
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
    }
    else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), 
            ncol(layout), heights = unit(heights, "null"), widths = unit(widths, "null"))))
        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                layout.pos.col = matchidx$col))
        }
    }
}