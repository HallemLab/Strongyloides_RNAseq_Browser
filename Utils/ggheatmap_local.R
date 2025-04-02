#' ggplot heatmap equivalent to heatmaply
#' 
#' This is a local version separated from the heatmaply package
#' Separation by Astra S. Bryant, PhD for the purpose of greater control over 
#' ggplot graphical output. Specifically, this local version differs from the
#' main pacakge version of ggheatmap in two major ways:
#' 1. New input showticklabels passed into the function by the user that
#' determines whether x and y tick labels are included in the plot. 
#' 2. Users can specify a main title for heatmap plot - done by 
#'    passing the main variable from the main ggheatmap function call through
#'    to the ggheatmap call.
#' 3. Margins are also hard coded around the plot to provide some addional spacing.
#' 
#' Note: the vast majority of the code is not difference from the main package
#' version - references to //TODO items reflect goals of the main package
#' authors.
#' 
#' This function produces a ggplot analogue of heatmaply figures 
#' using \link[egg]{ggarrange}. This function may not always support the same
#' set of features as , and exporting the heatmaply object with, for example,
#' \link[plotly]{orca} or \code{heatmaply(mtcars, file = "foo.png")}.
#' 
#' @param ... Passed to \link{heatmaply}
#' @param widths,heights Relative widths and heights of plots.
#' @param row_dend_left Logical argument controlling whether the row 
#'  dendrogram is placed on the left of the plot.
#' @examples
#' ggheatmap(mtcars)
#' @export
ggheatmap_local <- function(...,  showticklabels = c(TRUE, TRUE), widths = NULL, heights = NULL, row_dend_left = FALSE, main = NULL) {
    plots <- heatmaply(
        ..., 
        row_dend_left = row_dend_left, 
        return_ppxpy = TRUE,
        plot_method = "ggplot"
    )
    arrange_plots(
        plots,
        widths = widths,
        heights = heights,
        row_dend_left = row_dend_left,
        showticklabels = showticklabels,
        main = main
    )
    
}

default_dims <- function(px, pr) {
    if (!is.null(px)) {
        if (is.null(pr)) {
            widths <- c(0.8, 0.2)
        } else {
            widths <- c(0.7, 0.1, 0.2)
        }
    } else {
        if (is.null(pr)) {
            widths <- 1
        } else {
            widths <- c(0.9, 0.1)
        }
    }
    widths
}


## TODO: duplication with heatmap_subplot_from_ggplotly
arrange_plots <- function(
    plots, 
    widths = NULL, 
    heights = NULL, 
    row_dend_left = FALSE,
    showticklabels = showticklabels,
    main = NULL) {
   
    plots <- plots[!sapply(plots, is.null)]
    if (!row_dend_left) {
        plots$p <- plots$p + theme(legend.position = "left")
    }
    plots <- lapply(plots, function(x) {x + theme(plot.margin = unit(c(0, 0, 0, 0), "npc"))})
   
     plots$p <- plots$p + theme(axis.text.y = element_blank(), 
                                axis.ticks.y = element_blank())
     
     plots$py <- plots$py + theme(plot.margin = unit(c(0.5,0,0,0),"cm"))
     plots$px <- plots$px + theme(plot.margin = unit(c(0,0.5,0,0),"cm"))
    
    if (showticklabels[1] == FALSE) {
        plots <- lapply(plots, function(x) {x + theme(axis.text.x = element_blank())}
        )}
    
    if (showticklabels[2] == FALSE) {
    plots <- lapply(plots, function(x){ x + theme(axis.text.y = element_blank())}
    )}
    
    column_list <- list(plots$py, plots$pc, plots$p)
    ind_null_col <- sapply(column_list, is.null)
    
    row1_list <- list(plots$py, ggplot_empty(), ggplot_empty())
    row2_list <- list(plots$pc, ggplot_empty(), ggplot_empty())
    row3_list <- list(plots$p, plots$pr, plots$px)
    
    if (row_dend_left) {
        row3_list <- rev(row3_list)
        row2_list <- rev(row2_list)
        row1_list <- rev(row1_list)
    }
    plotlist <- c(
        row1_list,
        row2_list,
        row3_list
    )
    
    nrows <- sum(!ind_null_col)
    ind_remove_col <- rep(ind_null_col, each = length(plotlist) / 3)
    
    ind_null_row <- sapply(row3_list, is.null)
    ncols <- sum(!ind_null_row)
    ind_remove_row <- rep(ind_null_row, length.out = length(plotlist))
    plotlist <- plotlist[!(ind_remove_row | ind_remove_col)]
    
    egg::ggarrange(
        plots = plotlist,
        ncol = ncols,
        top = main,
        widths = widths %||% default_dims(plots$px, plots$pr),
        heights = heights %||% rev(default_dims(plots$py, plots$pc))
    )
}

ggplot_empty <- function() {
    ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "npc"))
}