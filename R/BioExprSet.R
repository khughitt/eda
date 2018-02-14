#' An S6 class representing an Expression dataset
#'
#' BioExprSet is a simple class for interacting with biological expression
#' data (e.g. from microarray or RNA-Seq experiments.) It can accept either
#' a matrix and optional column/row metadata dataframes, or a Bioconductor
#' ExpressionSet instance. Methods are provided for common transformations
#' and visualizations.
#'
#' @section Arguments:
#' \describe{
#'   \item{dat}{matrix|ExpressionSet An expression dataset with rows corresponding
#'       to genes, transcripts, probes, etc. and column corresponding to 
#'       individual samples.}
#'   \item{row_mdata}{Data frame with rows corresponding to the column names of 
#'       \code{dat}.}
#'   \item{row_mdata}{Data frame with rows corresponding to the row names of 
#'       \code{dat}.}
#' }
#'
#' @importFrom R6 R6Class
#' @name BioExprSet
#' @export
#'
NULL

BioExprSet <- R6::R6Class("BioExprSet",
    inherit = EDAMatrix,
    public = list(
        # BioExprSet constructor
        initialize = function(dat, col_mdata=NULL, row_mdata=NULL,
                              col_maxn=Inf, col_maxr=1.0,
                              row_maxn=Inf, row_maxr=1.0,
                              color_var=NULL, shape_var=NULL, label_var=NULL,
                              color_pal='Set1', ggplot_theme=ggplot2::theme_bw) { 
            # verify input data type and call parent constructor
            private$check_input(dat)

            # for ExpressionSets, retrieve relevant data and metadata
            if (class(dat) == 'ExpressionSet') {
                if (is.null(col_mdata)) {
                    col_mdata <- pData(dat)
                }
                if (is.null(row_mdata)) {
                    row_mdata <- fData(dat)
                }
                dat <- exprs(dat)
            }

            super$initialize(dat, col_mdata, row_mdata, col_maxn, col_maxr,
                             row_maxn, row_maxr, color_var, shape_var, label_var,
                             color_pal, ggplot_theme)
        },

        #' Prints an overview of the object instance
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= BioExprSet (%s)\n", class(self$dat[,1])))
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        },

        #' Performs a counts-per-million (CPM) transformation.
        #'
        #' @return A CPM-transformed version of the BioExprSet instance.
        cpm = function() {
            obj <- self$clone()
            obj$dat <- sweep(obj$dat, 2, colSums(obj$dat), '/') * 1E6
            obj
        },

        #' Log2 transforms data (adding 1 to ensure finite results).
        #'
        #' @return Log2-transformed version of the expression data.
        log2p = function() {
            self$log(2, offset=1)
        },

        #' Plots a histogram of sample library sizes.
        #'
        #' Generates a histogram of sample library sizes (sum of expression
        #' levels within each column/sample). 
        #'
        #' @param bins Number of bins to use for histogram
        #'
        #' @return A ggplot instance.
        plot_libsize_hist = function(...) {
            libsizes <- data.frame(libsize=colSums(self$dat))
            ggplot(aes(x=libsize), data=libsizes) +
                geom_histogram(..., fill='#CCCCCC', color='#333333') +
                private$ggplot_theme()
        },

        #' Plots bar graph of sample library sizes
        #'
        #' Generates a bargraph plot of sample library sizes (sum of expression
        #' levels within each column/sample). Each sample is shown as a separate
        #' bar in the plot.
        #'
        #' @return A ggplot instance.
        plot_libsize_bargraph = function(color_var=NULL) {
            # data frame with sample library sizes
            libsizes <- data.frame(libsize=colSums(self$dat))

            libsizes$color_var <- private$get_plot_color_column(color_var)

            plot_aes  <- private$get_plot_aes(color_var)
            plot_labs <- private$get_plot_legend_labels(color_var)

            ggplot(aes(x=rownames(libsizes), y=libsize), data=libsizes) +
                geom_bar(aes(fill=color_var), stat='identity') + 
                   plot_labs +
                   xlab("Samples") +
                   ylab("Total expression") +
                   private$ggplot_theme() +
                   theme(axis.text.x=element_text(angle=90),
                         legend.text=element_text(size=8))
        }
    ),
    private = list(
        #' Verifies that input data is an acceptable format.
        check_input = function(dat) {
            if(!is.matrix(dat) && (class(dat) != 'ExpressionSet')) {
                stop("Invalid input for BioExprSet: dat must be a matrix or ExpressionSet.")
            }
        }
    )
)

