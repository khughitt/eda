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
        #' @param bins Number of bins to use for histogram
        #'
        #' @return A ggplot instance.
        plot_libsizes = function(bins=50) {
            ggplot(aes(x=libsize), data=data.frame(libsize=colSums(self$dat))) +
                geom_histogram(bins=bins, fill='#CCCCCC', color='#333333') +
                private$ggplot_theme()
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

