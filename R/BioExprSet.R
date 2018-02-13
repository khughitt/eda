#' An S6 class representing an Expression dataset
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

        # class greeting
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

        # counts-per-million
        cpm = function() {
            obj <- self$clone()
            obj$dat <- sweep(obj$dat, 2, colSums(obj$dat), '/') * 1E6
            obj
        },

        log2p = function() {
            self$log(2, offset=1)
        },

        plot_libsizes = function(bins=50) {
            ggplot(aes(x=libsize), data=data.frame(libsize=colSums(self$dat))) +
                geom_histogram(bins=bins, fill='#CCCCCC', color='#333333') +
                private$ggplot_theme()
        }
    ),
    private = list(
        # verify that input data is of type matrix
        check_input = function(dat) {
            if(!is.matrix(dat) && (class(dat) != 'ExpressionSet')) {
                stop("Invalid input for BioExprSet: dat must be a matrix or ExpressionSet.")
            }
        }
    )
)

