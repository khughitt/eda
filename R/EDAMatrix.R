#' An S6 class representing a matrix dataset.
#'
#' EDAMatrix is a helper class for wrapping data matrices, with optional
#' support for row and column datadata. Methods are provided for common
#' exploratory data analysis summary statistics, transformations, and 
#' visualizations.
#'
#' @section Arguments:
#' \describe{
#'   \item{dat}{An m x n dataset.}
#'   \item{row_mdata}{Data frame with rows corresponding to the column names of 
#'       \code{dat}.}
#'   \item{row_mdata}{Data frame with rows corresponding to the row names of 
#'       \code{dat}.}
#' }
#'
#' @importFrom R6 R6Class
#' @name EDAMatrix
#' @export
#'
NULL

EDAMatrix <- R6::R6Class("EDAMatrix",
    inherit = EDADataSet,
    public = list(
        #' EDAMatrix constructor
        initialize = function(dat, col_mdata=NULL, row_mdata=NULL, title='',
                              col_maxn=Inf, col_maxr=1.0,
                              row_maxn=Inf, row_maxr=1.0,
                              color_var=NULL, shape_var=NULL, label_var=NULL,
                              color_pal='Set1', ggplot_theme=ggplot2::theme_bw) { 
            # verify input data type and call parent constructor
            private$check_input(dat)

            super$initialize(dat, col_mdata, row_mdata, title, col_maxn, col_maxr,
                             row_maxn, row_maxr, color_var, shape_var, label_var,
                             color_pal, ggplot_theme)
        },

        #' Prints an overview of the object instance
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= EDAMatrix (%s)\n", class(self$dat[,1])))
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        }
    ),
    private = list(
        #' Verifies that input data is of type matrix
        check_input = function(dat) {
            if(!is.matrix(dat)) {
                stop("Invalid input for EDAMatrix: dat must be a matrix.")
            }
        }
    )
)

