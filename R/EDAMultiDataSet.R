#' An S6 class representing collection of related datasets
#'
#' @section Usage:
#' ```
#' # TODO
#' ```
#'
#' @section Arguments:
#' - `dat`: Primary dataset (matrix, data frame, etc.)
#' - `row_data`: List of zero or more additional datasets which share some
#'     or all row identifiers with dat.
#' - `col_data`: List of zero or more additional datasets which share some
#'     or all column identifiers with dat.
#'
#' @section Fields:
#'  - `dat`: Primary dataset
#'  - `row_data`: List of additional data keyed on row identifiers
#'  - `col_data`: List of additional data keyed on column identifiers
#'
#' @section Methods:
#'  - `cross_cor(key1=1, key2=2, method='pearson')`: Computes cross-dataset
#'     correlation matrix between rows in two specified datasets.
#'  - `plot_cross_cor_heatmap(key1=1, key2=2, method='pearson', interactive=TRUE)`:
#'      Plots multidataset correlation heatmap.
#'  - `print()`: Prints an overview of the object instance.
#'
#' @section Examples:
#' ```
#' TODO
#' ```
#'
#' @importFrom R6 R6Class
#' @name EDAMultiDataSet
#' @export
#'
NULL

EDAMultiDataSet <- R6Class("EDAMultiDataSet",
    inherit = eda:::AbstractMultiDataSet,

    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # EDADataSet constructor
        initialize = function(dat, row_data=list(), col_data=list()) {
            super$initialize(dat, row_data = row_data, col_data = col_data)
        },

        # Computes cross-dataset correlation matrix
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        cross_cor = function(key1=1, key2=2, method='pearson') {
            super$cross_cor(key1, key2, method)
        },

        # Plots multidataset correlation heatmap
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        plot_cross_cor_heatmap = function(key1=1, key2=2, method='pearson', interactive=TRUE) {
            super$plot_cross_cor_heatmap(key1, key2, method, interactive)
        },

        # Prints an overview of the object instance
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= AbstractMultiDataSet (n=%d)\n", length(private$datasets)))
            cat("=\n")
            cat(sprintf("= dat: %s (%d x %d)\n", class(self$dat)[1], nrow(self$dat), ncol(self$dat)))
            if (length(private$row_data) > 1) {
                cat("=\n")
                cat("= Row data\n")
                cat("=\n")
                for (i in 2:length(private$row_data)) {
                    ds <- private$row_data[[i]]
                    cat(sprintf("= %02d. %s (%d x %d)\n", i, class(ds)[1], nrow(ds$dat), ncol(ds$dat)))
                }
            }
            if (length(private$col_data) > 1) {
                cat("=\n")
                cat("= Column data\n")
                cat("=\n")
                for (i in 2:length(private$col_data)) {
                    ds <- private$col_data[[i]]
                    cat(sprintf("= %02d. %s (%d x %d)\n", i, class(ds)[1], nrow(ds$dat), ncol(ds$dat)))
                }
            }
            cat("=\n")
            cat("=========================================\n")
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(),

    # ------------------------------------------------------------------------
    # active
    # ------------------------------------------------------------------------
    active = list(
        # Make datasets publically visible for EDAMultiDataSet instances
        row_data = function(value) {
            if (missing(value)) {
                private$row_data
            } else {
                private$row_data <- value
            }
        },
        col_data = function(value) {
            if (missing(value)) {
                private$col_data
            } else {
                private$col_data <- value
            }
        }
    )
)
