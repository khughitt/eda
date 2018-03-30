#' An S6 class representing collection of related datasets with a matrix
#' for the primary dataset.
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
#' @name EDAMultiMatrix
#' @export
#'
NULL

EDAMultiMatrix <- R6Class("EDAMultiMatrix",
    inherit = eda:::AbstractMatrixDataSet,

    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # EDADataSet constructor
        initialize = function(dat, row_data=list(), col_data=list(),
                              row_color=NULL, row_color_ds='dat',
                              row_shape=NULL, row_shape_ds='dat',
                              row_label=NULL, row_label_ds='dat',
                              col_color=NULL, col_color_ds='dat',
                              col_shape=NULL, col_shape_ds='dat',
                              col_label=NULL, col_label_ds='dat',
                              color_pal='Set1', title="", ggplot_theme=theme_bw) {
            super$initialize(dat, row_data = row_data, col_data = col_data,
                             row_color, row_color_ds, row_shape, row_shape_ds,
                             row_label, row_label_ds, col_color, col_color_ds,
                             col_shape, col_shape_ds, col_label, col_label_ds,
                             color_pal, title, ggplot_theme)
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
        # Make datasets publically visible for EDAMultiMatrix instances
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
