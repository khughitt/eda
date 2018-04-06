#' An R6 class representing collection of related datasets with a matrix
#' for the primary dataset.
#'
#' @section Usage:
#' ```
#' # TODO
#' ```
#'
#' @section Arguments:
#' - `dataset`: A list of datasets (matrices, data frames, etc.), each of
#'      which shared some column / row identifiers with the first entry in
#'      the list.
#'
#' @section Fields:
#'  - `datasets`: List of datasets
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
        # EDAMultiMatrix constructor
        initialize = function(datasets, color_pal='Set1', title="", ggplot_theme=theme_bw) {
            super$initialize(datasets, color_pal, title, ggplot_theme)
        },

        # Computes cross-dataset correlation matrix
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        cross_cor = function(key1=1, key2=2, method='pearson') {
            super$compute_cross_cor(key1, key2, method)
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
        datasets = function(value) {
            if (missing(value)) {
                lapply(self$edat, function(x) { x$dat })
            } else {
                # TO TEST (may need to make datasets read-only...)
                lapply(self$edat, function(x) { x$dat }) <- value
            }
        }
    )
)
