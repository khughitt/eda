#' An R6 class representing a high-throughput expression dataset along with
#' associated gene and sample metadata.
#'
#' @section Usage:
#' ```
#' # TODO
#' ```
#'
#' @section Arguments:
#' - `datasets`: A list of datasets (matrices, data frames, etc.), each of
#'      which shared some column / row identifiers with the first entry in
#'      the list.
#'
#' @section Fields:
#'  - `datasets`: List of datasets
#'
#' @section Methods:
#'  - `cross_cor(key1=1, key2=2, meas='pearson')`: Computes cross-dataset
#'     correlation matrix between rows in two specified datasets.
#'  - `plot_cross_cor_heatmap(key1=1, key2=2, meas='pearson', interactive=TRUE)`:
#'      Plots multidataset correlation heatmap.
#'  - `print()`: Prints an overview of the object instance.
#'
#' @section Examples:
#' ```
#' TODO
#' ```
#'
#' @importFrom R6 R6Class
#' @name BioExprSet
#' @export
#'
NULL

BioExprSet <- R6Class("BioExprSet",
    inherit = eda:::BioDataSet,

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
        # @param meas Correlation measure to use
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        cross_cor = function(key1=1, key2=2, meas='pearson', new_key=NULL, ...) {
            private$compute_cross_cor(key1, key2, meas, new_key, ...)
        },

        # Performs a counts-per-million (CPM) transformation.
        #
        # @return A CPM-transformed version of the BioExprSet instance.
        cpm = function(key=1) {
            obj <- private$clone_() 
            dat <- obj$edat[[key]]$dat 
            obj$update(key, sweep(dat, 2, colSums(dat), '/') * 1E6)
            obj
        },

        # Plots multidataset correlation heatmap
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param meas Correlation measure to use
        #
        plot_cross_cor_heatmap = function(key1=1, key2=2, meas='pearson', interactive=TRUE) {
            super$plot_cross_cor_heatmap(key1, key2, meas, interactive)
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
        datasets = function() {
            lapply(self$edat, function(x) { x$dat })
        }
    )
)
