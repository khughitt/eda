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
            # Unpack any ExpressionSet inputs
            num_esets <- 0

            for (i in seq_along(datasets)) {
                if (class(datasets[[i]])[1] == 'ExpressionSet') {
                    num_esets <- num_esets + 1

					# get eset
                    dat <- datasets[[i]]
                    datasets[[i]] <- NULL

					# in cases where multiple ExpressionSets are provided, make
					# sure keys don't overlap
                    key <- c('exprs', 'pdata', 'fdata')

					phenotype_xid <- 'samples'
					phenotype_yid <- 'phenotypes'
					feature_xid   <- 'features'
					feature_yid   <- 'feature_annotations'

                    if (num_esets > 1) {
                        keys          <- sprintf('%s_%d', keys, num_esets)
						phenotype_xid <- sprintf('%s_%d', phenotype_xid, num_esets)
						phenotype_yid <- sprintf('%s_%d', phenotype_yid, num_esets)
						feature_xid   <- sprintf('%s_%d', feature_xid, num_esets)
						feature_yid   <- sprintf('%s_%d', feature_yid, num_esets)
                    }

					# add expression data
                    data[[keys[1]]] = EDADat$new(exprs(dat), xid = feature_xid, yid = phenotype_xid)

					# add phenotype data 
                    if (!is.null(pData(dat))) {
                        datasets[[keys[2]]] = EDADat$new(pData(dat), xid = phenotype_xid, yid = phenotype_yid)
                    }
					# add feature data 
                    if (!is.null(fData(dat))) {
                        datasets[[keys[3]]] = EDADat$new(fData(dat), xid = feature_xid, yid = feature_yid)
					}
                }
            }
            super$initialize(datasets, color_pal, title, ggplot_theme)
        },

        # Computes cross-dataset correlation matrix
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        cross_cor = function(key1=1, key2=2, method='pearson', new_key=NULL, ...) {
            private$compute_cross_cor(key1, key2, method, new_key, ...)
        },

        # Performs a counts-per-million (CPM) transformation.
        #
        # @return A CPM-transformed version of the BioExprSet instance.
        cpm = function(key=1) {
            obj <- private$clone_() 
            dat <- obj$edat[[key]]$dat 
            obj$edat[[key]]$dat <- sweep(dat, 2, colSums(dat), '/') * 1E6
            obj
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
        datasets = function() {
			lapply(self$edat, function(x) { x$dat })
        }
    )
)
