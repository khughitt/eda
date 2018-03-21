#' An S6 class representing collection of related datasets
#'
#' @importFrom R6 R6Class
#' @name EDAMultiDataSet
#' @export
#'
NULL

EDAMultiDataSet <- R6Class("EDAMultiDataSet",
    inherit = AbstractMultiDataSet,

    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # EDADataSet constructor
        initialize = function(...) {
            super$initialize(...)
            private$check_inputs()
        },

        #' Computes cross-dataset correlation matrix
        #'
        #' @param key1 Numeric or character index of first dataset to use
        #' @param key2 Numeric or character index of second dataset to use
        #' @param method Correlation method to use (passed to `cor` function)
        #'
        #' @return Matrix of pairwise dataset1 - dataset2 correlations
        cross_cor = function(key1=1, key2=2, method='pearson') {
            super$cross_cor(key1, key2, method)
        },

        #' Plots multidataset correlation heatmap
        #'
        #' @param key1 Numeric or character index of first dataset to use
        #' @param key2 Numeric or character index of second dataset to use
        #' @param method Correlation method to use (passed to `cor` function)
        #'
        plot_cross_cor_heatmap = function(key1=1, key2=2, method='pearson', interactive=TRUE) {
            super$plot_cross_cor_heatmap(key1, key2, method, interactive)
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        #' Make sure input datasets are all EDADataSet instances and include
        #' the same column identifiers
        check_inputs = function() {
            # make sure inputs are all EDADataSets, and are indexed by the same columns
            ids <- sort(colnames(self$datasets[[1]]$dat))

            for (ds in self$datasets) {
                if (!'EDADataSet' %in% class(ds)) {
                    stop("Invalid input: each input dataset much be an EDADataSet instance.")
                }

                # TODO: Allow datasets with partially overlapping columns..
                if (!all(sort(colnames(ds$dat)) == ids)) {
                    stop("Input datasets must all include the same columns identifiers.")
                }
            }

        }
    ),

    # ------------------------------------------------------------------------
    # active
    # ------------------------------------------------------------------------
    active = list(
        #' Make datasets publically visible for EDAMultiDataSet instances
        datasets = function() {
            private$datasets
        }       
    )
)
