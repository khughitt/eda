#' An S6 class representing collection of related datasets
#'
#' @importFrom R6 R6Class
#' @name EDADataSet
#' @export
#'
NULL

EDAMultiDataSet <- R6Class("EDAMultiDataSet",
    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # params
        datasets = NULL,

        # EDADataset constructor
        initialize = function(...) {
            self$datasets <- list(...)
            private$check_inputs()
        },

        #' Plots multidataset correlation heatmap
        #'
        #' @param key1 Numeric or character index of first dataset to use
        #' @param key2 Numeric or character index of second dataset to use
        #' @param method Correlation method to use (passed to `cor` function)
        #'
        plot_cor_heatmap = function(key1=1, key2=2, method='pearson', interactive=TRUE) {
            # make sure datasets are ordered similarly
            # TODO: include row and col indices
            col_ind <- colnames(self$datasets[[key1]]$dat)
            cor_mat <- cor(t(rbind(self$datasets[[key1]]$dat, 
                                   self$datasets[[key2]]$dat[,col_ind])), 
                           method=method)

            # limit to cross-dataset correlations
            cor_mat <- cor_mat[1:nrow(self$datasets[[key1]]$dat), 
                               (nrow(self$datasets[[key1]]$dat) + 1):ncol(cor_mat)]

            # determine subsampling indices, if requested
            #indices <- private$get_indices(...)

            # list of parameters to pass to heatmaply
            params <- list(
                x=cor_mat,
                showticklabels=c(FALSE, FALSE),
                subplot_widths=c(0.65, 0.35),
                subplot_heights=c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$datasets[[key1]]$row_mdata)) {
                #params[['row_side_colors']] <- self$row_mdata[indices$row, 
                #                                              binary_vars, drop=FALSE]
                params[['row_side_colors']] <- self$datasets[[key1]]$row_mdata
                params[['subplot_widths']] <- c(0.15, 0.3, 0.55)
            }
            
            if (!is.null(self$datasets[[key2]]$row_mdata)) {
                #params[['col_side_colors']] <- self$col_mdata[indices$col, 
                #                                                !binary_vars, drop=FALSE]
                params[['col_side_colors']] <- self$datasets[[key2]]$row_mdata
                params[['subplot_heights']] <- c(0.55, 0.3, 0.15)
            }

            # add any additional function arguments
            #params <- c(params, private$get_custom_function_args(...))
            private$plot_heatmap(params, interactive)
        }
    ),
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

                if (!all(sort(colnames(ds$dat)) == ids)) {
                    stop("Input datasets must all include the same columns identifiers.")
                }
            }

        },

        #' Creates a static or interactive heatmap plot
        #'
        #' @param params A list of plotting parameters
        #' @param interactive Logical indicating whether an interactive heatmap
        #'     should be generated.
        plot_heatmap = function(params, interactive) {
            # interactive heatmap
            if (interactive) {
                return(do.call(heatmaply::heatmaply, params))
            }

            # static heatmap

            # convert colside and rowside colors to explicit color values,
            # if present
            if ('col_side_colors' %in% names(params)) {
                params$ColSideColors <- params$col_side_colors

                colors <- c('blue', 'yellow')

                for (col_name in colnames(params$ColSideColors)) {
                    col <- params$ColSideColors[,col_name]
                    params$ColSideColors[,col_name] <- colors[as.numeric(factor(col))]
                }
                params$ColSideColors <- as.matrix(params$ColSideColors)
            }

            if ('row_side_colors' %in% names(params)) {
                params$RowSideColors <- params$row_side_colors

                pal <- RColorBrewer::brewer.pal(8, 'Set1')

                for (col_name in colnames(params$RowSideColors)) {
                    col <- params$RowSideColors[,col_name]
                    colors <- colorRampPalette(pal)(min(1E4, length(unique(col))))
                    params$RowSideColors[,col_name] <- colors[as.numeric(factor(col))]
                }
                params$RowSideColors <- as.matrix(params$RowSideColors)
            }

            # remove irrelevant function arguments
            heatmaply_args <- c('showticklabels', 'subplot_widths', 'subplot_heights',
                                'col_side_colors', 'row_side_colors')
            params <- params[!names(params) %in% heatmaply_args]

            do.call(heatmap.plus::heatmap.plus, params)
        }
    )
)
