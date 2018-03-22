#' An S6 class representing a generic collection of datasets
#'
#' @importFrom R6 R6Class
#' @name AbstractMultiDataSet
#'
NULL

AbstractMultiDataSet <- R6Class("AbstractMultiDataSet",
    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # AbstractMultiDataSet constructor
        initialize = function(...) {
            private$datasets <- list(...)
        }
    ),
    private = list(
        # private params
        datasets = NULL,

        # Computes cross-dataset correlation matrix
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        cross_cor = function(key1=1, key2=2, method='pearson') {
            # make sure datasets are ordered similarly
            dat1 <- private$datasets[[key1]]
            dat2 <- private$datasets[[key2]]

            # for multidatasets, get underlying data
            if ('EDADataSet' %in% class(dat1)) {
                dat1 <- dat1$dat
                dat2 <- dat2$dat
            }

            # for single dataset metadata (feature) correlations, column 
            # metadata need to be transposed to be in the proper order
            if (key2 == 'col_mdata') {
                # transpose metadata
                dat2 <- t(dat2)

                # for metadata, we can also exclude factor fields with all
                # unique values (e.g. alternate identifiers)
                exclude <- apply(drug_mdata, 1, function(x) { is.factor(x) && length(unique(x)) == length(x) })

                if (sum(exclude) > 0) {
                    message(sprintf("Excluding %d unique factor fields", sum(exclude)))
                    dat2 <- dat2[!exclude,]
                }
            }

            # only compare matching columns
            col_ind <- intersect(colnames(dat1), colnames(dat2))

            dat1 <- dat1[,col_ind, drop=FALSE]
            dat2 <- dat2[,col_ind, drop=FALSE]

            # drop any rows with zero variance
            mask <- apply(dat1, 1, function(x) { length(table(x))} ) > 1
            if (sum(!mask) > 0) {
                message(sprintf("Excluding %d zero-variance rows from first dataset.", sum(!mask)))
                dat1 <- dat1[mask,, drop=FALSE]
            }

            mask <- apply(dat2, 1, function(x) { length(table(x))} ) > 1
            if (sum(!mask) > 0) {
                message(sprintf("Excluding %d zero-variance rows from second dataset.", sum(!mask)))
                dat2 <- dat2[mask,, drop=FALSE]
            }


            # linear model
            #
            # Based on code adapted from cbcbSEQ
            # (https://github.com/kokrah/cbcbSEQ/) originally written by
            # Kwame Okrah.
            if (method == 'lm') {
                # construct linear model to measure dependence of each projected
                # axis on column metadata
                cor_mat <- matrix(0, nrow=nrow(dat1), ncol=nrow(dat2))

                for (i in 1:nrow(dat2)) {
                    feature_cor <- function(y) {
                        round(summary(lm(y~dat2[i,]))$r.squared*100, 2)
                    }
                    cor_mat[,i] <- apply(dat1, 1, feature_cor)
                }
                colnames(cor_mat) <- rownames(dat2)
            } else {
                # Pearson correlation, etc.
                cor_mat <- cor(t(rbind(dat1, dat2)), method=method)

                # limit to cross-dataset correlations
                cor_mat <- cor_mat[1:nrow(dat1), (nrow(dat1) + 1):ncol(cor_mat)]
            }
            cor_mat
        },

        # Plots multidataset correlation heatmap
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        plot_cross_cor_heatmap = function(key1=1, key2=2, method='pearson', interactive=TRUE) {
            # compute cross correlations
            cor_mat <- self$cross_cor(key1, key2, method)

            # list of parameters to pass to heatmaply
            params <- list(
                x=cor_mat,
                showticklabels=c(FALSE, FALSE),
                subplot_widths=c(0.65, 0.35),
                subplot_heights=c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$datasets[[key1]]$row_mdata)) {
                mdata1 <- self$datasets[[key1]]$row_mdata
                mask1  <- sapply(mdata1, function(x) { max(table(x)) > 1 })
                mdata1 <- mdata1[,mask1, drop=FALSE]

                params[['row_side_colors']] <- mdata1 
                params[['subplot_widths']] <- c(0.15, 0.3, 0.55)
            }
            
            if (!is.null(self$datasets[[key2]]$row_mdata)) {
                mdata2 <- self$datasets[[key2]]$row_mdata
                mask2  <- sapply(mdata2, function(x) { max(table(x)) > 1 })
                mdata2 <- mdata2[,mask2, drop=FALSE]

                params[['col_side_colors']] <- mdata2 
                params[['subplot_heights']] <- c(0.55, 0.3, 0.15)
            }

            # add any additional function arguments
            private$plot_heatmap(params, interactive)
        },

        # Creates a static or interactive heatmap plot
        #
        # @param params A list of plotting parameters
        # @param interactive Logical indicating whether an interactive heatmap
        #     should be generated.
        plot_heatmap = function(params, interactive) {
            # interactive heatmap
            if (interactive) {
                return(do.call(heatmaply::heatmaply, params))
            }

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
