#' An S6 class representing a generic collection of datasets linked by
#' either column or row identifiers.
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
        initialize = function(dat, row_data=list(), col_data=list()) {
            private$row_data <- c(list('dat' = dat), row_data)
            private$col_data <- c(list('dat' = dat), col_data)
        }
    ),
    private = list(
        # private params
        row_data = NULL,
        col_data = NULL,

        #
        # Supported similarity measures
        #
        # Each of the below functions should accept either one or two matrices
        # or data frames, measures the similarity/dependence either within
        # (single dataset) or across (two datasets) columns in the dataset(s).
        #
        # The result is a matrix of similarity scores (e.g. correlation scores,
        # r^2, etc.) of dimensions:
        #
        # ncol(dat1) x ncol(dat1)  - single dataset
        # ncol(dat1) x nrow(dat2)  - two datasets
        #
        similarity_measures = list(
            #
            # linear model fit
            #
            'lm' = function(dat1, dat2=NULL, ...) {
                # single dataset
                if (is.null(dat2)) {
                    # construct linear model to measure pairwise dependence of
                    # dataset1 columns on dataset2 columns
                    cor_mat <- matrix(0, nrow = ncol(dat1), ncol = ncol(dat1))

                    # for each column in dataset1
                    for (i in 1:ncol(dat1)) {
                        # fit a linear model between the ith column of dataset1
                        # and the jth column of dataset1
                        feature_cor <- function(y) {
                            round(summary(lm(y ~ dat1[, i]))$r.squared * 100, 2)
                        }
                        # for each column in dataset1
                        cor_mat[, i] <- apply(dat1, 2, feature_cor)
                    }
                } else {
                    # two datasets

                    # construct linear model to measure pairwise dependence of
                    # dataset1 columns on dataset2 columns
                    cor_mat <- matrix(0, nrow = ncol(dat1), ncol = ncol(dat2))

                    # for each column in dataset2
                    for (i in 1:ncol(dat2)) {
                        # fit a linear model between the ith column of dat2 and the
                        # jth column of dat1
                        feature_cor <- function(y) {
                            round(summary(lm(y ~ dat2[, i]))$r.squared * 100, 2)
                        }
                        # for each column in dataset1
                        cor_mat[, i] <- apply(dat1, 2, feature_cor)
                    }
                }
                cor_mat
            },
            #
            # Mutual Information
            #
            'mi' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    # two datasets
                    cor_mat <- mpmi::cmi(cbind(dat1, dat2), ...)$bcmi
                    cor_mat[1:ncol(dat1), (ncol(dat1) + 1):nrow(cor_mat)]
                } else {
                    # single dataset
                    mpmi::cmi(dat1, ...)$bcmi
                }
            },
            #
            # Pearson correlation
            #
            'pearson' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    cor_mat <- cor(cbind(dat1, dat2), method = 'pearson', ...)

                    # limit to cross-dataset correlations
                    cor_mat[1:ncol(dat1), (ncol(dat1) + 1):nrow(cor_mat)]
                } else {
                    cor(dat1, method = 'pearson')
                }
            },

            #
            # Spearman correlation
            #
            'spearman' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    cor_mat <- cor(cbind(dat1, dat2), method = 'spearman', ...)

                    # limit to cross-dataset correlations
                    cor_mat[1:ncol(dat1), (ncol(dat1) + 1):nrow(cor_mat)]
                } else {
                    cor(dat1, method = 'spearman')
                }
            }
        ),

        check_input = function() {
            # TODO: Check to make sure row_data and col_data shared expected
            # keys with dat 
        },

        # Computes cross-dataset correlation matrix
        #
        # Operates on datasets sharing a common column name (i.e. members
        # of col_data).
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use. Supported options:
        #   - pearson  (Pearson correlation)
        #   - spearman (Spearman correlation)
        #   - lm       (Linear model)
        #   - mi       (Mututal information)
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        cross_cor = function(key1='dat', key2=2, method='pearson', ...) {
            # make sure datasets are ordered similarly
            dat1 <- private$col_data[[key1]]
            dat2 <- private$col_data[[key2]]

            # for multidatasets, get underlying data
            if ('EDADataSet' %in% class(dat1)) {
                dat1 <- dat1$dat
                dat2 <- dat2$dat
            }

            # TODO: Include checks?

            # for metadata, we can also exclude factor fields with all
            # unique values (e.g. alternate identifiers)
            #exclude <- apply(dat2, 1, function(x) {
            #    is.factor(x) && length(unique(x)) == length(x)
            #})

            #if (sum(exclude) > 0) {
            #    message(sprintf("Excluding %d unique factor fields", sum(exclude)))
            #    dat2 <- dat2[!exclude, ]

            # only compare matching columns
            col_ind <- intersect(colnames(dat1), colnames(dat2))

            dat1 <- dat1[, col_ind, drop = FALSE]
            dat2 <- dat2[, col_ind, drop = FALSE]

            # measure similarity between rows in datasets 1 and rows in
            # dataset 2; the similarity() method operates on columns so we
            # transposes the datasets first
            cor_mat <- private$similarity(t(dat1), t(dat2), method = method, ...)

            cor_mat
        },

        # Plots multidataset correlation heatmap
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        plot_cross_cor_heatmap = function(key1='dat', key2=2, method='pearson', interactive=TRUE) {
            # compute cross correlations
            cor_mat <- self$cross_cor(key1, key2, method)

            # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if ('row_mdata' %in% names(self$row_data)) {
                mdata1 <- self$row_data[['row_mdata']]
                mask1  <- sapply(mdata1, function(x) max(table(x)) > 1 )
                mdata1 <- mdata1[, mask1, drop = FALSE]

                params[['row_side_colors']] <- mdata1
                params[['subplot_widths']] <- c(0.15, 0.3, 0.55)
            }

            if ('col_mdata' %in% names(self$col_data)) {
                mdata1 <- self$col_data[['col_mdata']]
                mask2  <- sapply(mdata2, function(x) max(table(x)) > 1 )
                mdata2 <- mdata2[, mask2, drop = FALSE]

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
                    col <- params$ColSideColors[, col_name]
                    params$ColSideColors[, col_name] <- colors[as.numeric(factor(col))]
                }
                params$ColSideColors <- as.matrix(params$ColSideColors)
            }

            if ('row_side_colors' %in% names(params)) {
                params$RowSideColors <- params$row_side_colors

                pal <- RColorBrewer::brewer.pal(8, 'Set1')

                for (col_name in colnames(params$RowSideColors)) {
                    col <- params$RowSideColors[, col_name]
                    colors <- colorRampPalette(pal)(min(1E4, length(unique(col))))
                    params$RowSideColors[, col_name] <- colors[as.numeric(factor(col))]
                }
                params$RowSideColors <- as.matrix(params$RowSideColors)
            }

            # remove irrelevant function arguments
            heatmaply_args <- c('showticklabels', 'subplot_widths', 'subplot_heights',
                                'col_side_colors', 'row_side_colors')
            params <- params[!names(params) %in% heatmaply_args]

            do.call(heatmap.plus::heatmap.plus, params)
        },
        #
        # Measures similarity between columns within or across datasets
        #
        similarity = function(dat1, dat2=NULL, method='pearson', ...) {
            # check to make sure specified similarity measure is valid
            if (!method %in% names(private$similarity_measures)) {
                stop('Unsupported similarity method specified.')
            }

            # drop any columns with zero variance
            var_mask1 <- apply(dat1, 2, function(x) length(table(x)) ) > 1

            if (sum(!var_mask1) > 0) {
                message(sprintf("Excluding %d zero-variance entries from first dataset.", sum(!var_mask1)))
                dat1 <- dat1[, var_mask1, drop = FALSE]
            }

            if (!is.null(dat2)) {
                var_mask2 <- apply(dat2, 2, function(x) length(table(x)) ) > 1
                if (sum(!var_mask2) > 0) {
                    message(sprintf("Excluding %d zero-variance entries from second dataset.", sum(!var_mask2)))
                    dat2 <- dat2[, var_mask2, drop = FALSE]
                }
            }

            # construct similarity matrix
            cor_mat <- private$similarity_measures[[method]](dat1, dat2, ...)

            # fix row and column names
            rownames(cor_mat) <- colnames(dat1)

            if (!is.null(dat2)) {
                colnames(cor_mat) <- colnames(dat2)
            } else {
                colnames(cor_mat) <- colnames(dat1)
            }

            cor_mat
        }
    ),

    # ------------------------------------------------------------------------
    # active bindings
    # ------------------------------------------------------------------------
    active = list(
        # main dataset
        dat = function(value) {
            if (missing(value)) {
                private$col_data[['dat']]
            } else {
                private$row_data[['dat']] <- value
                private$col_data[['dat']] <- value
            }
        }
    )
)
