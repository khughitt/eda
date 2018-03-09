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
        initialize = function(dat, 
                              row_mdata=NULL, col_mdata=NULL, 
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
                              row_maxn=Inf, row_max_ratio=1.0, row_ind=NULL,
                              col_maxn=Inf, col_max_ratio=1.0, col_ind=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) { 
            # verify input data type and call parent constructor
            private$check_input(dat)

            super$initialize(dat, row_mdata, col_mdata, row_ids, col_ids,
                             row_mdata_ids, col_mdata_ids, row_color, row_shape,
                             row_labels, col_color, col_shape, col_labels,
                             row_maxn, row_max_ratio, row_ind,
                             col_maxn, col_max_ratio, col_ind,
                             color_pal, title, ggplot_theme)
        },

        #' Detects column outliers in the dataset
        #'
        #' Computes pairwise correlations between all columns in the dataset.
        #' Columns whose median pairwise correlation with all other columns
        #' is greater than \code{num_sd} standard deviations from the average
        #' median correlation are considered to be outliers.
        #'
        #' @param num_sd Number of standard deviations to use to determine
        #'      outliers.
        #'
        #' @return Character vector or column ids for columns with low
        #'     average pairwise correlations.
        detect_col_outliers = function(num_sd=2) {
            # TODO: include correlation in results?
            # TODO: Write alternative version for data frame datasets?
            cor_mat <- cor(self$dat)
            median_column_cors <- apply(cor_mat, 1, median)
            cutoff <- mean(median_column_cors) - (num_sd * sd(median_column_cors))
            colnames(self$dat)[median_column_cors < cutoff]
        },

        #' Detects row outliers in the dataset
        #'
        #' Computes pairwise correlations between all rows in the dataset.
        #' Rows whose median pairwise correlation with all other rows
        #' is greater than \code{num_sd} standard deviations from the average
        #' median correlation are considered to be outliers.
        #'
        #' @param num_sd Number of standard deviations to use to determine
        #'      outliers.
        #'
        #' @return Character vector or row ids for rows with low
        #'     average pairwise correlations.
        detect_row_outliers = function(num_sd=2) {
            cor_mat <- cor(t(self$dat))
            median_row_cors <- apply(cor_mat, 1, median)
            cutoff <- mean(median_row_cors) - num_sd * sd(median_row_cors)
            rownames(self$dat)[median_row_cors < cutoff]
        },

        #' Log-transforms data
        #'
        #' @param base Numeric logarithm base to use (default: e)
        #' @param offset Numeric offset to apply to data before taking the
        #'     logarithm (default: 0)
        #'
        #' @return A log-transformed version of the object.
        log = function(base=exp(1), offset=0) {
            obj <- private$clone_()
            obj$dat <- log(obj$dat + offset, base)
            obj
        },

        #' Log(x + 1) transforms data
        #'
        #' @return A log(x + 1) transformed version of the object.
        log1p = function() {
            obj <- private$clone_()
            obj$dat <- log(obj$dat + 1)
            obj
        },

        #' Removes column outliers from the dataset
        #'
        #' Computes pairwise correlations between all columns in the dataset.
        #' Columns whose median pairwise correlation with all other columns
        #' is greater than \code{num_sd} standard deviations from the average
        #' median correlation are removed.
        #'
        #' @param num_sd Number of standard deviations to use to determine
        #'      outliers.
        #'
        #' @return A filtered version of the original EDADataSet object.
        filter_col_outliers = function(num_sd=2) {
            obj <- private$clone_()
            cor_mat <- cor(obj$dat)
            median_column_cors <- apply(cor_mat, 1, median)
            cutoff <- mean(median_column_cors) - num_sd * sd(median_column_cors)
            obj$filter_cols(median_column_cors > cutoff)
        },

        #' Removes row outliers from the dataset
        #'
        #' Computes pairwise correlations between all rows in the dataset.
        #' rows whose median pairwise correlation with all other rows
        #' is greater than \code{num_sd} standard deviations from the average
        #' median correlation are removed.
        #'
        #' @param num_sd Number of standard deviations to use to determine
        #'      outliers.
        #'
        #' @return A filtered version of the original EDADataSet object.
        filter_row_outliers = function(num_sd=2) {
            obj <- private$clone_()
            cor_mat <- cor(t(obj$dat))
            median_row_cors <- apply(cor_mat, 1, median)
            cutoff <- mean(median_row_cors) - num_sd * sd(median_row_cors)
            obj$filter_rows(median_row_cors > cutoff)
        },

        #' Computes correlations between data principle components and column
        #' metadata entries (features).
        #' 
        #' Measures the predictive power of each feature (column in 
        #' \code{obj$col_mdata}) and the principle components (PCs) of the 
        #' dataset using a simple linear model.
        #'
        #' Based on code adapted from cbcbSEQ
        #' (https://github.com/kokrah/cbcbSEQ/) originally written by 
        #' Kwame Okrah.
        #'
        #' @param include Vector of strings indicating metadata columns which
        #' should be included in the analysis; takes priority over exclude
        #' if both are set.
        #' @param exclude Vector of strings indicating metadata columns which
        #' should be excluded from the analysis.
        #'
        #' @return Dataframe containing feature/PC correlations, as well as
        #'      information about the amount of variance explained by each PC.
        get_pca_feature_correlations = function(include=NULL, exclude=NULL) {
            # If already computed, return cached result
            #if (!is.null(private$cache[['pca_feature_cor']])) {
            #    return(private$cache[['pca_feature_cor']])
            #}

            # SVD
            s <- corpcor:::fast.svd(t(self$dat) - rowMeans(t(self$dat)))
            rownames(s$v) <- rownames(self$dat)

            # Create output dataframe
            pc_var <- round((s$d^2) / sum(s$d^2) * 100, 2)

            result <- data.frame(
                pc_var=pc_var,
                pc_var_cum=cumsum(pc_var) 
            )

            # determine which features to include
            if (!is.null(include)) {
                include <- include
            } else if (!is.null(exclude)) {
                include <- colnames(self$row_mdata)[!colnames(self$row_mdata) %in% exclude]
            } else {
                include <- colnames(self$row_mdata)
            }

            # measure feature correlations and add to result data frame
            result <- cbind(result, private$compute_feature_correlations(s$v, include))
            rownames(result) <- paste0("PC", 1:nrow(result)) 

            # cache result and return
            private$cache[['pca_feature_cor']] <- result

            result
        },

        #' Computes correlations between data principle components and column
        #' metadata entries (features).
        #' 
        #' Measures the predictive power of each feature (column in 
        #' \code{obj$row_mdata}) and the t-SNE projected axes of the 
        #' dataset using a simple linear model.
        #'
        #' Based on code adapted from cbcbSEQ
        #' (https://github.com/kokrah/cbcbSEQ/) originally written by 
        #' Kwame Okrah.
        #'
        #' @param include Vector of strings indicating metadata columns which
        #' should be included in the analysis.
        #' @param exclude Vector of strings indicating metadata columns which
        #' should be excluded from the analysis.
        #'
        #' @return Dataframe containing feature/t-SNE axes correlations.
        get_tsne_feature_correlations = function(include=NULL, exclude=NULL, ...) {
            # If already computed, return cached result
            #if (!is.null(private$cache[['tsne_feature_cor']])) {
            #    return(private$cache[['tsne_feature_cor']])
            #}

            # t-SNE
            tsne <- Rtsne::Rtsne(self$dat, ...)

            # determine which features to include
            if (!is.null(include)) {
                include <- include
            } else if (!is.null(exclude)) {
                include <- colnames(self$row_mdata)[!colnames(self$row_mdata) %in% exclude]
            } else {
                include <- colnames(self$row_mdata)
            }

            # measure feature correlations and add to result data frame
            result <- private$compute_feature_correlations(tsne$Y, include)
            rownames(result) <- paste0("Dim", 1:nrow(result)) 

            # cache result and return
            private$cache[['tsne_feature_cor']] <- result

            result
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
        },

        ######################################################################
        # plotting methods
        ######################################################################

        #' Correlation heatmap. 
        #'
        #' Generates a correlation heatmap depicting the pairwise column 
        #' correlations in the data.
        #'
        #' @param method String name of correlation method to use.
        #' @param ... Additional arguments
        #'
        #' @seealso \code{cor} for more information about supported correlation 
        #'      methods.
        plot_cor_heatmap = function(method='pearson', interactive=TRUE, ...) { 
            # determine subsampling indices, if requested
            indices <- private$get_indices(...)

            # generate correlation matrix
            cor_mat <- cor(t(self$dat[indices$row, indices$col]), method=method)

            # list of parameters to pass to heatmaply
            params <- list(
                x=cor_mat,
                showticklabels=c(FALSE, FALSE),
                subplot_widths=c(0.65, 0.35),
                subplot_heights=c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$row_mdata)) {
                # for heatmaps, show binary/logical variables on one side of the heatmap and
                # other variables on the other side
                lens <- apply(self$row_mdata, 2, function(x) { length(unique(x)) })
                binary_vars <- lens == 2

                # col colors (binary variables)
                if (sum(binary_vars) >= 1) {
                    params[['col_side_colors']] <- self$row_mdata[indices$row, 
                                                                  binary_vars, drop=FALSE]
                    params[['subplot_heights']] <- c(0.15, 0.3, 0.55)
                }

                # row colors (everything else)
                if (sum(!binary_vars) >= 1) {
                    params[['row_side_colors']] <- self$row_mdata[indices$row, 
                                                                  !binary_vars, drop=FALSE]
                    params[['subplot_widths']] <- c(0.55, 0.3, 0.15)
                }
            }

            # add any additional function arguments
            params <- c(params, private$get_custom_function_args(...))

            private$plot_heatmap(params, interactive)
        },

        #' Generates a two-dimensional PCA plot from the dataset 
        #'
        #' @param dat Data matrix to generate PCA plot for
        #' @param pcx integer PC number to plot along x-axis (default: 1)
        #' @param pcy integer PC number to plot along x-axis (default: 2)
        #'
        #' @return ggplot plot instance
        plot_pca = function(pcx=1, pcy=2, scale=FALSE,
                            color=NULL, shape=NULL, title=NULL,
                            text_labels=FALSE, ...) {
            # determine subsampling indices, if requested
            indices <- private$get_indices(...)

            # perform pca
            prcomp_results <- prcomp(self$dat[indices$row, indices$col], scale=scale)
            var_explained <- round(summary(prcomp_results)$importance[2,] * 100, 2)

            # create data frame for plotting
            dat <- data.frame(id=rownames(self$dat)[indices$row], 
                              pc1=prcomp_results$x[,pcx],
                              pc2=prcomp_results$x[,pcy])

            # get color/shape styles
            styles <- private$get_geom_point_styles(color, shape)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color[indices$row])
            }
            if (!is.null(styles$shape)) {
                dat <- cbind(dat, shape=styles$shape[indices$row])
            }

            xl <- sprintf("PC%d (%.2f%% variance)", pcx, var_explained[pcx])
            yl <- sprintf("PC%d (%.2f%% variance)", pcy, var_explained[pcy])

            # plot title
            if (is.null(title)) {
                title <- sprintf("PCA: %s", private$title)
            }

            # PC1 vs PC2
            plt <- ggplot(dat, aes(pc1, pc2)) +
                geom_point(stat="identity", styles$aes, size=0.5) +
                xlab(xl) + ylab(yl) +
                ggtitle(title) +
                private$ggplot_theme() +
                theme(axis.ticks=element_blank(), 
                      axis.text.x=element_text(angle=-90))

            # text labels
            if (text_labels) {
                plt <- plt + geom_text(aes(label=id), angle=45, size=0.5, vjust=2)
            }

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}
            plt
        },

		#' Plots PCA / feature correlations as a colored grid.
		#'
		#' @param num_pcs Number of principle components to include (default: 6)
        #' @param exclude Features (column metadata variables) to exclude from 
        #'     the analysis.
        #' @param low_color String indicating color to use for low correlation
        #'     values (default: green)
        #' @param high_color String indicating color to use for high correlation
        #'     values (default: red)
		#'
		#' @return ggplot plot instance
        plot_pca_feature_correlations = function(num_pcs=6, include=NULL, 
                                                 exclude=NULL, low_color='green',
											     high_color='red', ...) {
            # compute pca feature correlations or retrieved cached version
            if (!is.null(private$cache[['pca_feature_cor']])) {
                pca_cor <- private$cache[['pca_feature_cor']]
            } else {
                pca_cor <- self$get_pca_feature_correlations(include, exclude)
            }

            pca_cor <- pca_cor[1:num_pcs,]
            pca_long <- cbind(PC=sprintf("PC%d (%0.2f%%)", 
                                         1:nrow(pca_cor), pca_cor$pc_var), 
                              melt(pca_cor, id.vars=c('pc_var', 'pc_var_cum')))

            ggplot(pca_long, aes(x=PC, y=variable)) + 
                geom_tile(aes(fill=value)) + 
                geom_text(aes(label=value), size=2, show.legend=FALSE) +
                scale_fill_gradient(low=low_color, high=high_color) +
                private$ggplot_theme() +
                theme(axis.text.x=element_text(size=8, angle=45, vjust=1, hjust=1), 
                      axis.text.y=element_text(size=8)) + 
                xlab("Principle Components (% var explained)") +
                ylab("Features") +
                guides(fill=guide_legend("R^2"))
        },

		#' Plots t-SNE / feature correlations as a colored grid.
		#'
        #' @param exclude Features (column metadata variables) to exclude from 
        #'     the analysis.
        #' @param low_color String indicating color to use for low correlation
        #'     values (default: green)
        #' @param high_color String indicating color to use for high correlation
        #'     values (default: red)
		#'
		#' @return ggplot plot instance
		plot_tsne_feature_correlations = function(include=NULL, exclude=NULL,
												  low_color='green', 
												  high_color='red', ...) {
            # compute t-SNE feature correlations or retrieved cached version
            if (!is.null(private$cache[['tsne_feature_cor']])) {
                tsne_cor <- private$cache[['tsne_feature_cor']]
            } else {
                tsne_cor <- self$get_tsne_feature_correlations(include, exclude, ...)
            }

            tsne_long <- melt(tsne_cor)
            colnames(tsne_long) <- c('dim', 'variable', 'value')

            ggplot(tsne_long, aes(x=dim, y=variable)) + 
                geom_tile(aes(fill=value)) + 
                geom_text(aes(label=value), size=2, show.legend=FALSE) +
                scale_fill_gradient(low=low_color, high=high_color) +
                private$ggplot_theme() +
                theme(axis.text.x=element_text(size=8, angle=45, vjust=1, hjust=1), 
                      axis.text.y=element_text(size=8)) + 
                xlab("t-SNE dimension") +
                ylab("Features") +
                guides(fill=guide_legend("R^2"))
        },

        cluster_tsne = function(k=10, ...) {
            # perform t-sne and store results
            #if (is.null(private$cache[['tsne']])) {
            #    private$cache[['tsne']] <- Rtsne::Rtsne(t(self$dat), ...)
            #}
            #tsne_res <- as.data.frame(private$cache[['tsne']]$Y)

            # for clustering, we want to ensure that all points are assigned
            # to a cluster, so include all indices
            indices <- private$get_indices(row_max_ratio=1, col_max_ratio=1)

            tsne <- Rtsne::Rtsne(self$dat[indices$row, indices$col], ...)
            dat <- setNames(as.data.frame(tsne$Y), c('x', 'y'))

            # Cluster patients from t-sne results
            kmeans_clusters <- kmeans(dat, k)$cluster

            factor(paste0('cluster_', kmeans_clusters))
        },

        plot_tsne = function(color=NULL, shape=NULL, title=NULL,
                             text_labels=FALSE, ...) {
            # 2018/02/27: Disabling cache for now until it can be improved..
            # perform t-sne and store results
            #if (is.null(private$cache[['tsne']])) {
            #    private$cache[['tsne']] <- Rtsne::Rtsne(t(self$dat), ...)
            #}
            #dat <- as.data.frame(private$cache[['tsne']]$Y)

            # determine subsampling indices, if requested
            indices <- private$get_indices(...)

            tsne <- Rtsne::Rtsne(self$dat[indices$row, indices$col], ...)
            dat <- setNames(as.data.frame(tsne$Y), c('x', 'y'))

            # add column ids
            dat <- cbind(dat, id=rownames(self$dat))

            # get color/shape styles
            styles <- private$get_geom_point_styles(color, shape)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (!is.null(styles$shape)) {
                dat <- cbind(dat, shape=styles$shape)
            }

            # plot title
            if (is.null(title)) {
                title <- sprintf("t-SNE: %s", private$title)
            }

            # treatment response
            plt <- ggplot(dat, aes(x, y)) +
                   geom_point(styles$aes, stat="identity", size=0.5) +
                   ggtitle(title) +
                   private$ggplot_theme() +
                   theme(axis.ticks=element_blank(), 
                         axis.text.x=element_text(angle=-90),
                         legend.text=element_text(size=7))

            # text labels
            if (text_labels) {
                plt <- plt + geom_text(aes(label=id), angle=45, size=0.5, vjust=2)
            }

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}

            plt
        },

        #' Plot median pairwise variable correlations
        #'
        #' Plots the median correlation of each variable (column)
        #'
        #' @author V. Keith Hughitt, \email{keith.hughitt@nih.gov}
        #'
        #' @return None
        plot_var_correlations = function (color=NULL, title="", method='pearson',
                                          mar=c(12,6,4,6)) {
            # compute pairwise variable correlations
            median_pairwise_cor <- apply(cor(self$dat * 1.0, method=method), 1, median)

            quantiles <- quantile(median_pairwise_cor, probs=c(0.25, 0.75))
            iqr <- diff(quantiles)

            #outlimit
            cutoff <- quantiles[1] - 1.5 * iqr

            ylimit <- c(pmin(min(median_pairwise_cor), cutoff), 
                        max(median_pairwise_cor))

            # get color properties
			color_vector <- get_var_colors(color)
			label_vector <- get_var_labels()

            # variable labels
            if (!all(colnames(self$dat) == label_vector)) {
                var_labels <- sprintf("%s (%s)", colnames(self$dat), label_vector)
            } else {
                var_labels <- colnames(self$dat)
            }

            # render plot
            par(mar=mar)
            plot(median_pairwise_cor, xaxt="n", ylim=ylimit,
                 ylab="Median Pairwise Correlation", xlab="", main=title,
                 col=color_vector, pch=16, cex=2.2)
            axis(side=1, at=seq(along=median_pairwise_cor),
                 labels=var_labels, las=2)
            abline(h=cutoff, lty=2)
            abline(v=1:length(var_labels), lty=3, col="black")
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

