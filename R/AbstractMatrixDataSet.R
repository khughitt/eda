#' An R6 class representing collection of related datasets
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
#'  - `dat`: Primary dataset
#'  - `row_data`: List of additional data keyed on row identifiers
#'  - `col_data`: List of additional data keyed on column identifiers
#'
#' @section Methods:
#'  - `cross_cor(key1=1, key2=2, method='pearson')`: Computes cross-dataset
#'     correlation matrix between rows in two specified datasets.
#'  - `print()`: Prints an overview of the object instance.
#'
#' @section Examples:
#' ```
#' TODO
#' ```
#'
#' @importFrom R6 R6Class
#' @name AbstractMatrixDataSet
#' @export
#'
NULL

AbstractMatrixDataSet <- R6Class("AbstractMatrixDataSet",
    inherit = eda:::AbstractMultiDataSet,

    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # AbstractMatrixDataSet constructor
        initialize = function(datasets, color_pal='Set1', title="", ggplot_theme=theme_bw) {
            super$initialize(datasets, color_pal, title, ggplot_theme)
        },

        # Clusters dataset rows using k-means clustering in a t-SNE projected
        # space.
        #
        # k     Number of clusters to detect (default: 10)
        # ...   Additional arguments passed to Rtsne function
        #
        # return Vector of cluster assignments with length equal to the
        #     number of rows in the dataset.
        cluster_tsne = function(key=1, k=10, ...) {
            tsne <- Rtsne::Rtsne(self$edat[[key]]$dat, ...)
            dat <- setNames(as.data.frame(tsne$Y), c('x', 'y'))

            # Cluster patients from t-sne results
            kmeans_clusters <- kmeans(dat, k)$cluster

            factor(paste0('cluster_', kmeans_clusters))
        },

        # Detects column outliers in the dataset
        #
        # Computes pairwise correlations between all columns in the dataset.
        # Columns whose average pairwise correlation with all other columns
        # is greater than `num_sd` standard deviations from the average
        # average correlation are considered to be outliers.
        #
        # num_sd    Number of standard deviations to use to determine outliers.
        # ctend     Measure of central tendency to use (default: median)
        # method
        #
        # return Character vector or column ids for columns with low
        #     average pairwise correlations.
        detect_col_outliers = function(key=1, num_sd=2, ctend=median, method='pearson', ...) {
            # TODO: include correlation in results?
            # TODO: Write alternative version for data frame datasets?
            dat <- self$edat[[key]]$dat
            cor_mat <- private$similarity(dat, method=method, ...)

            avg_column_cors <- apply(cor_mat, 1, ctend)
            cutoff <- mean(avg_column_cors) - (num_sd * sd(avg_column_cors))
            colnames(dat)[avg_column_cors < cutoff]
        },

        # Detects row outliers in the dataset
        #
        # Computes pairwise correlations between all rows in the dataset.
        # Rows whose median pairwise correlation with all other rows
        # is greater than `num_sd` standard deviations from the average
        # median correlation are considered to be outliers.
        #
        # num_sd Number of standard deviations to use to determine
        #      outliers.
        #
        # return Character vector or row ids for rows with low
        #     average pairwise correlations.
        detect_row_outliers = function(key=1, num_sd=2, ctend=median, method='pearson', ...) {
            dat <- self$edat[[key]]$dat
            cor_mat <- private$similarity(t(dat), method=method, ...)

            avg_row_cors <- apply(cor_mat, 1, ctend)
            cutoff <- mean(avg_row_cors) - num_sd * sd(avg_row_cors)
            rownames(dat)[avg_row_cors < cutoff]
        },

        # Removes column outliers from the dataset
        #
        # Computes pairwise correlations between all columns in the dataset.
        # Columns whose median pairwise correlation with all other columns
        # is greater than `num_sd` standard deviations from the average
        # median correlation are removed.
        #
        # num_sd Number of standard deviations to use to determine
        #      outliers.
        #
        # return A filtered version of the original EDAMatrix object.
        filter_col_outliers = function(key=1, num_sd=2, ctend=median, method='pearson') {
            obj <- private$clone_()
            outliers <- obj$detect_col_outliers(key, num_sd, avg, method)
            obj$filter_cols(!colnames(obj$dat) %in% outliers)
        },

        # Removes row outliers from the dataset
        #
        # Computes pairwise correlations between all rows in the dataset.
        # rows whose median pairwise correlation with all other rows
        # is greater than `num_sd` standard deviations from the average
        # median correlation are removed.
        #
        # num_sd Number of standard deviations to use to determine
        #      outliers.
        #
        # return A filtered version of the original EDAMatrix object.
        filter_row_outliers = function(num_sd=2, ctend=median, method='pearson') {
            obj <- private$clone_()
            outliers <- obj$detect_row_outliers(num_sd, ctend, method)
            obj$filter_rows(!rownames(obj$dat) %in% outliers)
        },

        # Log-transforms data
        #
        # base Numeric logarithm base to use (default: e)
        # offset Numeric offset to apply to data before taking the
        #     logarithm (default: 0)
        #
        # return A log-transformed version of the object.
        log = function(key=1, base=exp(1), offset=0) {
            obj <- private$clone_()
            obj$edat[[key]]$dat <- log(obj$edat[[key]]$dat + offset, base)
            obj
        },

        # Log(x + 1) transforms data
        #
        # return A log(x + 1) transformed version of the object.
        log1p = function() {
            obj <- private$clone_()
            obj$edat[[key]]$dat <- log(obj$edat[[key]]$dat + 1)
            obj
        },

        # PCA
        #
        pca = function(key=1, ...) {
            obj <- private$clone_()
            obj$edat[[key]]$dat <- prcomp(obj$edat[[key]]$dat, ...)$x
            obj$col_mdata <- NULL
            obj
        },

        #
        # t-SNE
        #
        tsne = function(key=1, ...) {
            obj <- private$clone_()
            obj$edat[[key]]$dat <- Rtsne::Rtsne(obj$edat[[key]]$dat, ...)
            obj$col_mdata <- NULL
            obj
        },

        ######################################################################
        # plotting methods
        ######################################################################

        # Correlation heatmap.
        #
        # Generates a correlation heatmap depicting the pairwise column
        # correlations in the data.
        #
        # method String name of correlation method to use.
        # ... Additional arguments
        #
        # @seealso \code{cor} for more information about supported correlation
        #      methods.
        plot_cor_heatmap = function(key=1, method='pearson', interactive=TRUE, ...) {
            # generate correlation matrix
            cor_mat <- private$similarity(self$edat[[key]]$dat, method=method, ...)


            # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$col_mdata)) {
                # for heatmaps, show binary/logical variables on one side of the heatmap and
                # other variables on the other side
                lens <- apply(self$col_mdata, 2, function(x) {
                    length(unique(x))
                })
                binary_vars <- lens == 2

                # column colors (binary variables)
                if (sum(binary_vars) >= 1) {
                    params[['col_side_colors']] <- self$col_mdata[, binary_vars, drop = FALSE]
                    params[['subplot_heights']] <- c(0.15, 0.3, 0.55)
                }

                # column colors (everything else)
                if (sum(!binary_vars) >= 1) {
                    params[['col_side_colors']] <- self$col_mdata[, !binary_vars, drop = FALSE]
                    params[['subplot_widths']]  <- c(0.55, 0.3, 0.15)
                }
            }

            # add any additional function arguments
            params <- c(params, ...)

            private$construct_heatmap_plot(params, interactive)
        },

        # Generates a heatmap plot of the dataset
        #
        # ... Additional arguments
        plot_heatmap = function(key, interactive=TRUE, ...) {
            # list of parameters to pass to heatmaply
            params <- list(
                x               = self$edat[[key]]$dat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$row_mdata)) {
                params[['row_side_colors']] <- self$row_mdata
                params[['subplot_widths']]  <- c(0.15, 0.3, 0.55)
            }

            if (!is.null(self$col_mdata)) {
                params[['col_side_colors']] <- self$col_mdata
                params[['subplot_heights']] <- c(0.55, 0.3, 0.15)
            }

            # add any additional function arguments
            params <- c(params, ...)

            private$construct_heatmap_plot(params, interactive)
        },

        # Generates a two-dimensional PCA plot from the dataset
        #
        # pcx integer PC number to plot along x-axis (default: 1)
        # pcy integer PC number to plot along x-axis (default: 2)
        # scale Whether or not to scale variables prior to performing
        #     pca; passed to `prcomp` function.
        # color Column metadata field to use for coloring points.
        # shape Column metadata field to use to assign shapes to points
        # title Plot title.
        # text_labels Whether or not to include individual point labels
        #     plot (default: FALSE).
        # ...
        #
        # return ggplot plot instance
        plot_pca = function(key=1, pcx=1, pcy=2, scale=FALSE,
                            color_var=NULL, color_key=NULL,
                            shape_var=NULL, shape_key=NULL,
                            label_var=NULL, label_key=NULL,
                            title=NULL, text_labels=FALSE, ...) {
            dat <- self$edat[[key]]$dat

            # perform pca
            prcomp_results <- prcomp(dat, scale = scale)
            var_explained <- round(summary(prcomp_results)$importance[2, ] * 100, 2)

            # create data frame for plotting
            res <- data.frame(id = rownames(dat),
                              pc1 = prcomp_results$x[, pcx],
                              pc2 = prcomp_results$x[, pcy])

            # get color/shape styles
            styles <- private$get_geom_point_styles(key,
                                                    color_var, color_key,
                                                    shape_var, shape_key)

            if (!is.null(styles$color)) {
                res <- cbind(res, color_var = styles$color)
            }
            if (!is.null(styles$shape)) {
                res <- cbind(res, shape_var = styles$shape)
            }

            xl <- sprintf("PC%d (%.2f%% variance)", pcx, var_explained[pcx])
            yl <- sprintf("PC%d (%.2f%% variance)", pcy, var_explained[pcy])

            # plot title
            if (is.null(title)) {
                title <- sprintf("PCA: %s", private$title)
            }

            # PC1 vs PC2
            plt <- ggplot(res, aes(pc1, pc2)) +
                geom_point(stat = "identity", styles$aes, size = 0.5) +
                xlab(xl) + ylab(yl) +
                ggtitle(title) +
                private$ggplot_theme() +
                theme(axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = -90))

            # text labels
            if (text_labels) {
                plt <- plt + geom_text(aes(label = id), angle = 45,
                                       size = 0.5, vjust = 2)
            }

            # legend labels
            if (length(styles$labels) > 0) {
                plt <- plt + styles$labels
            }
            plt
        },

        # Generates a two-dimensional t-SNE plot from the dataset
        #
        # color Column metadata field to use for coloring points.
        # shape Column metadata field to use to assign shapes to points
        # title Plot title.
        # text_labels Whether or not to include individual point labels
        #     plot (default: FALSE).
        # ...
        #
        # return ggplot plot instance
        plot_tsne = function(key=1,
                             color_var=NULL, color_key=NULL,
                             shape_var=NULL, shape_key=NULL,
                             label_var=NULL, label_key=NULL,
                             title=NULL, text_labels=FALSE, ...) {
            dat <- self$edat[[key]]$dat

            # compute t-SNE projection
            tsne <- Rtsne::Rtsne(dat, ...)
            res <- setNames(as.data.frame(tsne$Y), c('x', 'y'))

            # add column ids
            res <- cbind(res, id = rownames(dat))

            # get color/shape styles
            styles <- private$get_geom_point_styles(key,
                                                    color_var, color_key,
                                                    shape_var, shape_key)

            if (!is.null(styles$color)) {
                res <- cbind(res, color_var = styles$color)
            }
            if (!is.null(styles$shape)) {
                res <- cbind(res, shape_var = styles$shape)
            }

            # plot title
            if (is.null(title)) {
                title <- sprintf("t-SNE: %s", private$title)
            }

            # treatment response
            plt <- ggplot(res, aes(x, y)) +
                   geom_point(styles$aes, stat = "identity", size = 0.5) +
                   ggtitle(title) +
                   private$ggplot_theme() +
                   theme(axis.ticks = element_blank(),
                         axis.text.x = element_text(angle = -90),
                         legend.text = element_text(size = 7))

            # text labels
            if (text_labels) {
                plt <- plt + geom_text(aes(label = id), angle = 45, size = 0.5,
                                       vjust = 2)
            }

            # legend labels
            if (length(styles$labels) > 0) {
                plt <- plt + styles$labels
            }

            plt
        },

        # Plot median pairwise column correlations
        #
        # Plots the median correlation of each variable (column). This is
        # useful for visually inspecting columns for possible outliers, when
        # the total number of columns is relatively small.
        #
        # @author V. Keith Hughitt, \email{keith.hughitt@nih.gov}
        #
        # return None
        plot_pairwise_column_cors = function(key=1,
                                             color_var=NULL, color_key=NULL,
                                             label_var=NULL, label_key=NULL,
                                             title="", method='pearson',
                                             mar=c(12, 6, 4, 6), ...) {
            dat <- self$edat[[key]]$dat

            # compute pairwise variable correlations
            cor_mat <- private$similarity(dat, method=method, ...)
            median_pairwise_cor <- apply(cor_mat, 1, median)

            quantiles <- quantile(median_pairwise_cor, probs = c(0.25, 0.75))
            iqr <- diff(quantiles)

            #outlimit
            cutoff <- quantiles[1] - 1.5 * iqr

            ylimit <- c(pmin(min(median_pairwise_cor), cutoff),
                        max(median_pairwise_cor))

            # get color properties
            color_vector <- private$get_var_colors(key, color_var, color_key)
            label_vector <- private$get_var_labels(key, label_var, label_key)

            # variable labels
            if (!all(colnames(dat) == label_vector)) {
                var_labels <- sprintf("%s (%s)", colnames(dat), label_vector)
            } else {
                var_labels <- colnames(dat)
            }

            # render plot
            par(mar = mar)
            plot(median_pairwise_cor, xaxt = "n", ylim = ylimit,
                 ylab = "Median Pairwise Correlation", xlab = "", main = title,
                 col = color_vector, pch = 16, cex = 2.2)
            axis(side = 1, at = seq(along = median_pairwise_cor),
                 labels = var_labels, las = 2)
            abline(h = cutoff, lty = 2)
            abline(v = 1:length(var_labels), lty = 3, col = "black")
        }
    )
)
