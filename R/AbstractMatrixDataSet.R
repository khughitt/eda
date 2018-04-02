#' An S6 class representing collection of related datasets
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
        # EDADataSet constructor
        initialize = function(datasets,
                              row_color=NULL, row_color_ds='dat',
                              row_shape=NULL, row_shape_ds='dat',
                              row_label=NULL, row_label_ds='dat',
                              col_color=NULL, col_color_ds='dat',
                              col_shape=NULL, col_shape_ds='dat',
                              col_label=NULL, col_label_ds='dat',
                              color_pal='Set1', title="", ggplot_theme=theme_bw) {
            super$initialize(datasets,
                             row_color, row_color_ds, row_shape, row_shape_ds,
                             row_label, row_label_ds, col_color, col_color_ds,
                             col_shape, col_shape_ds, col_label, col_label_ds,
                             color_pal, title, ggplot_theme)
        },

        # Clusters dataset rows using k-means clustering in a t-SNE projected
        # space.
        #
        # k     Number of clusters to detect (default: 10)
        # ...   Additional arguments passed to Rtsne function
        #
        # return Vector of cluster assignments with length equal to the
        #     number of rows in the dataset.
        cluster_tsne = function(k=10, ...) {
            tsne <- Rtsne::Rtsne(self$dat, ...)
            dat <- setNames(as.data.frame(tsne$Y), c('x', 'y'))

            # Cluster patients from t-sne results
            kmeans_clusters <- kmeans(dat, k)$cluster

            factor(paste0('cluster_', kmeans_clusters))
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
        detect_col_outliers = function(num_sd=2, ctend=median, method='pearson', ...) {
            # TODO: include correlation in results?
            # TODO: Write alternative version for data frame datasets?
            cor_mat <- private$similarity(self$dat, method=method, ...)

            avg_column_cors <- apply(cor_mat, 1, ctend)
            cutoff <- mean(avg_column_cors) - (num_sd * sd(avg_column_cors))
            colnames(self$dat)[avg_column_cors < cutoff]
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
        detect_row_outliers = function(num_sd=2, ctend=median, method='pearson', ...) {
            cor_mat <- private$similarity(t(self$dat), method=method, ...)

            avg_row_cors <- apply(cor_mat, 1, ctend)
            cutoff <- mean(avg_row_cors) - num_sd * sd(avg_row_cors)
            rownames(self$dat)[avg_row_cors < cutoff]
        },

        # Detects dependencies between column metadata entries (features) and
        # dataset rows.
        feature_cor = function(method='pearson') {
            if (is.null(self$col_mdata)) {
                stop("Error: missing column metadata.")
            }
            private$cross_cor('dat', 'col_mdata', method)
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
        filter_col_outliers = function(num_sd=2, ctend=median, method='pearson') {
            obj <- private$clone_()
            outliers <- obj$detect_col_outliers(num_sd, avg, method)
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
        log = function(base=exp(1), offset=0) {
            obj <- private$clone_()
            obj$dat <- log(obj$dat + offset, base)
            obj
        },

        # Log(x + 1) transforms data
        #
        # return A log(x + 1) transformed version of the object.
        log1p = function() {
            obj <- private$clone_()
            obj$dat <- log(obj$dat + 1)
            obj
        },

        # PCA
        #
        pca = function(...) {
            obj <- private$clone_()
            obj$dat <- prcomp(self$dat, ...)$x
            obj$col_mdata <- NULL
            obj
        },

        #
        pca_feature_cor = function(method='pearson', ...) {
            self$t$pca(...)$t$feature_cor(method)
        },

        #
        # t-SNE
        #
        tsne = function(...) {
            obj <- private$clone_()
            obj$dat <- Rtsne::Rtsne(obj$dat, ...)
            obj$col_mdata <- NULL
            obj
        },

        tsne_feature_cor = function(method='pearson', ...) {
            self$t$tsne(...)$t$feature_cor(...)
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
        plot_cor_heatmap = function(method='pearson', interactive=TRUE, ...) {
            # generate correlation matrix
            cor_mat <- private$similarity(self$dat, method=method, ...)

            # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$row_mdata)) {
                # for heatmaps, show binary/logical variables on one side of the heatmap and
                # other variables on the other side
                lens <- apply(self$row_mdata, 2, function(x) {
                    length(unique(x))
                })
                binary_vars <- lens == 2

                # col colors (binary variables)
                if (sum(binary_vars) >= 1) {
                    params[['col_side_colors']] <- self$row_mdata[, binary_vars, drop = FALSE]
                    params[['subplot_heights']] <- c(0.15, 0.3, 0.55)
                }

                # row colors (everything else)
                if (sum(!binary_vars) >= 1) {
                    params[['row_side_colors']] <- self$row_mdata[, !binary_vars, drop = FALSE]
                    params[['subplot_widths']]  <- c(0.55, 0.3, 0.15)
                }
            }

            # add any additional function arguments
            params <- c(params, ...)

            private$plot_heatmap(params, interactive)
        },

        # Generates a heatmap plot of the dataset
        #
        # ... Additional arguments
        plot_heatmap = function(interactive=TRUE, ...) {
            # list of parameters to pass to heatmaply
            params <- list(
                x               = self$dat,
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

            private$plot_heatmap(params, interactive)
        },

        # Creates a tile plot of projected data / feature correlations
        #
        # include Vector of strings indicating metadata columns which
        # should be included in the analysis.
        # exclude Features (column metadata variables) to exclude from
        #     the analysis.
        # color_scale Character vector containing colors to sue for
        #     low-correlation and high-correlation values
        #     (default: c('green', 'red')).
        #
        # return ggplot plot instance
        plot_feature_cor = function(method='pearson', color_scale=c('green', 'red')) {
            # compute feature correlations
            dat <- melt(private$feature_cor(method = method))
            colnames(dat) <- c('dim', 'variable', 'value')

            # Labels
            #if (method == 'pca') {
            #    xlab_text <- 'Principle Components'
            #} else if (method == 't-sne') {
            #    xlab_text <- "t-SNE dimension"
            #}
            xlab_text <- 'TODO...'

            # create plot
            ggplot(dat, aes(x = dim, y = variable)) +
                geom_tile(aes(fill = value)) +
                geom_text(aes(label = value), size = 2, show.legend = FALSE) +
                scale_fill_gradient(low = color_scale[1], high = color_scale[2]) +
                private$ggplot_theme() +
                theme(axis.text.x = element_text(size = 8, angle = 45,
                                                 vjust = 1, hjust = 1),
                      axis.text.y = element_text(size = 8)) +
                xlab(xlab_text) +
                ylab("Features") +
                guides(fill = guide_legend("R^2"))
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
        plot_pca = function(pcx=1, pcy=2, scale=FALSE,
                            color=NULL, shape=NULL, title=NULL,
                            text_labels=FALSE, ...) {
            # perform pca
            prcomp_results <- prcomp(self$dat, scale = scale)
            var_explained <- round(summary(prcomp_results)$importance[2, ] * 100, 2)

            # create data frame for plotting
            dat <- data.frame(id = rownames(self$dat),
                              pc1 = prcomp_results$x[, pcx],
                              pc2 = prcomp_results$x[, pcy])

            # get color/shape styles
            styles <- private$get_geom_point_styles(color, shape)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color = styles$color)
            }
            if (!is.null(styles$shape)) {
                dat <- cbind(dat, shape = styles$shape)
            }

            xl <- sprintf("PC%d (%.2f%% variance)", pcx, var_explained[pcx])
            yl <- sprintf("PC%d (%.2f%% variance)", pcy, var_explained[pcy])

            # plot title
            if (is.null(title)) {
                title <- sprintf("PCA: %s", private$title)
            }

            # PC1 vs PC2
            plt <- ggplot(dat, aes(pc1, pc2)) +
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
        plot_tsne = function(color=NULL, shape=NULL, title=NULL,
                             text_labels=FALSE, ...) {
            # compute t-SNE projection
            tsne <- Rtsne::Rtsne(self$dat, ...)
            dat <- setNames(as.data.frame(tsne$Y), c('x', 'y'))

            # add column ids
            dat <- cbind(dat, id = rownames(self$dat))

            # get color/shape styles
            styles <- private$get_geom_point_styles(color, shape)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color = styles$color)
            }
            if (!is.null(styles$shape)) {
                dat <- cbind(dat, shape = styles$shape)
            }

            # plot title
            if (is.null(title)) {
                title <- sprintf("t-SNE: %s", private$title)
            }

            # treatment response
            plt <- ggplot(dat, aes(x, y)) +
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
        plot_pairwise_column_cors = function(color=NULL, 
                                             label=NULL, title="",
                                             method='pearson',
                                             mar=c(12, 6, 4, 6),
                                             ...) {
            # compute pairwise variable correlations
            cor_mat <- private$similarity(self$dat, method=method, ...)
            median_pairwise_cor <- apply(cor_mat, 1, median)

            quantiles <- quantile(median_pairwise_cor, probs = c(0.25, 0.75))
            iqr <- diff(quantiles)

            #outlimit
            cutoff <- quantiles[1] - 1.5 * iqr

            ylimit <- c(pmin(min(median_pairwise_cor), cutoff),
                        max(median_pairwise_cor))

            # get color properties
            color_vector <- private$get_var_colors(color)
            label_vector <- private$get_var_labels(label)

            # variable labels
            if (!all(colnames(self$dat) == label_vector)) {
                var_labels <- sprintf("%s (%s)", colnames(self$dat), label_vector)
            } else {
                var_labels <- colnames(self$dat)
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
    )
)
