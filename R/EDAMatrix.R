#' An S6 class representing a matrix dataset.
#'
#' EDAMatrix is a helper class for wrapping data matrices, with optional
#' support for row and column datadata. Methods are provided for common
#' exploratory data analysis summary statistics, transformations, and
#' visualizations.
#'
#' @section Usage:
#' \preformatted{
#' dat <- as.matrix(iris[,1:4])
#' row_mdata <- iris[,5,drop=FALSE]
#'
#' edm <- EDAMatrix$new(dat, row_mdata=row_mdata, row_color='Species')
#'
#' edm
#' edm$summary()
#'
#' edm$plot_pca()
#' edm$log1p()$plot_cor_heatmap()
#' edm$subsample(100)$plot_tsne()
#' }
#'
#' @section Arguments:
#'  - **dat** An m x n dataset.
#'  - **row_mdata** A matrix or data frame with rows corresponding to the row 
#'      names of `dat`
#'  - **col_mdata** A matrix or data frame with rows corresponding to the 
#'      column names of `dat`
#'  - **row_ids** Column name or number containing row identifiers. If set to
#'      `rownames` (default), row names will be used as identifiers.
#'  - **col_ids** Column name or number containing column identifiers. If set to
#'      `colnames` (default), column names will be used as identifiers.
#'  - **row_mdata_ids** Column name or number containing row metadata row 
#'      identifiers. If set to `rownames` (default), row names will be used 
#'      as identifiers.
#'  - **col_mdata_ids** Column name or number containing col metadata row 
#'      identifiers. If set to `rownames` (default), row names will be used 
#'      as identifiers.
#'  - **row_color** Row metadata field to use for coloring rowwise plot elements.
#'  - **row_shape** Row metadata field to use for determine rowwise plot 
#'      element shape.
#'  - **row_labels** Row metadata field to use when labeling plot points or
#'      other elements.
#'  - **col_color** Column metadata field to use for coloring columnwise plot elements.
#'  - **col_shape** Column metadata field to use for determine columnwise plot 
#'      element shape.
#'  - **col_labels** Column metadata field to use when labeling plot points or
#'      other elements.
#'  - **color_pal** Color palette to use for relevant plotting methods 
#'      (default: `Set1`).
#'  - **title** Text to use as a title or subtitle for plots.
#'  - **ggplot_theme** Default theme to use for ggplot2 plots 
#'      (default: `theme_bw`).
#'
#' @section Fields:
#'  - **dat** Underlying data matrix
#'  - **row_mdata** Dataframe containing row metadata
#'  - **col_mdata** Dataframe containing column metadata
#'
#' @section Methods:
#'  - **clear_cache()** Clears EDAMatrix cache.
#'  - **clone()** Creates a copy of the EDAMatrix instance.
#'  - **cluster_tsne(k=10, ...)** Clusters rows in dataset using a combination
#'      of t-SNE and k-means clustering.
#' - **detect_col_outliers(num_sd=2, avg='median', sim_method='pearson')**
#'      Measures average pairwise similarities between all columns in the dataset.
#'      Outliers are considered to be those columns who mean similarity to
#'      all other columns is greater than `num_sd` standard deviations from the
#'      average of averages.
#' - **detect_row_outliers(num_sd=2, avg='median', sim_method='pearson')**
#'      Measures average pairwise similarities between all rows in the dataset.
#'      Outliers are considered to be those rows who mean similarity to
#'      all other rows is greater than `num_sd` standard deviations from the
#'      average of averages.
#'  - **feature_cor()** Detects dependencies between column metadata entries 
#'		(features) and dataset rows. 
#'  - **filter_col_outliers(num_sd=2, avg='median', sim_method='pearson')** 
#'		Removes column outliers from the dataset. See `detect_col_outliers()` 
#'		for details of outlier detection approach.
#'  - **filter_row_outliers(num_sd=2, avg='median', sim_method='pearson')** 
#'		Removes row outliers from the dataset. See `detect_row_outliers()` 
#'		for details of outlier detection approach.
#'  - **filter_cols(mask)** Accepts a logical vector of length `ncol(obj$dat)`
#'		and returns a new EDAMatrix instance with only the columns associated
#'      with `TRUE` values in the mask.
#'  - **filter_rows(mask)** Accepts a logical vector of length `nrow(obj$dat)`
#'		and returns a new EDAMatrix instance with only the rowsumns associated
#'      with `TRUE` values in the mask.
#'  - **impute(method='knn')** Imputes missing values in the dataset and stores
#'		the result _in-place_. Currently only k-Nearest Neighbors (kNN) 
#'		imputation is supported.
#'  - **log(base=exp(1), offset=0)** Log-transforms data.
#'  - **log1p()** Log(x + 1)-transforms data.
#'  - **pca(...)** Performs principle component analysis (PCA) on the dataset
#'		and returns a new EDAMatrix instance of the projected data points.
#'      Any additional arguements specified are passed to the `prcomp()` function.
#'  - **pca_feature_cor(method='pearson', ...)** Measures correlation between
#'		dataset features (column metadata fields) and dataset principle
#'      components.
#'  - **plot_cor_heatmap(method='pearson', interactive=TRUE, ...)** Plots a 
#'		correlation heatmap of the dataset.
#'  - **plot_densities(color=NULL, title="", ...)** Plots densities for each 
#'		column in the dataset.
#'  - **plot_feature_cor(method='pearson', color_scale=c('green', 'red'))**
#'		Creates a tile plot of projected data / feature correlations. See
#'		`feature_cor()` function.
#'  - **plot_heatmap(interactive=TRUE, ...)** Generates a heatmap plot of the 
#'		dataset
#'  - **plot_pairwise_column_cors(color=NULL, title="", method='pearson', mar=c(12,6,4,6))**
#'		Plot median pairwise column correlations for each variable (column)
#'		in the dataset.
#'  - **plot_pca(pcx=1, pcy=2, scale=FALSE, color=NULL, shape=NULL, title=NULL,
#'               text_labels=FALSE, ...)**
#'		Generates a two-dimensional PCA plot from the dataset.
#'  - **plot_tsne(color=NULL, shape=NULL, title=NULL, text_labels=FALSE, ...)**
#'		Generates a two-dimensional t-SNE plot from the dataset.
#'  - **print()** Prints an overview of the object instance.
#'  - **subsample(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL)**
#'		Subsamples dataset rows and/or columns.
#'  - **summary(markdown=FALSE, num_digits=2)** Summarizes overall 
#'		characteristics of a dataset.
#'  - **t()** Transposes dataset rows and columns.
#'  - **tsne(...)** Performs T-distributed stochastic neighbor embedding (t-SNE) 
#'		on the dataset and returns a new EDAMatrix instance of the projected 
#' 		data points. Any additional arguements specified are passed to the 
#'		`Rtsne()` function.
#'  - **tsne_feature_cor(method='pearson', ...)** Measures correlation between
#'		dataset features (column metadata fields) and dataset t-SNE projected
#'      axes.
#'
#' @importFrom R6 R6Class
#' @format \code{\link{R6Class}} object.
#' @export
#' @name EDAMatrix
#'
NULL

EDAMatrix <- R6::R6Class("EDAMatrix",
    inherit = EDADataSet,
    public = list(
        # EDAMatrix constructor
        initialize = function(dat,
                              row_mdata=NULL, col_mdata=NULL,
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) {
            # verify input data type and call parent constructor
            private$check_input(dat)

            super$initialize(dat, row_mdata, col_mdata, row_ids, col_ids,
                             row_mdata_ids, col_mdata_ids, row_color, row_shape,
                             row_labels, col_color, col_shape, col_labels,
                             color_pal, title, ggplot_theme)
        },

        # Clusters dataset rows using k-means clustering in a t-SNE projected
        # space.
        #
        # k Number of clusters to detect (default: 10)
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

        # Detects column outliers in the dataset
        #
        # Computes pairwise correlations between all columns in the dataset.
        # Columns whose average pairwise correlation with all other columns
        # is greater than `num_sd` standard deviations from the average
        # average correlation are considered to be outliers.
        #
        # num_sd Number of standard deviations to use to determine
        #      outliers.
        #
        # return Character vector or column ids for columns with low
        #     average pairwise correlations.
        detect_col_outliers = function(num_sd=2, avg=median, method='pearson') {
            # TODO: include correlation in results?
            # TODO: Write alternative version for data frame datasets?
            cor_mat <- cor(self$dat, method=method)
            avg_column_cors <- apply(cor_mat, 1, avg)
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
        detect_row_outliers = function(num_sd=2, avg=median, method='pearson') {
            cor_mat <- cor(t(self$dat), method=method)
            avg_row_cors <- apply(cor_mat, 1, avg)
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
        # return A filtered version of the original EDADataSet object.
        filter_col_outliers = function(num_sd=2, avg=median, method='pearson') {
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
        # return A filtered version of the original EDADataSet object.
        filter_row_outliers = function(num_sd=2) {
            obj <- private$clone_()
			outliers <- obj$detect_row_outliers(num_sd, avg, method)
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
        pca = function (...) {
            obj <- private$clone_()
            obj$dat <- prcomp(self$dat, ...)$x
            obj$col_mdata <- NULL
            obj
        },

        #
        pca_feature_cor = function(method='pearson', ...) {
            self$t$pca(...)$t$feature_cor(method)
        },

        # t-SNE
        #
        tsne = function (...) {
            obj <- private$clone_()
            obj$dat <- Rtsne::Rtsne(obj$dat, ...)
            obj$col_mdata <- NULL
            obj
        },

        tsne_feature_cor = function(method='pearson', ...) {
            self$t$tsne(...)$t$feature_cor(...)
        },

        # Prints an overview of the object instance
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
            cor_mat <- cor(t(self$dat), method=method)

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
                    params[['col_side_colors']] <- self$row_mdata[,binary_vars, drop=FALSE]
                    params[['subplot_heights']] <- c(0.15, 0.3, 0.55)
                }

                # row colors (everything else)
                if (sum(!binary_vars) >= 1) {
                    params[['row_side_colors']] <- self$row_mdata[,!binary_vars, drop=FALSE]
                    params[['subplot_widths']] <- c(0.55, 0.3, 0.15)
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
                x=self$dat,
                showticklabels=c(FALSE, FALSE),
                subplot_widths=c(0.65, 0.35),
                subplot_heights=c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$row_mdata)) {
                params[['row_side_colors']] <- self$row_mdata
                params[['subplot_widths']] <- c(0.15, 0.3, 0.55)
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
            dat <- melt(private$feature_cor(method=method))
            colnames(dat) <- c('dim', 'variable', 'value')

            # Labels
            #if (method == 'pca') {
            #    xlab_text <- 'Principle Components'
            #} else if (method == 't-sne') {
            #    xlab_text <- "t-SNE dimension"
            #}
            xlab_text <- 'TODO...'

            # create plot
            ggplot(dat, aes(x=dim, y=variable)) +
                geom_tile(aes(fill=value)) +
                geom_text(aes(label=value), size=2, show.legend=FALSE) +
                scale_fill_gradient(low=color_scale[1], high=color_scale[2]) +
                private$ggplot_theme() +
                theme(axis.text.x=element_text(size=8, angle=45, vjust=1, hjust=1),
                      axis.text.y=element_text(size=8)) +
                xlab(xlab_text) +
                ylab("Features") +
                guides(fill=guide_legend("R^2"))
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
            prcomp_results <- prcomp(self$dat, scale=scale)
            var_explained <- round(summary(prcomp_results)$importance[2,] * 100, 2)

            # create data frame for plotting
            dat <- data.frame(id=rownames(self$dat),
                              pc1=prcomp_results$x[,pcx],
                              pc2=prcomp_results$x[,pcy])

            # get color/shape styles
            styles <- private$get_geom_point_styles(color, shape)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (!is.null(styles$shape)) {
                dat <- cbind(dat, shape=styles$shape)
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

        # Plot median pairwise column correlations
        #
        # Plots the median correlation of each variable (column). This is
        # useful for visually inspecting columns for possible outliers, when
        # the total number of columns is relatively small.
        #
        # @author V. Keith Hughitt, \email{keith.hughitt@nih.gov}
        #
        # return None
        plot_pairwise_column_cors = function (color=NULL, title="",
                                              method='pearson',
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
        # Verifies that input data is of type matrix
        check_input = function(dat) {
            if(!is.matrix(dat)) {
                stop("Invalid input for EDAMatrix: dat must be a matrix.")
            }
        }
    )
)

