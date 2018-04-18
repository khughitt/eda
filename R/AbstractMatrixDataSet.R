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
#'  - `cross_cor(key1=1, key2=2, meas='pearson')`: Computes cross-dataset
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
        detect_col_outliers = function(key=1, num_sd=2, ctend=median, meas='pearson', ...) {
            # TODO: include correlation in results?
            # TODO: Write alternative version for data frame datasets?
            dat <- self$edat[[key]]$dat
            cor_mat <- private$similarity(dat, meas = meas, ...)

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
        detect_row_outliers = function(key=1, num_sd=2, ctend=median, meas='pearson', ...) {
            dat <- self$edat[[key]]$dat
            cor_mat <- private$similarity(t(dat), meas = meas, ...)

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
        filter_col_outliers = function(key=1, num_sd=2, ctend=median, meas='pearson') {
            obj <- private$clone_()
            outliers <- obj$detect_col_outliers(key, num_sd, avg, meas)
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
        filter_row_outliers = function(num_sd=2, ctend=median, meas='pearson') {
            obj <- private$clone_()
            outliers <- obj$detect_row_outliers(num_sd, ctend, meas)
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

        # NMF
        nmf = function(key=1, rank, ...) {
            res <- NMF::nmf(self$edat[[key]]$dat, rank = rank, ...)@fit@W

            # clone object and append result
            obj <- private$clone_()
            obj$edat[[key]]$dat <- res
            obj$edat[[key]]$ylab <- 'NMF factors'

            obj$remove_unlinked(key)
            obj
        },

        # PCA
        #
        # Methods:
        #
        # 1. prcomp - regular PCA
        #
        pca = function(key=1, method='prcomp', ...) {
            dat <- self$edat[[key]]$dat

            # compute PCA
            if (method == 'prcomp') {
                # regular PCA
                res <- prcomp(dat, ...)
                pca_dat <- res$x

                # update colnames to include variance explained
                var_explained <- round(summary(res)$importance[2, ] * 100, 2)

                colnames(pca_dat) <- sprintf("PC%d (%.2f%% variance)", 
                                             1:ncol(pca_dat),
                                             var_explained[1:ncol(pca_dat)])
            } else if (method == 'kpca') {
                # kernal PCA
                pca_dat <- kernlab::pcv(kernlab::kpca(dat, ...))
            } else if (method == 'nsprcomp') {
                # sparse PCA
                res <- nsprcomp::nsprcomp(dat, ...)

                pca_dat <- res$x

                # update colnames to include variance explained
                var_explained <- round(summary(res)$importance[2, ] * 100, 2)

                colnames(pca_dat) <- sprintf("PC%d (%.2f%% variance)", 
                                             1:ncol(pca_dat),
                                             var_explained[1:ncol(pca_dat)])
            } else if (method == 'robpca') {
                # robust PCA
                pca_dat <- rospca::robpca(t(dat), ...)$loadings
            } else if (method == 'rospca') {
                # robust sparse PCA
                pca_dat <- rospca::rospca(t(dat), ...)$loadings
            }

            # clone object and append result
            obj <- private$clone_()
            obj$edat[[key]]$dat <- pca_dat
            obj$edat[[key]]$ylab <- 'Priniple Components'

            obj$remove_unlinked(key)
            obj
        },

        #
        # t-SNE
        #
        tsne = function(key=1, ...) {
            dat <- Rtsne::Rtsne(self$edat[[key]]$dat, ...)$Y
            colnames(dat) <- sprintf("t-SNE Dim %d", 1:ncol(dat))

            obj <- private$clone_()
            obj$edat[[key]]$dat <- dat
            obj$edat[[key]]$ylab <- 't-SNE dimensions'

            obj$remove_unlinked(key)
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
        plot_cor_heatmap = function(key=1, meas='pearson', interactive=TRUE, ...) {
            # generate correlation matrix
            cor_mat <- private$similarity(self$edat[[key]]$dat, meas=meas, ...)

            # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # add any additional function arguments
            params <- c(params, ...)

            private$construct_heatmap_plot(params, interactive)
        },

        # Generates a heatmap plot of the dataset
        #
        # ... Additional arguments
        plot_heatmap = function(key=1,
                                row_label=NULL, row_edat=NULL,
                                col_label=NULL, col_edat=NULL,
                                interactive=TRUE, ...) {

            # determine row and column labels to use
            row_labels <- private$get_row_labels(key, label_var = row_label, label_key = row_edat)
            col_labels <- private$get_col_labels(key, label_var = col_label, label_key = col_edat)

            dat <- self$edat[[key]]$dat

            rownames(dat) <- row_labels
            colnames(dat) <- col_labels

            # list of parameters to pass to heatmaply
            params <- list(
                x               = dat,
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

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
        plot_pca = function(key=1, method='prcomp', pcx=1, pcy=2,
                            color_var=NULL, color_key=NULL,
                            shape_var=NULL, shape_key=NULL,
                            label_var=NULL, label_key=NULL,
                            title=NULL, text_labels=FALSE, ...) {
            # plot title
            if (is.null(title)) {
                key <- ifelse(is.numeric(key), names(self$edat)[key], key)

                if (private$title != '') {
                    title <- sprintf("PCA: %s (%s)", key, private$title)
                } else {
                    title <- sprintf("PCA: %s", key)
                }
            }

            self$pca(key = key, method = method, ...)$scatter_plot(
                key = 1, x = pcx, y = pcy,
                color_var = color_var, color_key = color_key, 
                shape_var = shape_var, shape_key = shape_key,
                label_var = label_var, label_key = label_key, 
                title = title, text_labels = text_labels)
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
        plot_tsne = function(key=1, dim1=1, dim2=2,
                             color_var=NULL, color_key=NULL,
                             shape_var=NULL, shape_key=NULL,
                             label_var=NULL, label_key=NULL,
                             title=NULL, text_labels=FALSE, ...) {
            # plot title
            if (is.null(title)) {
                key <- ifelse(is.numeric(key), names(self$edat)[key], key)

                if (private$title != '') {
                    title <- sprintf("t-SNE: %s (%s)", key, private$title)
                } else {
                    title <- sprintf("t-SNE: %s", key)
                }
            }

            self$tsne(key = key, ...)$scatter_plot(
                key = 1, x = dim1, y = dim2,
                color_var = color_var, color_key = color_key, 
                shape_var = shape_var, shape_key = shape_key,
                label_var = label_var, label_key = label_key, 
                title = title, text_labels = text_labels)
        },

        # creates a scatter plot of from two columns in a specified dataset
        scatter_plot = function(key=1, x=1, y=2, 
                                color_var=NULL, color_key=NULL,
                                shape_var=NULL, shape_key=NULL, 
                                label_var=NULL, label_key=NULL, 
                                title=NULL, text_labels=FALSE) {

            dat <- as.data.frame(self$edat[[key]]$dat[, c(x, y)])

            # get color/shape styles
            styles <- private$get_geom_point_styles(key,
                                                    color_var, color_key,
                                                    shape_var, shape_key)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color_var = styles$color)
            }
            if (!is.null(styles$shape)) {
                dat <- cbind(dat, shape_var = styles$shape)
            }

            # generate scatter plot
            xid <- sprintf("`%s`", colnames(dat)[1])
            yid <- sprintf("`%s`", colnames(dat)[2])

            plt <- ggplot(dat, aes_string(xid, yid)) +
                geom_point(stat = "identity", styles$aes, size = 0.5) +
                ggtitle(title) +
                private$ggplot_theme() +
                theme(axis.ticks = element_blank(),
                      axis.text.x = element_text(angle = -90))

            # text labels
            if (text_labels) {
                plt <- plt + geom_text(aes_q(label = rownames(dat)), angle = 45, size = 1, vjust = 2)
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
                                             title="", meas='pearson',
                                             mar=c(12, 6, 4, 6), ...) {
            dat <- self$edat[[key]]$dat

            # compute pairwise variable correlations
            cor_mat <- private$similarity(dat, meas = meas, ...)
            median_pairwise_cor <- apply(cor_mat, 1, median)

            quantiles <- quantile(median_pairwise_cor, probs = c(0.25, 0.75))
            iqr <- diff(quantiles)

            #outlimit
            cutoff <- quantiles[1] - 1.5 * iqr

            ylimit <- c(pmin(min(median_pairwise_cor), cutoff),
                        max(median_pairwise_cor))

            # determine colors and labels to use
            color_vector <- private$get_row_colors(key, color_var, color_key)
            label_vector <- private$get_row_labels(key, label_var, label_key)

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
        },

        # computes matrix of row z-scores
        zscores = function(key=1) {
            # compute z-scores
            dat <- t(scale(t(self$edat[[key]]$dat)))
            attr(dat, 'scaled:scale') <- NULL
            attr(dat, 'scaled:center') <- NULL

            # clone object and add new dataset
            obj <- private$clone_()

            # map numeric keys to string names
            key <- ifelse(is.numeric(key), names(self$edat)[key], key)

            # determine key to use for storing result
            new_key <- sprintf('%s_zscores', key)

            # add new matrix to front of edat list and return
            obj$add(EDADat$new(dat, xid = self$edat[[key]]$xid, yid = self$edat[[key]]$yid), new_key)

            obj
        }
    ),
)
