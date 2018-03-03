#' An S6 class representing a generic dataset.
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
#' @name EDADataSet
#' @export
#'
NULL

EDADataSet <- R6Class("EDADataSet",
    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # params
        dat = NULL,
        col_mdata = NULL,
        row_mdata = NULL,

        # EDADataset constructor
        initialize = function(dat, 
                              row_mdata=NULL, col_mdata=NULL, 
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
                              row_maxn=Inf, row_max_ratio=1.0, row_ind=NULL,
                              col_maxn=Inf, col_max_ratio=1.0, col_ind=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) { 

            # check to make sure row and columns identifiers are stored as
            # row and column names
            self$dat <- private$normalize_data_ids(dat, row_ids, col_ids)

            # store row and column metadata, if present
            self$col_mdata <- private$normalize_metadata_order(col_mdata, colnames(self$dat))
            self$row_mdata <- private$normalize_metadata_order(row_mdata, rownames(self$dat))
            
            # default variables to use for plot color, shape, and labels when
            # visualizing either columns or rows in the dataset
            private$row_color  <- row_color
            private$row_shape  <- row_shape
            private$row_labels <- row_labels

            private$col_color  <- col_color
            private$col_shape  <- col_shape
            private$col_labels <- col_labels
            private$color_pal  <- color_pal

            # determine subsampling indices, if requested
            private$row_ind <- private$get_subsample_indices(row_maxn, row_max_ratio, 
                                                             row_ind, nrow(dat))
            private$col_ind <- private$get_subsample_indices(col_maxn, col_max_ratio,
                                                             col_ind, ncol(dat))

            private$ggplot_theme <- ggplot_theme
            private$title  <- title
        },

        #' Clears any cached resuts and performs garbage collection to free 
        #' up memory.
        clear_cache = function() {
            private$cache <- list()
            invisible(gc())
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

        #' Applies a filter to rows of the dataset
        #'
        #' @param mask Logical vector of length equal to the number of rows in
        #'      the dataset.
        #'
        #' @return A filtered version of the original EDADataSet object.
        filter_rows = function(mask) {
            obj <- private$clone_()

            # update data and metadata matrices
            obj$dat <- obj$dat[mask,, drop=FALSE] 
            obj$row_mdata <- obj$row_mdata[mask,, drop=FALSE]

            # update subsampling indices
            rows_kept <- seq(nrow(self$dat))[mask]
            new_ind <- match(obj$get_row_indices(), rows_kept)
            new_ind <- new_ind[!is.na(new_ind)]
            obj$set_row_indices(new_ind)

            obj
        },

        #' Applies a filter to columns of the dataset
        #'
        #' @param mask Logical vector of length equal to the number of columns 
        #'     in the dataset.
        #'
        #' @return A filtered version of the original EDADataSet object.
        filter_cols = function(mask) {
            obj <- private$clone_()

            # update data and metadata matrices
            obj$dat <- obj$dat[,mask, drop=FALSE] 
            obj$col_mdata <- obj$col_mdata[mask,, drop=FALSE]

            # update subsampling indices
            cols_kept <- seq(ncol(self$dat))[mask]
            new_ind <- match(obj$get_col_indices(), cols_kept)
            new_ind <- new_ind[!is.na(new_ind)]
            obj$set_col_indices(new_ind)

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
        #' \code{obj$col_mdata}) and the t-SNE projected axes of the 
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

        #' Summarizes overall characteristics of a dataset
        #' 
        summary = function(markdown=FALSE, subsample=TRUE) {
            # collection summary info
            x <- self$dat

            # list to store summary info
            info <- list()

            # overall dataset
            info[['dat_range']] <- range(x, na.rm=TRUE)
            info[['dat_num_nas']] <- sum(is.na(x))
            info[['dat_num_zeros']] <- sum(x == 0)
            info[['dat_quartiles']] <- quantile(x, na.rm=TRUE)

            # rows
            info[['row_outliers']] <- self$detect_row_outliers()

            # columns
            info[['col_types']] <- table(sapply(x, class))
            info[['col_outliers']] <- self$detect_col_outliers()

            # row & column correlations
            if (subsample) {
                info[['col_cor_mat']] <- cor(x[private$row_ind, private$col_ind])
                info[['row_cor_mat']] <- cor(t(x[private$row_ind, private$col_ind]))
            } else {
                info[['col_cor_mat']] <- cor(x)
                info[['row_cor_mat']] <- cor(t(x))
            }
            
            # display using selected output format
            if (markdown) {
            } else {
                private$print_summary(info)
            }
        },

        ######################################################################
        # plotting methods
        ######################################################################

        #' Kernel density plot
        #'
        #' Plots densities for each row or column in the dataset. This is most
        #' useful when you are interested in similarties or differences in
        #' distributions across columns, for a relatively small number of
        #' columns.
        #'
        #' @param color Variable to color density curves by. If not is
        #' specified, uses variable specified at object construction time,
        #' or else, uses a separate color for each column.
        #'
        #' @return ggplot plot instance.
        #'
        plot_densities = function(color=NULL, title="", ...) {
            #if (target == 'rows') {
            #    private$transpose()
            #}

            # determine subsampling indices, if requested
            indices <- private$get_indices(...)

            dat <- setNames(melt(self$dat), c('row', 'column', 'val'))
            styles <- private$get_geom_density_styles(color)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Column densities: %s", private$title)
            }

            # construct density plot
            plt <- ggplot(dat, aes(x=val, group=column, color=column)) +
                geom_density(styles$aes) +
                ggtitle(title) +
                private$ggplot_theme()

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}

            # only show legend if there are relatively few groups
            if (length(unique(dat$column)) > 10) {
                plt <- plt + guides(color=FALSE)
            }

			plt
        },
        
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
                geom_point(stat="identity", styles$aes) +
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

        get_tsne_clusters = function(k=10, ...) {
            # perform t-sne and store results
            #if (is.null(private$cache[['tsne']])) {
            #    private$cache[['tsne']] <- Rtsne::Rtsne(t(self$dat), ...)
            #}
            #tsne_res <- as.data.frame(private$cache[['tsne']]$Y)

            # determine subsampling indices, if requested
            indices <- private$get_indices(...)

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
                   geom_point(styles$aes, stat="identity") +
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
        },

        #' Prints class greeting to the screen
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat("= EDADataSet\n")
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        },

        #' Gets the subsampling row indices for a dataset
        #'
        #' @return Numeric vector of row indices to use for intensive 
        #'     operations when subsampling is enabled.
        get_row_indices = function() {
            private$row_ind
        },

        #' Gets the subsampling column indices for a dataset
        #'
        #' @return Numeric vector of column indices to use for intensive 
        #'     operations when subsampling is enabled.
        get_col_indices = function() {
            private$col_ind
        },

        #' Sets the subsampling row indices for a dataset
        #'
        #' @param ind Numeric vector of row indices to use for intensive 
        #'     operations when subsampling is enabled.
        set_row_indices = function(ind) {
            private$row_ind <- ind 
        },

        #' Sets the subsampling column indices for a dataset
        #'
        #' @param ind Numeric vector of column indices to use for intensive 
        #'     operations when subsampling is enabled.
        set_col_indices = function(ind) {
            private$col_ind <- ind 
        }
    ),
    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # private params
        row_ind      = NULL,
        col_ind      = NULL,
        row_color    = NULL,
        row_shape    = NULL,
        row_labels   = NULL,
        col_color    = NULL,
        col_shape    = NULL,
        col_labels   = NULL,
        color_pal    = NULL,
        ggplot_theme = NULL,
        title        = NULL,
        cache        = list(),

        # private methods
        clone_ = function() {
            obj <- self$clone()
            obj$clear_cache()
            obj
        },

        #' Computes correlations between axes of a data projection (PCA or
        #' t-SNE) and column metadata.
        #' 
        #' @param mat Numeric projected data matrix
        #' @param include Vector of strings indicating metadata columns which
        #' should be included in the analysis.
        #'
        #' @return Dataframe of feature x project axes dependencies
        compute_feature_correlations = function(mat, include) {
            # drop any covariates with only a single level
            single_level <- apply(self$col_mdata, 2, function(x) {length(table(x))}) == 1
            features <- self$col_mdata[,!single_level, drop=FALSE]

            # drop any undesired features
            if (!is.null(include)) {
                features <- features %>%
                    select(one_of(include))
            }

            # construct linear model to measure dependence of each projected
            # axis on column metadata
            result <- matrix(0, nrow=ncol(mat), ncol=ncol(features))

            for (i in 1:ncol(features)) {
                feature_cor <- function(y) {
                    round(summary(lm(y~features[,i]))$r.squared*100, 2)
                }
                result[,i] <- apply(mat, 2, feature_cor)
            }
            colnames(result) <- colnames(features)
            result
        },

        #' Returns a list of arguments specific for a given function call.
        #'
        #' This function takes a list of function arguments and generates
        #' a new list without any shared arguments (row_maxn, etc.).
        #'
        #' @param ... Arguments passed to a given plotting, etc. function call.
        #'
        #' @return A list containing only function-specific arguments 
        get_custom_function_args = function(...) {
            args <- list(...)
      
            shared_args <- c('row_ind', 'row_maxn', 'row_max_ratio',
                             'col_ind', 'col_maxn', 'col_max_ratio')

            # return list of non-shared arguments
            args[!names(args) %in% shared_args]
        },

        get_subsample_indices = function(maxn, maxr, ind, n) {
            # if indices are explictly provided, use them
            if (!is.null(ind)) {
                return(ind)
            }

            # otherwise, if maxn or maxr, specified, use to randomly select
            # row/column indices to use for memory-/cpu-intensive operations
            maxn <- min(maxn, round(maxr * n))
            
            if (maxn < n) {
                sort(sample(n, maxn))
            } else {
                1:n
            }
        },

        #' Determines row and column indices to use for a function call
        #'
        #' Parses a list of function arguments and check whether user has
        #' specified any subsampling-related arguments (row_maxn, row_max_ratio,
        #' etc.) and, if so, select indices to use. Otherwise object-level
        #' indices are returned.
        #'
        #' @param ... Arguments passed to a given plotting, etc. function call.
        #'
        #' @return A list with 'row' and 'column' entries corresponding to the
        #'     specific numeric row and column indices to be used.
        get_indices = function(...) {
            args <- list(...)
            result <- list(row=NULL, col=NULL)

            # row indices
            if (!is.null(args$row_ind)) {
                result[['row']] <- args$row_ind
            } else if (!is.null(args$row_maxn)) {
                result[['row']] <- sample(nrow(self$dat), args$row_maxn) 
            } else if (!is.null(args$row_max_ratio) ){
                result[['row']] <- sample(nrow(self$dat), round(nrow(self$dat) * args$row_max_ratio))
            } else {
                result[['row']] <- private$row_ind
            }

            # col indices
            if (!is.null(args$col_ind)) {
                result[['col']] <- args$col_ind
            } else if (!is.null(args$col_maxn)) {
                result[['col']] <- sample(ncol(self$dat), args$col_maxn) 
            } else if (!is.null(args$col_max_ratio) ){
                result[['col']] <- sample(ncol(self$dat), round(ncol(self$dat) * args$col_max_ratio))
            } else {
                result[['col']] <- private$col_ind
            }

            result
        },

        #' Generates ggplot aesthetics for density plots 
        #'
        #' @param color Color variable as passed into plot function call
        #'
        #' @return List of style information
        get_geom_density_styles = function(color) {
            # list to store style properties
            res <- list(aes=aes(), labels=list())
            private$add_color_styles(res, color)
        },

        #' Generates ggplot aesthetics for histogram plots
        #'
        #' @param color Color variable as passed into plot function call
        #'
        #' @return List of style information
        get_geom_histogram_styles = function(color) {
            # list to store style properties
            res <- list(aes=aes(), labels=list())
            private$add_color_styles(res, color)
        },

        #' Generates ggplot aesthetics for a geom_point plot
        #'
        #' @param color Color variable as passed into plot function call
        #' @param shape Shape variable as passed into plot function call
        #'
        #' @return List of geom_point style information
        get_geom_point_styles = function(color, shape) {
            # list to store style properties
            res <- list(
                aes=aes(),
                labels=list()
            )
            res <- private$add_color_styles(res, color)
            res <- private$add_shape_styles(res, shape)

            res
        },

        #' Determines color-related style information to use for a plot
        #'
        #' @param styles List of color-related style info
        #' @param color Color variable passed into plot function call
        add_color_styles = function(styles, color) {
            color_info <- private$get_color_styles(color)

            styles[['color']] <- color_info[['color']]

			# update styles with color info
			if (length(color_info[['aes']]) > 0) {
				styles[['aes']] <- modifyList(color_info[['aes']], styles[['aes']])
				styles[['labels']] <- modifyList(color_info[['labels']], styles[['labels']])
			}
            styles
        },

        #' Determines shape-related style information to use for a plot
        #'
        #' @param styles List of shape-related style info
        #' @param shape shape variable passed into plot function call
        add_shape_styles = function(styles, shape) {
            shape_info <- private$get_shape_styles(shape)

            styles[['shape']] <- shape_info[['shape']]

			# update styles with shape info
			if (length(shape_info[['aes']]) > 0) {
				styles[['aes']] <- modifyList(shape_info[['aes']], styles[['aes']])
				styles[['labels']] <- modifyList(styles[['labels']], shape_info[['labels']])
			}
            styles
        },

        #' Returns a list of color-related style information
        #'
        #' @param color Color variable as passed into plot function call
        #' 
        #' @return list List of color-related properties
        get_color_styles = function(color) {
            res <- list(
				color  = NULL,
				aes    = aes(),
				labels = list()
			)

            # if specified as a function argument, override default color
            if (!is.null(color) && (color != FALSE)) {
                # color variable can either correspond to a column in the
                # dataset itself, or in the 
                res[['color']]  <- self$row_mdata[,color]
                res[['aes']]    <- aes(color=color)
                res[['labels']] <- labs(color=color)
            } else if (is.null(color) && !is.null(private$row_color)) {
                # otherwise, use object-level default value
                res[['color']] <- self$row_mdata[,private$row_color]
                res[['aes']] <- aes(color=color)
                res[['labels']] <- labs(color=private$row_color)
            }
            res
        },

        #' Returns a list of shape-related style information
        #'
        #' @param shape Shape variable as passed into plot function call
        #' 
        #' @return list List of shape-related properties
        get_shape_styles = function(shape) {
            res <- list(
				shape  = NULL,
				aes    = aes(),
				labels = list()
			)

            # if specified as a function argument, override default shape
            if (!is.null(shape) && (shape != FALSE)) {
                res[['shape']]  <- self$row_mdata[,shape]
                res[['aes']]    <- aes(shape=shape)
                res[['labels']] <- labs(shape=shape)
            } else if (is.null(shape) && !is.null(private$shape)) {
                # otherwise, use object-level default value
                res[['shape']] <- self$row_mdata[,private$shape]
                res[['aes']] <- aes(shape=shape)
                res[['labels']] <- labs(shape=private$shape)
            }
            
            res
        },

		#' Returns a vector of color codes associated with the specified
		#' variable.
		#'
		#' @param color Variable to use for assigning colors.
		#' @param color_pal Color palette to use
		#'
		#' @return Vector of colors with length equal to the number of columns
		#' 	       in the data.
        get_var_colors = function(color, color_pal) {
            # if no variable is specified, use default black for plots
            if (is.null(color)) {
				if (is.null(private$col_colo)) {
					return('black') 
				}
				color <- private$col_colo
            } else if (color == FALSE) {
				return('black')
			}

            # otherwise, assign colors based on the variable specified
            column_var <- as.numeric(factor(self$col_mdata[,color]))

            pal <- RColorBrewer::brewer.pal(9, color_pal)
            colors <- colorRampPalette(pal)(min(1E4, length(unique(column_var))))

            # return vector of column color assignments
            colors[as.integer(column_var)]
        },

		#' Returns a vector of text labels to use for a plot
		#'
		#' @param shape Variable to use for assigning shapes.
		#'
		#' @return Vector of labels with length equal to the number of columns
		#'         in the data.
        get_var_labels = function(label) {
            if (is.null(label)) {
                return(colnames(self$dat))
            }
            self$col_mdata[,label]
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
        },

        #' Normalizes handling of data row and column identifiers
        #'
        #' Checks dataset row and column identifiers and converts them to row 
        #' and column names, respectively, if they are not already stored there.
        #'
        #' @param dat Dataset
        #' @param row_ids Column id or number where row identifiers are stored.
        #' @param col_ids Row id or number where column identifiers are stored.
        #'
        #' @param Dataset with identifiers as rows and columns
        normalize_data_ids = function(dat, row_ids, col_ids) {
            # row ids
            if (row_ids != 'rownames') {
                # column number containing row ids specified
                if (is.numeric(row_ids)) {
                    rownames(dat) <- dat[,row_ids]
                    dat <- dat[,-row_ids]
                } else if (row_ids %in% colnames(dat)) {
                    # column name containing row ids specified
                    ind <- which(colnames(dat) == row_ids)
                    rownames(dat) <- dat[,ind]
                    dat <- dat[,-ind]
                }
            }

            # column ids
            if (col_ids != 'colnames') {
                # row number containing column ids
                if (is.numeric(col_ids)) {
                    colnames(dat) <- dat[col_ids,]
                    dat <- dat[-col_ids,]
                } else if (col_ids %in% colnames(dat)) {
                    # row name containing columns ids
                    ind <- which(rownames(dat) == col_ids)
                    colnames(dat) <- dat[ind,]
                    dat <- dat[-ind,]
                }
            }

            # if either row or column id's are missing, assign arbitrary numeric
            # identifiers; required for some plotting, etc. functionality.
            if (is.null(colnames(dat))) {
                colnames(dat) <- 1:ncol(dat)                  
            }
            if (is.null(rownames(dat))) {
                rownames(dat) <- 1:nrow(dat)                  
            }

            # return normalized dataset
            dat
        },

        # normalize dat / metadata order
        normalize_metadata_order = function(metadata, ids) {
            # If no metadata was provided, do nothing
            if (is.null(metadata)) {
                return(NULL)
            }

            # if already in the right order, return as-is
            if (all(rownames(metadata) == ids)) {
                return(metadata)
            }

            # otherwise match metadata row order to row/col names specified
            ind <- order(match(rownames(metadata), ids))

            # return result
            metadata[ind,, drop=FALSE]
        },

        #' Transpose the dataset and metadata in-place
        transpose = function() {
            self$dat <- t(self$dat) 

            row_mdata  <- self$row_mdata
            col_mdata  <- self$col_mdata
            row_color  <- private$row_color 
            row_shape  <- private$row_shape 
            row_labels <- private$row_labels 
            col_color  <- private$col_color 
            col_shape  <- private$col_shape 
            col_labels <- private$col_labels 
            row_ind    <- private$row_ind 
            col_ind    <- private$col_ind

            self$row_mdata     <- col_mdata
            self$col_mdata     <- row_mdata
            private$row_color  <- col_color 
            private$row_shape  <- col_shape 
            private$row_labels <- col_labels 
            private$col_color  <- row_color 
            private$col_shape  <- row_shape 
            private$col_labels <- row_labels 
            private$row_ind    <- col_ind 
            private$col_ind    <- row_ind
        }
    ),
    # ------------------------------------------------------------------------
    # active bindings
    # ------------------------------------------------------------------------
    active = list(
        t = function() {
            EDADataSet$new(t(self$dat), 
                           row_mdata=self$col_mdata,
                           col_mdata=self$row_mdata,
                           row_color=self$col_color, row_shape=self$col_shape, 
                           row_labels=self$col_labels, 
                           col_color=self$row_color, col_shape=self$row_shape, 
                           col_labels=self$row_labels, 
                           row_ind=self$col_ind, col_ind=self$row_ind,
                           color_pal=self$color_pal, ggplot_theme=private$ggplot_theme)
        }
    )
)

