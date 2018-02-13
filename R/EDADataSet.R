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

        # constructor
        initialize = function(dat, col_mdata=NULL, row_mdata=NULL, 
                              col_maxn=Inf, col_maxr=1.0,
                              row_maxn=Inf, row_maxr=1.0,
                              color_var=NULL, shape_var=NULL, label_var=NULL,
                              color_pal='Set1', ggplot_theme=ggplot2::theme_bw) { 
            # params
            self$dat <- dat
            self$col_mdata <- private$normalize_metadata_order(col_mdata, colnames(dat))
            self$row_mdata <- private$normalize_metadata_order(row_mdata, rownames(dat))
            
            # default variables to use for plot color, shape, and labels
            private$color_var <- color_var
            private$shape_var <- shape_var
            private$label_var <- label_var

            private$colors <- private$get_var_colors(color_var, color_pal)
            private$shapes <- private$get_var_shapes(shape_var)
            private$labels <- private$get_var_labels(label_var)

            private$ggplot_theme <- ggplot_theme

            # determine subsampling indices
            private$row_ind <- private$get_subsample_indices(row_maxn, row_maxr, nrow(dat))
            private$col_ind <- private$get_subsample_indices(col_maxn, col_maxr, ncol(dat))
        },

        clear_cache = function() {
            private$cache <- list()
            invisible(gc())
        },

        log = function(base=exp(1), offset=0) {
            obj <- private$clone_()
            obj$dat <- log(obj$dat + offset, base)
            obj
        },

        log1p = function() {
            obj <- private$clone_()
            obj$dat <- log(obj$dat + 1)
            obj
        },

        filter_rows = function(mask) {
            obj <- private$clone_()
            obj$dat <- obj$dat[mask,] 
            obj$row_mdata <- obj$row_mdata[mask,]
            obj
        },

        filter_cols = function(mask) {
            obj <- private$clone_()
            obj$dat <- obj$dat[,mask] 
            obj$col_mdata <- obj$col_mdata[mask,]
            obj
        },

        filter_col_outliers = function(num_sd=2) {
            obj <- private$clone_()
            cor_mat <- cor(obj$dat)
            median_column_cors <- apply(cor_mat, 1, median)
            cutoff <- mean(median_column_cors) - num_sd * sd(median_column_cors)
            obj$filter_cols(median_column_cors > cutoff)
        },

        filter_row_outliers = function(num_sd=2) {
            obj <- private$clone_()
            cor_mat <- cor(obj$dat)
            median_row_cors <- apply(cor_mat, 1, median)
            cutoff <- mean(median_row_cors) - num_sd * sd(median_row_cors)
            obj$filter_rows(median_row_cors > cutoff)
        },

        #' Measure the predictive power of each feature (col_mdata column)
        #' and the principle components of the dataset using a simple linear
        #' model.
        #'
        #' Based on code adapted from cbcbSEQ
        #' (https://github.com/kokrah/cbcbSEQ/) originally written by Kwame Okrah.
        #
        #' @author Kwame Okrah
        #' @author V. Keith Hughitt, \email{khughitt@umd.edu}
        #'
        ###############################################################################
        get_pca_feature_correlations = function(exclude=NULL) {
            # If already computed, return cached result
            if (!is.null(private$cache[['pca_feature_cor']])) {
                return(private$cache[['pca_feature_cor']])
            }

            # Drop any covariates with only a single level
            single_level <- apply(self$col_mdata, 2, function(x) {length(table(x))}) == 1

            features <- self$col_mdata[,!single_level]

            # drop any undesired features
            if (!is.null(exclude)) {
                features <- features %>%
                    select(-one_of(exclude))
            }

            # SVD
            s <- corpcor:::fast.svd(self$dat - rowMeans(self$dat))
            rownames(s$v) <- colnames(self$dat)

            # Create output dataframe
            pc_var <- round((s$d^2) / sum(s$d^2) * 100, 2)

            result <- data.frame(
                pc_var=pc_var,
                pc_var_cum=cumsum(pc_var) 
            )

            # Measure feature correlations
            for (i in 1:ncol(features)) {
                feature_cor <- function(y) {
                    round(summary(lm(y~features[,i]))$r.squared*100, 2)
                }
                result <- cbind(result, apply(s$v, 2, feature_cor))
            }
            colnames(result) <- c('pc_var', 'pc_var_cum', colnames(features))
            rownames(result) <- paste0("PC", 1:nrow(result)) 

            # cache result and return
            private$cache[['pca_feature_cor']] <- result

            result
        },

        #
        # plotting methods
        #
        plot_heatmap = function(...) { 
            # determine subsampling indices, if requested
            indices <- private$get_indices(...)

            # generate correlation matrix
            cor_mat <- cor(self$dat[indices$row, indices$col])

            # for heatmaps, show binary/logical variables on one side of the heatmap and
            # other variables on the other side; first column (patient_id) is excluded.
            binary_vars <- apply(self$col_mdata, 2, function(x) { length(unique(x)) }) == 2

            # additional arguments to heatmaply
            params <- list(
                x=cor_mat,
                showticklabels=c(FALSE, FALSE),
                subplot_widths=c(0.65, 0.35),
                subplot_heights=c(0.35, 0.65)
            )

            # row and column annotations
            if (ncol(self$col_mdata) == 1) {
                metadata_col <- data.frame(self$col_mdata[indices$col])
                colnames(metadata_col) <- colnames(self$col_mdata)

                params[['row_side_colors']] <- metadata_col
            } else {
                if (sum(binary_vars) > 0) {
                    params[['col_side_colors']] <- self$col_mdata[indices$col, binary_vars]
                }

                if (sum(!binary_vars) > 0) {
                    params[['row_side_colors']] <- self$col_mdata[indices$col,!binary_vars]
                }
            }

            if ('row_side_colors' %in% names(params)) {
                params[['subplot_widths']] <- c(0.55, 0.3, 0.15)
            }
            if ('col_side_colors' %in% names(params)) {
                params[['subplot_heights']] <- c(0.15, 0.3, 0.55)
            }

            do.call(heatmaply, params)
        },

        #
        # plot_pca
        #
        # Helper function for generating PCA plots
        #
        # dat  matrix Data matrix to generate PCA plot for
        # pc_x integer PC # to plot along x-axis (default: 1)
        # pc_y integer PC # to plot along x-axis (default: 2)
        #
        plot_pca = function(pc_x=1, pc_y=2, scale=FALSE, color_var=NULL, 
                            shape_var=NULL, plot_title='') {
            prcomp_results <- prcomp(t(self$dat), scale=scale)
            var_explained <- round(summary(prcomp_results)$importance[2,] * 100, 2)

            xl <- sprintf("PC%d (%.2f%% variance)", pc_x, var_explained[pc_x])
            yl <- sprintf("PC%d (%.2f%% variance)", pc_y, var_explained[pc_y])

            df <- data.frame(id=colnames(self$dat), 
                             pc1=prcomp_results$x[,pc_x],
                             pc2=prcomp_results$x[,pc_y])

            # color (optional)
            if (!is.null(color_var)) {
                df <- cbind(df, color_var=self$col_mdata[,color_var])
                plot_aes    <- aes(color=color_var)
                plot_labels <- labs(col=color_var)
            } else if (!is.null(private$color_var)) {
                df <- cbind(df, color_var=self$col_mdata[,private$color_var])
                plot_aes    <- aes(color=color_var)
                plot_labels <- labs(col=private$color_var)
            } else {
                plot_aes    <- aes() 
                plot_labels <- list()
            }

            # shape (optional)
            if (!is.null(shape_var)) {
                df <- cbind(df, shape_var=self$col_mdata[,shape_var])
                plot_aes    <- modifyList(plot_aes, aes(color=shape_var))
                plot_labels <- modifyList(plot_labels, labs(col=shape_var))
            } else if (!is.null(private$shape_var)) {
                df <- cbind(df, shape_var=self$col_mdata[,private$shape_var])
                plot_aes    <- modifyList(plot_aes, aes(color=shape_var))
                plot_labels <- modifyList(plot_labels, labs(col=private$shape_var))
            }

            # PC1 vs PC2
            ggplot2::ggplot(df, aes(pc1, pc2)) +
                geom_point(stat="identity", size=1, plot_aes) +
                geom_text(aes(label=id), angle=45, size=1, vjust=2) +
                xlab(xl) + ylab(yl) +
                ggtitle(plot_title) +
                private$ggplot_theme() +
                plot_labels +
                theme(axis.ticks=element_blank(), 
                      axis.text.x=element_text(angle=-90))
        },

        plot_pca_feature_correlations = function(num_pcs=6, ...) {
            # compute pca feature correlations or retrieved cached version
            if (!is.null(private$cache[['pca_feature_cor']])) {
                pca_cor <- private$cache[['pca_feature_cor']]
            } else {
                pca_cor <- self$get_pca_feature_correlations(...)
            }

            pca_cor <- pca_cor[1:num_pcs,]
            pca_long <- cbind(PC=sprintf("PC%d (%0.2f%%)", 
                                         1:nrow(pca_cor), pca_cor$pc_var), 
                              reshape2::melt(pca_cor, id.vars=c('pc_var', 'pc_var_cum')))

            ggplot(pca_long, aes(x=PC, y=variable)) + 
                geom_tile(aes(fill=value)) + 
                geom_text(aes(label=value), size=2, show.legend=FALSE) +
                scale_fill_gradient(low="green", high="red") +
                private$ggplot_theme() +
                theme(axis.text.x=element_text(size=8, angle=45, vjust=1, hjust=1), 
                      axis.text.y=element_text(size=8)) + 
                xlab("Principle Components (% var explained)") +
                ylab("Features") +
                guides(fill=guide_legend("R^2"))
        },

        get_tsne_clusters = function(k=10, ...) {
            # perform t-sne and store results
            if (is.null(private$cache[['tsne']])) {
                private$cache[['tsne']] <- Rtsne::Rtsne(t(self$dat), ...)
            }
            tsne_res <- as.data.frame(private$cache[['tsne']]$Y)
            colnames(tsne_res) <- c('x', 'y')

            # Cluster patients from t-sne results
            kmeans_clusters <- kmeans(tsne_res, k)$cluster
            factor(paste0('cluster_', kmeans_clusters))
        },

        plot_tsne = function(plot_title='t-SNE', color_var=NULL, shape_var=NULL, ...) {
            # perform t-sne and store results
            if (is.null(private$cache[['tsne']])) {
                private$cache[['tsne']] <- Rtsne::Rtsne(t(self$dat), ...)
            }

            tsne_res <- as.data.frame(private$cache[['tsne']]$Y)
            colnames(tsne_res) <- c('x', 'y')

            # add color and shape info
            tsne_res <- cbind(tsne_res,
                              color_var=private$get_plot_color_column(color_var),
                              shape_var=private$get_plot_color_column(shape_var))

            plot_aes <- private$get_plot_aes(color_var, shape_var)
            plot_labs <- private$get_plot_legend_labels(color_var, shape_var)

            # treatment response
            ggplot2::ggplot(tsne_res, aes(x, y)) +
                   geom_point(plot_aes, stat="identity", size=1) +
                   ggtitle(plot_title) +
                   plot_labs +
                   private$ggplot_theme() +
                   theme(axis.ticks=element_blank(), 
                         axis.text.x=element_text(angle=-90),
                         legend.text=element_text(size=8))
        },

        #' Plot median pairwise variable correlations
        #'
        #' Plots the median correlation of each variable (column)
        #'
        #' @author V. Keith Hughitt, \email{keith.hughitt@nih.gov}
        #'
        #' @return None
        plot_var_correlations = function (main="", method='pearson',
                                          mar=c(12,6,4,6)) {
            # compute pairwise variable correlations
            median_pairwise_cor <- apply(cor(self$dat * 1.0, method=method), 1, median)

            quantiles <- quantile(median_pairwise_cor, probs=c(0.25, 0.75))
            iqr <- diff(quantiles)

            #outlimit
            cutoff <- quantiles[1] - 1.5 * iqr

            ylimit <- c(pmin(min(median_pairwise_cor), cutoff), 
                        max(median_pairwise_cor))

            # variable labels
            if (!all(colnames(self$dat) == private$labels)) {
                var_labels <- sprintf("%s (%s)", colnames(self$dat), private$labels)
            } else {
                var_labels <- colnames(self$dat)
            }

            # render plot
            par(mar=mar)
            plot(median_pairwise_cor, xaxt="n", ylim=ylimit,
                 ylab="Median Pairwise Correlation", xlab="", main=main,
                 col=private$colors, pch=16, cex=2.2)
            axis(side=1, at=seq(along=median_pairwise_cor),
                 labels=var_labels, las=2)
            abline(h=cutoff, lty=2)
            abline(v=1:length(var_labels), lty=3, col="black")
        },

        # class greeting
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat("= EDADataSet\n")
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        }
    ),
    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # private params
        row_ind      = NULL,
        col_ind      = NULL,
        color_var    = NULL,
        shape_var    = NULL,
        label_var    = NULL,
        colors       = NULL,
        shapes       = NULL,
        labels       = NULL,
        ggplot_theme = NULL,
        cache        = list(),

        # private methods
        clone_ = function() {
            obj <- self$clone()
            obj$clear_cache()
            obj
        },

        get_subsample_indices = function(maxn, maxr, n) {
            maxn <- min(maxn, round(maxr * n))
            
            if (maxn < n) {
                sort(sample(n, maxn))
            } else {
                1:n
            }
        },

        get_indices = function(...) {
            args <- list(...)
            result <- list(row=NULL, col=NULL)

            # row indices
            if (!is.null(args$row_maxn)) {
                result[['row']] <- sample(nrow(self$dat), args$row_maxn) 
            } else if (!is.null(args$row_maxr) ){
                result[['row']] <- sample(nrow(self$dat), round(nrow(self$dat) * args$row_maxr))
            } else {
                result[['row']] <- private$row_ind
            }

            # col indices
            if (!is.null(args$col_maxn)) {
                result[['col']] <- sample(ncol(self$dat), args$col_maxn) 
            } else if (!is.null(args$col_maxr) ){
                result[['col']] <- sample(ncol(self$dat), round(ncol(self$dat) * args$col_maxr))
            } else {
                result[['col']] <- private$col_ind
            }

            result
        },

        get_plot_color_column = function(color_var) {
            if (!is.null(color_var)) {
                factor(as.character(self$col_mdata[,color_var]))
            } else if (!is.null(private$color_var)) {
                factor(as.character(self$col_mdata[,private$color_var]))
            }
        },

        get_plot_shape_column = function(shape_var) {
            if (!is.null(shape_var)) {
                factor(as.character(self$col_mdata[,shape_var]))
            } else if (!is.null(private$shape_var)) {
                factor(as.character(self$col_mdata[,private$shape_var]))
            }
        },

        get_plot_color_aes = function(color_var=NULL) {
            if (!is.null(color_var)) {
                aes(color=color_var)
            } else if (!is.null(private$color_var)) {
                aes(color=color_var)
            } else {
                aes() 
            }
        },

        get_plot_shape_aes = function(shape_var=NULL) {
            if (!is.null(shape_var)) {
                aes(shape=shape_var)
            } else if (!is.null(private$shape_var)) {
                aes(shape=shape_var)
            } else {
                aes() 
            }
        },

        # determine ggplot2 color and shape aesthetics to use
        get_plot_aes = function(color_var=NULL, shape_var=NULL) {
            color_aes <- private$get_plot_color_aes(color_var)
            shape_aes <- private$get_plot_shape_aes(shape_var)

            # combine color and shape aesthetics and return
            modifyList(color_aes, shape_aes)
        },

        get_plot_color_legend_label = function(color_var=NULL) {
            if (!is.null(color_var)) {
                labs(col=color_var)
            } else if (!is.null(private$color_var)) {
                labs(col=private$color_var)
            } else {
                list()
            }
        },

        get_plot_shape_legend_label = function(shape_var=NULL) {
            if (!is.null(shape_var)) {
                labs(shape=shape_var)
            } else if (!is.null(private$shape_var)) {
                labs(shape=private$shape_var)
            } else {
                list()
            }
        },

        # determine ggplot2 color and shape aesthetics to use
        get_plot_legend_labels = function(color_var=NULL, shape_var=NULL) {
            color_label <- private$get_plot_color_legend_label(color_var)
            shape_label <- private$get_plot_shape_legend_label(shape_var)

            modifyList(color_label, shape_label)
        },

        # expects a categorical column variable
        get_var_colors = function(color_var, color_pal) {
            # if no variable is specified, use default black for plots
            if (is.null(color_var)) {
                return('black') 
            }

            # get column variable to use for assigning colors
            column_var <- as.numeric(factor(self$col_mdata[,color_var]))
    
            # otherwise, assign colors based on the variable specified
            pal <- RColorBrewer::brewer.pal(9, color_pal)
            colors <- colorRampPalette(pal)(min(1E4, length(unique(column_var))))

            # return vector of column color assignments
            colors[as.integer(column_var)]
        },

        get_var_shapes = function(shape_var) {
            if (is.null(shape_var)) {
                return(NULL)
            }
            as.numeric(factor(self$col_mdata[,shape_var]))
        },

        get_var_labels = function(label_var) {
            if (is.null(label_var)) {
                return(colnames(self$dat))
            }
            self$col_mdata[,label_var]
        },

        # normalize dat / metadata order
        normalize_metadata_order = function(metadata, to_match) {
            # If no metadata was provided, do nothing
            if (is.null(metadata)) {
                return(NULL)
            }

            # if already in the right order, return as-is
            if (all(rownames(metadata) == to_match)) {
                return(metadata)
            }

            # otherwise match metadata row order to row/col names specified
            ind <- order(match(rownames(metadata), to_match))

            # multi-column metadata
            if (ncol(metadata) > 1) {
                return(metadata[ind,])
            }

            # single-column metadata (need to recreate matrix after reordering)
            rnames <- rownames(metadata)
            cnames <- colnames(metadata)

            metadata <- matrix(metadata[ind,])
            rownames(metadata) <- rnames[ind]
            colnames(metadata) <- cnames

            metadata
        }
    ),
    # ------------------------------------------------------------------------
    # active bindings
    # ------------------------------------------------------------------------
    active = list(
        t = function() {
            EDADataSet$new(t(self$dat), self$col_mdata, self$row_mdata)
        }
    )
)

