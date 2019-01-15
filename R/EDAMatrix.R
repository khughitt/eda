#' An R6 class representing a matrix dataset.
#'
#' EDAMatrix is a helper class for wrapping data matrices, with optional
#' support for row and column datadata. Methods are provided for common
#' exploratory data analysis summary statistics, transformations, and
#' visualizations.
#'
#' @section Usage:
#' ```
#' edm <- EDAMatrix$new(mat, row_mdata=row_mdata_df, row_color='some_var')
#' edm$summary()
#'
#' edm$plot_pca()
#'
#' edm$t$subsample(100)$plot_heatmap()
#' ```
#'
#' @section Arguments:
#' - `dat`: An m x n dataset.
#' - `row_mdata`: A matrix or data frame with rows corresponding to the row
#'      names of `dat`
#' - `col_mdata`: A matrix or data frame with rows corresponding to the
#'      column names of `dat`
#' - `row_names`: Column name or number containing row identifiers. If set to
#'      `rownames` (default), row names will be used as identifiers.
#' - `col_names`: Column name or number containing column identifiers. If set to
#'      `colnames` (default), column names will be used as identifiers.
#' - `row_mdata_row_names`: Column name or number containing row metadata row
#'      identifiers. If set to `rownames` (default), row names will be used
#'      as identifiers.
#' - `col_mdata_row_names`: Column name or number containing col metadata row
#'      identifiers. If set to `rownames` (default), row names will be used
#'      as identifiers.
#' - `row_color`: Row metadata field to use for coloring rowwise plot elements.
#' - `row_shape`: Row metadata field to use for determine rowwise plot
#'      element shape.
#' - `row_label`: Row metadata field to use when labeling plot points or
#'      other elements.
#' - `col_color`: Column metadata field to use for coloring columnwise plot elements.
#' - `col_shape`: Column metadata field to use for determine columnwise plot
#'      element shape.
#' - `col_label`: Column metadata field to use when labeling plot points or
#'      other elements.
#' - `color_pal`: Color palette to use for relevant plotting methods
#'      (default: `Set1`).
#' - `title`: Text to use as a title or subtitle for plots.
#' - `ggplot_theme`: Default theme to use for ggplot2 plots
#'      (default: `theme_bw`).
#'
#' @section Fields:
#'  - `dat`: Underlying data matrix
#'  - `row_mdata`: Dataframe containing row metadata
#'  - `col_mdata`: Dataframe containing column metadata
#'
#' @section Methods:
#'  - `clear_cache()`: Clears EDAMatrix cache.
#'  - `clone()`: Creates a copy of the EDAMatrix instance.
#'  - `cluster_tsne(num_clusters=10, ...)`: Clusters rows in dataset using a combination
#'      of t-SNE and k-means clustering.
#' - `detect_col_outliers(num_sd=2, ctend='median', meas='pearson')`:
#'      Measures average pairwise similarities between all columns in the dataset.
#'      Outliers are considered to be those columns who mean similarity to
#'      all other columns is greater than `num_sd` standard deviations from the
#'      average of averages.
#' - `detect_row_outliers(num_sd=2, ctend='median', meas='pearson')`:
#'      Measures average pairwise similarities between all rows in the dataset.
#'      Outliers are considered to be those rows who mean similarity to
#'      all other rows is greater than `num_sd` standard deviations from the
#'      average of averages.
#'  - `feature_cor()`: Detects dependencies between column metadata entries
#'        (features) and dataset rows.
#'  - `filter_col_outliers(num_sd=2, ctend='median', meas='pearson')`:
#'        Removes column outliers from the dataset. See `detect_col_outliers()`
#'        for details of outlier detection approach.
#'  - `filter_row_outliers(num_sd=2, ctend='median', meas='pearson')`:
#'        Removes row outliers from the dataset. See `detect_row_outliers()`
#'        for details of outlier detection approach.
#'  - `filter_cols(mask)`: Accepts a logical vector of length `ncol(obj$dat)`
#'        and returns a new EDAMatrix instance with only the columns associated
#'      with `TRUE` values in the mask.
#'  - `filter_rows(mask)`: Accepts a logical vector of length `nrow(obj$dat)`
#'        and returns a new EDAMatrix instance with only the rowsumns associated
#'      with `TRUE` values in the mask.
#'  - `impute(method='knn')`: Imputes missing values in the dataset and stores
#'        the result _in-place_. Currently only k-Nearest Neighbors (kNN)
#'        imputation is supported.
#'  - `log(base=exp(1), offset=0)`: Log-transforms data.
#'  - `log1p()`: Log(x + 1)-transforms data.
#'  - `pca(...)`: Performs principle component analysis (PCA) on the dataset
#'        and returns a new EDAMatrix instance of the projected data points.
#'      Any additional arguements specified are passed to the `prcomp()` function.
#'  - `pca_feature_cor(meas='pearson', ...)`: Measures correlation between
#'        dataset features (column metadata fields) and dataset principle
#'      components.
#'  - `plot_cor_heatmap(meas='pearson', interactive=TRUE, ...)`: Plots a
#'        correlation heatmap of the dataset.
#'  - `plot_densities(color=NULL, title="", ...)`: Plots densities for each
#'        column in the dataset.
#'  - `plot_feature_cor(meas='pearson', color_scale=c('green', 'red')`:
#'        Creates a tile plot of projected data / feature correlations. See
#'        `feature_cor()` function.
#'  - `plot_heatmap(interactive=TRUE, ...)`: Generates a heatmap plot of the
#'        dataset
#'  - `plot_pairwise_column_cors(color=NULL, title="", meas='pearson', mar=c(12,6,4,6))`:
#'        Plot median pairwise column correlations for each variable (column)
#'        in the dataset.
#'  - `plot_pca(pcx=1, pcy=2, scale=FALSE, color=NULL, shape=NULL, title=NULL,
#'               text_labels=FALSE, ...)`:
#'        Generates a two-dimensional PCA plot from the dataset.
#'  - `plot_tsne(color=NULL, shape=NULL, title=NULL, text_labels=FALSE, ...)`:
#'        Generates a two-dimensional t-SNE plot from the dataset.
#'  - `print()`: Prints an overview of the object instance.
#'  - `subsample(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL)`:
#'        Subsamples dataset rows and/or columns.
#'  - `summary(markdown=FALSE, num_digits=2)`: Summarizes overall
#'        characteristics of a dataset.
#'  - `t()`: Transposes dataset rows and columns.
#'  - `tsne(...)`: Performs T-distributed stochastic neighbor embedding (t-SNE)
#'        on the dataset and returns a new EDAMatrix instance of the projected
#'         data points. Any additional arguements specified are passed to the
#'        `Rtsne()` function.
#'  - `tsne_feature_cor(meas='pearson', ...)`: Measures correlation between
#'        dataset features (column metadata fields) and dataset t-SNE projected
#'      axes.
#'
#' @section Examples:
#' ```
#' library('eda')
#'
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
#' ```
#'
#' @importFrom R6 R6Class
#' @export
#' @name EDAMatrix
#'
NULL

EDAMatrix <- R6::R6Class("EDAMatrix",
    inherit = eda:::AbstractMatrixDataSet,
    public = list(
        # EDAMatrix constructor
        initialize = function(dat,
                              row_mdata=NULL, col_mdata=NULL,
                              row_names='rownames', col_names='colnames',
                              row_mdata_row_names='rownames',
                              col_mdata_row_names='rownames',
                              row_color=NULL, row_shape=NULL, row_label=NULL,
                              col_color=NULL, col_shape=NULL, col_label=NULL,
                              color_pal='Set1', title='', ggplot_theme=theme_bw) {
            # verify input data type and call parent constructor
            private$check_input(dat)

            # associate styles with row / column metadata, if present
            row_edat <- NULL
            col_edat <- NULL

            if (!is.null(row_mdata)) {
                row_edat <- 'row_mdata'
            }
            if (!is.null(col_mdata)) {
                col_edat <- 'col_mdata'
            }

            # create EDADat instances
            edats <- list('dat' = EDADat$new(dat, xid = 'x', yid = 'y',
                                             row_names = row_names, col_names = col_names,
                                             row_color = row_color, row_shape = row_shape,
                                             row_label = row_label, row_edat = row_edat,
                                             col_color = col_color, col_shape = col_shape,
                                             col_label = col_label, col_edat = col_edat))

            # add row and column metadata, if provided
            if (!is.null(row_mdata)) {
                # create row metadata EDADat instance
                row_edat <- EDADat$new(row_mdata, row_names = row_mdata_row_names)

                # determine data orientation and assign x and y axis id's
                if (length(intersect(rownames(edats[['dat']]$dat), rownames(row_edat$dat))) > 0) {
                    row_edat$xid <- 'x'
                    row_edat$yid <- 'row metadata'
                } else {
                    row_edat$xid <- 'row metadata'
                    row_edat$yid <- 'x'
                }
                edats[['row_mdata']] <- row_edat 
            }
            if (!is.null(col_mdata)) {
                # create column metadata EDADat instance
                col_edat <- EDADat$new(col_mdata, row_names = col_mdata_row_names)

                # determine data orientation and assign x and y axis id's
                if (length(intersect(colnames(edats[['dat']]$dat), rownames(col_edat$dat))) > 0) {
                    col_edat$xid <- 'y'
                    col_edat$yid <- 'column metadata'
                } else {
                    col_edat$xid <- 'column metadata'
                    col_edat$yid <- 'y'
                }
                edats[['col_mdata']] <- col_edat 
            }

            super$initialize(edats, color_pal, title, ggplot_theme)
        },

        cluster_tsne = function(num_clusters=10, ...) {
            super$cluster_tsne(key='dat', num_clusters=k, ...)
        },

        detect_col_outliers = function(num_sd=2, ctend=median, meas='pearson', ...) {
            super$detect_col_outliers(key='dat', num_sd=num_sd, ctend=ctend, meas = meas, ...)
        },

        detect_row_outliers = function(num_sd=2, ctend=median, meas='pearson', ...) {
            super$detect_row_outliers(key='dat', num_sd=num_sd, ctend=ctend, meas = meas, ...)
        },

        # Detects dependencies between column metadata entries (features) and
        # dataset rows
        #
        # Note: If metedata is not all-numeric, than a similarity measure which
        # supports categorical data (currently only 'lm') must be chosen.
        feature_cor = function(meas='lm') {
            if (is.null(self$col_mdata)) {
                stop("Error: missing column metadata.")
            }
            private$compute_cross_cor('dat', 'col_mdata', meas)
        },

        filter_col_outliers = function(num_sd=2, ctend=median, meas='pearson') {
            super$filter_col_outliers(key = 'dat', num_sd = num_sd, ctend = ctend, meas = meas)
        },

        filter_rows = function(mask) {
            super$filter_rows(key='dat', mask=mask)
        },

        filter_cols = function(mask) {
            super$filter_cols(key='dat', mask=mask)
        },

        filter_row_outliers = function(num_sd=2, ctend=median, meas='pearson') {
            super$filter_row_outliers(key = 'dat', num_sd = num_sd, ctend = ctend, meas = meas)
        },

        impute = function(method='knn') {
            super$impute(key='dat', method=method)
        },

        log = function(base=exp(1), offset=0) {
            super$log(key = 'dat', base = base, offset = offset)
        },

        log1p = function() {
            self$log(base = 1, offset = 1)
        },

        log2p = function() {
            self$log(base = 2, offset = 1)
        },

        pca_feature_cor = function(num_dims=10, meas='pearson', ...) {
            x <- self$t()$pca(...)
            x$dat <- x$dat[, 1:min(ncol(x$dat), num_dims)]
            x$t()$feature_cor(meas)
        },

        tsne_feature_cor = function(num_dims=10, meas='pearson', ...) {
            x <- self$t()$tsne(...)
            x$dat <- x$dat[, 1:min(ncol(x$dat), num_dims)]
            x$t()$feature_cor(meas)
        },

        plot_pairwise_column_cors = function(color=NULL, label=NULL, title="",
                                             meas='pearson',
                                             mar=c(12, 6, 4, 6), ...) {
            super$plot_pairwise_column_cors(key = 'dat', color = color,
                                            label = label, title = title,
                                            meas = meas, mar = mar, ...)
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
        plot_feature_cor = function(meas='lm', color_scale=c('green', 'red')) {
            # compute feature correlations
            cor_mat <- self$feature_cor(meas = meas)$edat[[1]]$dat
            dat <- melt(cor_mat)
            colnames(dat) <- c('dim', 'variable', 'value')

            # Labels
            #if (method == 'pca') {
            #    xlab_text <- 'Principle Components'
            #} else if (method == 't-sne') {
            #    xlab_text <- "t-SNE dimension"
            #}
            xlab_text <- ''

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

        # Correlation heatmap.
        #
        # Generates a correlation heatmap depicting the pairwise column
        # correlations in the data.
        #
        # meas String name of correlation measure to use.
        # ... Additional arguments
        #
        # @seealso \code{cor} for more information about supported correlation
        #      methods.
        plot_cor_heatmap = function(meas='pearson', interactive=TRUE, cor_args = list(), heatmap_args = list()) {
            # generate correlation matrix
            #cor_mat <- private$similarity(self$edat[['dat']]$dat, meas = meas, ...)
            cor_args <- c(list(self$edat[['dat']]$dat, meas = meas), cor_args)
            cor_mat <- do.call(private$similarity, cor_args)

            # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if (!is.null(self$col_mdata)) {
                # get metadata with variables stored as columns
                if (self$edat[['dat']]$yid == self$edat[['col_mdata']]$yid) {
                    mdata <- self$edat[['col_mdata']]$tdat
                } else {
                    mdata <- self$edat[['col_mdata']]$dat
                }

                # for heatmaps, show binary/logical variables on one side of the heatmap and
                # other variables on the other side
                lens <- apply(mdata, 2, function(x) {
                    length(unique(x))
                })
                binary_vars <- lens == 2

                # column colors (binary variables)
                if (sum(binary_vars) >= 1) {
                    params[['col_side_colors']] <- mdata[rownames(cor_mat), 
                                                         binary_vars, drop = FALSE]
                    params[['subplot_heights']] <- c(0.15, 0.3, 0.55)
                }

                # row colors (everything else)
                if (sum(!binary_vars) >= 1) {
                    params[['row_side_colors']] <- mdata[rownames(cor_mat), 
                                                         !binary_vars, drop = FALSE]
                    params[['subplot_widths']]  <- c(0.55, 0.3, 0.15)
                }
            }

            # add any additional function arguments
            params <- c(params, heatmap_args)

            private$construct_heatmap_plot(params, interactive)
        },

        # Generates a heatmap plot of the dataset
        #
        # ... Additional arguments
        plot_heatmap = function(interactive=TRUE, ...) {
            # list of parameters to pass to heatmaply
            params <- list(
                x               = self$edat[['dat']]$dat,
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
            params <- c(params, list(...))

            private$construct_heatmap_plot(params, interactive)
        },

        plot_densities = function(color=NULL, title="", ...) {
            super$plot_densities(key='dat', color=color, title=title, ...)
        },

        plot_pca = function(method='prcomp', pcx=1, pcy=2,
                            color=NULL, shape=NULL, label=NULL,
                            title=NULL, text_labels=NULL, ...) {
            super$plot_pca(key='dat', method=method, pcx=pcx, pcy=pcy,
                           color_var=color, shape_var=shape, label_var=label,
                           title=title, text_labels=text_labels, ...)
        },

        plot_tsne = function(dim1=1, dim2=2,
                             color=NULL, shape=NULL, label=NULL, title=NULL,
                             text_labels=NULL, ...) {
            super$plot_tsne(key='dat', dim1=dim1, dim2=dim2, color_var=color,
                            shape_var=shape, label_var=label, title=title,
                            text_labels=text_labels, ...)
        },

        # Prints class greeting to the screen
        print = function() {
            rm <- ifelse(!is.null(self$row_mdata), '(m)', '')
            cm <- ifelse(!is.null(self$col_mdata), '(m)', '')

            cat("=========================================\n")
            cat("=\n")
            cat("= EDAMatrix\n")
            cat("=\n")
            cat(sprintf("=   rows   : %d %s\n", nrow(self$dat), rm))
            cat(sprintf("=   columns: %d %s\n", ncol(self$dat), cm))
            cat("=\n")
            cat("=========================================\n")
        },

        # Subsample dataset
        subsample = function(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL) {
            # clone and subsample dataset
            obj <- private$clone_()
            obj$edat[['dat']]$subsample(row_n, col_n, row_ratio, col_ratio)

            # update metadata
            if ('row_mdata' %in% names(obj$edat)) {
                row_ind <- rownames(obj$row_mdata) %in% rownames(obj$dat)
                obj$row_mdata <- obj$row_mdata[row_ind, ]
            }
            if ('col_mdata' %in% names(obj$edat)) {
                col_ind <- colnames(obj$col_mdata) %in% colnames(obj$dat)
                obj$col_mdata <- obj$col_mdata[, col_ind]
            }

            obj
        },

        summary = function(markdown=FALSE, num_digits=2) {
            super$summary(key='dat', markdown=markdown, num_digits=num_digits)
        },

        transpose = function() {
            # transpose datasets in-place
            super$transpose()

            # swap keys for row/column metadata
            row_mdat <- self$edat[['row_mdata']]
            self$edat[['row_mdata']] <- self$edat[['col_mdata']]
            self$edat[['col_mdata']] <- row_mdat

            if (!is.null(self$edat[['dat']]$row_edat)) {
            	self$edat[['dat']]$row_edat <- 'row_mdata'
            }
            if (!is.null(self$edat[['dat']]$col_edat)) {
            	self$edat[['dat']]$col_edat <- 'col_mdata'
            }
        },

        # transpose (out-of-place)
        t = function() {
            obj <- private$clone_()
            obj$transpose()
            obj
        }
    ),
    private = list(
        # Verifies that input data is of type matrix
        check_input = function(dat) {
            if (!is.matrix(dat)) {
                stop("Invalid input for EDAMatrix: dat must be a matrix.")
            }
        },

        # Prints dataset summary to screen
        print_summary = function(x) {
            cat("=========================================\n")
            cat("\n")
            cat(sprintf(" %s\n", class(self)[1]))
            cat("\n")

            cat(" OVERVIEW\n")
            cat(" --------\n")
            cat("\n")
            cat(sprintf(" Rows      : %d\n", nrow(self$dat)))
            cat(sprintf(" Columns   : %d\n", ncol(self$dat)))
            cat(sprintf(" Min value : %s\n", x$dat_range[1]))
            cat(sprintf(" Max value : %s\n", x$dat_range[2]))
            cat(sprintf(" # NA's    : %d\n", x$dat_num_nas))
            cat(sprintf(" # 0's     : %d\n", x$dat_num_zeros))

            cat("\n")
            cat(' Quartiles:\n')
            quartiles <- setNames(as.data.frame(x$dat_quartiles), '')
            rownames(quartiles) <- paste0("  ", rownames(quartiles))
            print(quartiles)
            cat("\n")

            cat(" Column types:\n")
            col_types <- setNames(as.data.frame(x$col_types), c('', ''))
            col_types[, 1] <- paste0(" ", col_types[, 1])
            rownames(col_types) <- ""
            print(col_types)
            cat("\n")

            cat(" COLUMNS\n")
            cat(" -------\n")
            cat("\n")

            cat(sprintf(" Mean range   : %s - %s\n", x$col_means[1], x$col_means[2]))
            cat(sprintf(" Median range : %s - %s\n", x$col_medians[1], x$col_medians[2]))
            cat(sprintf(" Stdev range  : %s - %s\n", x$col_std_devs[1], x$col_std_devs[2]))
            cat(sprintf(" Cor range    : %s - %s\n",
                        min(x$col_cor_mat, na.rm = TRUE),
                        max(x$col_cor_mat, na.rm = TRUE)))
            cat("\n")

            cat(" Outliers:\n\n ")
            cat(paste0(sprintf("%2d", 1:length(x$col_outliers)), ". ", x$col_outliers, "\n"))
            cat("\n")

            cat(" ROWS\n")
            cat(" ----\n")
            cat("\n")

            cat(sprintf(" Mean range   : %s - %s\n", x$row_means[1], x$row_means[2]))
            cat(sprintf(" Median range : %s - %s\n", x$row_medians[1], x$row_medians[2]))
            cat(sprintf(" Stdev range  : %s - %s\n", x$row_std_devs[1], x$row_std_devs[2]))
            cat(sprintf(" Cor range    : %s - %s\n",
                        min(x$row_cor_mat, na.rm = TRUE),
                        max(x$row_cor_mat, na.rm = TRUE)))

            cat("\n")
            cat(" Outliers:\n\n ")
            cat(paste0(sprintf("%2d", 1:length(x$row_outliers)), ". ", x$row_outliers, "\n"))
            cat("\n")
            cat("=========================================\n")
        }

        # Determines metadata columns to include for heatmap, etc. functions
        #
        # @param include Vector of strings indicating metadata columns which
        # should be included in the analysis.
        # @param exclude Vector of strings indicating metadata columns which
        # should be excluded from the analysis.
        #
        # @return Character vector of fields to include in analysis.
        #select_features = function(include, exclude) {
        #    # determine fields to include based on user arguments
        #    if (!is.null(include)) {
        #        include <- include
        #    } else if (!is.null(exclude)) {
        #        include <- colnames(self$row_mdata)[!colnames(self$row_mdata) %in% exclude]
        #    } else {
        #        include <- colnames(self$row_mdata)
        #    }

        #    # always exclude fields with all unique values (e.g. identifiers)
        #    mask <- sapply(self$row_mdata[, include], function(x) {
        #        max(table(x)) > 1
        #    })
        #    include[mask]
        #}
    ),
    # ------------------------------------------------------------------------
    # active bindings
    # ------------------------------------------------------------------------
    active = list(
        dat = function(value) {
            if (missing(value)) {
                self$edat[['dat']]$dat
            } else {
                self$edat[['dat']]$dat <- value
            }
        },
        row_mdata = function(value) {
            if (missing(value)) {
                if ('row_mdata' %in% names(self$edat)) {
                    self$edat[['row_mdata']]$dat
                } else {
                    NULL
                }
            } else {
                if ('row_mdata' %in% names(self$edat)) {
                    self$edat[['row_mdata']]$dat <- value
                } else {
                    self$edat[['row_mdata']] <- EDADat$new(value, xid = 'x', yid = 'row metadata')
                }
            }
        },
        col_mdata = function(value) {
            if (missing(value)) {
                if ('col_mdata' %in% names(self$edat)) {
                    self$edat[['col_mdata']]$dat
                } else {
                    NULL
                }
            } else {
                if ('col_mdata' %in% names(self$edat)) {
                    self$edat[['col_mdata']]$dat <- value
                } else {
                    self$edat[['col_mdata']] <- EDADat$new(value, xid = 'y', yid = 'column metadata')
                }
            }
        }
    )
)
