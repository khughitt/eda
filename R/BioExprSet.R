#' An R6 class representing an Expression dataset
#'
#' BioExprSet is a simple class for interacting with biological expression
#' data (e.g. from microarray or RNA-Seq experiments.) It can accept either
#' a matrix and optional column/row metadata dataframes, or a Bioconductor
#' ExpressionSet instance. Methods are provided for common transformations
#' and visualizations.
#'
#' @section Usage:
#' ```
#' bset <- BioExprSet$new(count_matrix, row_mdata=gene_annotations,
#'                        col_mdata=sample_annotations, row_color='some_var')
#' bset$summary()
#'
#' bset$cpm()$plot_pca()
#' bset$t$subsample(1500)$plot_cor_heatmap()
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
#' - `row_mdata_rownames`: Column name or number containing row metadata row
#'      identifiers. If set to `rownames` (default), row names will be used
#'      as identifiers.
#' - `col_mdata_ids`: Column name or number containing col metadata row
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
#' - `clear_cache()`: Clears BioExprSet cache.
#' - `clone()`: Creates a copy of the BioExprSet instance.
#' - `cluster_tsne(k=10, ...)`: Clusters rows in dataset using a combination
#'      of t-SNE and k-means clustering.
#' - `cpm()`: Performs counts-per-million (CPM) transformation.
#' - `detect_col_outliers(num_sd=2, avg='median', method='pearson')`:
#'      Measures average pairwise similarities between all columns in the dataset.
#'      Outliers are considered to be those columns who mean similarity to
#'      all other columns is greater than `num_sd` standard deviations from the
#'      average of averages.
#' - `detect_row_outliers(num_sd=2, avg='median', method='pearson')`:
#'      Measures average pairwise similarities between all rows in the dataset.
#'      Outliers are considered to be those rows who mean similarity to
#'      all other rows is greater than `num_sd` standard deviations from the
#'      average of averages.
#'  - `feature_cor()`: Detects dependencies between column metadata entries
#'		(features) and dataset rows.
#'  - `filter_col_outliers(num_sd=2, avg='median', method='pearson')`:
#'		Removes column outliers from the dataset. See `detect_col_outliers()`
#'		for details of outlier detection approach.
#'  - `filter_row_outliers(num_sd=2, avg='median', method='pearson')`:
#'		Removes row outliers from the dataset. See `detect_row_outliers()`
#'		for details of outlier detection approach.
#'  - `filter_cols(mask)`: Accepts a logical vector of length `ncol(obj$dat)`
#'		and returns a new BioExprSet instance with only the columns associated
#'      with `TRUE` values in the mask.
#'  - `filter_rows(mask)`: Accepts a logical vector of length `nrow(obj$dat)`
#'		and returns a new BioExprSet instance with only the rowsumns associated
#'      with `TRUE` values in the mask.
#'  - `impute(method='knn')`: Imputes missing values in the dataset and stores
#'		the result _in-place_. Currently only k-Nearest Neighbors (kNN)
#'		imputation is supported.
#'  - `log(base=exp(1), offset=0)`: Log-transforms data.
#'  - `log1p()`: Logn(x + 1)-transforms data.
#'  - `log2p()`: Log2(x + 1)-transforms data.
#'  - `pca(...)`: Performs principle component analysis (PCA) on the dataset
#'		and returns a new BioExprSet instance of the projected data points.
#'      Any additional arguements specified are passed to the `prcomp()` function.
#'  - `pca_feature_cor(method='pearson', ...)`: Measures correlation between
#'		dataset features (column metadata fields) and dataset principle
#'      components.
#'  - `plot_cor_heatmap(method='pearson', interactive=TRUE, ...)`: Plots a
#'		correlation heatmap of the dataset.
#'  - `plot_densities(color=NULL, title="", ...)`: Plots densities for each
#'		column in the dataset. (For histogram and bar plots, see the
#'		`plot_libsizes` method.)
#'  - `plot_feature_cor(method='pearson', color_scale=c('green', 'red')`:
#'		Creates a tile plot of projected data / feature correlations. See
#'		`feature_cor()` function.
#'  - `plot_heatmap(interactive=TRUE, ...)`: Generates a heatmap plot of the
#'		dataset
#'  - `plot_libsizes(color=NULL, title=NULL, geom='bar')`: Creates a histogram
#'     or bar plot of sample library sizes.
#'  - `plot_pairwise_column_cors(color=NULL, title="", method='pearson', mar=c(12,6,4,6))`:
#'		Plot median pairwise column correlations for each variable (column)
#'		in the dataset.
#'  - `plot_pca(pcx=1, pcy=2, scale=FALSE, color=NULL, shape=NULL, title=NULL,
#'               text_labels=FALSE, ...)`:
#'		Generates a two-dimensional PCA plot from the dataset.
#'  - `plot_tsne(color=NULL, shape=NULL, title=NULL, text_labels=FALSE, ...)`:
#'		Generates a two-dimensional t-SNE plot from the dataset.
#'  - `print()`: Prints an overview of the object instance.
#'  - `subsample(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL)`:
#'		Subsamples dataset rows and/or columns.
#'  - `summary(markdown=FALSE, num_digits=2)`: Summarizes overall
#'		characteristics of a dataset.
#'  - `t()`: Transposes dataset rows and columns.
#'  - `tsne(...)`: Performs T-distributed stochastic neighbor embedding (t-SNE)
#'		on the dataset and returns a new BioExprSet instance of the projected
#' 		data points. Any additional arguements specified are passed to the
#'		`Rtsne()` function.
#'  - `tsne_feature_cor(method='pearson', ...)`: Measures correlation between
#'		dataset features (column metadata fields) and dataset t-SNE projected
#'      axes.
#'
#' @section Examples:
#' ```
#' library('eda')
#'
#' dat <- as.matrix(iris[,1:4])
#' row_mdata <- iris[,5,drop=FALSE]
#'
#' bset <- BioExprSet$new(dat, row_mdata=row_mdata, row_color='Species')
#'
#' bset
#' bset$summary()
#'
#' bset$plot_pca()
#' bset$log1p()$plot_cor_heatmap()
#' bset$subsample(100)$plot_tsne()
#' ```
#'
#' @importFrom R6 R6Class
#' @export
#' @name BioExprSet
#'
NULL

BioExprSet <- R6::R6Class("BioExprSet",
    inherit = eda:::BioDataSet,
    public = list(
        # BioExprSet constructor
        initialize = function(dat,
                              row_mdata=NULL, col_mdata=NULL,
                              row_names='rownames', col_names='colnames',
                              row_mdata_row_names='rownames', 
                              col_mdata_row_names='rownames',
                              row_color=NULL, row_shape=NULL, row_label=NULL,
                              col_color=NULL, col_shape=NULL, col_label=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) {
            # verify input data type and call parent constructor
            private$check_input(dat)

            # for ExpressionSets, retrieve relevant data and metadata
            if (class(dat) == 'ExpressionSet') {
                if (is.null(col_mdata)) {
                    col_mdata <- pData(dat)
                }
                if (is.null(row_mdata)) {
                    row_mdata <- fData(dat)
                }
                dat <- exprs(dat)
            }

            # create EDADat instances
            edats <- list('dat' = EDADat$new(dat, 
                                             xid = 'genes', yid = 'samples',
                                             row_names = row_names, col_names = col_names, 
                                             row_color = row_color, row_shape = row_shape,
                                             row_label = row_label, row_edat = 'row_mdata',
                                             col_color = col_color, col_shape = col_shape, 
                                             col_label = col_label, col_edat = 'col_mdata'))

            # add row and column metadata, if provided
            if (!is.null(row_mdata)) {
                edats[['row_mdata']] <- EDADat$new(row_mdata,
                                                   xid = 'genes', yid = 'gene metadata',
                                                   row_names = row_mdata_row_names)
            }
            if (!is.null(col_mdata)) {
                edats[['col_mdata']] <- EDADat$new(col_mdata, 
                                                   xid = 'samples', yid = 'sample metadata',
                                                   row_names = col_mdata_row_names)
            }

            # call parent constructor
            super$initialize(edats, color_pal, title, ggplot_theme)
        },

        cluster_tsne = function(k=10, ...) {
            super$cluster_tsne(key='dat', k=k, ...)
        },

        # Performs a counts-per-million (CPM) transformation.
        #
        # @return A CPM-transformed version of the BioExprSet instance.
        cpm = function() {
            obj <- private$clone_() 
            obj$dat <- sweep(obj$dat, 2, colSums(obj$dat), '/') * 1E6
            obj
        },

        #diff_expr = function() {
        #},

        # Detects dependencies between column metadata entries (features) and
        # dataset rows
        #
        # Note: If metedata is not all-numeric, than a similarity method which
        # supports categorical data (currently only 'lm') must be chosen.
        feature_cor = function(method='lm', include=NULL, exclude=NULL, ...) {
            if (is.null(self$col_mdata)) {
                stop("Error: missing column metadata.")
            }

            # determine which features to include
            features <- colnames(self$col_mdata)

            if (!is.null(include)) {
                features <- include
            }
            if (!is.null(exclude)) {
                features <- features[!features %in% exclude]
            }
            mask <- colnames(self$col_mdata) %in% features

            # for metadata, we can also exclude factor fields with all
            # unique values (e.g. alternate identifiers)
            #exclude <- apply(dat2, 1, function(x) {
            #    is.factor(x) && length(unique(x)) == length(x)
            #})

            #if (sum(exclude) > 0) {
            #    message(sprintf("Excluding %d unique factor fields", sum(exclude)))
            #    dat2 <- dat2[!exclude, ]
            if (sum(!mask) > 0) {
                super$filter_cols(key='col_mdata', mask)$cross_cor('dat', 'col_mdata', method, ...)
            } else {
                self$cross_cor('dat', 'col_mdata', method, ...)
            }
        },

        pca_feature_cor = function(num_dims=10, method='lm', include=NULL, exclude=NULL, ...) {
            self$t()$pca(num_dims = num_dims, ...)$t()$feature_cor(method, include, exclude)
        },

        plot_pca_feature_cor = function(num_dims=10, method='lm', include=NULL, exclude=NULL, top_n=NULL, ...) {
            self$t()$pca(num_dims = num_dims, ...)$t()$plot_feature_cor(method, include, exclude, top_n = top_n)
        },

        tsne_feature_cor = function(num_dims=10, method='lm', include=NULL, exclude=NULL, ...) {
            self$t()$tsne(num_dims = num_dims, ...)$t()$feature_cor(method, include, exclude)
        },

        plot_tsne_feature_cor = function(num_dims=10, method='lm', include=NULL, exclude=NULL, top_n=NULL, ...) {
            self$t()$tsne(num_dims = num_dims, ...)$t()$plot_feature_cor(method, include, exclude, top_n = top_n)
        },

        filter_rows = function(mask) {
            super$filter_rows(key='dat', mask=mask)
        },

        filter_cols = function(mask) {
            super$filter_cols(key='dat', mask=mask)
        },

        filter_row_outliers = function(num_sd=2, ctend=median, method='pearson') {
            super$filter_row_outliers(key='dat', num_sd=num_sd, ctend=ctend, method=method)
        },

        filter_col_outliers = function(num_sd=2, ctend=median, method='pearson') {
            super$filter_col_outliers(key='dat', num_sd=num_sd, ctend=ctend, method=method)
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

        pca = function(num_dims=NULL, ...) {
            super$pca(key='dat', num_dims = num_dims, ...)
        },

        tsne = function(num_dims=NULL, ...) {
            super$tsne(key='dat', num_dims = num_dims, ...)
        },

        plot_feature_cor = function(method='pearson', include=NULL, exclude=NULL, top_n=NULL, color_scale=c('green', 'red'), ...) {
            # compute feature correlations
            cor_mat <- self$feature_cor(method = method, include = include,
                                        exclude = exclude, ...)

            # if requested, limit to top N features with the highest average correlation
            # to data rows
            if (!is.null(top_n)) {
                top_features <- names(head(sort(colSums(cor_mat), decreasing=TRUE), top_n))
                cor_mat <- cor_mat[, colnames(cor_mat) %in% top_features]
            }
            dat <- melt(cor_mat)
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

        plot_pca = function(pcx=1, pcy=2, scale=FALSE,
                            color=NULL, shape=NULL, label=NULL, 
                            title=NULL, text_labels=FALSE, ...) {
            super$plot_pca(key='dat', pcx=pcx, pcy=pcy, scale=scale,
                           color_var=color, shape_var=shape, label_var=label,
                           title=title, text_labels=text_labels, ...)
        },

        plot_tsne = function(color=NULL, shape=NULL, label=NULL, title=NULL,
                             text_labels=FALSE, ...) {
            super$plot_tsne(key='dat', color_var=color, shape_var=NULL,
                            label_var=label, title=NULL,
                            text_labels=text_labels, ...)
        },

        # Prints an overview of the object instance
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= BioExprSet (%s)\n", class(self$dat[,1])))
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        },
        transpose = function() {
            # transpose datasets in-place
            super$transpose()

            # swap keys for row/column metadata
            row_mdat <- self$edat[['row_mdata']] 
            self$edat[['row_mdata']] <- self$edat[['col_mdata']]
            self$edat[['col_mdata']] <- row_mdat
        }
    ),
    private = list(
        # Verifies that input data is an acceptable format.
        check_input = function(dat) {
            if (!is.matrix(dat) && (class(dat) != 'ExpressionSet')) {
                stop("Invalid input for BioExprSet: dat must be a matrix or ExpressionSet.")
            }
        }
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
                    self$edat[['row_mdata']] <- EDADat$new(value)
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
                    self$edat[['col_mdata']] <- EDADat$new(value)
                }
            }
        }
    )
)

