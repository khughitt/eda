#' An S6 class representing an Expression dataset
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
#' - `row_ids`: Column name or number containing row identifiers. If set to
#'      `rownames` (default), row names will be used as identifiers.
#' - `col_ids`: Column name or number containing column identifiers. If set to
#'      `colnames` (default), column names will be used as identifiers.
#' - `row_mdata_ids`: Column name or number containing row metadata row 
#'      identifiers. If set to `rownames` (default), row names will be used 
#'      as identifiers.
#' - `col_mdata_ids`: Column name or number containing col metadata row 
#'      identifiers. If set to `rownames` (default), row names will be used 
#'      as identifiers.
#' - `row_color`: Row metadata field to use for coloring rowwise plot elements.
#' - `row_shape`: Row metadata field to use for determine rowwise plot 
#'      element shape.
#' - `row_labels`: Row metadata field to use when labeling plot points or
#'      other elements.
#' - `col_color`: Column metadata field to use for coloring columnwise plot elements.
#' - `col_shape`: Column metadata field to use for determine columnwise plot 
#'      element shape.
#' - `col_labels`: Column metadata field to use when labeling plot points or
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
#' - `detect_col_outliers(num_sd=2, avg='median', sim_method='pearson')`:
#'      Measures average pairwise similarities between all columns in the dataset.
#'      Outliers are considered to be those columns who mean similarity to
#'      all other columns is greater than `num_sd` standard deviations from the
#'      average of averages.
#' - `detect_row_outliers(num_sd=2, avg='median', sim_method='pearson')`:
#'      Measures average pairwise similarities between all rows in the dataset.
#'      Outliers are considered to be those rows who mean similarity to
#'      all other rows is greater than `num_sd` standard deviations from the
#'      average of averages.
#'  - `feature_cor()`: Detects dependencies between column metadata entries 
#'		(features) and dataset rows. 
#'  - `filter_col_outliers(num_sd=2, avg='median', sim_method='pearson')`: 
#'		Removes column outliers from the dataset. See `detect_col_outliers()` 
#'		for details of outlier detection approach.
#'  - `filter_row_outliers(num_sd=2, avg='median', sim_method='pearson')`: 
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
#'		column in the dataset.
#'  - `plot_feature_cor(method='pearson', color_scale=c('green', 'red')`:
#'		Creates a tile plot of projected data / feature correlations. See
#'		`feature_cor()` function.
#'  - `plot_heatmap(interactive=TRUE, ...)`: Generates a heatmap plot of the 
#'		dataset
#'  - `plot_libsize_hist(color=NULL, title=NULL)`: Plots a histogram of sample 
#'     library sizes.
#'  - `plot_libsize_bargraph(color=NULL, title=NULL)`: Plots a bargraph of 
#'     sample library sizes.
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
    inherit = EDAMatrix,
    public = list(
        # BioExprSet constructor
        initialize = function(dat, 
                              row_mdata=NULL, col_mdata=NULL, 
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
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

            super$initialize(dat, row_mdata, col_mdata, row_ids, col_ids,
                             row_mdata_ids, col_mdata_ids, row_color, row_shape,
                             row_labels, col_color, col_shape, col_labels,
                             color_pal, title, ggplot_theme)
        },

        # Performs a counts-per-million (CPM) transformation.
        #
        # @return A CPM-transformed version of the BioExprSet instance.
        cpm = function() {
            obj <- self$clone()
            obj$dat <- sweep(obj$dat, 2, colSums(obj$dat), '/') * 1E6
            obj
        },

        diff_expr = function() {
        
        },

        # Log2 transforms data (adding 1 to ensure finite results).
        #
        # @return Log2-transformed version of the expression data.
        log2p = function() {
            self$log(2, offset=1)
        },

        # Plots a histogram of sample library sizes.
        #
        # Generates a histogram of sample library sizes (sum of expression
        # levels within each column/sample). 
        #
        # @param color Column metadata field to use for coloring points.
        # @param title Title to use for plot.
        #
        # @return A ggplot instance.
        plot_libsize_hist = function(color=NULL, title=NULL) {
            dat <- data.frame(libsize=colSums(self$dat))

            styles <- private$get_geom_histogram_styles(color)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Library sizes: %s", private$title)
            }

            plt <- ggplot(aes(x=libsize), data=dat) +
                geom_histogram(styles$aes) +
                private$ggplot_theme()

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}
            plt
        },

        # Plots bar graph of sample library sizes
        #
        # Generates a bargraph plot of sample library sizes (sum of expression
        # levels within each column/sample). Each sample is shown as a separate
        # bar in the plot.
        #
        # @param color Column metadata field to use for coloring points.
        # @param title Title to use for plot.
        #
        # @return A ggplot instance.
        plot_libsize_bargraph = function(color=NULL, title=NULL) {
            # data frame with sample library sizes
            dat <- data.frame(libsize=colSums(self$dat))

            styles <- private$get_geom_histogram_styles(color)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Library sizes: %s", private$title)
            }

            plt <- ggplot(aes(x=rownames(libsizes), y=libsize), data=dat) +
                geom_bar(styles$aes, stat='identity') + 
                   xlab("Samples") +
                   ylab("Total expression") +
                   private$ggplot_theme() +
                   theme(axis.text.x=element_text(angle=90),
                         legend.text=element_text(size=8))

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}
            plt
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
        }
    ),
    private = list(
        # Verifies that input data is an acceptable format.
        check_input = function(dat) {
            if(!is.matrix(dat) && (class(dat) != 'ExpressionSet')) {
                stop("Invalid input for BioExprSet: dat must be a matrix or ExpressionSet.")
            }
        }
    )
)

