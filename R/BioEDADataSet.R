#' An S6 class representing an biological dataset
#'
#' BioEDADataSet is an abstract class for interacting with biological datasets,
#' typically those collected from high-throughput experiments. It is not meant
#' to be used directly, but rather should is subclassed by other datatype-
#' specific classes such as BioExprDataSet.
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
#' - `clear_cache()`: Clears BioEDADataSet cache.
#' - `clone()`: Creates a copy of the BioEDADataSet instance.
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
#'		and returns a new BioEDADataSet instance with only the columns associated
#'      with `TRUE` values in the mask.
#'  - `filter_rows(mask)`: Accepts a logical vector of length `nrow(obj$dat)`
#'		and returns a new BioEDADataSet instance with only the rowsumns associated
#'      with `TRUE` values in the mask.
#'  - `impute(method='knn')`: Imputes missing values in the dataset and stores
#'		the result _in-place_. Currently only k-Nearest Neighbors (kNN)
#'		imputation is supported.
#'  - `log(base=exp(1), offset=0)`: Log-transforms data.
#'  - `log1p()`: Logn(x + 1)-transforms data.
#'  - `log2p()`: Log2(x + 1)-transforms data.
#'  - `pca(...)`: Performs principle component analysis (PCA) on the dataset
#'		and returns a new BioEDADataSet instance of the projected data points.
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
#'		on the dataset and returns a new BioEDADataSet instance of the projected
#' 		data points. Any additional arguements specified are passed to the
#'		`Rtsne()` function.
#'  - `tsne_feature_cor(method='pearson', ...)`: Measures correlation between
#'		dataset features (column metadata fields) and dataset t-SNE projected
#'      axes.
#'
#' @importFrom R6 R6Class
#' @export
#' @name BioEDADataSet
#'
NULL

BioEDADataSet <- R6::R6Class("BioEDADataSet",
    inherit = EDAMatrix,
    public = list(
        # BioEDADataSet constructor
        initialize = function(dat,
                              row_mdata=NULL, col_mdata=NULL,
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) {

            super$initialize(dat, row_mdata, col_mdata, row_ids, col_ids,
                             row_mdata_ids, col_mdata_ids, row_color, row_shape,
                             row_labels, col_color, col_shape, col_labels,
                             color_pal, title, ggplot_theme)
        },

        # Computes pathway- or annotation-level statistics for a given dataset and
        # annotation mapping.
        #
        # @param annot A n x 2 gene/pathway mapping where each row is a pathway,gene
        #     pair; should contain columns named 'pathway' and 'gene'.
        # @param stat The pathway-level statistic to compute. Can either be a
        #   function (e.g. sum, median, or var), or a string indicating a specific statistic
        #   to compute.
        #
        # Currently supported statistics include:
        #
        #   - `num_above_cutoff`    number of genes w/ values above some cutoff
        #   - `num_below_cutoff`    number of genes w/ values below some cutoff
        #   - `ratio_above_cutoff`  ratio of genes w/ values above some cutoff
        #   - `ratio_below_cutoff`  ratio of genes w/ values below some cutoff
        #   - `num_nonzero`         number of genes with values not equal to zero
        #   - `num_zero`            number of genes with values equal to zero
        #   - `ratio_nonzero`       ratio of genes with values not equal to zero
        #   - `ratio_zero`          ratio of genes with values equal to zero
        #
        compute_pathway_stats = function(annot, stat=median, ...) {
            # output data frame
            res <- data.frame()

            # determine statistic to use
            if (!is.function(stat)) {
                if (stat %in% names(private$pathway_stats)) {
                    stat <- private$pathway_stats[[stat]]
                } else {
                    stop("Invalid pathway statistic specified.")
                }
            }

            # iterate over annotations / pathways
            for (pathway in unique(annot$pathway)) {
                # get list of genes in the pathway
                genes <- annot$gene[annot$pathway == pathway]

                # get data values for relevant genes
                path_dat <- self$dat[rownames(self$dat) %in% genes,, drop = FALSE]

                # compute statistic for each column (cell line) and append to result
                res <- rbind(res, apply(path_dat, 2, stat, ...))
            }

            # fix column and row names and return result
            rownames(res) <- unique(annot$pathway)
            colnames(res) <- colnames(self$dat)

            as.matrix(res)
        },

        # Log2 transforms data (adding 1 to ensure finite results).
        #
        # @return Log2-transformed version of the expression data.
        log2p = function() {
            self$log(2, offset = 1)
        },

        # Creates a histogram or bar plot of sample library sizes.
        #
        # Generates a bargraph or histogram plot of sample library sizes
        # (sum of expression levels within each column/sample). For the bar
        # graph plot, each sample is shown as a separate bar in the plot.
        # For larger datasets, a histogram of library sizes may be more useful
        # to display.
        #
        # @param color Column metadata field to use for coloring points.
        # @param title Title to use for plot.
        # @param geom ggplot geometry to use (bar|histogram)
        #
        # @return A ggplot instance.
        plot_libsizes = function(color=NULL, title=NULL, geom='hist') {
            # create a data frame of library sizes
            dat <- data.frame(libsize = colSums(self$dat))

            # determine plot styles to use
            if (geom == 'hist') {
                styles <- private$get_geom_histogram_styles(color)
            } else {
                styles <- private$get_geom_bar_styles(color)
            }

            # update styles and title
            if (!is.null(styles$color)) {
                dat <- cbind(dat, color = styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Library sizes: %s", private$title)
            }

            # create plot
            if (geom == 'hist') {
                # histogram
                plt <- ggplot(aes(x = libsize), data = dat) +
                    geom_histogram(styles$aes) +
                    private$ggplot_theme()
            } else {
                # bar plot
                plt <- ggplot(aes(x = rownames(libsizes), y = libsize), data = dat) +
                    geom_bar(styles$aes, stat = 'identity') +
                    xlab("Samples") +
                    ylab("Total expression") +
                    private$ggplot_theme() +
                    theme(axis.text.x = element_text(angle = 90),
                          legend.text = element_text(size = 8))
            }

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
            cat(sprintf("= BioEDADataSet (%s)\n", class(self$dat[,1])))
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        }
    ),
    private = list(
        # Helper functions for pathway-level statistics; used by the
        # `compute_pathway_stats` method.
        pathway_stats = list(
            'num_nonzero' = function(x) {
                sum(x != 0)
            },
            'num_zero' = function(x) {
                sum(x == 0)
            },
            'num_above_cutoff' = function(x, cutoff=0) {
                sum(x > cutoff)
            },
            'num_below_cutoff' = function(x, cutoff=Inf) {
                sum(x < cutoff)
            },
            'ratio_nonzero' = function(x) {
                sum(x != 0) / length(x)
            },
            'ratio_zero' = function(x) {
                sum(x == 0) / length(x)
            },
            'ratio_above_cutoff' = function(x, cutoff=0) {
                sum(x > cutoff) / length(x)
            },
            'ratio_below_cutoff' = function(x, cutoff=Inf) {
                sum(x < cutoff) / length(x)
            }
        ),

        # Verifies that input data is an acceptable format.
        check_input = function(dat) {
            if (!is.matrix(dat) && (class(dat) != 'ExpressionSet')) {
                stop("Invalid input for BioEDADataSet: dat must be a matrix or ExpressionSet.")
            }
        }
    )
)

