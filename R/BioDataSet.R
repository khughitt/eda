#' An S6 class representing collection of related biological datasets
#'
#' BioEDADataSet is a class for interacting with one or more biological
#' datasets, typically those collected from high-throughput experiments.
#'
#' @section Arguments:
#' - `dataset`: A list of datasets (matrices, data frames, etc.), each of
#'      which shared some column / row identifiers with the first entry in
#'      the list.
#'
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
#'  - `edat`: List of EDADat dataset instances
#'  - `annotations`: A list of gene, etc. annotations from external sources.
#'
#' @section Methods:
#' - `clear_cache()`: Clears BioEDADataSet cache.
#' - `clone()`: Creates a copy of the BioEDADataSet instance.
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
#' @name BioDataSet
#' @export
#'
NULL

BioDataSet <- R6Class("BioDataSet",
    inherit = eda:::EDAMultiMatrix,

    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # List of loaded annotations
        annotations = list(),

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

        # Computes pathway- or annotation-level statistics for a given dataset and
        # annotation mapping.
        #
        # @param key Name of dataset to analyze.
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
        annotation_stats = function(key, annotation, stat=median, ...) {
            # check for valid dataset key
            if (!key %in% c(1:length(self$edat), names(self$edat))) {
                stop(sprintf("Invalid dataset specified: %s", key))
            }

            # determine statistic to use
            if (!is.function(stat)) {
                if (stat %in% names(private$annotation_stat_fxns)) {
                    stat <- private$annotation_stat_fxns[[stat]]
                } else {
                    stop("Invalid statistic specified.")
                }
            }

            # check to see if already computed
            stat_name <- as.character(substitute(stat))

            # output data frame
            res <- data.frame()

            # check to see if annotation has been loaded, and if not, load it
            # TODO
            mapping <- self$load_annotations(annotation)

            # annotations are parsed into n x 2 dataframes consisting of
            # with each row containing an (annotation, gene id) pair.
            ANNOT_IND <- 1
            ITEM_IND  <- 2

            # dataset to compute statistics on
            dat <- self$edat[[key]]$dat

            message("Computing annotation statistics...")

            # iterate over annotations
            for (annot in unique(mapping[, ANNOT_IND])) {
                # get list of genes, etc. associated with the annotation
                annot_items <- mapping[mapping[, ANNOT_IND] == annot, ITEM_IND]

                # get data values for relevant genes
                dat_subset <- dat[rownames(dat) %in% annot_items,, drop = FALSE]

                # compute statistic for each column (cell line) and append to result
                res <- rbind(res, apply(dat_subset, 2, stat, ...))
            }

            # fix column and row names and return result
            rownames(res) <- unique(mapping[, ANNOT_IND])
            colnames(res) <- colnames(dat)

            # clone BioDataSet instance and replace main dataset with 
            obj <- private$clone_()
            obj$edat[[key]] <- NULL

            # if replacing the main dataset, also remove any additional datasets
            # that may be linked by row ids (for now, it is assumed that
            # function is passed expression, etc. data in the usual orientation)
            if ((key == 1) || (!is.null(names(obj$edat)) && names(obj$edat)[1] == key)) {
                for (i in seq_along(obj$edat)) {
                    if (obj$edat[[i]]$orientation == 'rows') {
                        obj$edat[[i]] <- NULL
                    }    
                }
            }

            # set annotation_stat result as main dataset and return
            obj$edat[[1]] <- res

            # key form: <old key>_<annotation>_<stat>
            stat_name <- as.character(substitute(stat))
            res_key <- paste(c(key, annotation, stat_name), collapse='_')
            names(obj$edat)[1] <- res_key

            obj
        },
        
        # load annotations from a supported source
        load_annotations = function(source='cpdb', keytype='ensembl') {
            # check to make sure annotation type is valid
            if (!source %in% c('cpdb')) {
                stop("Unsupported annotation source specified.") 
            }

            # If annotations have already been retrieved, simply return them
            if (source %in% names(self$annotations)) {
                return(self$annotations[[source]])
            }

            message(sprintf("Loading %s annotations...", source))

            # CPDB
            if (source == 'cpdb') {
                private$get_cpdb_annotations(keytype)
            }
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
            dat <- data.frame(libsize = colSums(self$edat[[1]]$dat))

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
            cls <- class(self)[1]

            # create a vector of dataset keys to display
            if (!is.null(names(self$edat))) {
                keys <- names(self$edat)

                # default to numeric key if no character key exists
                keys[is.na(keys)] <- seq_along(keys)[is.na(keys)]
            }

            # determine lengtht to use for justification
            key_format <- sprintf("%%%ds", max(nchar(keys)) + 1)

            # entry output string
            entry_template <- sprintf("= %s : %%s (%%d x %%d)\n", key_format)

            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= %s (n=%d)\n", cls, length(self$edat)))
            cat("=\n")
            for (i in seq_along(self$edat)) {
                ds <- self$edat[[i]]$dat
                
                # print dataset entry
                cat(sprintf(entry_template, keys[i], class(ds), nrow(ds), ncol(ds)))
            }
            cat("=\n")
            if (length(self$annotations) > 0) {
                cat("= Annotations:\n")
                cat("=\n")
                for (source in names(self$annotations)) {
                    annot <- self$annotations[[source]]
                    cat(sprintf("= - %s (%d annotations, %d genes)\n", source, 
                                length(unique(annot[,1])),
                                length(annot[,2])))
                }
                cat("=\n")
            }
            cat("=========================================\n")
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

        # Plots multidataset correlation heatmap
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        plot_cross_cor_heatmap = function(key1=1, key2=2, method='pearson', interactive=TRUE) {
            super$plot_cross_cor_heatmap(key1, key2, method, interactive)
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # Helper functions for pathway-level statistics; used by the
        # `annotation_stats` method.
        annotation_stat_fxns = list(
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
        },

        # load ConsensusPathDB annotations
        get_cpdb_annotations = function(keytype) {
            # Check to make sure key type is valid
            if (!keytype %in% c('ensembl', 'entrez-gene', 'hgnc-symbol')) {
                stop("Invalid key type specified.")
            }

            # Load CPDB pathways
            url <- sprintf('http://cpdb.molgen.mpg.de/CPDB/getPathwayGenes?idtype=%s', keytype)
            cpdb <- read.delim(url, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

            # There are a few pathways with multiple entries, each with a
            # different list of genes; for now, arbitrarily choose one of the
            # mappings.
            cpdb <- cpdb[!duplicated(cpdb$external_id), ]

            header_key <- sprintf('%s_ids', sub('-', '_', keytype))

            # convert to an n x 2 mapping of pathway, gene pairs
            cpdb_pathway_list <- apply(cpdb, 1, function(x) {
                cbind(x['pathway'], unlist(strsplit(x[header_key], ',')))
            })
            cpdb <- as.data.frame(do.call('rbind', cpdb_pathway_list))

            colnames(cpdb) <- c('pathway', 'gene')
            rownames(cpdb) <- NULL

            # For now, remove all genes from mapping that are not in our dataset
            #cpdb <- cpdb[cpdb$gene %in% rownames(self$dat),]

            # exclude any pathways with only a single gene (not that interesting..)
            mask <- cpdb$pathway %in% names(which(table(cpdb$pathway) > 1))
            cpdb <- cpdb[mask, ]

            # discard unused factor levels
            cpdb$pathway <- factor(cpdb$pathway)

            # store and return mapping
            self$annotations[['cpdb']] <- cpdb

            cpdb
        }
    )
)
