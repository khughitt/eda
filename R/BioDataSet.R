#' An R6 class representing collection of related biological datasets
#'
#' BioEDADataSet is a class for interacting with one or more biological
#' datasets, typically those collected from high-throughput experiments.
#'
#' @section Arguments:
#' - `dataset`: A list of datasets (matrices, data frames, etc.), each of
#'    which shared some column / row identifiers with the first entry in
#'    the list.
#'
#' - `row_color`: Row metadata field to use for coloring rowwise plot elements.
#' - `row_shape`: Row metadata field to use for determine rowwise plot
#'    element shape.
#' - `row_label`: Row metadata field to use when labeling plot points or
#'    other elements.
#' - `col_color`: Column metadata field to use for coloring columnwise plot elements.
#' - `col_shape`: Column metadata field to use for determine columnwise plot
#'    element shape.
#' - `col_label`: Column metadata field to use when labeling plot points or
#'    other elements.
#' - `color_pal`: Color palette to use for relevant plotting methods
#'    (default: `Set1`).
#' - `title`: Text to use as a title or subtitle for plots.
#' - `ggplot_theme`: Default theme to use for ggplot2 plots
#'    (default: `theme_bw`).
#'
#' @section Fields:
#'  - `edat`: List of EDADat dataset instances
#'  - `annotations`: A list of gene, etc. annotations from external sources.
#'
#' @section Methods:
#' - `clear_cache()`: Clears BioEDADataSet cache.
#' - `clone()`: Creates a copy of the BioEDADataSet instance.
#' - `cluster_tsne(num_clusters=10, ...)`: Clusters rows in dataset using a combination
#'    of t-SNE and k-means clustering.
#' - `cpm()`: Performs counts-per-million (CPM) transformation.
#' - `detect_col_outliers(num_sd=2, avg='median', meas='pearson')`:
#'    Measures average pairwise similarities between all columns in the dataset.
#'    Outliers are considered to be those columns who mean similarity to
#'    all other columns is greater than `num_sd` standard deviations from the
#'    average of averages.
#' - `detect_row_outliers(num_sd=2, avg='median', meas='pearson')`:
#'    Measures average pairwise similarities between all rows in the dataset.
#'    Outliers are considered to be those rows who mean similarity to
#'    all other rows is greater than `num_sd` standard deviations from the
#'    average of averages.
#'  - `feature_cor()`: Detects dependencies between column metadata entries
#'    (features) and dataset rows.
#'  - `filter_col_outliers(num_sd=2, avg='median', meas='pearson')`:
#'    Removes column outliers from the dataset. See `detect_col_outliers()`
#'    for details of outlier detection approach.
#'  - `filter_row_outliers(num_sd=2, avg='median', meas='pearson')`:
#'    Removes row outliers from the dataset. See `detect_row_outliers()`
#'    for details of outlier detection approach.
#'  - `filter_cols(mask)`: Accepts a logical vector of length `ncol(obj$dat)`
#'    and returns a new BioEDADataSet instance with only the columns associated
#'    with `TRUE` values in the mask.
#'  - `filter_rows(mask)`: Accepts a logical vector of length `nrow(obj$dat)`
#'    and returns a new BioEDADataSet instance with only the rowsumns associated
#'    with `TRUE` values in the mask.
#'  - `impute(method='knn')`: Imputes missing values in the dataset and stores
#'    the result _in-place_. Currently only k-Nearest Neighbors (kNN)
#'    imputation is supported.
#'  - `log(base=exp(1), offset=0)`: Log-transforms data.
#'  - `log1p()`: Logn(x + 1)-transforms data.
#'  - `log2p()`: Log2(x + 1)-transforms data.
#'  - `pca(...)`: Performs principle component analysis (PCA) on the dataset
#'    and returns a new BioEDADataSet instance of the projected data points.
#'    Any additional arguements specified are passed to the `prcomp()` function.
#'  - `pca_feature_cor(meas='pearson', ...)`: Measures correlation between
#'    dataset features (column metadata fields) and dataset principle
#'    components.
#'  - `plot_cor_heatmap(meas='pearson', interactive=TRUE, ...)`: Plots a
#'    correlation heatmap of the dataset.
#'  - `plot_densities(color=NULL, title="", ...)`: Plots densities for each
#'    column in the dataset.
#'  - `plot_feature_cor(meas='pearson', color_scale=c('green', 'red')`:
#'    Creates a tile plot of projected data / feature correlations. See
#'    `feature_cor()` function.
#'  - `plot_heatmap(x=1, interactive=TRUE, ...)`: Generates a heatmap plot of the
#'    dataset
#'  - `plot_pairwise_column_cors(color=NULL, title="", meas='pearson', mar=c(12,6,4,6))`:
#'    Plot median pairwise column correlations for each variable (column)
#'    in the dataset.
#'  - `plot_pca(pcx=1, pcy=2, scale=FALSE, color=NULL, shape=NULL, title=NULL,
#'         text_labels=FALSE, ...)`:
#'    Generates a two-dimensional PCA plot from the dataset.
#'  - `plot_tsne(color=NULL, shape=NULL, title=NULL, text_labels=FALSE, ...)`:
#'    Generates a two-dimensional t-SNE plot from the dataset.
#'  - `print()`: Prints an overview of the object instance.
#'  - `subsample(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL)`:
#'    Subsamples dataset rows and/or columns.
#'  - `summary(markdown=FALSE, num_digits=2)`: Summarizes overall
#'    characteristics of a dataset.
#'  - `t()`: Transposes dataset rows and columns.
#'  - `tsne(...)`: Performs T-distributed stochastic neighbor embedding (t-SNE)
#'    on the dataset and returns a new BioEDADataSet instance of the projected
#'     data points. Any additional arguements specified are passed to the
#'    `Rtsne()` function.
#'  - `tsne_feature_cor(meas='pearson', ...)`: Measures correlation between
#'    dataset features (column metadata fields) and dataset t-SNE projected
#'    axes.
#'  - `cross_cor(key1=1, key2=2, meas='pearson')`: Computes cross-dataset
#'   correlation matrix between rows in two specified datasets.
#'  - `plot_cross_cor_heatmap(key1=1, key2=2, meas='pearson', interactive=TRUE)`:
#'    Plots multidataset correlation heatmap.
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

    # BioDataSet constructor
    initialize = function(datasets, color_pal='Set1', title="", ggplot_theme=theme_bw) {
      # supported bioconductor data types
      bioc_types <- c('RangedSummarizedExperiment', 'ExpressionSet')

      # create a list of counters to keep track of the number of each object type encountered
      bioc_counters <- as.list(setNames(rep(1, length(bioc_types)), bioc_types))

      # if a single Bioconductor object instance is passed in, wrap it in a list to
      # be consistent
      if (class(datasets) %in% bioc_types) {
        datasets <- list(datasets) 
      }

      for (i in seq_along(datasets)) {
        # get a single dataset entry
        dat <- datasets[[i]]

        # skip anything that isn't a supported bioc object instance
        if (!any(class(dat) %in% c('RangedSummarizedExperiment', 'ExpressionSet'))) {
          next
        }

        # remove original object from input list
        datasets[[i]] <- NULL

        # base axis labels (rowData xid, rowData yid, colData xid, colData yid)
        axis_ids <- c(feat_xid = 'features', feat_yid = 'feature_annotations', 
                      pheno_xid = 'samples', pheno_yid = 'phenotypes')

        # unpack bioconductor object and determine base EDADat keys to use 
        if ('ExpressionSet' %in% class(dat)) {
          edat_keys   <- c(assay = 'exprs', col_dat = 'pdata', row_dat = 'fdata')
          bioc_dat     <- Biobase::exprs(dat)
          bioc_row_dat <- Biobase::fData(dat)
          bioc_col_dat <- Biobase::pData(dat)
        } else if ('RangedSummarizedExperiment' %in% class(dat)) {
          edat_keys    <- c(assay = 'assays', col_dat = 'coldata', row_dat = 'rowdata')
          # TODO: extend support for RSE objects with multiple datasets
          bioc_dat     <- SummarizedExperiment::assays(dat)[[1]]
          bioc_row_dat <- as.data.frame(SummarizedExperiment::rowData(dat))
          bioc_col_dat <- as.data.frame(SummarizedExperiment::colData(dat))
        }

        # if multiple instances of the same object are passed in, append a numeric suffix
        # to differentiate them
        if (bioc_counters[[class(dat)]] > 1) {
          edat_keys <- sprintf('%s_%d', edat_keys, bioc_counters[[class(dat)]])
          axis_ids  <- sprintf('%s_%d', axis_ids, bioc_counters[[class(dat)]])
        }

        # add assay data
        datasets[[edat_keys['assay']]] = EDADat$new(bioc_dat, 
                                                    xid = axis_ids['feat_xid'],
                                                    yid = axis_ids['pheno_xid'])

        # add feature data 
        if (!is.null(bioc_row_dat)) {
          datasets[[edat_keys['row_dat']]] = EDADat$new(bioc_row_dat, 
                                                        xid = axis_ids['feat_xid'], 
                                                        yid = axis_ids['feat_yid'])
        }

        # add phenotype data 
        if (!is.null(bioc_col_dat)) {
          datasets[[edat_keys['col_dat']]] = EDADat$new(bioc_col_dat, 
                                                        xid = axis_ids['pheno_xid'], 
                                                        yid = axis_ids['pheno_yid'])
        }
        # increment bioc obj counter
        bioc_counters[[class(dat)]] <- bioc_counters[[class(dat)]] + 1
      }


      # call parent constructor
      super$initialize(datasets, color_pal, title, ggplot_theme)
    },

    # Computes pathway- or annotation-level statistics for a given dataset and
    # annotation mapping.
    #
    # @param key Name of dataset to analyze.
    # @param annot 
    # @param stat The annotation-level statistic to compute. Can either be a
    #   function (e.g. sum, median, or var), or a string indicating a specific statistic
    #   to compute.
    #
    # Currently supported statistics include:
    #
    #   - `num_above_cutoff`  number of genes w/ values above some cutoff
    #   - `num_below_cutoff`  number of genes w/ values below some cutoff
    #   - `ratio_above_cutoff`  ratio of genes w/ values above some cutoff
    #   - `ratio_below_cutoff`  ratio of genes w/ values below some cutoff
    #   - `num_nonzero`     number of genes with values not equal to zero
    #   - `num_zero`      number of genes with values equal to zero
    #   - `ratio_nonzero`     ratio of genes with values not equal to zero
    #   - `ratio_zero`      ratio of genes with values equal to zero
    #
    aapply = function(key, fun, annot_key, annot_keytype='ensgene', result_key=NULL, fun_args=list(), ...) {
      # annotations are parsed into n x 2 dataframes consisting of
      # with each row containing an (annotation, gene id) pair.
      ANNOT_IND <- 1
      ITEM_IND  <- 2

      # check for valid dataset key
      private$check_key(key)

      # check for valid function
      if (!is.function(fun)) {
        valid_funcs <- c('gsva', names(private$stat_fxns))

        if (!fun %in% valid_funcs) {
          stop("Invalid function specified.")
        }
      }

      # convert numeric keys
      key <- ifelse(is.numeric(key), names(self$edat)[key], key)

      # key form: <old key>_<annot_name>_<fun>
      if (is.null(result_key)) {
        # if a function is provided, require manually specifying result key
        if (is.function(fun)) {
          stop("When a function is specified for 'fun', 'result_key' must also be provided")
        }
        result_key <- paste(c(key, annot_key, fun), collapse='_')
      }

      # check to see if result has already been computed, and if so, move to front and return
      if (result_key %in% names(self$edat)) {
        # move dataset key to top of list and return self-reference
        self$edat <- self$edat[c(result_key, names(self$edat)[names(self$edat) != result_key])]
        return(self)
      }

      # get annotation mapping 
      mapping <- self$load_annotations(annot_key, keytype = self$edat[[key]]$xid)

      # dataset to compute statistics on
      dat <- self$edat[[key]]$dat

      # exclude mapping pairs for items not in target dataset
      mapping <- mapping[mapping[, ITEM_IND] %in% rownames(dat), ]

      # exclude any annotations with less than two entries in the target dataset
      gset_counts <- table(mapping[, ANNOT_IND])
      to_exclude <- names(gset_counts)[gset_counts == 1]

      mapping <- mapping[!mapping[, ANNOT_IND] %in% to_exclude, ]

      message("Computing annotation statistics...")

      #
      # GSVA-based aggregation
      #
      if (!is.function(fun) && fun == 'gsva') {
        # convert gene set mapping dataframe to a list
        gset_list <- split(mapping$gene, mapping$annotation)

        # get any addition function arguments, and set verbose to FALSE
        gsva_args <- list(expr = dat, gset.idx.list = gset_list, verbose = FALSE)
        gsva_args <- modifyList(gsva_args, fun_args)

        # call GSVA
        res <- do.call(GSVA::gsva, gsva_args)
      } else {
        #
        # All other aggegration methods
        #
        if (!is.function(fun)) {
          # built-in aggregation functions
          fun <- private$stat_fxns[[fun]]
        }

        # output data frame
        res <- data.frame()

        # iterate over annotations
        for (annot_key in unique(mapping[, ANNOT_IND])) {
          # get list of genes, etc. associated with the annotation
          annot_items <- mapping[mapping[, ANNOT_IND] == annot_key, ITEM_IND]

          # get data values for relevant genes
          dat_subset <- dat[rownames(dat) %in% annot_items,, drop = FALSE]

          # combine arguments to apply and the aggregation function into a single list
          iter_args <- c(list(X = dat_subset, MARGIN = 2, FUN = fun), fun_args)

          # compute statistic for each column (often, samples) and append to result
          res <- rbind(res, do.call(apply, iter_args))
        }
        rownames(res) <- unique(mapping[, ANNOT_IND])
        colnames(res) <- colnames(dat)
      }

      # fix rownames and remove any duplicates after converting to valid rownames
      rn <- gsub('\\.+', '\\.', make.names(rownames(res)))
      mask <- !duplicated(rn)

      # fix row names and return result; annotation names
      # may change as a result, but this will help to avoid downstream
      # confusion when saving to CSV, etc.
      res <- res[mask, ]
      rownames(res) <- rn[mask]

      # clone BioDataSet instance and add new edat
      obj <- private$clone_()
      obj$add(result_key, EDADat$new(res, xid=annot_key, yid=self$edat[[key]]$yid))

      obj
    },

    # returns annotation mapping using specified gene identifiers
    load_annotations = function(annot, keytype='ensgene', annot_file=NULL, annot_keytype='ensgene',
                                annot_filetype='gmt', annot_field=1, gene_field=2, 
                                exclude_annotations=NULL) {
      #
      # TODO: decide on best place for exclude annotations logic...
      #

      # check if annotations have been previously loaded
      if (annot %in% names(self$annotations)) {
        # if annotations have been loaded for the requested gene identifiers, return mapping
        if (keytype %in% names(self$annotations[[annot]])) {
          return(self$annotations[[annot]][[keytype]])
        } 

        # if annotations have been loaded for a different key type, map identifiers and return
        from <- names(self$annotations[[annot]])[1]
        self$annotations[[annot]][[keytype]] <- private$map_gene_ids(self$annotations[[annot]][[from]], 
                                                                     from, keytype)
        return(self$annotations[[annot]][[keytype]])
      } 

      # 
      if (is.null(annot_file)) {
        stop("Missing required parameter: annot_file")
      }

      # make sure valid filepath provided (TODO: extend support to URLs?)
      if (!file.exists(annot_file)) {
        stop(sprintf("Invalid annotation filepath specified: %s", annot_file))
      }

      message(sprintf("Loading %s annotations...", annot))

      # load annotations from file
      if (annot_filetype == 'gmt') {
        mapping <- private$load_gmt(annot_file) 
      } else if (annot_filetype == 'tsv') {
        mapping <- private$load_tsv(annot_file, annot_field, gene_field, exclude_annotations)
      }

      # store returned annotation mapping
      if (!annot %in% names(self$annotations)) {
        self$annotations[[annot]] <- list()
      }
      self$annotations[[annot]][[annot_keytype]] <- mapping

      # if dataset keytype differs from annotation keytype, map ids
      if (keytype != annot_keytype) {
        self$annotations[[annot]][[keytype]] <- private$map_gene_ids(mapping, annot_keytype, keytype)
      }

      return(self$annotations[[annot]][[keytype]])
    },

    # Log2 transforms data (adding 1 to ensure finite results).
    #
    # @return Log2-transformed version of the expression data.
    log2p = function(key = 1) {
      self$log(key = key, 2, offset = 1)
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
    plot_libsizes = function(key=1, color=NULL, title=NULL, geom='hist') {
      # create a data frame of library sizes
      dat <- data.frame(libsize = colSums(self$edat[[key]]$dat))

      # determine plot styles to use
      if (geom == 'hist') {
        styles <- private$get_geom_histogram_styles(key, color)
      } else {
        styles <- private$get_geom_bar_styles(key, color)
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

      # determine length to use for justification
      key_format <- sprintf("%%%ds", max(nchar(keys)) + 1)

      # entry output string
      entry_template <- sprintf("= %s : %%s (%%d x %%d) %%s\n", key_format)

      cat("=============================================\n")
      cat("=\n")
      cat(sprintf("= %s (n=%d)\n", cls, length(self$edat)))
      cat("=\n")
      for (i in seq_along(self$edat)) {
        ds <- self$edat[[i]]$dat

        # print dataset entry
        missing_flag <- ifelse(sum(is.na(ds)) > 0, '[M]', '')
        cat(sprintf(entry_template, keys[i], class(ds), nrow(ds), ncol(ds), missing_flag))
      }

      cat("=\n")

      if (length(self$annotations) > 0) {
        cat("= Annotations:\n")
        cat("=\n")
        for (annot_name in names(self$annotations)) {
          for (annot_keytype in names(self$annotations[[annot_name]])) {
            annot <- self$annotations[[annot_name]][[annot_keytype]]

            cat(sprintf("= - %s (%d annotations, %d genes)\n", annot_name,
                  length(unique(annot[, 1])),
                  length(unique(annot[, 2]))))
          }
        }
        cat("=\n")
      }
      cat("=============================================\n")
    },

    # Plots multidataset correlation heatmap
    #
    # @param key1 Numeric or character index of first dataset to use
    # @param key2 Numeric or character index of second dataset to use
    # @param meas Correlation meas to use
    #
    plot_cross_cor_heatmap = function(key1=1, key2=2, meas='pearson', new_key=NULL, interactive=TRUE) {
      super$plot_cross_cor_heatmap(key1, key2, meas, new_key, interactive)
    }
  ),

  # ------------------------------------------------------------------------
  # private
  # ------------------------------------------------------------------------
  private = list(
    # Verifies that input data is an acceptable format.
    check_input = function(dat) {
      if (!is.matrix(dat) && (class(dat) != 'ExpressionSet')) {
        stop("Invalid input for BioEDADataSet: dat must be a matrix or ExpressionSet.")
      }
    },

    # load gene set gmt file
    load_gmt = function(annot_file=NULL) {
      # check for valid annotation filepath
      if (is.null(annot_file)) {
        stop("Filepath to annotation gmt file must be provided.")
      } else if (!file.exists(annot_file)) {
        stop("Invalid filepath for annotations specified.")
      }

      mapping <- private$parse_gmt(annot_file)
      
      colnames(mapping) <- c('annotation', 'gene')
      rownames(mapping) <- NULL

      # exclude any annotations with only a single gene (not interesting...)
      mask <- mapping$annotation %in% names(which(table(mapping$annotation) > 1))
      mapping <- mapping[mask, ]

      # discard unused factor levels and return result
      mapping$annotation <- factor(mapping$annotation)

      mapping
    },

    load_tsv = function(infile, annot_field, gene_field, exclude_annotations=NULL) {
      # load annotation tab-delimited file
      mapping <- read.delim(infile, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
      mapping <- mapping[, c(annot_field, gene_field)]

      # for some annotation sources (e.g. CPDB), there are multiple entries for the same annotation
      # name, each with a different list of genes; for now, arbitrarily choose one of the mappings
      # to use
      mapping <- mapping[!duplicated(mapping[, annot_field]), ]

      # exclude specific annotations, if requested
      if (!is.null(exclude_annotations)) {
        mapping <- mapping[!mapping[, annot_field] %in% exclude_annotations, ]
      }

      # convert to an n x 2 mapping of (annotation, gene) pairs
      mapping_list <- apply(mapping, 1, function(x) {
        cbind(x[annot_field], unlist(strsplit(x[gene_field], ',')))
      })
      mapping <- as.data.frame(do.call('rbind', mapping_list))

      colnames(mapping) <- c('annotation', 'gene')
      rownames(mapping) <- NULL

      # exclude any annotations with only a single gene (not that interesting..)
      mask <- mapping$annotation %in% names(which(table(mapping$annotation) > 1))
      mapping <- mapping[mask, ]

      # discard unused factor levels
      mapping$annotation <- factor(mapping$annotation)

      mapping
    },

    # maps gene identifiers for a specified annotation mapping
    map_gene_ids = function(annot_mapping, from, to) {
      # annotation mapping indices
      GID_IND  <- 2

      # map identifier names, if needed
      if (from == 'entrez-gene') {
        from <- 'entrez'
      } else if (from == 'hgnc-symbol') {
        from <- 'symbol'
      }

      if (to == 'entrez-gene') {
        to <- 'entrez'
      } else if (to == 'hgnc-symbol') {
        to <- 'symbol'
      }

      # mapping gene identifiers and return result
      ind <- match(annot_mapping[, GID_IND], annotables::grch37[, from, drop = TRUE])
      annot_mapping$gene <- annotables::grch37[, to, drop = TRUE][ind]

      annot_mapping
    },

    parse_gmt = function(infile) {
      # determine maximum number columns
      max_cols <- max(count.fields(infile))

      # read in table, filling empty cells with NA's
      gmt <- read.delim(infile, sep = '\t', header = FALSE, 
                        col.names = paste0("gene_", seq_len(max_cols)), 
                        fill = TRUE, stringsAsFactors = FALSE)

      # fix column names
      colnames(gmt)[1:2] <- c('annotation', 'source')
      
      # convert to an n x 2 (annotation, gene) mapping
      res <- do.call('rbind', apply(gmt, 1, function(entry) {
        # gene gene id's starting in third column
        gids <- entry[3:length(entry)]
        gids <- gids[!is.na(gids)]

        data.frame(annotation = entry[1], gene = gids, row.names=NULL)
      }))

      # TODO: Don't assume entrez ids!!!
      res$gene <- as.numeric(as.character(res$gene))

      # return dataframe
      res
    }
  )
)
