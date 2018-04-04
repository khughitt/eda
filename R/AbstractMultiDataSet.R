#' An R6 class representing a generic collection of datasets linked by
#' either column or row identifiers.
#'
#' Includes fucntionality and methods that are relevant to all EDA child
#' subclasses.
#'
#' @importFrom R6 R6Class
#' @name AbstractMultiDataSet
#'
NULL

AbstractMultiDataSet <- R6Class("AbstractMultiDataSet",
    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        edat = NULL,

        # AbstractMultiDataSet constructor
        initialize = function(datasets,
                              row_color=NULL, row_color_ds='dat',
                              row_shape=NULL, row_shape_ds='dat',
                              row_label=NULL, row_label_ds='dat',
                              col_color=NULL, col_color_ds='dat',
                              col_shape=NULL, col_shape_ds='dat',
                              col_label=NULL, col_label_ds='dat',
                              color_pal='Set1', title="", ggplot_theme=theme_bw) {

            # drop any empty datasets
            datasets <- datasets[!sapply(datasets, is.null)]
            
            # assign names to datasets list entries if not already present
            if (is.null(names(datasets))) {
                names(datasets) <- paste0('dat', seq_along(datasets))
            } else if (sum(names(datasets) == '') > 0) {
                num_missing <- sum(names(datasets) == '')
                names(datasets)[names(datasets) == ''] <- paste0('dat', 1:num_missing)
            }

            # store datasets as a list of EDADat objects
            for (id in names(datasets)) {
                if (!class(datasets[[id]])[1] == 'EDADat') {
                    # TODO: Guess proper orientation?
                    datasets[[id]] <- EDADat$new(datasets[[id]])
                }
            }
            self$edat <- datasets

            # default variables to use for plot color, shape, and labels when
            # visualizing either columns or rows in the dataset
            private$row_color <- row_color
            private$row_shape <- row_shape
            private$row_label <- row_label

            private$col_color <- col_color
            private$col_shape <- col_shape
            private$col_label <- col_label

            private$row_color_ds <- row_color_ds
            private$row_shape_ds <- row_shape_ds
            private$row_label_ds <- row_label_ds
            private$col_color_ds <- col_color_ds
            private$col_shape_ds <- col_shape_ds
            private$col_label_ds <- col_label_ds

            private$color_pal    <- color_pal
            private$ggplot_theme <- ggplot_theme
            private$title        <- title
        },

        # Clears any cached resuts and performs garbage collection to free
        # up memory.
        clear_cache = function() {
            private$cache <- list()
            invisible(gc())
        },

        # Measure similarity between columns
        cor = function(key=1, method='pearson', ...) {
            private$similarity(self$edat[[key]]$dat, method=method, ...)
        },

        # Applies a filter to rows of the dataset
        #
        # @param mask Logical vector of length equal to the number of rows in
        #      the dataset.
        #
        # @return A filtered version of the original EDADataSet object.
        filter_rows = function(key=1, mask) {
            obj <- private$clone_()
            obj$edat[[key]]$dat <- obj$edat[[key]]$dat[mask,, drop = FALSE] 
            obj
        },

        # Applies a filter to columns of the dataset
        #
        # @param mask Logical vector of length equal to the number of columns
        #     in the dataset.
        #
        # @return A filtered version of the original EDADataSet object.
        filter_cols = function(key=1, mask) {
            obj <- private$clone_()
            obj$edat[[key]]$dat <- obj$edat[[key]]$dat[, mask, drop = FALSE] 
            obj
        },

        # Imputes missing values in the dataset
        #
        # Imputes missing values in the dataset using a specified method
        # and stores the result in-place. Currently only support k-Nearest
        # Neighbors (kNN) method.
        #
        # Note: When using the `knn` method, it may be neccessary to set R
        # the environmental variable `R_MAX_NUM_DLLS` to some value larger
        # than its default of 100 (150 should be sufficient), prior to
        # launching R. This can be set in the ~/.Renviron.
        #
        # @param method Character array specifying imputation method to use
        #     (options: knn)
        impute = function(key=1, method='knn') {
            # Keep track of original dataset class
            dat <- self$edat[[key]]$dat
            cls <- class(dat)

            message(sprintf("Imputing %d missing values...", sum(is.na(dat))))

            # kNN
            if (method == 'knn') {
                imputed <- VIM::kNN(dat)[, 1:ncol(dat)]

                if (cls == 'matrix') {
                    imputed <- as.matrix(imputed)
                }
                rownames(imputed) <- rownames(dat)
            }
            message("Done.")
            self$edat[[key]]$dat <- imputed
        },

        # Subsample dataset
        #
        # Randomly subsamples dataset and returns a new EDADataSet instance.
        # If both `xx_n` and `xx_ratio` arguments are specified, the `xx_n`
        # setting will be used.
        # For removing specific rows and columns, see the `filter_rows` and
        # `filter_cols` functions.
        #
        # @param row_n  Number of rows to randomly select
        # @param col_n  Number of columns to randomly select
        # @param row_ratio Ratio of rows to randomly select
        # @param col_ratio Ratio of columns to randomly select
        #
        # @return An EDADataSet instance
        subsample = function(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL) {
            # clone and subsample dataset
            obj <- private$clone_()
            obj$edat[[1]]$subsample(row_n, col_n, row_ratio, col_ratio)
            obj
        },

        # Summarizes overall characteristics of a dataset
        #
        # @param markdown Logical indicating whether summary output should be
        #   formatted using markdown. (NOTE: Not yet implemented...)
        # @param num_digits Number of digits to display when printing summary
        #   output (default: 2).
        summary = function(key=1, markdown=FALSE, num_digits=2) {
            # collection summary info
            x <- self$edat[[key]]$dat

            # list to store summary info
            info <- list()

            # include trailing zeros when rounding
            round_ <- function(y, digits) {
                formatC(round(y, digits), digits, format = "f")
            }

            # overall dataset
            info[['dat_range']]     <- round_(range(x, na.rm =  TRUE), num_digits)
            info[['dat_num_nas']]   <- sum(is.na(x))
            info[['dat_num_zeros']] <- sum(x == 0)
            info[['dat_quartiles']] <- round_(quantile(x, na.rm = TRUE), num_digits)

            # column types
            info[['col_types']] <- table(apply(x, 2, class))

            # rows
            info[['row_means']]    <- round_(range(apply(x, 1, mean, na.rm = TRUE)), num_digits)
            info[['row_medians']]  <- round_(range(apply(x, 1, median, na.rm = TRUE)), num_digits)
            info[['row_std_devs']] <- round_(range(apply(x, 1, sd, na.rm = TRUE)), num_digits)
            info[['row_outliers']] <- self$detect_row_outliers()

            # columns
            info[['col_means']]    <- round_(range(apply(x, 2, mean, na.rm = TRUE)), num_digits)
            info[['col_medians']]  <- round_(range(apply(x, 2, median, na.rm = TRUE)), num_digits)
            info[['col_std_devs']] <- round_(range(apply(x, 2, sd, na.rm = TRUE)), num_digits)
            info[['col_outliers']] <- self$detect_col_outliers()

            # clip outlier lists if more than a few
            if (length(info[['row_outliers']]) > 5) {
                info[['row_outliers']] <- c(head(info[['row_outliers']], 5), '...')
            }
            if (length(info[['col_outliers']]) > 5) {
                info[['col_outliers']] <- c(head(info[['col_outliers']], 5), '...')
            }

            # row & column correlations
            info[['col_cor_mat']] <- round_(cor(x), num_digits)
            info[['row_cor_mat']] <- round_(cor(t(x)), num_digits)

            diag(info[['col_cor_mat']]) <- NA
            diag(info[['row_cor_mat']]) <- NA

            # display using selected output format
            if (markdown) {
                # TODO
            } else {
                private$print_summary(info)
            }
        },

        ######################################################################
        # plotting methods
        ######################################################################

        # Plot column densities
        #
        # Plots densities for each column in the dataset. This is most
        # useful when you are interested in similarties or differences in
        # distributions across columns for a relatively small number of
        # columns.
        #
        # @param color Variable to color density curves by. If not is
        # specified, uses variable specified at object construction time,
        # or else, uses a separate color for each column.
        #
        # @return ggplot plot instance.
        plot_densities = function(key=1, color=NULL, title="", ...) {
            dat <- setNames(melt(self$edat[[key]]$dat), c('row', 'column', 'val'))
            styles <- private$get_geom_density_styles(color)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color = styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Column densities: %s", private$title)
            }

            # construct density plot
            plt <- ggplot(dat, aes(x = val, group = column, color = column)) +
                geom_density(styles$aes) +
                ggtitle(title) +
                private$ggplot_theme()

            # legend labels
            if (length(styles$labels) > 0) {
                plt <- plt + styles$labels
            }

            # only show legend if there are relatively few groups
            if (length(unique(dat$column)) > 10) {
                plt <- plt + guides(color = FALSE)
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
            cat(sprintf("= %s (n = %d)\n", cls, length(self$edat)))
            cat("=\n")
            for (i in seq_along(self$edat)) {
                # print dataset entry
                ds <- self$edat[[i]]$dat
                cat(sprintf(entry_template, keys[i], class(ds), nrow(ds), ncol(ds)))
            }
            cat("=\n")
            cat("=========================================\n")
        },

        # transpose (out-of-place)
        t = function() {
            obj <- private$clone_()
            obj$transpose()
            obj
        },

        # transpose (in-place)
        transpose = function() {
            # transpose data
            for (id in names(self$edat)) {
               self$edat[[id]]$transpose()
            }

            # swap row and column style parameters
            row_color    <- private$row_color
            row_shape    <- private$row_shape
            row_label    <- private$row_label
            row_color_ds <- private$row_color_ds
            row_shape_ds <- private$row_shape_ds
            row_label_ds <- private$row_label_ds

            private$row_color    <- private$col_color
            private$row_shape    <- private$col_shape
            private$row_label    <- private$col_label
            private$row_color_ds <- private$col_shape_ds
            private$row_shape_ds <- private$col_shape_ds
            private$row_label_ds <- private$col_label_ds

            private$col_color    <- row_color
            private$col_shape    <- row_shape
            private$col_shape_ds <- row_shape_ds
            private$col_label_ds <- row_label_ds
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # private params
        row_color    = NULL,
        row_shape    = NULL,
        row_label    = NULL,
        col_color    = NULL,
        col_shape    = NULL,
        col_label    = NULL,

        row_color_ds = NULL,
        row_shape_ds = NULL,
        row_label_ds = NULL,
        col_color_ds = NULL,
        col_shape_ds = NULL,
        col_label_ds = NULL,

        color_pal    = NULL,
        ggplot_theme = NULL,
        title        = NULL,

        cache        = list(),

        # clone object
        clone_ = function() {
            obj <- self$clone(deep = TRUE)
            obj$clear_cache()
            obj
        },

        # deep_clone function necessary to ensure edat list of R6 classes
        # are copied during cloning
        deep_clone = function(name, value) {
            if (name == "edat") {
                copied <- list()

                for (key in names(value)) {
                    copied[[key]] <- value[[key]]$clone()
                }
                copied
            } else {
                value
            }
        },

        # Generates ggplot aesthetics for bar plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_bar_styles = function(color, key=NULL) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(res, color, key)
        },

        # Generates ggplot aesthetics for density plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_density_styles = function(color, key=NULL) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(res, color, key)
        },

        # Generates ggplot aesthetics for histogram plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_histogram_styles = function(color, key=NULL) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(res, color, key)
        },

        # Generates ggplot aesthetics for a geom_point plot
        #
        # @param color Color variable as passed into plot function call
        # @param shape Shape variable as passed into plot function call
        #
        # @return List of geom_point style information
        get_geom_point_styles = function(color, shape,
                                         color_key=NULL,
                                         shape_key=NULL) {
            # list to store style properties
            res <- list(
                aes = aes(),
                labels = list()
            )
            res <- private$add_color_styles(res, color, color_key)
            res <- private$add_shape_styles(res, shape, shape_key)
            res
        },

        # Determines color-related style information to use for a plot
        #
        # @param styles List of color-related style info
        # @param color Color variable passed into plot function call
        add_color_styles = function(styles, color, key=NULL) {
            color_info <- private$get_color_styles(color, key)

            styles[['color']] <- color_info[['color']]

            # update styles with color info
            if (length(color_info[['aes']]) > 0) {
                styles[['aes']] <- modifyList(color_info[['aes']], styles[['aes']])
                styles[['labels']] <- modifyList(color_info[['labels']], styles[['labels']])
            }
            styles
        },

        # Determines shape-related style information to use for a plot
        #
        # @param styles List of shape-related style info
        # @param shape shape variable passed into plot function call
        add_shape_styles = function(styles, shape, key=NULL) {
            shape_info <- private$get_shape_styles(shape, key)

            styles[['shape']] <- shape_info[['shape']]

            # update styles with shape info
            if (length(shape_info[['aes']]) > 0) {
                styles[['aes']] <- modifyList(shape_info[['aes']], styles[['aes']])
                styles[['labels']] <- modifyList(styles[['labels']], shape_info[['labels']])
            }
            styles
        },

        # Returns a list of color-related style information
        #
        # @param key Name of dataset containing field to use for color
        # @param color Color variable as passed into plot function call
        #
        # @return list List of color-related properties
        get_color_styles = function(color, key) {
            res <- list(
                color  = NULL,
                aes    = aes(),
                labels = list()
            )

            # dataset containing color field
            if (is.null(key)) {
                key <- private$row_color_ds
            }

            # if specified as a function argument, override default color
            if (!is.null(color) && (color != FALSE)) {
                # color variable can either correspond to a column in the
                # dataset itself, or in the
                res[['color']]  <- self$edat[[key]]$dat[, color]
                res[['aes']]    <- aes(color = color)
                res[['labels']] <- labs(color = color)
            } else if (is.null(color) && !is.null(private$row_color)) {
                # otherwise, use object-level default value
                res[['color']]  <- self$edat[[key]]$dat[, private$row_color]
                res[['aes']]    <- aes(color = color)
                res[['labels']] <- labs(color = private$row_color)
            }
            res
        },

        # Returns a list of shape-related style information
        #
        # @param shape Shape variable as passed into plot function call
        #
        # @return list List of shape-related properties
        get_shape_styles = function(shape, key=NULL) {
            res <- list(
                shape  = NULL,
                aes    = aes(),
                labels = list()
            )

            # dataset containing color field
            if (is.null(key)) {
                key <- private$row_shape_ds
            }

            # if specified as a function argument, override default shape
            if (!is.null(shape) && (shape != FALSE)) {
                res[['shape']]  <- self$edat[[key]]$dat[, shape]
                res[['aes']]    <- aes(shape = shape)
                res[['labels']] <- labs(shape = shape)
            } else if (is.null(shape) && !is.null(private$shape)) {
                # otherwise, use object-level default value
                res[['shape']]  <- self$edat[[key]]$dat[, private$shape]
                res[['aes']]    <- aes(shape = shape)
                res[['labels']] <- labs(shape = private$shape)
            }

            res
        },

        # Returns a vector of color codes associated with the specified
        # variable.
        #
        # @param color Variable to use for assigning colors.
        # @param color_pal Color palette to use
        #
        # @return Vector of colors with length equal to the number of columns
        #            in the data.
        get_var_colors = function(color, key, color_pal) {
            # if no variable is specified, use default black for plots
            if (is.null(color)) {
                if (is.null(private$col_color)) {
                    return('black')
                }
                color <- private$col_color
            } else if (color == FALSE) {
                return('black')
            }

            # dataset containing color field
            if (is.null(key)) {
                key <- private$row_color_ds
            }

            # otherwise, assign colors based on the variable specified
            column_var <- as.numeric(factor(self$edat[[key]]$dat[, color]))

            pal <- RColorBrewer::brewer.pal(9, color_pal)
            colors <- colorRampPalette(pal)(min(1E4, length(unique(column_var))))

            # return vector of column color assignments
            colors[as.integer(column_var)]
        },

        # Returns a vector of text labels to use for a plot
        #
        # @param shape Variable to use for assigning shapes.
        #
        # @return Vector of labels with length equal to the number of columns
        #         in the data.
        get_var_labels = function(label, key=NULL) {
            if (is.null(label)) {
                return(colnames(self$dat))
            }

            # dataset containing color field
            if (is.null(key)) {
                key <- private$row_label_ds
            }

            self$edat[[key]]$dat[, label]
        },

        #
        # Supported similarity measures
        #
        # Each of the below functions should accept either one or two matrices
        # or data frames, measures the similarity/dependence either within
        # (single dataset) or across (two datasets) columns in the dataset(s).
        #
        # The result is a matrix of similarity scores (e.g. correlation scores,
        # r^2, etc.) of dimensions:
        #
        # ncol(dat1) x ncol(dat1)  - single dataset
        # ncol(dat1) x nrow(dat2)  - two datasets
        #
        similarity_measures = list(
            #
            # linear model fit
            #
            'lm' = function(dat1, dat2=NULL, ...) {
                # single dataset
                if (is.null(dat2)) {
                    # construct linear model to measure pairwise dependence of
                    # dataset1 columns on dataset2 columns
                    cor_mat <- matrix(0, nrow = ncol(dat1), ncol = ncol(dat1))

                    # for each column in dataset1
                    for (i in 1:ncol(dat1)) {
                        # fit a linear model between the ith column of dataset1
                        # and the jth column of dataset1
                        feature_cor <- function(y) {
                            round(summary(lm(y ~ dat1[, i]))$r.squared * 100, 2)
                        }
                        # for each column in dataset1
                        cor_mat[, i] <- apply(dat1, 2, feature_cor)
                    }
                } else {
                    # two datasets

                    # construct linear model to measure pairwise dependence of
                    # dataset1 columns on dataset2 columns
                    cor_mat <- matrix(0, nrow = ncol(dat1), ncol = ncol(dat2))

                    # for each column in dataset2
                    for (i in 1:ncol(dat2)) {
                        # fit a linear model between the ith column of dat2 and the
                        # jth column of dat1
                        feature_cor <- function(y) {
                            round(summary(lm(y ~ dat2[, i]))$r.squared * 100, 2)
                        }
                        # for each column in dataset1
                        cor_mat[, i] <- apply(dat1, 2, feature_cor)
                    }
                }
                cor_mat
            },
            #
            # Mutual Information
            #
            'mi' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    # two datasets
                    cor_mat <- mpmi::cmi(cbind(dat1, dat2), ...)$bcmi
                    cor_mat[1:ncol(dat1), (ncol(dat1) + 1):nrow(cor_mat)]
                } else {
                    # single dataset
                    mpmi::cmi(dat1, ...)$bcmi
                }
            },
            #
            # Pearson correlation
            #
            'pearson' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    cor_mat <- cor(cbind(dat1, dat2), method = 'pearson', ...)

                    # limit to cross-dataset correlations
                    cor_mat[1:ncol(dat1), (ncol(dat1) + 1):nrow(cor_mat)]
                } else {
                    cor(dat1, method = 'pearson')
                }
            },

            #
            # Spearman correlation
            #
            'spearman' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    cor_mat <- cor(cbind(dat1, dat2), method = 'spearman', ...)

                    # limit to cross-dataset correlations
                    cor_mat[1:ncol(dat1), (ncol(dat1) + 1):nrow(cor_mat)]
                } else {
                    cor(dat1, method = 'spearman')
                }
            }
        ),

        check_input = function() {
            # TODO: Check to make sure all datasets overalp in keys with dat
        },

        # Computes cross-dataset correlation matrix
        #
        # Measures similarity between rows in two datasets which share the
        # same set of column ids.
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use. Supported options:
        #   - pearson  (Pearson correlation)
        #   - spearman (Spearman correlation)
        #   - lm       (Linear model)
        #   - mi       (Mututal information)
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        compute_cross_cor = function(key1='dat', key2=2, method='pearson', ...) {
            # make sure both datasets are both row-oriented
            if (self$edat[[key1]]$orientation != 'rows' || self$edat[[key2]]$orientation != 'rows') {
                stop("Specified datasets should both be row-oriented")
            }

            # TODO: make sure datasets are ordered similarly

            # for metadata, we can also exclude factor fields with all
            # unique values (e.g. alternate identifiers)
            #exclude <- apply(dat2, 1, function(x) {
            #    is.factor(x) && length(unique(x)) == length(x)
            #})

            #if (sum(exclude) > 0) {
            #    message(sprintf("Excluding %d unique factor fields", sum(exclude)))
            #    dat2 <- dat2[!exclude, ]

            # the similarity() method operates on columns so for datasets
            # that share common row id's, we transpose the datasets first
            dat1 <- self$edat[[key1]]$tdat
            dat2 <- self$edat[[key2]]$tdat

            # only operate on shared columns (rows after transposition)
            row_ind <- intersect(rownames(dat1), rownames(dat2))

            dat1 <- dat1[match(rownames(dat1), row_ind),, drop = FALSE]
            dat2 <- dat2[match(rownames(dat2), row_ind),, drop = FALSE]

            # measure similarity between rows in datasets 1 and rows in
            # dataset 2
            cor_mat <- private$similarity(dat1, dat2, method = method, ...)

            cor_mat
        },

        # Plots multidataset correlation heatmap
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param method Correlation method to use (passed to `cor` function)
        #
        plot_cross_cor_heatmap = function(key1='dat', key2=2, method='pearson', interactive=TRUE) {
            # compute cross correlations
            cor_mat <- private$compute_cross_cor(key1, key2, method)

           # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            # TODO: Move this logic to subclass (parent shouldn't know about
            # row / col mdat)
            if ('row_mdata' %in% names(self$edat)) {
                mdata1 <- self$edat[['row_mdata']]$dat
                mask1  <- sapply(mdata1, function(x) max(table(x)) > 1 )
                mdata1 <- mdata1[, mask1, drop = FALSE]

                params[['row_side_colors']] <- mdata1
                params[['subplot_widths']] <- c(0.15, 0.3, 0.55)
            }

            if ('col_mdata' %in% names(self$edat)) {
                mdata1 <- self$edat[['col_mdata']]$dat
                mask2  <- sapply(mdata2, function(x) max(table(x)) > 1 )
                mdata2 <- mdata2[, mask2, drop = FALSE]

                params[['col_side_colors']] <- mdata2
                params[['subplot_heights']] <- c(0.55, 0.3, 0.15)
            }

            # add any additional function arguments
            private$construct_heatmap_plot(params, interactive)
        },

        # Creates a static or interactive heatmap plot
        #
        # @param params A list of plotting parameters
        # @param interactive Logical indicating whether an interactive heatmap
        #     should be generated.
        construct_heatmap_plot = function(params, interactive) {
            # interactive heatmap
            if (interactive) {
                return(do.call(heatmaply::heatmaply, params))
            }

            # convert colside and rowside colors to explicit color values,
            # if present
            if ('col_side_colors' %in% names(params)) {
                params$ColSideColors <- params$col_side_colors

                colors <- c('blue', 'yellow')

                for (col_name in colnames(params$ColSideColors)) {
                    col <- params$ColSideColors[, col_name]
                    params$ColSideColors[, col_name] <- colors[as.numeric(factor(col))]
                }
                params$ColSideColors <- as.matrix(params$ColSideColors)
            }

            if ('row_side_colors' %in% names(params)) {
                params$RowSideColors <- params$row_side_colors

                pal <- RColorBrewer::brewer.pal(8, 'Set1')

                for (col_name in colnames(params$RowSideColors)) {
                    col <- params$RowSideColors[, col_name]
                    colors <- colorRampPalette(pal)(min(1E4, length(unique(col))))
                    params$RowSideColors[, col_name] <- colors[as.numeric(factor(col))]
                }
                params$RowSideColors <- as.matrix(params$RowSideColors)
            }

            # remove irrelevant function arguments
            heatmaply_args <- c('showticklabels', 'subplot_widths', 'subplot_heights',
                                'col_side_colors', 'row_side_colors')
            params <- params[!names(params) %in% heatmaply_args]

            do.call(heatmap.plus::heatmap.plus, params)
        },
        #
        # Measures similarity between columns within or across datasets
        #
        similarity = function(dat1, dat2=NULL, method='pearson', ...) {
            # check to make sure specified similarity measure is valid
            if (!method %in% names(private$similarity_measures)) {
                stop('Unsupported similarity method specified.')
            }

            # drop any columns with zero variance
            var_mask1 <- apply(dat1, 2, function(x) length(table(x)) ) > 1

            if (sum(!var_mask1) > 0) {
                message(sprintf("Excluding %d zero-variance entries from first dataset.", sum(!var_mask1)))
                dat1 <- dat1[, var_mask1, drop = FALSE]
            }

            if (!is.null(dat2)) {
                var_mask2 <- apply(dat2, 2, function(x) length(table(x)) ) > 1
                if (sum(!var_mask2) > 0) {
                    message(sprintf("Excluding %d zero-variance entries from second dataset.", sum(!var_mask2)))
                    dat2 <- dat2[, var_mask2, drop = FALSE]
                }
            }

            # construct similarity matrix
            cor_mat <- private$similarity_measures[[method]](dat1, dat2, ...)

            # fix row and column names
            rownames(cor_mat) <- colnames(dat1)

            if (!is.null(dat2)) {
                colnames(cor_mat) <- colnames(dat2)
            } else {
                colnames(cor_mat) <- colnames(dat1)
            }

            cor_mat
        }
    )
)
