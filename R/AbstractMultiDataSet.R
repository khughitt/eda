#' An R6 class representing a generic collection of datasets linked by
#' either column or row identifiers.
#'
#' Includes fucntionality and methods that are relevant to all EDA child
#' subclasses.
#'
#' @import ggplot2
#' @import viridis
#' @importFrom R6 R6Class
#' @importFrom reshape2 melt
#' @importFrom NMF aheatmap nmf
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
        initialize = function(datasets, color_pal='Set1', title="", ggplot_theme=theme_bw) {

            # if a single dataset was provided, wrap as a list for consistency
            if (class(datasets) != 'list') {
                datasets <- list(dat = datasets)
            }

            # drop any empty datasets
            datasets <- datasets[!sapply(datasets, is.null)]

            # assign names to datasets if not provided
            if (is.null(names(datasets))) {
                names(datasets) <- paste0('dat', seq_along(datasets))
            } else if (sum(names(datasets) == '') > 0) {
                num_missing <- sum(names(datasets) == '')
                names(datasets)[names(datasets) == ''] <- paste0('dat', 1:num_missing)
            }

            # store datasets as a list of EDADat objects
            for (id in names(datasets)) {
                if (!class(datasets[[id]])[1] == 'EDADat') {
                    # by default, do not assume that any datasets have
                    # shared axes; this must be explicitly specified by the user
                    datasets[[id]] <- EDADat$new(datasets[[id]],
                                                 xid=paste0(id, '_x'),
                                                 yid=paste0(id, '_y'))
                }
            }
            self$edat <- datasets

            # check to make sure valid edat keys specified
            keys <- names(self$edat)

            for (key in keys) {
                if (!is.null(self$edat[[key]]$row_edat) && (!self$edat[[key]]$row_edat %in% keys)) {
                    stop(sprintf("Invalid row edat key specified for dataset: %s", key))
                }
                if (!is.null(self$edat[[key]]$col_edat) && (!self$edat[[key]]$col_edat %in% keys)) {
                    stop(sprintf("Invalid column edat key specified for dataset: %s", key))
                }
            }

            # check to make sure datasets which are specified as sharing an axis
            # have at least some row/column id's in common
            for (k1 in keys) {
                e1 <- self$edat[[k1]]
                for (k2 in keys) {
                    e2 <- self$edat[[k2]]
                    # check for shared axis identifiers
                    # TODO: clean up / move to separate function
                    if ((e1$xid == e2$xid && length(intersect(rownames(e1$dat), rownames(e2$dat))) == 0) ||
                        (e1$xid == e2$yid && length(intersect(rownames(e1$dat), colnames(e2$dat))) == 0) ||
                        (e1$yid == e2$xid && length(intersect(colnames(e1$dat), rownames(e2$dat))) == 0) ||
                        (e1$yid == e2$yid && length(intersect(colnames(e1$dat), colnames(e2$dat))) == 0)) {
                        stop(sprintf("Datasets '%s' and '%s' specified as having a common axis, but no shared id's found.", k1, k2))
                    }
                } 
            }

            # share styles
            private$color_pal    <- color_pal
            private$ggplot_theme <- ggplot_theme
            private$title        <- title
        },

        # Adds a new dataset to front of datasets list, in-place
        add = function (key, edat) {
            edat_names <- c(key, names(self$edat))
            self$edat <- c(edat, self$edat)
            names(self$edat) <- edat_names
        },

        # update edat data and move to front of the edat list
        update = function (key, dat) {
            edat <- self$edat[[key]]
            edat$dat <- dat

            # convert numeric key to string
            key <- ifelse(is.numeric(key), names(self$edat)[key], key)

            self$edat[[key]] <- NULL
            self$edat <- c(edat, self$edat)
            names(self$edat)[1] <- key

            # remove any unlinked datasets
            private$remove_unlinked(key)
        },

        # Clears any cached resuts and performs garbage collection to free
        # up memory.
        clear_cache = function() {
            private$cache <- list()
            invisible(gc())
        },

        # cluster dataset and return new object instance which including the 
        # results
        cluster = function(key=1, method='kmeans', ...) {
            # check to make sure specified clustering method is valid
            if (!method %in% names(private$cluster_methods)) {
                stop('Unsupported clustering method specified.')
            }

            # convert key to string name, if not already provided as such
            if (is.numeric(key)) {
                key <- names(self$edat)[key]
            }

            # key for cluster result dataset
            cluster_key <- sprintf("clusters-%s", key)

            # perform clustering
            clusters <- private$cluster_methods[[method]](self$edat[[key]]$dat, ...)

            # convert result to an n x 1 mapping matrix
            mat <- as.matrix(clusters)
            rownames(mat) <- rownames(self$edat[[key]])

            # if cluster table doesn't exist, create it
            if (!cluster_key %in% names(self$edat)) {
                # crate a new EDADat instance
                self$edat[[cluster_key]] <- EDADat$new(
                    mat,
                    xid = self$edat[[key]]$xid,
                    yid = cluster_key
                )
            } else {
                # otherwise, add new column to clusters table
                self$edat[[cluster_key]]$dat <- cbind(self$edat[[cluster_key]]$dat, mat)
            }

            clusters
        },

        # cluster dataset and apply function to resulting clusters
        capply = function(key, method='kmeans', fun=median, fun_args=list(), 
                          edat_suffix='auto', ...) {
            # check for valid dataset key
            private$check_key(key)

            # determine function to use
            if (!is.function(fun)) {
                if (fun %in% names(private$stat_fxns)) {
                    fun <- private$stat_fxns[[fun]]
                } else {
                    stop("Invalid function specified.")
                }
            }

            # output data frame
            res <- data.frame()

            # check to see if clustering has been performed, and if not,
            # compute clustering
            clusters <- self$cluster(key = key, method = method, ...)

            # dataset to compute statistics on
            dat <- self$edat[[key]]$dat

            message("Computing cluster statistics...")

            # iterate over clusters
            cluster_ids <- sort(unique(clusters))

            for (cluster in cluster_ids) {
                # data for cluster members
                dat_subset <- dat[clusters == cluster,, drop = FALSE]

                # compute statistic for each column and append to result
                stats <- do.call(apply, c(list(X=dat_subset, MARGIN=2, FUN=fun), fun_args))
                res <- rbind(res, stats)
            }

            # fix column and row names and return result
            rownames(res) <- cluster_ids
            colnames(res) <- colnames(dat)

            # convert numeric keys
            key <- ifelse(is.numeric(key), names(self$edat)[key], key)

            # key form: <old key>_<cluster_method>_<edat_suffix>
            if (edat_suffix == 'auto') {
                edat_suffix <- stringi::stri_rand_strings(n=1, length=6, pattern="[A-Za-z0-9]")
            }
            res_key <- paste(c(key, method, edat_suffix), collapse='_')

            # clone dataset instance and add new edat
            obj <- private$clone_()
            obj$add(res_key, EDADat$new(res, xid=method, yid=self$edat[[key]]$yid))

            obj
        },

        # Measure similarity between columns
        cor = function(key=1, meas='pearson', ...) {
            private$similarity(self$edat[[key]]$dat, meas=meas, ...)
        },

        # Applies a filter to rows of the dataset
        #
        # @param mask Logical vector of length equal to the number of rows in
        #      the dataset.
        #
        # @return A filtered version of the original EDADataSet object.
        filter_rows = function(key=1, mask) {
            obj <- private$clone_()
            obj$update(key, obj$edat[[key]]$dat[mask,, drop = FALSE])
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
            obj$update(key, obj$edat[[key]]$dat[, mask, drop = FALSE])
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
        # @param key    Dataset key
        # @param row_n  Number of rows to randomly select
        # @param col_n  Number of columns to randomly select
        # @param row_ratio Ratio of rows to randomly select
        # @param col_ratio Ratio of columns to randomly select
        #
        # @return An EDADataSet instance
        subsample = function(key=1, row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL) {
            # clone and subsample dataset
            obj <- private$clone_()
            obj$edat[[key]]$subsample(row_n, col_n, row_ratio, col_ratio)

            # move to front of edat list
            edat <- obj$edat[[key]]
            obj$edat[[key]] <- NULL
            obj$edat <- c(edat, obj$edat)
            names(obj$edat)[1] <- key

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

        # returns a new object instance with only the top N rows / columns
        # of specified dataset
        bottom_n = function(key=1, num_rows=NULL, num_cols=NULL, fxn=sum) {
            dat <- self$edat[[key]]$dat

            # masks to keep track of which rows / columns to keep
            row_mask <- rep(TRUE, nrow(dat))
            col_mask <- rep(TRUE, ncol(dat))

            # get top n rows
            if (!is.null(num_rows)) {
                row_scores <- apply(dat, 1, fxn)
                row_mask <- rank(row_scores) <= num_rows
            }

            # get top n columns
            if (!is.null(num_cols)) {
                col_scores <- apply(dat, 2, fxn)
                col_mask <- rank(col_scores) <= num_cols
            }

            obj <- private$clone_()
            obj$update(key, dat[row_mask, col_mask])
            obj
        },

        # returns a new object instance with only the top N rows / columns
        # of specified dataset
        top_n = function(key=1, num_rows=NULL, num_cols=NULL, fxn=sum) {
            dat <- self$edat[[key]]$dat

            # masks to keep track of which rows / columns to keep
            row_mask <- rep(TRUE, nrow(dat))
            col_mask <- rep(TRUE, ncol(dat))

            # get top n rows
            if (!is.null(num_rows)) {
                row_scores <- apply(dat, 1, fxn)
                row_mask <- rank(-row_scores) <= num_rows
            }

            # get top n columns
            if (!is.null(num_cols)) {
                col_scores <- apply(dat, 2, fxn)
                col_mask <- rank(-col_scores) <= num_cols
            }

            obj <- private$clone_()
            obj$update(key, dat[row_mask, col_mask])
            obj
        },

        # filters dataset to include both the rows and columns associated with
        # each of the top and/or bottom N values in the dataset
        extreme_n = function(key=1, top=NULL, bottom=NULL, unique=FALSE) {
            dat <- self$edat[[key]]$dat

            # check limits
            min_dim <- min(dim(dat))

            if (top > min_dim || bottom > min_dim) {
                stop("Limits must be at least as small as the smallest data dimension.")
            }

            # vector of rows / columns to keep
            row_ind <- c()
            col_ind <- c()

            # indices of top N values
            if (!is.null(top)) {
                # require <top> unique rows/columns
                if (unique) {
                    # iterate over max values until <top> unique rows and columns
                    # have been added
                    while (length(row_ind) < top || length(col_ind) < top) {
                        # get row and column indices for next highest value;
                        # in cases of ties, just use the first match
                        ind <- head(which(dat == max(dat, na.rm=TRUE), arr.ind=TRUE), 1)

                        rowi <- ind[, 1]
                        coli <- ind[, 2]

                        # set column / row to NA
                        dat[rowi, ] <- NA
                        dat[, coli] <- NA

                        if (length(row_ind) < top) {
                            row_ind <- c(row_ind, rowi)
                        }
                        if (length(col_ind) < top) {
                            col_ind <- c(col_ind, coli)
                        }
                    }
                } else {
                    # non-unique rows / columns
                    ind <- which(dat >= sort(dat, decreasing=TRUE)[top], arr.ind=TRUE)
                    row_ind <- ind[, 'row']
                    col_ind <- ind[, 'col']
                }
            }

            # indices of bottom N values
            if (!is.null(bottom)) {
                if (unique) {
                    # iterate until either all rows/cols have been added, or
                    # the requested number of bottom hits have been included
                    while (length(row_ind) < min(nrow(dat), (top + bottom)) || 
                           length(col_ind) < min(ncol(dat), (top + bottom))) {
                        ind <- head(which(dat == min(dat, na.rm=TRUE), arr.ind=TRUE), 1)

                        rowi <- ind[, 1]
                        coli <- ind[, 2]

                        # set column / row to NA
                        dat[rowi, ] <- NA
                        dat[, coli] <- NA

                        if (length(row_ind) < (top + bottom)) {
                            row_ind <- c(row_ind, rowi)
                        }
                        if (length(col_ind) < (top + bottom)) {
                            col_ind <- c(col_ind, coli)
                        }
                    }
                } else {
                    # non-unique rows / columns
                    ind <- which(dat <= sort(dat, decreasing=FALSE)[bottom], arr.ind=TRUE)
                    row_ind <- unique(c(row_ind, ind[, 'row']))
                    col_ind <- unique(c(col_ind, ind[, 'col']))
                }
            }

            # replace original dataset and move to front of edat list
            obj <- private$clone_()
            obj$update(key, self$edat[[key]]$dat[row_ind, col_ind])

            obj
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
        plot_densities = function(key=1, color_var=NULL, color_key=NULL, title="", ...) {
            dat <- setNames(melt(self$edat[[key]]$dat), c('row', 'column', 'val'))
            styles <- private$get_geom_density_styles(key, color_var, color_key)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color_var = styles$color)
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

        # replaces a dataset with a new one, in-place
        replace = function(old_key, new_key, dat, xid=NULL, yid=NULL) {
            # Assume original dataset row and column ids, unless specified
            orig_xid <- self$edat[[old_key]]$xid
            orig_yid <- self$edat[[old_key]]$yid

            # Remove any additional datasets which use the same
            # identifiers as the target dataset
            for (ds in names(self$edat)) {
                if (self$edat[[ds]]$xid == orig_xid || self$edat[[ds]]$yid == orig_yid) {
                    self$edat[[ds]] <- NULL
                }
            }

            # determine new row and column ids
            xid <- ifelse(is.null(xid), orig_xid, xid)
            yid <- ifelse(is.null(yid), orig_yid, yid)

            # replace original dataset with annotation stat result matrix
            self$edat[[old_key]] <- NULL
            self$edat[[new_key]] <- EDADat$new(dat, xid=xid, yid=yid)
        },

        # transpose (out-of-place)
        t = function(key=NULL) {
            obj <- private$clone_()
            obj$transpose(key)
            obj
        },

        # transpose (in-place)
        transpose = function(key=NULL) {
            # if key provided, transpose specific dataset
            if (!is.null(key)) {
                self$edat[[key]]$transpose()
            } else {
                # otherwise, transpose all datasets
                for (id in names(self$edat)) {
                self$edat[[id]]$transpose()
                }
            }
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # private properties
        color_pal    = NULL,
        ggplot_theme = NULL,
        title        = NULL,

        cache        = list(),

        # Helper functions for computing various statistics; used by the
        # `AbstractMultiDataSet$capply` and `BioDataSet$aapply`
        # method.
        stat_fxns = list(
            'num_nonzero' = function(x) {
                sum(x != 0)
            },
            'num_zero' = function(x) {
                sum(x == 0)
            },
            'num_above_cutoff' = function(x, cutoff=0) {
                sum(x > cutoff)
            },
            'num_below_cutoff' = function(x, cutoff=0) {
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
            'ratio_below_cutoff' = function(x, cutoff=0) {
                sum(x < cutoff) / length(x)
            }
        ),

        # check for valid dataset key
        check_key = function(key) {
            # check for valid dataset key
            if (!key %in% c(1:length(self$edat), names(self$edat))) {
                stop(sprintf("Invalid dataset specified: %s", key))
            }
        },

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
        get_geom_bar_styles = function(key, color_var, color_key) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(key, res, color_var, color_key)
        },

        # Generates ggplot aesthetics for density plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_density_styles = function(key, color_var, color_key) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(key, res, color_var, color_key)
        },

        # Generates ggplot aesthetics for histogram plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_histogram_styles = function(key, color_var, color_key) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(key, res, color_var, color_key)
        },

        # Generates ggplot aesthetics for a geom_point plot
        #
        # @param color Color variable as passed into plot function call
        # @param shape Shape variable as passed into plot function call
        #
        # @return List of geom_point style information
        get_geom_point_styles = function(key,
                                         color_var, color_key,
                                         shape_var, shape_key) {
            # list to store style properties
            res <- list(
                aes = aes(),
                labels = list()
            )
            res <- private$add_color_styles(key, res, color_var, color_key)
            res <- private$add_shape_styles(key, res, shape_var, shape_key)
            res
        },

        # Determines color-related style information to use for a plot
        #
        # @param styles List of color-related style info
        # @param color Color variable passed into plot function call
        add_color_styles = function(key, styles, color_var, color_key) {
            color_vec <- private$get_color_vector(key, color_var, color_key) 

            # if no color variable specified, or it is disabled, return styles as-is
            if (is.null(color_vec)) {
                return(styles)
            }

            styles[['color']] <- color_vec

            # update styles with color info
            styles[['aes']]    <- modifyList(aes(color = color_var),  styles[['aes']])
            styles[['labels']] <- modifyList(labs(color = color_var), styles[['labels']])

            styles
        },

        # returns a vector of color assignments to use for plotting
        get_color_vector = function(key, color_var, color_key, target='rows') {
            # determine color variable and data source to use
            if (target == 'rows') {
                if (is.null(color_var)) {
                    color_var <- self$edat[[key]]$row_color
                }
                if (is.null(color_key)) {
                    color_key <- self$edat[[key]]$row_edat
                }
            } else {
                if (is.null(color_var)) {
                    color_var <- self$edat[[key]]$col_color
                }
                if (is.null(color_key)) {
                    color_key <- self$edat[[key]]$col_edat
                }
            }

            # if no color variable specified, or disabled, stop here
            if (is.null(color_var) || color_var == FALSE) {
                return(NULL)
            }

            # convert numeric keys
            key <- ifelse(is.numeric(key), names(self$edat)[key], key)

            # check to make sure valid style key specified
            if (!color_key %in% names(self$edat)) {
                stop(sprintf("Invalid style source key specified: %s", color_key))
            }

            # otherwise, retrieve 1d vector to use for color assignment
            if (target == 'rows') {
                matched_axis <- self$edat[[key]]$xid
            } else {
                matched_axis <- self$edat[[key]]$yid
            }
            vals <- self$edat[[color_key]]$get(matched_axis, color_var, other_axis=TRUE)

            # drop elements not found in target dataset and store result
            matched_ids <- self$edat[[color_key]]$get(matched_axis, other_axis=TRUE)

            # return color vector
            if (target == 'rows') {
                axis_names <- rownames(self$edat[[key]]$dat)
            } else {
                axis_names <- colnames(self$edat[[key]]$dat)
            }
            vals[match(axis_names, matched_ids)]
        },

        # Determines shape-related style information to use for a plot
        #
        # @param styles List of shape-related style info
        # @param shape shape variable passed into plot function call
        add_shape_styles = function(key, styles, shape_var, shape_key) {
            if (is.null(shape_var)) {
                shape_var <- self$edat[[key]]$row_shape
            }
            if (is.null(shape_key)) {
                shape_key <- self$edat[[key]]$row_edat
            }

            # if no shape variable specified, or disabled, return styles as-is
            if (is.null(shape_var) || shape_var == FALSE) {
                return(styles)
            }

            # otherwise, retrieve 1d vector to use for shape assignment
            matched_axis <- self$edat[[key]]$xid
            vals <- self$edat[[shape_key]]$get(matched_axis, shape_var, other_axis=TRUE)

            # drop elements not found in target dataset and store result
            matched_ids <- self$edat[[shape_key]]$get(matched_axis, other_axis=TRUE)

            styles[['shape']] <- vals[match(rownames(self$edat[[key]]$dat), matched_ids)]

            # update styles with shape info
            styles[['aes']]    <- modifyList(aes(shape = shape_var),  styles[['aes']])
            styles[['labels']] <- modifyList(labs(shape = shape_var), styles[['labels']])

            styles
        },

        # Returns a vector of color codes associated with the specified
        # variable.
        #
        # @param color Variable to use for assigning colors.
        # @param color_pal Color palette to use
        #
        # @return Vector of colors with length equal to the number of columns
        #            in the data.
        get_row_colors = function(key, color_var=NULL, color_key=NULL) {
            # determine color dataset and variable name
            if (is.null(color_var)) {
                color_var <- self$edat[[key]]$row_color
            }
            if (is.null(color_key)) {
                color_key <- self$edat[[key]]$row_edat
            }

            # if no variable is specified, use default black for plots
            if (is.null(color_var) || color == FALSE) {
                return('black')
            }

            # if variable specified, but no edat, use self
            if (is.null(color_key)) {
                color_key <- key
            }

            # otherwise, assign colors based on the variable specified
            color_vector <- self$edat[[color_key]]$get(self$edat[[key]]$xid, color_var, other_axis=TRUE)
            color_vector <- as.numeric(factor(color_vector))

            pal <- RColorBrewer::brewer.pal(9, private$color_pal)
            colors <- colorRampPalette(pal)(min(1E4, length(unique(color_vector))))

            # return vector of column color assignments
            colors[as.integer(color_vector)]
        },

        # Returns a vector of text labels to use when plotting row elements of
        # a dataset.
        #
        # @return Vector of labels with length equal to the number of columns
        #         in the data.
        get_col_labels = function(key, label_var=NULL, label_key=NULL) {
            # target dataset edat
            edat <- self$edat[[key]]

            # determine key and variable name associated with labels, if specified
            if (is.null(label_var)) {
                label_var <- edat$col_label 
            }
            if (is.null(label_key)) {
                label_key <- edat$col_edat
            }

            # colnames of target dataset
            cnames <- colnames(edat$dat)

            # if none provided, default to using colnames
            if (is.null(label_var) || label_var == FALSE) {
                return(cnames)
            }

            # if variable specified, but no edat, use self
            if (is.null(label_key)) {
                label_key <- key
            }

            # otherwise, generate a label mapping and determine labels to use
            mapping <- data.frame(
                id    = self$edat[[label_key]]$get(edat$yid, other_axis = TRUE),
                label = self$edat[[label_key]]$get(edat$yid, label_var, other_axis = TRUE)
            )

            # return labels in proper order
            mapping$label[match(cnames, mapping$id)]
        },


        # Returns a vector of color codes associated with the specified
        # variable.
        #
        # @param color Variable to use for assigning colors.
        # @param color_pal Color palette to use
        #
        # @return Vector of colors with length equal to the number of columns
        #            in the data.
        get_col_colors = function(key, color_var=NULL, color_key=NULL) {
            # determine color dataset and variable name
            if (is.null(color_var)) {
                color_var <- self$edat[[key]]$col_color
            }
            if (is.null(color_key)) {
                color_key <- self$edat[[key]]$col_edat
            }

            # if no variable is specified, use default black for plots
            if (is.null(color_var) || color == FALSE) {
                return('black')
            }

            # if variable specified, but no edat, use self
            if (is.null(color_key)) {
                color_key <- key
            }

            # otherwise, assign colors based on the variable specified
            color_vector <- self$edat[[color_key]]$get(self$edat[[key]]$yid, color_var, other_axis=TRUE)
            color_vector <- as.numeric(factor(color_vector))

            pal <- RColorBrewer::brewer.pal(9, private$color_pal)
            colors <- colorRampPalette(pal)(min(1E4, length(unique(color_vector))))

            # return vector of column color assignments
            colors[as.integer(color_vector)]
        },

        # Returns a vector of text labels to use when plotting row elements of
        # a dataset.
        #
        # @return Vector of labels with length equal to the number of columns
        #         in the data.
        get_row_labels = function(key, label_var=NULL, label_key=NULL) {
            # target dataset edat
            edat <- self$edat[[key]]

            # determine key and variable name associated with labels, if specified
            if (is.null(label_var)) {
                label_var <- edat$row_label 
            }
            if (is.null(label_key)) {
                label_key <- edat$row_edat
            }

            # rownames of target dataset
            rnames <- rownames(edat$dat)


            # if none provided, default to using rownames
            if (is.null(label_var) || label_var == FALSE) {
                return(rnames)
            }

            # if variable specified, but no edat, use self
            if (is.null(label_key)) {
                label_key <- key
            }

            # otherwise, generate a label mapping and determine labels to use
            mapping <- data.frame(
                id    = self$edat[[label_key]]$get(edat$xid, other_axis = TRUE),
                label = self$edat[[label_key]]$get(edat$xid, label_var, other_axis = TRUE)
            )

            # return labels in proper order
            mapping$label[match(rnames, mapping$id)]
        },

        #
        # Supported cluster methods
        #
        cluster_methods = list(
            'kmeans' = function(dat, k, ...) {
                factor(as.numeric(kmeans(dat, k, ...)$cluster))
            },

            'hclust' = function(dat, cor_method='cor', 
                                hclust_method='average', 
                                cutree_method='cutree', k=NULL, h=NULL, ...) {
                # hierachical clustering
                # TODO: make similarity, etc. functions currently stored as
                # private lists acessible to one another (store directly as
                # methods in private, or make available somewhere else in eda
                # package namespace)
                #cor_mat <- private$similarity(dat, meas=cor_method)
                cor_mat <- tryCatch({
                    get(cor_method)(t(dat))
                }, warning = function(w) {
                    message("Warning encountered during correlation matrix generation:")
                    stop(w$message)
                })
                
                hc <- flashClust::flashClust(as.dist(1 - abs(cor_mat)), method=hclust_method)

                # tree cut
                cutree_fxn <- get(cutree_method)
                factor(as.numeric(cutree_fxn(hc, k=k, h=h, ...)))
            },

            # k-means clustering of t-SNE projected data
            #
            # k     Number of clusters to detect (default: 10)
            # ...   Additional arguments passed to Rtsne function
            #
            # return Vector of cluster assignments with length equal to the
            #     number of rows in the dataset.
            'tsne-kmeans' = function(dat, k=10, ...) {
                tsne <- Rtsne::Rtsne(dat, ...)
                res <- setNames(as.data.frame(tsne$Y), c('x', 'y'))
                kmeans(res, k)$cluster
            }
        ),

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
            # linear model fit r^2
            #
            # Computes a pairwise matrix of linear model r^2 values; useful,
            # for example, when measuring relationship between a numeric matrix
            # and a categorical one.
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
                            # y is the input from apply call; similar effect
                            # to having a nested for-loop.
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
            # Mutual Information (jackknife bias-corrected)
            # https://cran.r-project.org/web/packages/mpmi/vignettes/Vignette.pdf
            #
            'cmi' = function(dat1, dat2=NULL, ...) {
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
            # Maximal Information Coefficient (MIC)
            # http://www.exploredata.net/
            #
            'mic' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    # two datasets
                    mic_mat <- minerva::mine(cbind(dat1, dat2), ...)$MIC
                    mic_mat[1:ncol(dat1), (ncol(dat1) + 1):nrow(mic_mat)]
                } else {
                    # single dataset
                    minerva::mine(dat1, ...)$MIC
                }
            },

            #
            # Mutual Information (infotheo discretized estimate)
            # https://cran.r-project.org/web/packages/infotheo/index.html
            #
            'mutinformation' = function(dat1, dat2=NULL, ...) {
                # discretize continuous data
                if (!is.integer(dat1)) {
                    dat1 <- infotheo::discretize(dat1, ...)
                }

                # two datasets
                if (!is.null(dat2)) {
                    # result matrix
                    mi_mat <- matrix(nrow=ncol(dat1), ncol=ncol(dat2))

                    if (!is.integer(dat2)) {
                        dat2 <- infotheo::discretize(dat2, ...)
                    }

                    # iterate over pairs of columns and compute mutual information
                    # for each pair
                    for (i in 1:ncol(dat1)) {
                        for (j in 1:ncol(dat2)) {
                            mi_mat[i, j] <- infotheo::mutinformation(dat1[, i], dat2[, j])
                        }
                    }
                } else {
                    # single dataset

                    # result matrix
                    mi_mat <- matrix(nrow=ncol(dat1), ncol=ncol(dat1))

                    # iterate over pairs of columns and compute mutual information
                    # for each pair
                    for (i in 1:ncol(dat1)) {
                        for (j in 1:ncol(dat1)) {
                            if (is.na(mi_mat[i, j])) {
                                mi_mat[i, j] <- infotheo::mutinformation(dat1[, i], dat1[, j])
                                mi_mat[j, i] <- mi_mat[i, j]
                            }
                        }
                    }
                }
                mi_mat
            },

            # Kendall correlation
            #
            'kendall' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    cor(dat1, dat2, method = 'kendall', ...)
                } else {
                    cor(dat1, method = 'kendall', ...)
                }
            },

            # Pearson correlation
            #
            'pearson' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    # data is passed in with shared rows, so columns can
                    # be correlation directly
                    cor(dat1, dat2, ...)
                } else {
                    cor(dat1, ...)
                }
            },

            #
            # Spearman correlation
            #
            'spearman' = function(dat1, dat2=NULL, ...) {
                if (!is.null(dat2)) {
                    cor(dat1, dat2, method = 'spearman', ...)
                } else {
                    cor(dat1, method = 'spearman', ...)
                }
            }
        ),

        check_input = function() {
            # TODO: Check to make sure all datasets overalp in keys with dat
        },

        # Computes cross-dataset correlation matrix
        #
        # Measures similarity between rows or columns in two datasets which
        # share the same set of column or row ids.
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param meas Correlation measure to use. Supported options:
        #   - kendall  (Kendall correlation)
        #   - pearson  (Pearson correlation)
        #   - spearman (Spearman correlation)
        #   - lm       (Linear model)
        #   - cmi      (Mututal information)
        #
        # @return Matrix of pairwise dataset1 - dataset2 correlations
        compute_cross_cor = function(key1=1, key2=2, meas='pearson', new_key=NULL, ...) {
            # shortcuts to edats
            e1 <- self$edat[[key1]]
            e2 <- self$edat[[key2]]

            # make sure datasets share some common ids
            if (length(unique(c(e1$xid, e2$xid, e1$yid, e2$yid))) == 4) {
                stop("Specified datasets must share a common axis.")
            }

            # determine orientation to use for comparison; similarity is
            # measured across columns, so data will be rearranged such that
            # shared ids are along rows (then column 1 from dat1 can be
            # compared with column 1 from dat2, etc.)
            dat1_shared_axis <- ifelse(e1$xid == e2$xid || e1$xid == e2$yid, 'rows', 'columns')
            dat2_shared_axis <- ifelse(e2$xid == e1$xid || e2$xid == e1$yid, 'rows', 'columns')

            # get row-oriented datasets and associated style information
            if (dat1_shared_axis == 'rows') {
                # if dat1 is already oriented with shared row ids, then it
                # is in the expected orientation
                dat1 <- e1$dat

                # each row in the resulting correlation matrix will correspond
                # to a column in dat1
                xid <- e1$yid

                # similarly, the styles and labels for the rows in the cor
                # matrix correspond to that for the columns in dat1
                xcolor <- e1$col_color
                xshape <- e1$col_shape
                xlabel <- e1$col_label
                xedat  <- e1$col_edat 
            } else {
                # if dat1 shared its column names with dat2, transpose to
                # match expected orientation
                dat1 <- e1$tdat

                # in this case, each row in the resulting dataset will correspond
                # to a single row in dat1
                xid <- e1$xid

                xcolor <- e1$row_color
                xshape <- e1$row_shape
                xlabel <- e1$row_label
                xedat  <- e1$row_edat 
            }

            # similar logic to above, but for dat2 which will appear as column
            # entries in the resulting correlation matrix
            if (dat2_shared_axis == 'rows') {
                dat2 <- e2$dat
                yid <- e1$yid

                ycolor <- e2$col_color
                yshape <- e2$col_shape
                ylabel <- e2$col_label
                yedat  <- e2$col_edat 
            } else {
                dat2 <- e2$tdat
                yid  <- e2$xid

                ycolor <- e2$row_color
                yshape <- e2$row_shape
                ylabel <- e2$row_label
                yedat  <- e2$row_edat 
            }

            # normalize order of shared axis entries (ordered as rows now) and
            # only operate on shared entries
            row_ind <- intersect(rownames(dat1), rownames(dat2))

            dat1 <- dat1[order(match(rownames(dat1), row_ind)),, drop = FALSE]
            dat2 <- dat2[order(match(rownames(dat2), row_ind)),, drop = FALSE]

            # measure similarity between rows in datasets 1 and rows in dataset 2
            cor_mat <- private$similarity(dat1, dat2, meas = meas, ...)

            # clone object and add new dataset
            obj <- private$clone_()

            # map numeric keys to string names
            key1 <- ifelse(is.numeric(key1), names(self$edat)[key1], key1)
            key2 <- ifelse(is.numeric(key2), names(self$edat)[key2], key2)

            # determine key to use for storing result
            new_key <- ifelse(is.null(new_key), sprintf('%s_%s_%s', key1, key2, meas), new_key)

            # add new matrix to front of edat list and return
            obj$add(new_key,
                    EDADat$new(cor_mat, xid = xid, yid = yid, 
                               row_color = xcolor, row_shape = xshape, 
                               row_label = xlabel, row_edat  = xedat,
                               col_color = ycolor, col_shape = yshape, 
                               col_label = ylabel, col_edat  = yedat))

            obj
        },

        # Plots multidataset correlation heatmap
        #
        # @param key1 Numeric or character index of first dataset to use
        # @param key2 Numeric or character index of second dataset to use
        # @param meas Correlation measure to use (passed to `cor` function)
        #
        plot_cross_cor_heatmap = function(key1=1, key2=2, meas='pearson', show_tick_labels=c(TRUE, TRUE), interactive=TRUE) {
            # compute cross correlations
            cor_mat <- private$compute_cross_cor(key1, key2, meas)

           # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = show_tick_labels,
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

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

            # rename column and row annotation matrix parameters, if present
            if ('col_side_colors' %in% names(params)) {
                params$annCol <- params$col_side_colors
            }

            if ('row_side_colors' %in% names(params)) {
                params$annRow <- params$row_side_colors
            }

            # remove irrelevant heatmaply function arguments
            heatmaply_args <- c('showticklabels', 'subplot_widths',
                                'subplot_heights', 'col_side_colors',
                                'row_side_colors')
            params <- params[!names(params) %in% heatmaply_args]

            # default to viridis color scheme
            if (!'color' %in% names(params)) {
                params$color <- viridis(100)
            }

            # 2018/04/17 - adding border just around heatmap not working at
            # this time and results in warnings being generated
            #params$border <- list(matrix = TRUE)

            do.call(NMF::aheatmap, params)
        },

        # Removes any datasets which are not linked to the target dataset by
        # either shared row or column ids
        remove_unlinked = function(key) {
            # iterate over datasets
            edat_keys <- names(self$edat)

            if (is.numeric(key)) {
                edat_keys <- edat_keys[-key]
            } else {
                edat_keys <- edat_keys[names(edat_keys) != key]
            }

            for (eid in edat_keys) {
                # create a vector of axis ids for all other datasets
                axis_ids <- c()

                for (edat in self$edat[names(self$edat) != eid]) {
                    axis_ids <- c(axis_ids, edat$xid, edat$yid)
                }

                # if a dataset does not overlap in ids with any of the other
                # datasets (e.g. after a projection which replaces the original
                # dataset), remove it.
                if (length(intersect(axis_ids, c(self$edat[[eid]]$xid, self$edat[[eid]]$yid))) == 0) {
                    self$edat[[eid]] <- NULL
                }
            }
        },

        #
        # Measures similarity between columns within or across datasets
        #
        similarity = function(dat1, dat2=NULL, meas='pearson', ...) {
            # check to make sure specified similarity measure is valid
            if (!meas %in% names(private$similarity_measures)) {
                stop('Unsupported similarity measure specified.')
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
            cor_mat <- private$similarity_measures[[meas]](dat1, dat2, ...)

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
