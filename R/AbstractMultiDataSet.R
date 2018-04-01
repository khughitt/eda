#' An S6 class representing a generic collection of datasets linked by
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
        # AbstractMultiDataSet constructor
        initialize = function(dat, row_data=list(), col_data=list(),
                              row_color=NULL, row_color_ds='dat',
                              row_shape=NULL, row_shape_ds='dat',
                              row_label=NULL, row_label_ds='dat',
                              col_color=NULL, col_color_ds='dat',
                              col_shape=NULL, col_shape_ds='dat',
                              col_label=NULL, col_label_ds='dat',
                              color_pal='Set1', title="", ggplot_theme=theme_bw) {

            #dat <- private$normalize_data_ids(dat, row_ids, col_ids)

            # determine orientation of any addition row and column datasets,
            # relative to dat

            # TODO: revisit handling of row / column id's for multiple datasets
            # with MultiDataset classes... (skipping for now..)
            row_mdata <- private$normalize_row_order(row_data, rownames(dat))
            col_mdata <- private$normalize_col_order(col_data, colnames(dat))

            # store datasets
            private$row_data <- c(list('dat' = dat), row_data)
            private$col_data <- c(list('dat' = dat), col_data)

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
        cor = function(method='pearson', ...) {
            private$similarity(self$dat, method=method, ...)
        },

        # Applies a filter to rows of the dataset
        #
        # @param mask Logical vector of length equal to the number of rows in
        #      the dataset.
        #
        # @return A filtered version of the original EDADataSet object.
        filter_rows = function(mask) {
            obj <- private$clone_()

            # update data and metadata matrices
            obj$dat <- obj$dat[mask,, drop = FALSE]

            # TODO: Create setter method to allow metadata to be updated?
            # or, implement downstream...
            #if ('metadata' %in% names(private$row_data)) {
                #obj$metadata <- obj$metadata[mask,, drop = FALSE]
            #}

            obj
        },

        # Applies a filter to columns of the dataset
        #
        # @param mask Logical vector of length equal to the number of columns
        #     in the dataset.
        #
        # @return A filtered version of the original EDADataSet object.
        filter_cols = function(mask) {
            obj <- private$clone_()

            # update data and metadata matrices
            obj$dat <- obj$dat[, mask, drop = FALSE]

            # TODO (see above)
            #if (!is.null(obj$metadata)) {
            #    obj$metadata <- obj$metadata[mask,, drop = FALSE]
            #}

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
        impute = function(method='knn') {
            # Keep track of original dataset class
            cls <- class(self$dat)

            message(sprintf("Imputing %d missing values...", sum(is.na(self$dat))))

            # kNN
            if (method == 'knn') {
                imputed <- VIM::kNN(self$dat)[, 1:ncol(self$dat)]
                rownames(imputed) <- rownames(self$dat)

                if (cls == 'matrix') {
                    imputed <- as.matrix(imputed)
                }
            }

            message("Done.")

            self$dat <- imputed
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
            row_ind <- 1:nrow(self$dat)
            col_ind <- 1:ncol(self$dat)

            # subsample rows
            if (!is.null(row_n)) {
                row_ind <- sample(row_ind, row_n)
            } else if (!is.null(row_ratio)) {
                row_ind <- sample(row_ind, round(row_ratio * nrow(self$dat)))
            }

            # subsample columns
            if (!is.null(col_n)) {
                col_ind <- sample(col_ind, col_n)
            } else if (!is.null(col_ratio)) {
                col_ind <- sample(col_ind, round(col_ratio * ncol(self$dat)))
            }

            # clone and subsample dataset
            obj <- private$clone_()

            # update data and metadata matrices
            obj$dat <- obj$dat[row_ind, col_ind, drop = FALSE]

            # TODO: Move to child method... (no row_mdata here...)
            #if (!is.null(obj$metadata)) {
            #    obj$metadata <- obj$metadata[row_ind,, drop = FALSE]
            #}

            #if (!is.null(obj$metadata)) {
            #    obj$metadata <- obj$metadata[col_ind,, drop = FALSE]
            #}

            obj
        },

        # Summarizes overall characteristics of a dataset
        #
        # @param markdown Logical indicating whether summary output should be
        #   formatted using markdown. (NOTE: Not yet implemented...)
        # @param num_digits Number of digits to display when printing summary
        #   output (default: 2).
        summary = function(markdown=FALSE, num_digits=2) {
            # collection summary info
            x <- self$dat

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
        plot_densities = function(color=NULL, title="", ...) {
            dat <- setNames(melt(self$dat), c('row', 'column', 'val'))
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
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= AbstractMultiDataSet (n=%d)\n", length(private$datasets)))
            cat("=\n")
            cat(sprintf("= dat: %s (%d x %d)\n", class(self$dat)[1], nrow(self$dat), ncol(self$dat)))
            if (length(private$row_data) > 1) {
                cat("=\n")
                cat("= Row data\n")
                cat("=\n")
                for (i in 2:length(private$row_data)) {
                    ds <- private$row_data[[i]]
                    cat(sprintf("= %02d. %s (%d x %d)\n", i, class(ds)[1], nrow(ds$dat), ncol(ds$dat)))
                }
            }
            if (length(private$col_data) > 1) {
                cat("=\n")
                cat("= Column data\n")
                cat("=\n")
                for (i in 2:length(private$col_data)) {
                    ds <- private$col_data[[i]]
                    cat(sprintf("= %02d. %s (%d x %d)\n", i, class(ds)[1], nrow(ds$dat), ncol(ds$dat)))
                }
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
            # transpose main dataset, preserving underlying class
            dat_cls <- get(sprintf("as.%s", class(self$dat)))
            self$dat <- dat_cls(t(self$dat))

            # transpose any remaining row- and column-oriented datasets
            rdat <- private$row_data
            cdat <- private$col_data

            private$row_data <- cdat
            private$col_data <- rdat

            # TODO: Find a better way to transpose dataframes, preserving col types and
            # factor levels
            #> sapply(x, class)
            #obs_prop1 obs_prop2                                                                                       
            #"factor"  "factor"                                                                                       
            #> sapply(y, class)                                                                                        
            #obs_prop1 obs_prop2                                                                                       
            #"integer"  "factor" 

            if (length(private$row_data) > 1) {
                for (rdat in names(private$row_data)[2:length(private$row_data)]) {
                    dat_cls <- get(sprintf("as.%s", class(private$row_data[[rdat]])))
                    private$row_data[[rdat]] <- dat_cls(t(private$row_data[[rdat]]))
                }
            }

            if (length(private$col_data) > 1) {
                for (cdat in names(private$col_data)[2:length(private$col_data)]) {
                    dat_cls <- get(sprintf("as.%s", class(private$col_data[[cdat]])))
                    private$col_data[[cdat]] <- dat_cls(t(private$col_data[[cdat]]))
                }
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
            private$col_label    <- row_label    
            private$col_shape_ds <- row_color_ds 
            private$col_shape_ds <- row_shape_ds 
            private$col_label_ds <- row_label_ds 

            # transpose row datasets, if they exist
            #if (length(private$row_data) > 1) {
            #    for (key in names(private$row_data)[2:length(private$row_data)]) {
            #        dat <- private$row_data[[key]]

            #        # Preserve underlying dataset class (dataframe or matrix)
            #        dat_cls <- get(sprintf("as.%s", class(dat)))
            #        col_data[[key]] <- dat_cls(t(dat))

            #        rownames(col_data[[key]]) <- colnames(dat)
            #        colnames(col_data[[key]]) <- rownames(dat)
            #    }
            #}

            ## transpose column datasets, if they exist
            #if (length(private$col_data) > 1) {
            #    for (key in names(private$col_data)[2:length(private$col_data)]) {
            #        dat <- private$col_data[[key]]

            #        # Preserve underlying dataset class (dataframe or matrix)
            #        dat_cls <- get(sprintf("as.%s", class(dat)))
            #        row_data[[key]] <- dat_cls(t(dat))

            #        rownames(row_data[[key]]) <- colnames(dat)
            #        colnames(row_data[[key]]) <- rownames(dat)
            #    }
            #}

            ## create transposed object instance
            #cls$new(tdat, row_data = row_data, col_data = col_data,
            #        row_color = private$col_color, row_color_ds = private$col_color_ds,
            #        row_shape = private$col_shape, row_shape_ds = private$col_shape_ds,
            #        row_label = private$col_label, row_label_ds = private$col_label_ds,
            #        col_color = private$row_color, col_color_ds = private$row_color_ds,
            #        col_shape = private$row_shape, col_shape_ds = private$row_shape_ds,
            #        col_label = private$row_label, col_label_ds = private$row_label_ds,
            #        color_pal = private$color_pal,
            #        title     = private$title,
            #        ggplot_theme = private$ggplot_theme)
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # private params
        row_data = NULL,
        col_data = NULL,

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

        # Generates ggplot aesthetics for bar plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_bar_styles = function(color, dataset_key=NULL) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(res, color, dataset_key)
        },

        # Generates ggplot aesthetics for density plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_density_styles = function(color, dataset_key=NULL) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(res, color, dataset_key)
        },

        # Generates ggplot aesthetics for histogram plots
        #
        # @param color Color variable as passed into plot function call
        #
        # @return List of style information
        get_geom_histogram_styles = function(color, dataset_key=NULL) {
            # list to store style properties
            res <- list(aes = aes(), labels = list())
            private$add_color_styles(res, color, dataset_key)
        },

        # Generates ggplot aesthetics for a geom_point plot
        #
        # @param color Color variable as passed into plot function call
        # @param shape Shape variable as passed into plot function call
        #
        # @return List of geom_point style information
        get_geom_point_styles = function(color, shape,
                                         color_dataset_key=NULL,
                                         shape_dataset_key=NULL) {
            # list to store style properties
            res <- list(
                aes = aes(),
                labels = list()
            )
            res <- private$add_color_styles(res, color, color_dataset_key)
            res <- private$add_shape_styles(res, shape, shape_dataset_key)
            res
        },

        # Determines color-related style information to use for a plot
        #
        # @param styles List of color-related style info
        # @param color Color variable passed into plot function call
        add_color_styles = function(styles, color, dataset_key=NULL) {
            color_info <- private$get_color_styles(color, dataset_key)

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
        add_shape_styles = function(styles, shape, dataset_key=NULL) {
            shape_info <- private$get_shape_styles(shape, dataset_key)

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
        # @param color Color variable as passed into plot function call
        # @param dataset_key Name of row_data dataset containing field to use for color
        #
        # @return list List of color-related properties
        get_color_styles = function(color, dataset_key=NULL) {
            res <- list(
                color  = NULL,
                aes    = aes(),
                labels = list()
            )

            # dataset containing color field
            if (is.null(dataset_key)) {
                dataset_key <- private$row_color_ds
            }

            # if specified as a function argument, override default color
            if (!is.null(color) && (color != FALSE)) {
                # color variable can either correspond to a column in the
                # dataset itself, or in the
                res[['color']]  <- private$row_data[[dataset_key]][, color]
                res[['aes']]    <- aes(color = color)
                res[['labels']] <- labs(color = color)
            } else if (is.null(color) && !is.null(private$row_color)) {
                # otherwise, use object-level default value
                res[['color']]  <- private$row_data[[dataset_key]][, private$row_color]
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
        get_shape_styles = function(shape, dataset_key=NULL) {
            res <- list(
                shape  = NULL,
                aes    = aes(),
                labels = list()
            )

            # dataset containing color field
            if (is.null(dataset_key)) {
                dataset_key <- private$row_shape_ds
            }

            # if specified as a function argument, override default shape
            if (!is.null(shape) && (shape != FALSE)) {
                res[['shape']]  <- private$row_data[[dataset_key]][, shape]
                res[['aes']]    <- aes(shape = shape)
                res[['labels']] <- labs(shape = shape)
            } else if (is.null(shape) && !is.null(private$shape)) {
                # otherwise, use object-level default value
                res[['shape']]  <- private$row_data[[dataset_key]][, private$shape]
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
        get_var_colors = function(color, dataset_key, color_pal) {
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
            if (is.null(dataset_key)) {
                dataset_key <- private$row_color_ds
            }

            # otherwise, assign colors based on the variable specified
            column_var <- as.numeric(factor(self$row_data[[dataset_key]][, color]))

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
        get_var_labels = function(label, dataset_key=NULL) {
            if (is.null(label)) {
                return(colnames(self$dat))
            }

            # dataset containing color field
            if (is.null(dataset_key)) {
                dataset_key <- private$row_label_ds
            }

            private$row_data[[dataset_key]][, label]
        },

        # Normalizes handling of data row and column identifiers
        #
        # Checks dataset row and column identifiers and converts them to row
        # and column names, respectively, if they are not already stored there.
        #
        # @param dat Dataset
        # @param row_ids Column id or number where row identifiers are stored.
        # @param col_ids Row id or number where column identifiers are stored.
        #
        # @param Dataset with identifiers as rows and columns
        normalize_data_ids = function(dat, row_ids, col_ids) {
            # row ids
            if (row_ids != 'rownames') {
                # column number containing row ids specified
                if (is.numeric(row_ids)) {
                    rownames(dat) <- dat[, row_ids]
                    dat <- dat[, -row_ids]
                } else if (row_ids %in% colnames(dat)) {
                    # column name containing row ids specified
                    ind <- which(colnames(dat) == row_ids)
                    rownames(dat) <- dat[, ind]
                    dat <- dat[, -ind]
                }
            }

            # column ids
            if (col_ids != 'colnames') {
                # row number containing column ids
                if (is.numeric(col_ids)) {
                    colnames(dat) <- dat[col_ids, ]
                    dat <- dat[-col_ids, ]
                } else if (col_ids %in% colnames(dat)) {
                    # row name containing columns ids
                    ind <- which(rownames(dat) == col_ids)
                    colnames(dat) <- dat[ind, ]
                    dat <- dat[-ind, ]
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

        # normalize row-oriented dataset row order
        normalize_row_order = function(row_data, ids) {
            # iterate over datasets
            for (rdat in names(row_data)) {
                # if order is not the same as main dataset, reorder
                if (!(all(rownames(row_data[[rdat]]) == ids))) {
                    # match row dataset row order to row names specified
                    ind <- order(match(rownames(row_data[[rdat]]), ids))

                    # return result
                    row_data[[rdat]] <- row_data[[rdat]][ind,, drop = FALSE]
                }
            }
        },

        # normalize column-oriented dataset column order
        normalize_col_order = function(col_data, ids) {
            # iterate over datasets
            for (cdat in names(col_data)) {
                # if order is not the same as main dataset, reorder
                if (!(all(colnames(col_data[[cdat]]) == ids))) {
                    # match column dataset column order to column names specified
                    ind <- order(match(colnames(col_data[[cdat]]), ids))

                    # return result
                    col_data[[cdat]] <- col_data[[cdat]][ind,, drop = FALSE]
                }
            }
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
            # TODO: Check to make sure row_data and col_data shared expected
            # keys with dat
        },

        # Computes cross-dataset correlation matrix
        #
        # Operates on datasets sharing a common column name (i.e. members
        # of col_data).
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
        cross_cor = function(key1='dat', key2=2, method='pearson', ...) {
            # make sure datasets are ordered similarly
            dat1 <- private$col_data[[key1]]
            dat2 <- private$col_data[[key2]]

            # for multidatasets, get underlying data
            if ('EDADataSet' %in% class(dat1)) {
                dat1 <- dat1$dat
                dat2 <- dat2$dat
            }

            # TODO: Include checks?

            # for metadata, we can also exclude factor fields with all
            # unique values (e.g. alternate identifiers)
            #exclude <- apply(dat2, 1, function(x) {
            #    is.factor(x) && length(unique(x)) == length(x)
            #})

            #if (sum(exclude) > 0) {
            #    message(sprintf("Excluding %d unique factor fields", sum(exclude)))
            #    dat2 <- dat2[!exclude, ]

            # only compare matching columns
            col_ind <- intersect(colnames(dat1), colnames(dat2))

            dat1 <- dat1[, col_ind, drop = FALSE]
            dat2 <- dat2[, col_ind, drop = FALSE]

            # measure similarity between rows in datasets 1 and rows in
            # dataset 2; the similarity() method operates on columns so we
            # transposes the datasets first
            cor_mat <- private$similarity(t(dat1), t(dat2), method = method, ...)

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
            cor_mat <- self$cross_cor(key1, key2, method)

            # list of parameters to pass to heatmaply
            params <- list(
                x               = cor_mat,
                showticklabels  = c(FALSE, FALSE),
                subplot_widths  = c(0.65, 0.35),
                subplot_heights = c(0.35, 0.65)
            )

            # if metadata is availble, display along side of heatmap
            if ('metadata' %in% names(private$row_data)) {
                mdata1 <- private$row_data[['metadata']]
                mask1  <- sapply(mdata1, function(x) max(table(x)) > 1 )
                mdata1 <- mdata1[, mask1, drop = FALSE]

                params[['row_side_colors']] <- mdata1
                params[['subplot_widths']] <- c(0.15, 0.3, 0.55)
            }

            if ('metadata' %in% names(self$col_data)) {
                mdata1 <- self$col_data[['metadata']]
                mask2  <- sapply(mdata2, function(x) max(table(x)) > 1 )
                mdata2 <- mdata2[, mask2, drop = FALSE]

                params[['col_side_colors']] <- mdata2
                params[['subplot_heights']] <- c(0.55, 0.3, 0.15)
            }

            # add any additional function arguments
            private$plot_heatmap(params, interactive)
        },

        # Creates a static or interactive heatmap plot
        #
        # @param params A list of plotting parameters
        # @param interactive Logical indicating whether an interactive heatmap
        #     should be generated.
        plot_heatmap = function(params, interactive) {
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
    ),

    # ------------------------------------------------------------------------
    # active bindings
    # ------------------------------------------------------------------------
    active = list(
        # main dataset
        dat = function(value) {
            if (missing(value)) {
                private$col_data[['dat']]
            } else {
                private$row_data[['dat']] <- value
                private$col_data[['dat']] <- value
            }
        }
    )
)
