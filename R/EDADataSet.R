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

        # EDADataset constructor
        initialize = function(dat, 
                              row_mdata=NULL, col_mdata=NULL, 
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
                              row_maxn=Inf, row_max_ratio=1.0, row_ind=NULL,
                              col_maxn=Inf, col_max_ratio=1.0, col_ind=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) { 

            # check to make sure row and columns identifiers are stored as
            # row and column names
            self$dat <- private$normalize_data_ids(dat, row_ids, col_ids)

            # store row and column metadata, if present
            self$col_mdata <- private$normalize_metadata_order(col_mdata, colnames(self$dat))
            self$row_mdata <- private$normalize_metadata_order(row_mdata, rownames(self$dat))
            
            # default variables to use for plot color, shape, and labels when
            # visualizing either columns or rows in the dataset
            private$row_color  <- row_color
            private$row_shape  <- row_shape
            private$row_labels <- row_labels

            private$col_color  <- col_color
            private$col_shape  <- col_shape
            private$col_labels <- col_labels
            private$color_pal  <- color_pal

            # determine subsampling indices, if requested
            private$row_ind <- private$get_subsample_indices(row_maxn, row_max_ratio, 
                                                             row_ind, nrow(dat))
            private$col_ind <- private$get_subsample_indices(col_maxn, col_max_ratio,
                                                             col_ind, ncol(dat))

            private$ggplot_theme <- ggplot_theme
            private$title  <- title
        },

        #' Clears any cached resuts and performs garbage collection to free 
        #' up memory.
        clear_cache = function() {
            private$cache <- list()
            invisible(gc())
        },

        #' Applies a filter to rows of the dataset
        #'
        #' @param mask Logical vector of length equal to the number of rows in
        #'      the dataset.
        #'
        #' @return A filtered version of the original EDADataSet object.
        filter_rows = function(mask) {
            obj <- private$clone_()

            # update data and metadata matrices
            obj$dat <- obj$dat[mask,, drop=FALSE] 
            obj$row_mdata <- obj$row_mdata[mask,, drop=FALSE]

            # update subsampling indices
            rows_kept <- seq(nrow(self$dat))[mask]
            new_ind <- match(obj$get_row_indices(), rows_kept)
            new_ind <- new_ind[!is.na(new_ind)]
            obj$set_row_indices(new_ind)

            obj
        },

        #' Applies a filter to columns of the dataset
        #'
        #' @param mask Logical vector of length equal to the number of columns 
        #'     in the dataset.
        #'
        #' @return A filtered version of the original EDADataSet object.
        filter_cols = function(mask) {
            obj <- private$clone_()

            # update data and metadata matrices
            obj$dat <- obj$dat[,mask, drop=FALSE] 
            obj$col_mdata <- obj$col_mdata[mask,, drop=FALSE]

            # update subsampling indices
            cols_kept <- seq(ncol(self$dat))[mask]
            new_ind <- match(obj$get_col_indices(), cols_kept)
            new_ind <- new_ind[!is.na(new_ind)]
            obj$set_col_indices(new_ind)

            obj
        },

        #' Imputes missing values in the dataset
        #'
        #' Imputes missing values in the dataset using a specified method
        #' and stores the result in-place. Currently only support k-Nearest
        #' Neighbors (kNN) method.
        #'
        #' Note: When using the `knn` method, it may be neccessary to set R 
        #' the environmental variable `R_MAX_NUM_DLLS` to some value larger
        #' than its default of 100 (150 should be sufficient), prior to
        #' launching R. This can be set in the ~/.Renviron.
        #'
        #' @param method Character array specifying imputation method to use 
        #'     (knn)
        impute = function(method='knn') {
            # Keep track of original dataset class
            cls <- class(self$dat)

            message("Imputing %d missing values...", sum(is.na(self$dat)))
           
            # kNN
            if (method == 'knn') {
                imputed <- VIM::kNN(self$dat)[,1:ncol(self$dat)]
                rownames(imputed) <- rownames(self$dat)
                
                if (cls == 'matrix') {
                    imputed <- as.matrix(imputed)
                }
            }

            message("Done.")

            self$dat <- imputed
        },

        #' Summarizes overall characteristics of a dataset
        #' 
        summary = function(markdown=FALSE, subsample=TRUE, num_digits=2) {
            # collection summary info
            x <- self$dat

            # list to store summary info
            info <- list()

            # include trailing zeros when rounding
            round_ <- function(y, digits) {
                formatC(round(y, digits), digits, format="f")
            }

            # overall dataset
            info[['dat_range']]     <- round_(range(x, na.rm=TRUE), num_digits)
            info[['dat_num_nas']]   <- sum(is.na(x))
            info[['dat_num_zeros']] <- sum(x == 0)
            info[['dat_quartiles']] <- round_(quantile(x, na.rm=TRUE), num_digits)

            # column types
            info[['col_types']] <- table(apply(x, 2, class))

            # rows
            info[['row_means']]    <- round_(range(apply(x, 1, mean, na.rm=TRUE)), num_digits)
            info[['row_medians']]  <- round_(range(apply(x, 1, median, na.rm=TRUE)), num_digits)
            info[['row_std_devs']] <- round_(range(apply(x, 1, sd, na.rm=TRUE)), num_digits)
            info[['row_outliers']] <- self$detect_row_outliers()

            # columns
            info[['col_means']]    <- round_(range(apply(x, 2, mean, na.rm=TRUE)), num_digits)
            info[['col_medians']]  <- round_(range(apply(x, 2, median, na.rm=TRUE)), num_digits)
            info[['col_std_devs']] <- round_(range(apply(x, 2, sd, na.rm=TRUE)), num_digits)
            info[['col_outliers']] <- self$detect_col_outliers()

            # clip outlier lists if more than a few
            if (length(info[['row_outliers']]) > 5) {
                info[['row_outliers']] <- c(head(info[['row_outliers']], 5), '...')
            }
            if (length(info[['col_outliers']]) > 5) {
                info[['col_outliers']] <- c(head(info[['col_outliers']], 5), '...')
            }

            # row & column correlations
            if (subsample) {
                info[['col_cor_mat']] <- round_(cor(x[private$row_ind, private$col_ind]), num_digits)
                info[['row_cor_mat']] <- round_(cor(t(x[private$row_ind, private$col_ind])), num_digits)
            } else {
                info[['col_cor_mat']] <- round_(cor(x), num_digits)
                info[['row_cor_mat']] <- round_(cor(t(x)), num_digits)
            }
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

        #' Kernel density plot
        #'
        #' Plots densities for each row or column in the dataset. This is most
        #' useful when you are interested in similarties or differences in
        #' distributions across columns, for a relatively small number of
        #' columns.
        #'
        #' @param color Variable to color density curves by. If not is
        #' specified, uses variable specified at object construction time,
        #' or else, uses a separate color for each column.
        #'
        #' @return ggplot plot instance.
        #'
        plot_densities = function(color=NULL, title="", ...) {
            #if (target == 'rows') {
            #    private$transpose()
            #}

            # determine subsampling indices, if requested
            indices <- private$get_indices(...)

            dat <- setNames(melt(self$dat), c('row', 'column', 'val'))
            styles <- private$get_geom_density_styles(color)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Column densities: %s", private$title)
            }

            # construct density plot
            plt <- ggplot(dat, aes(x=val, group=column, color=column)) +
                geom_density(styles$aes) +
                ggtitle(title) +
                private$ggplot_theme()

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}

            # only show legend if there are relatively few groups
            if (length(unique(dat$column)) > 10) {
                plt <- plt + guides(color=FALSE)
            }

			plt
        },
        
        #' Prints class greeting to the screen
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat("= EDADataSet\n")
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        },

        #' Gets the subsampling row indices for a dataset
        #'
        #' @return Numeric vector of row indices to use for intensive 
        #'     operations when subsampling is enabled.
        get_row_indices = function() {
            private$row_ind
        },

        #' Gets the subsampling column indices for a dataset
        #'
        #' @return Numeric vector of column indices to use for intensive 
        #'     operations when subsampling is enabled.
        get_col_indices = function() {
            private$col_ind
        },

        #' Sets the subsampling row indices for a dataset
        #'
        #' @param ind Numeric vector of row indices to use for intensive 
        #'     operations when subsampling is enabled.
        set_row_indices = function(ind) {
            private$row_ind <- ind 
        },

        #' Sets the subsampling column indices for a dataset
        #'
        #' @param ind Numeric vector of column indices to use for intensive 
        #'     operations when subsampling is enabled.
        set_col_indices = function(ind) {
            private$col_ind <- ind 
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # private params
        row_ind      = NULL,
        col_ind      = NULL,
        row_color    = NULL,
        row_shape    = NULL,
        row_labels   = NULL,
        col_color    = NULL,
        col_shape    = NULL,
        col_labels   = NULL,
        color_pal    = NULL,
        ggplot_theme = NULL,
        title        = NULL,
        cache        = list(),

        # private methods
        clone_ = function() {
            obj <- self$clone()
            obj$clear_cache()
            obj
        },

        #' Computes correlations between axes of a data projection (PCA or
        #' t-SNE) and column metadata.
        #' 
        #' @param mat Numeric projected data matrix
        #' @param include Vector of strings indicating metadata columns which
        #' should be included in the analysis.
        #'
        #' @return Dataframe of feature x project axes dependencies
        compute_feature_correlations = function(mat, include) {
            # drop any covariates with only a single level
            single_level <- apply(self$row_mdata, 2, function(x) {length(table(x))}) == 1
            features <- self$row_mdata[,!single_level, drop=FALSE]

            # drop any undesired features
            if (!is.null(include)) {
                features <- features %>%
                    select(one_of(include))
            }

            # construct linear model to measure dependence of each projected
            # axis on column metadata
            result <- matrix(0, nrow=ncol(mat), ncol=ncol(features))

            for (i in 1:ncol(features)) {
                feature_cor <- function(y) {
                    round(summary(lm(y~features[,i]))$r.squared*100, 2)
                }
                result[,i] <- apply(mat, 2, feature_cor)
            }
            colnames(result) <- colnames(features)
            result
        },

        #' Returns a list of arguments specific for a given function call.
        #'
        #' This function takes a list of function arguments and generates
        #' a new list without any shared arguments (row_maxn, etc.).
        #'
        #' @param ... Arguments passed to a given plotting, etc. function call.
        #'
        #' @return A list containing only function-specific arguments 
        get_custom_function_args = function(...) {
            args <- list(...)
      
            shared_args <- c('row_ind', 'row_maxn', 'row_max_ratio',
                             'col_ind', 'col_maxn', 'col_max_ratio')

            # return list of non-shared arguments
            args[!names(args) %in% shared_args]
        },

        get_subsample_indices = function(maxn, maxr, ind, n) {
            # if indices are explictly provided, use them
            if (!is.null(ind)) {
                return(ind)
            }

            # otherwise, if maxn or maxr, specified, use to randomly select
            # row/column indices to use for memory-/cpu-intensive operations
            maxn <- min(maxn, round(maxr * n))
            
            if (maxn < n) {
                sort(sample(n, maxn))
            } else {
                1:n
            }
        },

        #' Determines row and column indices to use for a function call
        #'
        #' Parses a list of function arguments and check whether user has
        #' specified any subsampling-related arguments (row_maxn, row_max_ratio,
        #' etc.) and, if so, select indices to use. Otherwise object-level
        #' indices are returned.
        #'
        #' @param ... Arguments passed to a given plotting, etc. function call.
        #'
        #' @return A list with 'row' and 'column' entries corresponding to the
        #'     specific numeric row and column indices to be used.
        get_indices = function(...) {
            args <- list(...)
            result <- list(row=NULL, col=NULL)

            # row indices
            if (!is.null(args$row_ind)) {
                result[['row']] <- args$row_ind
            } else if (!is.null(args$row_maxn)) {
                result[['row']] <- sample(nrow(self$dat), args$row_maxn) 
            } else if (!is.null(args$row_max_ratio) ){
                result[['row']] <- sample(nrow(self$dat), round(nrow(self$dat) * args$row_max_ratio))
            } else {
                result[['row']] <- private$row_ind
            }

            # col indices
            if (!is.null(args$col_ind)) {
                result[['col']] <- args$col_ind
            } else if (!is.null(args$col_maxn)) {
                result[['col']] <- sample(ncol(self$dat), args$col_maxn) 
            } else if (!is.null(args$col_max_ratio) ){
                result[['col']] <- sample(ncol(self$dat), round(ncol(self$dat) * args$col_max_ratio))
            } else {
                result[['col']] <- private$col_ind
            }

            result
        },

        #' Generates ggplot aesthetics for density plots 
        #'
        #' @param color Color variable as passed into plot function call
        #'
        #' @return List of style information
        get_geom_density_styles = function(color) {
            # list to store style properties
            res <- list(aes=aes(), labels=list())
            private$add_color_styles(res, color)
        },

        #' Generates ggplot aesthetics for histogram plots
        #'
        #' @param color Color variable as passed into plot function call
        #'
        #' @return List of style information
        get_geom_histogram_styles = function(color) {
            # list to store style properties
            res <- list(aes=aes(), labels=list())
            private$add_color_styles(res, color)
        },

        #' Generates ggplot aesthetics for a geom_point plot
        #'
        #' @param color Color variable as passed into plot function call
        #' @param shape Shape variable as passed into plot function call
        #'
        #' @return List of geom_point style information
        get_geom_point_styles = function(color, shape) {
            # list to store style properties
            res <- list(
                aes=aes(),
                labels=list()
            )
            res <- private$add_color_styles(res, color)
            res <- private$add_shape_styles(res, shape)

            res
        },

        #' Determines color-related style information to use for a plot
        #'
        #' @param styles List of color-related style info
        #' @param color Color variable passed into plot function call
        add_color_styles = function(styles, color) {
            color_info <- private$get_color_styles(color)

            styles[['color']] <- color_info[['color']]

			# update styles with color info
			if (length(color_info[['aes']]) > 0) {
				styles[['aes']] <- modifyList(color_info[['aes']], styles[['aes']])
				styles[['labels']] <- modifyList(color_info[['labels']], styles[['labels']])
			}
            styles
        },

        #' Determines shape-related style information to use for a plot
        #'
        #' @param styles List of shape-related style info
        #' @param shape shape variable passed into plot function call
        add_shape_styles = function(styles, shape) {
            shape_info <- private$get_shape_styles(shape)

            styles[['shape']] <- shape_info[['shape']]

			# update styles with shape info
			if (length(shape_info[['aes']]) > 0) {
				styles[['aes']] <- modifyList(shape_info[['aes']], styles[['aes']])
				styles[['labels']] <- modifyList(styles[['labels']], shape_info[['labels']])
			}
            styles
        },

        #' Returns a list of color-related style information
        #'
        #' @param color Color variable as passed into plot function call
        #' 
        #' @return list List of color-related properties
        get_color_styles = function(color) {
            res <- list(
				color  = NULL,
				aes    = aes(),
				labels = list()
			)

            # if specified as a function argument, override default color
            if (!is.null(color) && (color != FALSE)) {
                # color variable can either correspond to a column in the
                # dataset itself, or in the 
                res[['color']]  <- self$row_mdata[,color]
                res[['aes']]    <- aes(color=color)
                res[['labels']] <- labs(color=color)
            } else if (is.null(color) && !is.null(private$row_color)) {
                # otherwise, use object-level default value
                res[['color']] <- self$row_mdata[,private$row_color]
                res[['aes']] <- aes(color=color)
                res[['labels']] <- labs(color=private$row_color)
            }
            res
        },

        #' Returns a list of shape-related style information
        #'
        #' @param shape Shape variable as passed into plot function call
        #' 
        #' @return list List of shape-related properties
        get_shape_styles = function(shape) {
            res <- list(
				shape  = NULL,
				aes    = aes(),
				labels = list()
			)

            # if specified as a function argument, override default shape
            if (!is.null(shape) && (shape != FALSE)) {
                res[['shape']]  <- self$row_mdata[,shape]
                res[['aes']]    <- aes(shape=shape)
                res[['labels']] <- labs(shape=shape)
            } else if (is.null(shape) && !is.null(private$shape)) {
                # otherwise, use object-level default value
                res[['shape']] <- self$row_mdata[,private$shape]
                res[['aes']] <- aes(shape=shape)
                res[['labels']] <- labs(shape=private$shape)
            }
            
            res
        },

		#' Returns a vector of color codes associated with the specified
		#' variable.
		#'
		#' @param color Variable to use for assigning colors.
		#' @param color_pal Color palette to use
		#'
		#' @return Vector of colors with length equal to the number of columns
		#' 	       in the data.
        get_var_colors = function(color, color_pal) {
            # if no variable is specified, use default black for plots
            if (is.null(color)) {
				if (is.null(private$col_colo)) {
					return('black') 
				}
				color <- private$col_colo
            } else if (color == FALSE) {
				return('black')
			}

            # otherwise, assign colors based on the variable specified
            column_var <- as.numeric(factor(self$row_mdata[,color]))

            pal <- RColorBrewer::brewer.pal(9, color_pal)
            colors <- colorRampPalette(pal)(min(1E4, length(unique(column_var))))

            # return vector of column color assignments
            colors[as.integer(column_var)]
        },

		#' Returns a vector of text labels to use for a plot
		#'
		#' @param shape Variable to use for assigning shapes.
		#'
		#' @return Vector of labels with length equal to the number of columns
		#'         in the data.
        get_var_labels = function(label) {
            if (is.null(label)) {
                return(colnames(self$dat))
            }
            self$row_mdata[,label]
        },

        #' Creates a static or interactive heatmap plot
        #'
        #' @param params A list of plotting parameters
        #' @param interactive Logical indicating whether an interactive heatmap
        #'     should be generated.
        plot_heatmap = function(params, interactive) {
            # interactive heatmap
            if (interactive) {
                return(do.call(heatmaply::heatmaply, params))
            }

            # static heatmap

            # convert colside and rowside colors to explicit color values,
            # if present
            if ('col_side_colors' %in% names(params)) {
                params$ColSideColors <- params$col_side_colors

                colors <- c('blue', 'yellow')

                for (col_name in colnames(params$ColSideColors)) {
                    col <- params$ColSideColors[,col_name]
                    params$ColSideColors[,col_name] <- colors[as.numeric(factor(col))]
                }
                params$ColSideColors <- as.matrix(params$ColSideColors)
            }

            if ('row_side_colors' %in% names(params)) {
                params$RowSideColors <- params$row_side_colors

                pal <- RColorBrewer::brewer.pal(8, 'Set1')

                for (col_name in colnames(params$RowSideColors)) {
                    col <- params$RowSideColors[,col_name]
                    colors <- colorRampPalette(pal)(min(1E4, length(unique(col))))
                    params$RowSideColors[,col_name] <- colors[as.numeric(factor(col))]
                }
                params$RowSideColors <- as.matrix(params$RowSideColors)
            }

            # remove irrelevant function arguments
            heatmaply_args <- c('showticklabels', 'subplot_widths', 'subplot_heights',
                                'col_side_colors', 'row_side_colors')
            params <- params[!names(params) %in% heatmaply_args]


            do.call(heatmap.plus::heatmap.plus, params)
        },

        #' Normalizes handling of data row and column identifiers
        #'
        #' Checks dataset row and column identifiers and converts them to row 
        #' and column names, respectively, if they are not already stored there.
        #'
        #' @param dat Dataset
        #' @param row_ids Column id or number where row identifiers are stored.
        #' @param col_ids Row id or number where column identifiers are stored.
        #'
        #' @param Dataset with identifiers as rows and columns
        normalize_data_ids = function(dat, row_ids, col_ids) {
            # row ids
            if (row_ids != 'rownames') {
                # column number containing row ids specified
                if (is.numeric(row_ids)) {
                    rownames(dat) <- dat[,row_ids]
                    dat <- dat[,-row_ids]
                } else if (row_ids %in% colnames(dat)) {
                    # column name containing row ids specified
                    ind <- which(colnames(dat) == row_ids)
                    rownames(dat) <- dat[,ind]
                    dat <- dat[,-ind]
                }
            }

            # column ids
            if (col_ids != 'colnames') {
                # row number containing column ids
                if (is.numeric(col_ids)) {
                    colnames(dat) <- dat[col_ids,]
                    dat <- dat[-col_ids,]
                } else if (col_ids %in% colnames(dat)) {
                    # row name containing columns ids
                    ind <- which(rownames(dat) == col_ids)
                    colnames(dat) <- dat[ind,]
                    dat <- dat[-ind,]
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

        # normalize dat / metadata order
        normalize_metadata_order = function(metadata, ids) {
            # If no metadata was provided, do nothing
            if (is.null(metadata)) {
                return(NULL)
            }

            # if already in the right order, return as-is
            if (all(rownames(metadata) == ids)) {
                return(metadata)
            }

            # otherwise match metadata row order to row/col names specified
            ind <- order(match(rownames(metadata), ids))

            # return result
            metadata[ind,, drop=FALSE]
        },

        #' Prints dataset summary to screen
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
            col_types[,1] <- paste0(" ", col_types[,1])
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
                        min(x$col_cor_mat, na.rm=TRUE),
                        max(x$col_cor_mat, na.rm=TRUE)))
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
                        min(x$row_cor_mat, na.rm=TRUE),
                        max(x$row_cor_mat, na.rm=TRUE)))

            cat("\n")
            cat(" Outliers:\n\n ")
            cat(paste0(sprintf("%2d", 1:length(x$row_outliers)), ". ", x$row_outliers, "\n"))
            cat("\n")
            cat("=========================================\n")
        },

        #' Transpose the dataset and metadata in-place
        transpose = function() {
            self$dat <- t(self$dat) 

            row_mdata  <- self$row_mdata
            col_mdata  <- self$col_mdata
            row_color  <- private$row_color 
            row_shape  <- private$row_shape 
            row_labels <- private$row_labels 
            col_color  <- private$col_color 
            col_shape  <- private$col_shape 
            col_labels <- private$col_labels 
            row_ind    <- private$row_ind 
            col_ind    <- private$col_ind

            self$row_mdata     <- col_mdata
            self$col_mdata     <- row_mdata
            private$row_color  <- col_color 
            private$row_shape  <- col_shape 
            private$row_labels <- col_labels 
            private$col_color  <- row_color 
            private$col_shape  <- row_shape 
            private$col_labels <- row_labels 
            private$row_ind    <- col_ind 
            private$col_ind    <- row_ind
        }
    ),
    # ------------------------------------------------------------------------
    # active bindings
    # ------------------------------------------------------------------------
    active = list(
        t = function() {
            cls <- get(class(self)[1])
            cls$new(t(self$dat), 
                           row_mdata=self$col_mdata,
                           col_mdata=self$row_mdata,
                           row_color=self$col_color, row_shape=self$col_shape, 
                           row_labels=self$col_labels, 
                           col_color=self$row_color, col_shape=self$row_shape, 
                           col_labels=self$row_labels, 
                           row_ind=self$col_ind, col_ind=self$row_ind,
                           color_pal=self$color_pal, ggplot_theme=private$ggplot_theme)
        }
    )
)

