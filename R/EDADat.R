#' An R6 class representing a single dataset.
#'
#' EDADat is a container class representing a single dataset (matrix or data
#' frame) along with some basic information about its indexing, orientation,
#' and representation.
#'
#' The purpose of this class is to provide a way for each atomic dataset to
#' carry its own relevant metadata, as needed by the other major `eda` classes.
#' This class is different from all other `eda` classes in that it does not
#' inherit from `AbstractMultiDataSet`, and, aside from being passed into other
#' `eda` class constructors, it is not intended to be used directly by the user.
#'
#' @section Arguments:
#' - `dat`: A data frame or matrix
#' - `key`: Character string or number indicating the row or column containg
#'      the dataset primary keys. If 'rownames', or 'colnames', row or column
#'      names will be used, respectively. If NULL, defaults to 'rownames' or
#'      'colnames', depending on orientation.
#' - `transposed`: Logical indicating whether the dataset orientation is
#'     transposed, relative to the primary dataset. For example, if a primary
#'     and secondary dataset both share the same column identifiers, but the
#'     second dataset is oriented with those id's as rows, then transposed
#'     should be set to TRUE to indicate such.
#' `xlab`: Character string containing x-axis label
#' `ylab`: Character string containing y-axis label
#'
#' @section Fields:
#'  - `dat`: Formatted data frame or matrix
#'
#' @importFrom R6 R6Class
#' @name EDADat
#' @export
#'
NULL

EDADat <- R6Class("EDADat",
    # ------------------------------------------------------------------------
    # public
    # ------------------------------------------------------------------------
    public = list(
        # TODO: add parameter to keep track of original data transposition;
        # have date show data in expected transposition state... 
        dat         = NULL,
        orientation = NULL,

        # EDADat constructor
        initialize = function(dat, orientation='rows', key=NULL,
                              transposed=FALSE, xlab=NULL, ylab=NULL) {
            # determine where primary keys are stored
            if (is.null(key)) {
                key <- ifelse(orientation == 'rows', 'rownames', 'colnames')
            }

            private$transposed <- transposed
            private$col_types  <- NULL

            # make sure a valid orientation / key combination is specified
            private$check_input(orientation, key)

            #
            # public
            #
            self$dat            <- private$format_data(dat, orientation, key)
            self$orientation    <- orientation

            #
            # private
            #
            private$key         <- key
            private$xlab        <- xlab
            private$ylab        <- ylab

            if (is.null(private$col_types)) {
                private$col_types <- sapply(self$dat, class)
            }

            # if data is provided in transposed form relative to main dataset,
            # transpose it so that datasets are in the same orientation
            if (transposed) {
                self$transpose()
            }
        },

        # subsamples dataset rows and/or columns in-place
        subsample = function(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL) {
            # indices to sample from
            row_ind <- 1:nrow(self$dat)
            col_ind <- 1:ncol(self$dat)

            # subsample rows
            if (!is.null(row_n)) {
                row_ind <- sample(row_ind, row_n)
            } else if (!is.null(row_ratio)) {
                row_ind <- sample(row_ind, round(row_ratio * length(row_ind)))
            }

            # subsample columns
            if (!is.null(col_n)) {
                col_ind <- sample(col_ind, col_n)
            } else if (!is.null(col_ratio)) {
                col_ind <- sample(col_ind, round(col_ratio * length(col_ind)))
            }

            # update data and metadata matrices
            self$dat <- self$dat[row_ind, col_ind, drop = FALSE]
        },

        # transpose data in-place
        transpose = function() {
            # Work-around for improved data frame support (2018/04/03)
            if (is.data.frame(self$dat)) {
                # store column types
                private$transposed = !private$transposed
            }

            # transpose data
            if (is.matrix(self$dat)) {
                self$dat <- t(self$dat)
            } else if (is.data.frame(self$dat)) {
                rn <- rownames(self$dat)
                cn <- colnames(self$dat)
                self$dat <- data.table::transpose(self$dat)
                rownames(self$dat) <- cn
                colnames(self$dat) <- rn
            }

            # Work-around: fix data frame column types
            if (is.data.frame(self$dat) && !private$transposed) {
                for(i in 1:ncol(self$dat)) {
                    cast_fxn <- get(paste0('as.', private$col_types[i]))
                    self$dat[,i] <- cast_fxn(self$dat[,i])
                }  
            }

            # update data orientation status
            self$orientation    <- ifelse(self$orientation == 'rows', 'columns', 'rows')

            # swap keynames if relevant
            if (private$key == 'rownames') {
                private$key = 'colnames'
            } else if (private$key == 'colnames') {
                private$key = 'rownames'
            }

            # swap axes labels
            xlab <- private$xlab
            private$xlab <- private$ylab
            private$ylab <- xlab
        }
    ),

    private = list(
        # Parameters
        key         = NULL,
        xlab        = NULL,
        ylab        = NULL,
        transposed  = NULL,
        col_types   = NULL,

        # validate input
        check_input = function(orientation, key) {
            # check key/orientation combo
            if ((orientation == 'rows' && key == 'colnames') ||
                (orientation == 'columns' && key == 'rownames')) {
                stop('Invalid combination of orientation and key specified.')
            }
        },

        # normalize row-oriented dataset row order
        #normalize_row_order = function(row_data, ids) {
        #    # iterate over datasets
        #    for (rdat in names(row_data)) {
        #        # if order is not the same as main dataset, reorder
        #        if (!(all(rownames(row_data[[rdat]]) == ids))) {
        #            # match row dataset row order to row names specified
        #            ind <- order(match(rownames(row_data[[rdat]]), ids))

        #            # return result
        #            row_data[[rdat]] <- row_data[[rdat]][ind,, drop = FALSE]
        #        }
        #    }
        #},

        # normalize column-oriented dataset column order
        #normalize_col_order = function(col_data, ids) {
        #    # iterate over datasets
        #    for (cdat in names(col_data)) {
        #        # if order is not the same as main dataset, reorder
        #        if (!(all(colnames(col_data[[cdat]]) == ids))) {
        #            # match column dataset column order to column names specified
        #            ind <- order(match(colnames(col_data[[cdat]]), ids))

        #            # return result
        #            col_data[[cdat]] <- col_data[[cdat]][ind,, drop = FALSE]
        #        }
        #    }
        #},
        get_key_index = function(dat, key, names_fxn=colnames) {
            # column number containing row ids specified, or
            # row number containing column ids
            if (is.numeric(key)) {
                key
            } else if (key %in% names_fxn(dat)) {
                # column name containing row ids specified, or
                # row name containing column ids
                which(names_fxn(dat) == key)
            }
        },

        # returns data with row / column names formatted as expected
        format_data = function(dat, orientation, key) {
            # row-oriented data
            if (orientation == 'rows') {
                if (key != 'rownames') {
                    ind <- private$get_key_index(dat, key, colnames)
                    rownames(dat) <- dat[, ind]
                    dat <- dat[, -ind]
                }
            } else {
                # column-oriented data
                if (key != 'colnames') {
                    ind <- private$get_key_index(dat, key, rownames)
                    colnames(dat) <- dat[ind, ]
                    dat <- dat[-ind, ]
                }
            }

            # ensure that formatted data has both row and column names;
            # used by some functions.
            if (is.null(rownames(dat))) {
                rownames(dat) <- paste0('row_', 1:nrow(dat))
            }
            if (is.null(colnames(dat))) {
                colnames(dat) <- paste0('col_', 1:ncol(dat))
            }
            dat
        }
    ),

    active = list(
        # return data in transposed order
        tdat = function() {
            self$transpose()
            dat <- self$dat
            self$transpose()
            dat
        }
    )
)
