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
        # public properties
        xid  = NULL,
        yid  = NULL,
        xlab = NULL,
        ylab = NULL,

        # EDADat constructor
        initialize = function(dat, xid='x', yid='y',
                              row_names='rownames', col_names='colnames',
                              xlab=NULL, ylab=NULL) {
            # properties
            self$xid  <- xid
            self$yid  <- yid
            self$xlab <- xlab
            self$ylab <- ylab

            # store data
            private$data <- private$format_data(dat, row_names, col_names)
        },

        # subsamples dataset rows and/or columns in-place
        subsample = function(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL) {
            # get underlying data
            dat <- private$data

            # indices to sample from
            if (is.data.frame(private$data) && private$transposed) {
                # for transposed data frames, swap row and column indices
                row_ind <- 1:ncol(dat)
                col_ind <- 1:nrow(dat)
            } else {
                # otherwise operate as-is
                row_ind <- 1:nrow(dat)
                col_ind <- 1:ncol(dat)
            }

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

            # update data
            if (is.data.frame(private$data) && private$transposed) {
                # transposed data frames
                private$data <- dat[col_ind, row_ind, drop = FALSE]
            } else {
                # everything else
                private$data <- dat[row_ind, col_ind, drop = FALSE]
            }
        },

        # transpose data;
        # matrices are transposed in-places while data frames are left alone,
        # but have a transposition flag toggled
        transpose = function() {
            # for data frames, keep track of transposition status
            if (is.data.frame(private$data)) {
                private$transposed = !private$transposed
            }

            # transpose matrix data
            if (is.matrix(private$data)) {
                private$data <- t(private$data)
            }

            # swap axis ids and labels
            xid <- self$xid
            self$xid <- self$yid
            self$yid <- xid

            xlab <- self$xlab
            self$xlab <- self$ylab
            self$ylab <- xlab
        }
    ),

    # ------------------------------------------------------------------------
    # private
    # ------------------------------------------------------------------------
    private = list(
        # Parameters
        data        = NULL,
        transposed  = FALSE,

        get_names_index = function(dat, key, names_fxn=colnames) {
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
        format_data = function(dat, row_names, col_names) {
            # normalize rownames
            if (row_names != 'rownames') {
                ind <- private$get_names_index(dat, row_names, colnames)
                rownames(dat) <- dat[, ind]
                dat <- dat[, -ind]
            }
            # normalize colnames
            if (col_names != 'colnames') {
                ind <- private$get_names_index(dat, col_names, rownames)
                colnames(dat) <- dat[ind, ]
                dat <- dat[-ind, ]
            }

            # ensure that formatted data has both row and column names;
            # required by some functions.
            if (is.null(rownames(dat))) {
                rownames(dat) <- paste0('row_', 1:nrow(dat))
            }
            if (is.null(colnames(dat))) {
                colnames(dat) <- paste0('col_', 1:ncol(dat))
            }
            dat
        }
    ),

    # ------------------------------------------------------------------------
    # active
    # ------------------------------------------------------------------------
    active = list(
        # return data in expected orientation
        dat = function(value) {
            # return data
            if (missing(value)) {
                # for transposed data frames, transpose on the fly 
                if (is.data.frame(private$data) && private$transposed) {
                    rn <- rownames(private$data)
                    cn <- colnames(private$data)
                    dat <- data.table::transpose(private$data)
                    rownames(dat) <- cn
                    colnames(dat) <- rn
                    return(dat)
                } else {
                    # otherwise, return as-is
                    return(private$data)
                }
            } else {
                # store updated matrix / data frame
                private$data <- value

                # for data frames, reset transposed flag
                if (is.data.frame(value)) {
                    private$transposed <- FALSE
                }
            }
        },
        # return data in transposed order
        tdat = function() {
            self$transpose()
            dat <- self$dat
            self$transpose()
            dat
        }
    )
)