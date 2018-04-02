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
#'  - `dat`: Original dataset
#'  - `fdat`: Dataset formatted in the manner expected by `eda` classes.
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
        dat         = NULL,
        orientation = NULL,

        # EDADat constructor
        initialize = function(dat, orientation='rows', key=NULL,
                              transposed=FALSE, xlab=NULL, ylab=NULL) {
            # original data frame or matrix
            self$dat <- dat
            self$orientation <- orientation

            private$key         <- key
            private$transposed  <- transposed
            private$xlab        <- xlab
            private$ylab        <- ylab

            if (is.null(key)) {
                private$key <- ifelse(orientation == 'rows', 'rownames', 'colnames')
            }

            private$check_input()

            # determine data orientation

        },

        # subsamples dataset rows and/or columns in-place
        subsample = function(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL) {
            dat <- self$dat

            # indices to sample from
            row_ind <- 1:nrow(dat)
            col_ind <- 1:ncol(dat)

            # exclude primary key, if present
            if (self$orientation == 'rows') {
                if (private$key != 'rownames') {
                    # for row-oriented data, exclude column containing id's
                    col_ind <- col_ind[-private$get_row_key_idx()]
                }
            } else {
                # column-oriented data
                if (private$key != 'colnames') {
                    # for column-oriented data, exclude row containing id's
                    row_ind <- row_ind[-private$get_col_key_idx()]
                }
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

            # update data and metadata matrices
            obj$dat <- obj$dat[row_ind, col_ind, drop = FALSE]
        },

        transpose = function() {
            private$transposed <- !private$transposed
        }
    ),

    private = list(
        # Parameters
        key         = NULL,
        transposed  = NULL,
        xlab        = NULL,
        ylab        = NULL,

        # validate input
        check_input = function() {
            # check key/orientation combo
            if ((self$orientation == 'rows' && private$key == 'colnames') ||
                (self$orientation == 'columns' && private$key == 'rownames')) {
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

        get_row_key_index = function () {
            # column number containing row ids specified
            if (is.numeric(private$key)) {
                private$key
            } else if (private$key %in% colnames(self$dat)) {
                # column name containing row ids specified
                which(colnames(self$dat) == private$key)
            }
        },

        get_col_key_index = function () {
            # column number containing column ids specified
            if (is.numeric(private$key)) {
                private$key
            } else if (private$key %in% rownames(self$dat)) {
                # row name containing column ids specified
                which(rownames(self$dat) == private$key)
            }
        },

        # returns data with row / column names formatted as expected
        get_formatted_data = function () {
            dat <- self$dat

            # row-oriented data
            if (self$orientation == 'rows') {
                if (private$key != 'rownames') {
                    ind <- private$get_row_key_index()
                    rownames(dat) <- dat[, ind]
                    dat <- dat[, -ind]
                }
            } else {
                # column-oriented data
                if (private$key != 'colnames') {
                    ind <- private$get_col_key_index()
                    colnames(dat) <- dat[ind, ]
                    dat <- dat[-ind, ]
                }
            }

            dat
        }
    ),

    active = list(
        # return formatted dataset
        fdat = function() {
            dat <- private$get_formatted_data()

            # return dataset in proper orientation
            if (private$transposed) {
                if (is.matrix(dat)) {
                    dat <- t(dat)
                } else if (is.data.frame(dat)) {
                    rn <- rownames(dat)
                    cn <- colnames(dat)
                    dat <- data.table::transpose(dat)
                    rownames(dat) <- cn
                    colnames(dat) <- rn
                }
            }
            dat
        },

        # return formatted data in transposed order
        tdat = function() {
            dat <- private$get_formatted_data()

            # if dataset is not already transposed, transpose it
            if (!private$transposed) {
                if (is.matrix(dat)) {
                    dat <- t(dat)
                } else if (is.data.frame(dat)) {
                    rn <- rownames(dat)
                    cn <- colnames(dat)
                    dat <- data.table::transpose(dat)
                    rownames(dat) <- cn
                    colnames(dat) <- rn
                }
            }
            dat
        }
    )
)
