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

        # style properties
        row_color    = NULL,
        row_shape    = NULL,
        row_label    = NULL,
        row_edat     = NULL,
        col_color    = NULL,
        col_shape    = NULL,
        col_label    = NULL,
        col_edat     = NULL,

        # EDADat constructor
        initialize = function(dat, xid=NULL, yid=NULL,
                              row_names='rownames', col_names='colnames',
                              row_color=NULL, row_shape=NULL, row_label=NULL, row_edat=NULL,
                              col_color=NULL, col_shape=NULL, col_label=NULL, col_edat=NULL,
                              xlab=NULL, ylab=NULL) {

            # assign random axis identifiers if none specified
            if (is.null(xid)) {
                xid <- private$get_hash()
            }
            if (is.null(yid)) {
                yid <- private$get_hash()
            }

            # public properties
            self$xid  <- xid
            self$yid  <- yid
            self$xlab <- xlab
            self$ylab <- ylab

            self$row_color  <- row_color   
            self$row_shape  <- row_shape   
            self$row_label  <- row_label   
            self$row_edat   <- row_edat
            self$col_color  <- col_color   
            self$col_shape  <- col_shape   
            self$col_label  <- col_label   
            self$col_edat   <- col_edat

            # store data
            private$data <- private$format_data(dat, row_names, col_names)
        },

        # returns the vector for the specified axis and row/column name
        # "other_axis" indicates that the row/column to be retrieved is
        # on the opposite axis as to that which was specified by the "axis"
        # parameter.
        # If not axis name is specified, row or column names are returned
        get = function(axis, name=NULL, other_axis=FALSE) {
            if (!axis %in% c(self$xid, self$yid)) {
                stop("Invalid axis ID specified.")
            }
            # matching axis is x (rows)
            if ((axis == self$xid  && !other_axis) || 
                (axis == self$yid && other_axis) || 
                (is.data.frame(private$data) && private$transposed)) {

                # for row requests on transposed data frames, retrieve column
                # or rownames from original data to preserve type
                if (is.data.frame(private$data) && private$transposed && (axis != self$yid)) {
                    if (is.null(name)) {
                        rownames(private$data)
                    } else {
                        private$data[, name] 
                    }
                } else {
                    # otherwise return row or column names
                    # unlist ensures that 1d data frames are converted to vectors;
                    # as.vector drops any associated names
                    if (is.null(name)) {
                        colnames(private$data)
                    } else {
                        as.vector(unlist(private$data[name, ]))
                    }
                }
            } else {
                # return rownames or a column in matrix / data frame
                if (is.null(name)) {
                    rownames(private$data)
                } else {
                    as.vector(unlist(private$data[, name]))
                }
            }
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

            # preserve row and column order
            row_ind <- sort(row_ind)
            col_ind <- sort(col_ind)

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
            
            # swap row and column style elements
            row_color <- self$row_color
            row_shape <- self$row_shape
            row_label <- self$row_label
            row_edat  <- self$row_edat

            self$row_color <- self$col_color
            self$row_shape <- self$col_shape
            self$row_label <- self$col_label
            self$row_edat  <- self$col_edat

            self$col_color <- row_color
            self$col_shape <- row_shape
            self$col_label <- row_label
            self$col_edat  <- row_edat
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
            # convert tibbles to normal data frames
            if (class(dat)[1] == 'tbl_df') {
                dat <- as.data.frame(dat)
            }

            # normalize rownames
            if (row_names != 'rownames') {
                ind <- private$get_names_index(dat, row_names, colnames)
                if (sum(duplicated(dat[, ind])) > 0) {
                    stop("Row identifiers must be unique.") 
                }
                rownames(dat) <- dat[, ind]
                dat <- dat[, -ind]
            }
            # normalize colnames
            if (col_names != 'colnames') {
                ind <- private$get_names_index(dat, col_names, rownames)
                if (sum(duplicated(dat[ind, ])) > 0) {
                    stop("Column identifiers must be unique.")
                }
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
        },

        get_hash = function(len=6) {
            paste0(sample(c(0:9, letters, toupper(letters)), len), collapse='')
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
