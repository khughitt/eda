#' An S6 class representing an Expression dataset
#'
#' BioExprSet is a simple class for interacting with biological expression
#' data (e.g. from microarray or RNA-Seq experiments.) It can accept either
#' a matrix and optional column/row metadata dataframes, or a Bioconductor
#' ExpressionSet instance. Methods are provided for common transformations
#' and visualizations.
#'
#' @section Arguments:
#' \describe{
#'   \item{dat}{matrix|ExpressionSet An expression dataset with rows corresponding
#'       to genes, transcripts, probes, etc. and column corresponding to 
#'       individual samples.}
#'   \item{row_mdata}{Data frame with rows corresponding to the column names of 
#'       \code{dat}.}
#'   \item{row_mdata}{Data frame with rows corresponding to the row names of 
#'       \code{dat}.}
#' }
#'
#' @importFrom R6 R6Class
#' @name BioExprSet
#' @export
#'
NULL

BioExprSet <- R6::R6Class("BioExprSet",
    inherit = EDAMatrix,
    public = list(
        # BioExprSet constructor
        initialize = function(dat, 
                              row_mdata=NULL, col_mdata=NULL, 
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
                              row_maxn=Inf, row_max_ratio=1.0, row_ind=NULL,
                              col_maxn=Inf, col_max_ratio=1.0, col_ind=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) { 
            # verify input data type and call parent constructor
            private$check_input(dat)

            # for ExpressionSets, retrieve relevant data and metadata
            if (class(dat) == 'ExpressionSet') {
                if (is.null(col_mdata)) {
                    col_mdata <- pData(dat)
                }
                if (is.null(row_mdata)) {
                    row_mdata <- fData(dat)
                }
                dat <- exprs(dat)
            }

            super$initialize(dat, row_mdata, col_mdata, row_ids, col_ids,
                             row_mdata_ids, col_mdata_ids, row_color, row_shape,
                             row_labels, col_color, col_shape, col_labels,
                             row_maxn, row_max_ratio, row_ind,
                             col_maxn, col_max_ratio, col_ind,
                             color_pal, title, ggplot_theme)
        },

        #' Performs a counts-per-million (CPM) transformation.
        #'
        #' @return A CPM-transformed version of the BioExprSet instance.
        cpm = function() {
            obj <- self$clone()
            obj$dat <- sweep(obj$dat, 2, colSums(obj$dat), '/') * 1E6
            obj
        },

        diff_expr = function() {
        
        },

        #' Log2 transforms data (adding 1 to ensure finite results).
        #'
        #' @return Log2-transformed version of the expression data.
        log2p = function() {
            self$log(2, offset=1)
        },

        #' Plots a histogram of sample library sizes.
        #'
        #' Generates a histogram of sample library sizes (sum of expression
        #' levels within each column/sample). 
        #'
        #' @param bins Number of bins to use for histogram
        #'
        #' @return A ggplot instance.
        plot_libsize_hist = function(color=NULL, title=NULL) {
            dat <- data.frame(libsize=colSums(self$dat))

            styles <- private$get_geom_histogram_styles(color)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Library sizes: %s", private$title)
            }

            plt <- ggplot(aes(x=libsize), data=dat) +
                geom_histogram(styles$aes) +
                private$ggplot_theme()

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}
            plt
        },

        #' Plots bar graph of sample library sizes
        #'
        #' Generates a bargraph plot of sample library sizes (sum of expression
        #' levels within each column/sample). Each sample is shown as a separate
        #' bar in the plot.
        #'
        #' @return A ggplot instance.
        plot_libsize_bargraph = function(color=NULL, title=NULL) {
            # data frame with sample library sizes
            dat <- data.frame(libsize=colSums(self$dat))

            styles <- private$get_geom_histogram_styles(color)

            if (!is.null(styles$color)) {
                dat <- cbind(dat, color=styles$color)
            }
            if (is.null(title)) {
                title <- sprintf("Library sizes: %s", private$title)
            }

            plt <- ggplot(aes(x=rownames(libsizes), y=libsize), data=dat) +
                geom_bar(styles$aes, stat='identity') + 
                   xlab("Samples") +
                   ylab("Total expression") +
                   private$ggplot_theme() +
                   theme(axis.text.x=element_text(angle=90),
                         legend.text=element_text(size=8))

			# legend labels
			if (length(styles$labels) > 0) {
				plt <- plt + styles$labels
			}
            plt
        },

        #' Prints an overview of the object instance
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= BioExprSet (%s)\n", class(self$dat[,1])))
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        }
    ),
    private = list(
        #' Verifies that input data is an acceptable format.
        check_input = function(dat) {
            if(!is.matrix(dat) && (class(dat) != 'ExpressionSet')) {
                stop("Invalid input for BioExprSet: dat must be a matrix or ExpressionSet.")
            }
        }
    )
)

