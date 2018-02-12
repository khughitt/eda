#'
#' EDAMatrix
#'
EDAMatrix <- R6Class("EDAMatrix",
    inherit = EDADataSet,
    public = list(
        # EDAMatrix constructor
        initialize = function(dat, row_metadata=NULL, col_metadata=NULL,
                              row_maxn=Inf, row_maxr=1.0,
                              col_maxn=Inf, col_maxr=1.0,
                              color_var=NULL, shape_var=NULL, label_var=NULL,
                              color_pal='Set1') { 
            # verify input data type and call parent constructor
            private$check_input(dat)

            super$initialize(dat, row_metadata, col_metadata, row_maxn, row_maxr,
                             col_maxn, col_maxr, color_var, shape_var, label_var,
                             color_pal)
        },

        # class greeting
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= EDAMatrix (%s)\n", class(self$dat[,1])))
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        }
    ),
    private = list(
        # verify that input data is of type matrix
        check_input = function(dat) {
            if(!is.matrix(dat)) {
                stop("Invalid input for EDAMatrix: dat must be a matrix.")
            }
        }
    )
)

