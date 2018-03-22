#' An S6 class representing a dataframe dataset.
#'
#' EDADataFrame is a helper class for wrapping data frames, with optional
#' support for row and column datadata. Methods are provided for common
#' exploratory data analysis summary statistics, transformations, and 
#' visualizations.
#'
#' @section Usage:
#' ```
#' edf <- EDADataFrame$new(df, row_mdata=row_mdata_df, row_color='some_var')
#' edf$summary()
#' ``` 
#'
#' @section Arguments:
#' - `dat`: An m x n dataset.
#' - `row_mdata`: A matrix or data frame with rows corresponding to the row 
#'      names of `dat`
#' - `col_mdata`: A matrix or data frame with rows corresponding to the 
#'      column names of `dat`
#' - `row_ids`: Column name or number containing row identifiers. If set to
#'      `rownames` (default), row names will be used as identifiers.
#' - `col_ids`: Column name or number containing column identifiers. If set to
#'      `colnames` (default), column names will be used as identifiers.
#' - `row_mdata_ids`: Column name or number containing row metadata row 
#'      identifiers. If set to `rownames` (default), row names will be used 
#'      as identifiers.
#' - `col_mdata_ids`: Column name or number containing col metadata row 
#'      identifiers. If set to `rownames` (default), row names will be used 
#'      as identifiers.
#' - `row_color`: Row metadata field to use for coloring rowwise plot elements.
#' - `row_shape`: Row metadata field to use for determine rowwise plot 
#'      element shape.
#' - `row_labels`: Row metadata field to use when labeling plot points or
#'      other elements.
#' - `col_color`: Column metadata field to use for coloring columnwise plot elements.
#' - `col_shape`: Column metadata field to use for determine columnwise plot 
#'      element shape.
#' - `col_labels`: Column metadata field to use when labeling plot points or
#'      other elements.
#' - `color_pal`: Color palette to use for relevant plotting methods 
#'      (default: `Set1`).
#' - `title`: Text to use as a title or subtitle for plots.
#' - `ggplot_theme`: Default theme to use for ggplot2 plots 
#'      (default: `theme_bw`).
#'
#' @section Fields:
#'  - `dat`: Underlying data frame
#'  - `row_mdata`: Dataframe containing row metadata
#'  - `col_mdata`: Dataframe containing column metadata
#'
#' @section Methods:
#'  - `clear_cache()`: Clears EDADataFrame cache.
#'  - `clone()`: Creates a copy of the EDADataFrame instance.
#' - `detect_col_outliers(num_sd=2, avg='median', sim_method='pearson')`:
#'      Measures average pairwise similarities between all columns in the dataset.
#'      Outliers are considered to be those columns who mean similarity to
#'      all other columns is greater than `num_sd` standard deviations from the
#'      average of averages.
#' - `detect_row_outliers(num_sd=2, avg='median', sim_method='pearson')`:
#'      Measures average pairwise similarities between all rows in the dataset.
#'      Outliers are considered to be those rows who mean similarity to
#'      all other rows is greater than `num_sd` standard deviations from the
#'      average of averages.
#'  - `filter_col_outliers(num_sd=2, avg='median', sim_method='pearson')`: 
#'		Removes column outliers from the dataset. See `detect_col_outliers()` 
#'		for details of outlier detection approach.
#'  - `filter_row_outliers(num_sd=2, avg='median', sim_method='pearson')`: 
#'		Removes row outliers from the dataset. See `detect_row_outliers()` 
#'		for details of outlier detection approach.
#'  - `filter_cols(mask)`: Accepts a logical vector of length `ncol(obj$dat)`
#'		and returns a new EDADataFrame instance with only the columns associated
#'      with `TRUE` values in the mask.
#'  - `filter_rows(mask)`: Accepts a logical vector of length `nrow(obj$dat)`
#'		and returns a new EDADataFrame instance with only the rowsumns associated
#'      with `TRUE` values in the mask.
#'  - `impute(method='knn')`: Imputes missing values in the dataset and stores
#'		the result _in-place_. Currently only k-Nearest Neighbors (kNN) 
#'		imputation is supported.
#'  - `log(base=exp(1), offset=0)`: Log-transforms data.
#'  - `log1p()`: Log(x + 1)-transforms data.
#'  - `plot_densities(color=NULL, title="", ...)`: Plots densities for each 
#'		column in the dataset.
#'  - `print()`: Prints an overview of the object instance.
#'  - `subsample(row_n=NULL, col_n=NULL, row_ratio=NULL, col_ratio=NULL)`:
#'		Subsamples dataset rows and/or columns.
#'  - `summary(markdown=FALSE, num_digits=2)`: Summarizes overall 
#'		characteristics of a dataset.
#'  - `t()`: Transposes dataset rows and columns.
#'
#' @section Examples:
#' ```
#' library('eda')
#'
#' dat <- iris[,1:4]
#' row_mdata <- iris[,5,drop=FALSE]
#'
#' edf <- EDADataFrame$new(dat, row_mdata=row_mdata, row_color='Species')
#'
#' edf
#' edf$summary()
#' ``` 
#'
#' @importFrom R6 R6Class
#' @name EDADataFrame
#' @export
#'
NULL

EDADataFrame <- R6::R6Class("EDADataFrame",
    inherit = eda:::EDADataSet,
    public = list(
        #' EDADataFrame constructor
        initialize = function(dat, 
                              row_mdata=NULL, col_mdata=NULL, 
                              row_ids='rownames', col_ids='colnames',
                              row_mdata_ids='rownames', col_mdata_ids='rownames',
                              row_color=NULL, row_shape=NULL, row_labels=NULL,
                              col_color=NULL, col_shape=NULL, col_labels=NULL,
                              color_pal='Set1', title="", ggplot_theme=theme_bw) { 
            # verify input data type and call parent constructor
            private$check_input(dat)

            super$initialize(dat, row_mdata, col_mdata, row_ids, col_ids,
                             row_mdata_ids, col_mdata_ids, row_color, row_shape,
                             row_labels, col_color, col_shape, col_labels,
                             color_pal, title, ggplot_theme)
        },

        #' Prints an overview of the object instance
        print = function() {
            cat("=========================================\n")
            cat("=\n")
            cat(sprintf("= EDADataFrame (%s)\n", class(self$dat[,1])))
            cat("=\n")
            cat(sprintf("=   rows   : %d\n", nrow(self$dat)))
            cat(sprintf("=   columns: %d\n", ncol(self$dat)))
            cat("=\n")
            cat("=========================================\n")
        }
    ),
    private = list(
        #' Verifies that input data is of type dataframe
        check_input = function(dat) {
            if(!is.data.frame(dat)) {
                stop("Invalid input for EDADataFrame: dat must be a dataframe")
            }
        }
    )
)

