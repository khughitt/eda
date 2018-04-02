#
# EDAMatrix Unit Tests
#
context("EDAMatrix")

##############################
# Setup
##############################

# ensure reproducibility
set.seed(1)

# generate a test dataset
num_rows <- 10
num_cols <- 5

mat <- matrix(rnorm(num_rows * num_cols), ncol = num_cols)
rownames(mat) <- sprintf('obs%02d', 1:num_rows)
colnames(mat) <- sprintf('var%02d', 1:num_cols)

# row and column metadata
row_mdat <- data.frame(obs_prop1 = sample(10, num_rows, replace = TRUE),
                       obs_prop2 = factor(sample(letters[1:5], num_rows, replace = TRUE)))
rownames(row_mdat) <- rownames(mat)

col_mdat <- data.frame(var_prop1 = sample(5, num_cols, replace = TRUE),
                       var_prop2 = factor(sample(c(T, F), num_cols, replace = TRUE)),
                       var_prop3 = rnorm(num_cols))
col_mdat <- data.table::transpose(col_mdat)
colnames(col_mdat) <- colnames(mat)
rownames(col_mdat) <- paste0('var_prop', 1:3)

# transposed versions
trow_mdat <- data.table::transpose(row_mdat)
rownames(trow_mdat) <- colnames(row_mdat)
colnames(trow_mdat) <- rownames(row_mdat)

tcol_mdat <- data.table::transpose(col_mdat)
colnames(tcol_mdat) <- rownames(col_mdat)
rownames(tcol_mdat) <- colnames(col_mdat)

# create EDAMatrix
edm <- EDAMatrix$new(mat, row_mdata = row_mdat, col_mdata = col_mdat)

##############################
# Tests
##############################

# EDAMatrix construction
test_that("initialization works", {
    expect_equal(edm$dat, mat)
    expect_equal(edm$row_mdata, row_mdat)
    expect_equal(edm$col_mdata, col_mdat)
})

test_that("getters work", {
    expect_equal(edm$get('dat'), edm$fget('dat'))
})

# transpose
test_that("transposition works", {
    expect_equal(edm$t()$t()$dat, edm$dat)
    expect_equal(edm$t()$t()$row_mdata, edm$row_mdata)
    expect_equal(edm$t()$t()$col_mdata, edm$col_mdata)
    expect_equal(edm$t()$dat, t(edm$dat))
    expect_equal(edm$t()$col_mdata, trow_mdat)
    expect_equal(edm$t()$row_mdata, tcol_mdat)
    expect_equal(class(edm$t())[1], 'EDAMatrix')
    expect_equal(class(edm$t()$dat), class(edm$dat))
    expect_equal(rownames(edm$t()$dat), colnames(edm$dat))
    expect_equal(colnames(edm$t()$dat), rownames(edm$dat))
    expect_equal(rownames(edm$t()$row_mdata), colnames(edm$col_mdata))
    expect_equal(colnames(edm$t()$row_mdata), rownames(edm$col_mdata))
})

# impute
test_that("imputation works", {
    # matrix with missing data
    na_mat <- mat
    na_mat[sample(50, 10)] <- NA

    # data frame with missing data
    na_dat <- head(iris, 20)
    na_dat[1, 1] <- NA
    na_dat[8, 3] <- NA
    na_dat[12,5] <- NA

    # test each available imputation method
    IMPUTE_METHODS = c('knn')

    for (method in IMPUTE_METHODS) {
        # test matrix imputation
        edm <- EDAMatrix$new(na_mat, row_mdata = row_mdat, col_mdata = col_mdat)
        edm$impute(method=method)

        expect_equal(class(edm$dat), 'matrix')
        expect_equal(rownames(edm$dat), rownames(na_mat))
        expect_equal(colnames(edm$dat), colnames(na_mat))
        expect_equal(sum(is.na(edm$dat)), 0)

        # test data frame imputation
        #edf <- EDADataFrame$new(na_mat, row_mdata = row_mdat, col_mdata = col_mdat)
        #edf$impute(method=method)

        #expect_equal(class(edf$dat), 'matrix')
        #expect_equal(rownames(edf$dat), rownames(na_mat))
        #expect_equal(colnames(edf$dat), colnames(na_mat))
        #expect_equal(sum(is.na(edf$dat)), 0)
    }
})
