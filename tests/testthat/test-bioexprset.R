#
# BioExprSet Unit Tests
#
context("BioExprSet")

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

col_mdat <- as.data.frame(t(data.frame(var_prop1 = sample(5, num_cols, replace = TRUE),
                                       var_prop2 = factor(sample(c(T,F), num_cols, replace = TRUE)),
                                       var_prop3 = rnorm(num_cols))))
colnames(col_mdat) <- colnames(mat)

# create EDAMatrix
bset <- BioExprSet$new(mat, row_mdata = row_mdat, col_mdata = col_mdat)

##############################
# Tests
##############################

test_that("pathway statistics works", {
    expect_error(bset$annotation_stats(NULL, stat='nonexistent_stat'), 
                 "Invalid statistic specified.", fixed = TRUE)
})

