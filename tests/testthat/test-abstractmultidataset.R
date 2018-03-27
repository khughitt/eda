#
# AbstractMultiDataSet Unit Tests
#
context("AbstractMultiDataSet")

##############################
# Setup
##############################

# ensure reproducibility
set.seed(1)

#
# generate a single dataset with rows corresponding to:
#
#   y =  x
#   y =  1
#   y = -x
#
# this will then be used for both "dat1" and "dat2", allowing us to easily
# check the expected outputs for various similarity measures.
#
num_cols <- 5

dat1 <- rbind(
    1:num_cols,
    rep(1, num_cols),
    num_cols:1
)
dat2 <- dat1

rownames(dat1) <- sprintf('obs%02d', 1:3)
rownames(dat2) <- sprintf('obs%02d', 4:6)
colnames(dat1) <- sprintf('var%02d', 1:5)
colnames(dat2) <- sprintf('var%02d', 1:5)

# next, we will create instances of two subclasses of AbstractMultiDataset,
# in order to test various inherited methods.
sds <- EDAMatrix$new(rbind(dat1, dat2))
mds <- EDAMultiDataSet$new(EDAMatrix$new(dat1), EDAMatrix$new(dat2))

# expected correlation result (zero-variance middle rows removed)

##############################
# Tests
##############################

# EDAMultiDataSet construction
test_that("initialization works", {
    expect_equal(mds$datasets[[1]]$dat, dat1)
    expect_equal(mds$datasets[[2]]$dat, dat2)
})

# Cross-correlation
test_that("Correlation measures work", {
    # Cross-correlation sub-indices
    row_ind <- 1:2
    col_ind <- 3:4

    # Input datasets with zero-variance entries removed
    d1 <- dat1[-2, ]
    d2 <- dat2[-2, ]

    # combined dataset
    combined_dat <- t(rbind(d1, d2))

    # Pearson correlation matrix
    cor_mat <- cor(combined_dat)

    # Mutual information matrix
    mut_mat <- mpmi::cmi(combined_dat)$bcmi
    rownames(mut_mat) <- colnames(combined_dat)
    colnames(mut_mat) <- colnames(combined_dat)

    # Linear model (perfect fits, so r^2 are all 100's)
    lm_mat <- matrix(100, 4, 4)
    rownames(lm_mat) <- colnames(combined_dat)
    colnames(lm_mat) <- colnames(combined_dat)

    # cor() operates on columns, so we transpose before comparing
    expect_equal(sds$t()$cor(method = 'pearson'), cor_mat)
    expect_equal(sds$t()$cor(method = 'mi'),      mut_mat)
    expect_equal(expect_warning(sds$t()$cor(method = 'lm')), lm_mat)

    # cross_cor() compares rows from two datasets
    expect_equal(mds$cross_cor(method = 'pearson'), cor_mat[row_ind, col_ind])
    expect_equal(mds$cross_cor(method = 'mi'),      mut_mat[row_ind, col_ind])
    expect_equal(expect_warning(mds$cross_cor(method = 'lm')), lm_mat[row_ind, col_ind])
})
