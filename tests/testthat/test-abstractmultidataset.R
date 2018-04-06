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
dats <- rbind(dat1, dat2)

# create a third dataset for evaluate behavior related to datasets that share
# either row or column ids
dat3 <- matrix(rnorm(15), nrow=5)
rownames(dat3) <- letters[1:5]
colnames(dat3) <- sprintf('col%02d', 1:3)

sds <- EDAMatrix$new(dats)
mds <- EDAMultiMatrix$new(list(a=dat1, b=dat2, c=dat3))

# expected correlation result (zero-variance middle rows removed)

##############################
# Tests
##############################

# AbstractMultiDataSet construction
test_that("initialization works", {
    expect_equal(mds$datasets[[1]], dat1)
    expect_equal(mds$datasets[[2]], dat2)
})

# Filtering
test_that("Filtering works", {
    expect_equal(sds$filter_rows(rep(c(T, F), 3))$dat, dats[rep(c(T, F), 3),])
    expect_equal(mds$filter_cols(1, c(F, T, F, T, F))$datasets[[1]],
                 dat1[, c(F, T, F, T, F)])
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

# Plot styles
test_that("Handling of plot styles works", {
    # EDADat$get()
    expect_equal(mds$edat[['a']]$get('x', 'obs01'), as.vector(dat1[1, ]))
    expect_equal(mds$edat[['a']]$get('y', 'var01'), as.vector(dat1[, 1]))
    expect_equal(mds$edat[['a']]$get('x', 'var01', other_axis=TRUE), as.vector(dat1[, 1]))
    expect_equal(mds$edat[['a']]$get('y', 'obs01', other_axis=TRUE), as.vector(dat1[1, ]))
})

# Sub-sampling
test_that("Sub-sampling works", {
    # 10 x 10 matrix
    mat <- matrix(rep(0, 100), 10)

    # EDAMatrix$subsample()
    edm <- EDAMatrix$new(mat)
    expect_equal(dim(edm$subsample(3, 5)$dat), c(3,5))

    # EDAMultiMatrix$subsample()
    emm <- EDAMultiMatrix$new(list(foo=mat))
    expect_equal(dim(edm$subsample(row_ratio=0.6, col_ratio=0.4)$dat), c(6, 4))
})
