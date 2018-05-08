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


# |      | var01| var02| var03| var04| var05|
# |:-----|-----:|-----:|-----:|-----:|-----:|
# |obs01 |     1|     2|     3|     4|     5|
# |obs02 |     1|     1|     1|     1|     1|
# |obs03 |     5|     4|     3|     2|     1|
dat1 <- rbind(
    1:num_cols,
    rep(1, num_cols),
    num_cols:1
)
rownames(dat1) <- sprintf('obs%02d', 1:3)
colnames(dat1) <- sprintf('var%02d', 1:5)

# |      | var01| var02| var03| var04| var05|                                                               
# |:-----|-----:|-----:|-----:|-----:|-----:|                                                               
# |obs04 |     1|     2|     3|     4|     5|                                                               
# |obs05 |     1|     1|     1|     1|     1|                                                               
# |obs06 |     5|     4|     3|     2|     1| 
dat2 <- dat1
rownames(dat2) <- sprintf('obs%02d', 4:6)
colnames(dat2) <- sprintf('var%02d', 1:5)

# randomize order of columns for one of the datasets
# to ensure proper handling of mixed column order
dat2_col_order <- sample(5)
dat2 <- dat2[, dat2_col_order]


# dat3 includes column labels for dat1/2
#
# |      |      col01|      col02|      col03|var_labels  |
# |:-----|----------:|----------:|----------:|:-----------|
# |var02 |  1.2724293| -0.0057672| -0.2894616|var_label02 |
# |var05 |  0.4146414|  2.4046534| -0.2992151|var_label05 |
# |var04 | -1.5399500|  0.7635935| -0.4115108|var_label04 |
# |var03 | -0.9285670| -0.7990092|  0.2522234|var_label03 |
# |var01 | -0.2947204| -1.1476570| -0.8919211|var_label01 |
dat3 <- data.frame(
    col01=rnorm(5),
    col02=rnorm(5),
    col03=rnorm(5),
    row.names=sample(colnames(dat1))
)
dat3 <- cbind(dat3, var_labels=sub('var', 'var_label', rownames(dat3)))

# dat4 includes row labels for dat2
#
# |      |obs_var1 | obs_var2|obs_labels |
# |:-----|:--------|--------:|:----------|
# |obs05 |a        |        3|obs05_lab  |
# |obs06 |b        |        1|obs06_lab  |
# |obs04 |c        |        2|obs04_lab  |
dat4 <- data.frame(
    obs_var1   = letters[1:3],
    obs_var2   = sample(3),
    obs_labels = c('obs05_lab', 'obs06_lab', 'obs04_lab'),
    row.names  = c('obs05', 'obs06', 'obs04')
)

# next, we will create instances of two subclasses of AbstractMultiDataset,
# in order to test various inherited methods.
dats <- rbind(dat1, dat2[, colnames(dat1)])

sds <- EDAMatrix$new(dats)

mds <- EDAMultiMatrix$new(list(a=dat1,
                               b=EDADat$new(dat2, xid = 'b_x', yid='a_y',
                                            row_label = 'obs_labels', row_edat = 'd',
                                            col_label = 'var_labels', col_edat = 'c'),
                               c=EDADat$new(dat3, xid = 'a_y'),
                               d=EDADat$new(dat4, xid = 'b_x')))

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

    # Input datasets with zero-variance rows removed
    d1 <- dat1[-2, ]
    d2 <- dat2[-2, ]

    # combined dataset
    combined_dat <- t(rbind(d1, d2[, colnames(dat1)]))

    #
    # Similarity matrices
    #

    # Pearson correlation matrix
    cor_mat <- cor(combined_dat)

    # Mutual information matrix
    mut_mat <- mpmi::cmi(combined_dat)$bcmi
    rownames(mut_mat) <- colnames(combined_dat)
    colnames(mut_mat) <- colnames(combined_dat)

    # Linear model matrix (perfect fits, so r^2 are all 100's)
    lm_mat <- matrix(100, 4, 4)
    rownames(lm_mat) <- colnames(combined_dat)
    colnames(lm_mat) <- colnames(combined_dat)

    # cor() operates on columns, so we transpose before comparing
    expect_equal(sds$t()$cor(meas = 'pearson'), cor_mat)
    expect_equal(sds$t()$cor(meas = 'cmi'),     mut_mat)
    expect_equal(expect_warning(sds$t()$cor(meas = 'lm')), lm_mat)

    #
    # Cross correlation
    #

    # Pearson correlation
    res <- mds$cross_cor(meas = 'pearson')
    expect_equal(res$edat[['a_b_pearson']]$dat, cor_mat[row_ind, col_ind])

    # Mutual information
    res <- mds$cross_cor(meas = 'cmi')
    expect_equal(res$edat[['a_b_cmi']]$dat, mut_mat[row_ind, col_ind])

    # Linear model
    # Gives an expected warning because of the perfect fit in fake data
    res <- expect_warning(mds$cross_cor(meas = 'lm'))
    expect_equal(res$edat[['a_b_lm']]$dat, lm_mat[row_ind, col_ind])

    # Check handling of axes and transposed data

    # check axis ids for cross cor matrix; should have axes corresponding
    # to the non-shared axes in original two datasets
    expect_equal(res$edat[['a_b_lm']]$xid, 'a_x')
    expect_equal(res$edat[['a_b_lm']]$yid, 'b_x')

    # should get the same result, even if relative dataset orientation differs
    res <- mds$t('a')$cross_cor(meas = 'pearson')

    expect_equal(res$edat[['a_b_pearson']]$dat, cor_mat[row_ind, col_ind])
    expect_equal(res$edat[['a_b_pearson']]$xid, 'a_x')
    expect_equal(res$edat[['a_b_pearson']]$yid, 'b_x')
})

# Plot styles
test_that("Handling of plot styles works", {
    expect_equal(mds$edat[['a']]$get('a_x', 'obs01'), as.vector(dat1[1, ]))
    expect_equal(mds$edat[['a']]$get('a_y', 'var01'), as.vector(dat1[, 1]))
    expect_equal(mds$edat[['a']]$get('a_x', 'var01', other_axis=TRUE), as.vector(dat1[, 1]))
    expect_equal(mds$edat[['a']]$get('a_y', 'obs01', other_axis=TRUE), as.vector(dat1[1, ]))

    # Row / Column labels (have to check private methods)
    private <- mds$.__enclos_env__$private

    expect_equal(private$get_row_labels('b', label_var=FALSE),    rownames(dat2))
    expect_equal(private$get_row_labels('b', label_var='obs_labels'), sort(dat4$obs_labels))
    expect_equal(private$get_col_labels('b', label_var=FALSE),    colnames(dat2))
    expect_equal(private$get_col_labels('b', label_var='var_labels'), sort(dat3$var_labels)[dat2_col_order])
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
