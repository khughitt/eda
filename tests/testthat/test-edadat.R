#
# EDADat Unit Tests
#
context("EDADat")

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
dat <- data.frame(
    id=c('alt1', 'alt2', 'alt3'),
    x=rnorm(3),
    y=c('a', 'b', 'c')
)
rownames(dat) <- paste0('id', 1:3)

# transposed version of data
tdat <- data.table::transpose(dat)
rownames(tdat) <- colnames(dat)
colnames(tdat) <- rownames(dat)

##############################
# Tests
##############################

# test various construction options
test_that("initialization works", {
    expect_equal(EDADat$new(dat)$dat, dat)
    expect_equal(EDADat$new(dat)$fdat, dat)
    expect_equal(EDADat$new(dat, transposed=TRUE)$fdat, tdat)

    # rownames = 1
    expected <- dat[, -1]
    rownames(expected) <- dat[, 1]
    expect_equal(EDADat$new(dat, key=1)$fdat, expected)

    # rownames = 'id'
    expect_equal(EDADat$new(dat, key='id')$fdat, expected)

    # colnames = 1
    expected <- dat[-1, ]
    colnames(expected) <- dat[1, ]
    expect_equal(EDADat$new(dat, orientation='columns', key=1)$fdat, expected)

    # colnames = 'id1'
    expect_equal(EDADat$new(dat, orientation='columns', key='id1')$fdat, expected)
})

