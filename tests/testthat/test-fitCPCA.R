context("Test fitting of contrastive PCA")
library(BiocParallel)
library(Matrix)
library(DelayedArray)

# initialization of inputs for testing function
target <- toy_df[, -31]
background <- background_df
center <- TRUE
scale <- TRUE
contrasts <- c(0, 1, 10, 100)
c_contrasts <- contrastiveCov(target, background, contrasts, TRUE, TRUE)
n_eigen <- 2
n_medoids <- 2
eigdecomp_tol <- 1e-10
eigdecomp_iter <- 1000

# data frame
fit <- fitCPCA(
  target, center, scale, c_contrasts, contrasts,
  n_eigen, n_medoids, eigdecomp_tol, eigdecomp_iter
)

# sparse matrix
dgC_target <- as(as.matrix(toy_df[, -31]), "dgCMatrix")
dgC_background <- as(as.matrix(background_df), "dgCMatrix")
c_contrasts <- contrastiveCov(dgC_target, dgC_background, contrasts, TRUE, TRUE)
dgC_fit <- fitCPCA(
  dgC_target, center, scale, c_contrasts, contrasts,
  n_eigen, n_medoids, eigdecomp_tol, eigdecomp_iter
)

# delayed array
dm_target <- DelayedArray(toy_df[, -31])
dm_background <- DelayedArray(background_df)
c_contrasts <- contrastiveCov(dm_target, dm_background, contrasts, TRUE, TRUE)
dm_fit <- fitCPCA(
  dm_target, center, scale, c_contrasts, contrasts,
  n_eigen, n_medoids, eigdecomp_tol, eigdecomp_iter
)

test_that(paste(
  "The length of the ouput list is 4 (rotations, projections,",
  "contrasts and penalties)"
), {
  expect_equal(length(fit), 4)
  expect_equal(length(dgC_fit), 4)
  expect_equal(length(dm_fit), 4)
})

test_that("The number of medoids output equals the initial specified option", {
  expect_equal(length(fit$contrast), 2)
  expect_equal(length(dgC_fit$contrast), 2)
  expect_equal(length(dm_fit$contrast), 2)
})

test_that("Dimensions of the rotation matrices are ncol(target) by n_eigen", {
  expect_equal(nrow(fit$rotation[[1]]), ncol(target))
  expect_equal(ncol(fit$rotation[[1]]), n_eigen)
  expect_equal(nrow(fit$rotation[[2]]), ncol(target))
  expect_equal(ncol(fit$rotation[[2]]), n_eigen)
  expect_equal(nrow(dgC_fit$rotation[[1]]), ncol(target))
  expect_equal(ncol(dgC_fit$rotation[[1]]), n_eigen)
  expect_equal(nrow(dgC_fit$rotation[[2]]), ncol(target))
  expect_equal(ncol(dgC_fit$rotation[[2]]), n_eigen)
  expect_equal(nrow(dm_fit$rotation[[1]]), ncol(target))
  expect_equal(ncol(dm_fit$rotation[[1]]), n_eigen)
  expect_equal(nrow(dm_fit$rotation[[2]]), ncol(target))
  expect_equal(ncol(dm_fit$rotation[[2]]), n_eigen)
})

test_that("The projection matrices are nrow(target) by n_eigen", {
  expect_equal(nrow(fit$x[[1]]), nrow(target))
  expect_equal(ncol(fit$x[[1]]), n_eigen)
  expect_equal(nrow(fit$x[[2]]), nrow(target))
  expect_equal(ncol(fit$x[[2]]), n_eigen)
  expect_equal(nrow(dgC_fit$x[[1]]), nrow(target))
  expect_equal(ncol(dgC_fit$x[[1]]), n_eigen)
  expect_equal(nrow(dgC_fit$x[[2]]), nrow(target))
  expect_equal(ncol(dgC_fit$x[[2]]), n_eigen)
  expect_equal(nrow(dm_fit$x[[1]]), nrow(target))
  expect_equal(ncol(dm_fit$x[[1]]), n_eigen)
  expect_equal(nrow(dm_fit$x[[2]]), nrow(target))
  expect_equal(ncol(dm_fit$x[[2]]), n_eigen)
})

test_that("The penalties output vector is a zero vector", {
  expect_equal(fit$penalty, rep(0, n_medoids))
  expect_equal(dgC_fit$penalty, rep(0, n_medoids))
  expect_equal(dm_fit$penalty, rep(0, n_medoids))
})

test_that("The parallelized analog matches the sequential variant exactly", {
  # fit parallelized variant
  register(SerialParam())
  bpfit <- bpFitCPCA(
    target, center, scale, c_contrasts, contrasts,
    n_eigen, n_medoids, eigdecomp_tol, eigdecomp_iter
  )
  expect_equal(fit, bpfit)
  expect_equal(dgC_fit, bpfit)
  expect_equivalent(dm_fit, bpfit)
})

test_that("Sparse and DelayedMatrix versions identical to dataframe version", {
  
  # multiple contrasts
  expect_equivalent(dm_fit, fit)
  expect_equivalent(dm_fit, dgC_fit)
  
  # single contrasts
  fit <- fitCPCA(
    target, center, scale, list(c_contrasts[[2]]), contrasts[2],
    n_eigen, n_medoids, eigdecomp_tol, eigdecomp_iter
  )
  dgC_fit <- fitCPCA(
    dgC_target, center, scale, list(c_contrasts[[2]]), contrasts[2],
    n_eigen, n_medoids, eigdecomp_tol, eigdecomp_iter
  )
  dm_fit <- fitCPCA(
    dm_target, center, scale, list(c_contrasts[[2]]), contrasts[2],
    n_eigen, n_medoids, eigdecomp_tol, eigdecomp_iter
  )
  expect_equivalent(dm_fit, fit)
  expect_equivalent(dm_fit, dgC_fit)
})
