context("Test fitting of contrastive PCA")
library(BiocParallel)

# initialization of inputs for testing function
target <- toy_df[, -31]
background <- background_df
center <- TRUE
scale <- TRUE
contrasts <- c(0, 1, 10, 100)
c_contrasts <- contrastiveCov(target, background, contrasts, TRUE, TRUE)
n_eigen <- 2
n_medoids <- 2
fit <- fitCPCA(
  target, center, scale, c_contrasts, contrasts,
  n_eigen, n_medoids
)

test_that(paste(
  "The length of the ouput list is 4 (rotations, projections,",
  "contrasts and penalties)"
), {
  expect_equal(length(fit), 4)
})

test_that("The number of medoids output equals the initial specified option", {
  expect_equal(length(fit$contrast), 2)
})

test_that("Dimensions of the rotation matrices are ncol(target) by n_eigen", {
  expect_equal(nrow(fit$rotation[[1]]), ncol(target))
  expect_equal(ncol(fit$rotation[[1]]), n_eigen)
  expect_equal(nrow(fit$rotation[[2]]), ncol(target))
  expect_equal(ncol(fit$rotation[[2]]), n_eigen)
})

test_that("The projection matrices are nrow(target) by n_eigen", {
  expect_equal(nrow(fit$x[[1]]), nrow(target))
  expect_equal(ncol(fit$x[[1]]), n_eigen)
  expect_equal(nrow(fit$x[[2]]), nrow(target))
  expect_equal(ncol(fit$x[[2]]), n_eigen)
})

test_that("The penalties output vector is a zero vector", {
  expect_equal(fit$penalty, rep(0, n_medoids))
})

test_that("The parallelized analog matches the sequential variant exactly", {
  # fit sequential variant
  fit <- fitCPCA(
    target, center, scale, c_contrasts, contrasts,
    n_eigen, n_medoids
  )
  # fit parallelized variant
  register(SerialParam())
  bpfit <- bpFitCPCA(
    target, center, scale, c_contrasts, contrasts,
    n_eigen, n_medoids
  )
  expect_equal(fit, bpfit)
})