context("Test fitting of sparse contrastive PCA")
library(BiocParallel)

test_that("scPCA outputs a list of length 6", {
  cPCA_res <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = 0,
    n_centers = 4
  )
  expect_equal(length(cPCA_res), 6)

  scPCA_res <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = seq(0.1, 1, length.out = 3),
    n_centers = 4, alg = "iterative"
  )
  expect_equal(length(scPCA_res), 6)
  
  
  set.seed(978)
  cv_cPCA_res <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = 0,
    n_centers = 4,
    cv = 2
  )
  expect_equal(length(cv_cPCA_res), 6)
  
  cv_scPCA_res <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = seq(0.1, 1, length.out = 3),
    n_centers = 4,
    alg = "iterative",
    cv = 2
  )
  expect_equal(length(cv_scPCA_res), 6)
})

test_that(paste(
  "scPCA outputs a list of length 6, where each element is a",
  "sublist of length n_medoids when n_centers is set to 1"
), {
  cPCA_res <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 10)),
    penalties = 0,
    n_centers = 1,
    n_medoids = 5
  )

  expect_equal(length(cPCA_res$contrast), 5)
  expect_equal(length(cPCA_res$penalty), 5)
  expect_equal(length(cPCA_res$rotation), 5)
  expect_equal(length(cPCA_res$x), 5)
})

test_that("The parallelized analog matches the sequential variant exactly", {
  set.seed(653)
  cPCA_res <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = 0,
    n_centers = 4, parallel = FALSE
  )
  register(SerialParam())
  set.seed(653)
  cPCA_bpres <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = 0,
    n_centers = 4, parallel = TRUE
  )
  expect_equal(cPCA_res, cPCA_bpres)
  
  set.seed(123)
  scPCA_res <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = seq(0.1, 1, length.out = 3),
    n_centers = 4, parallel = FALSE
  )
  set.seed(123)
  scPCA_bpres <- scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = seq(0.1, 1, length.out = 3),
    n_centers = 4, parallel = TRUE
  )
  expect_equal(scPCA_res, scPCA_bpres)
})
