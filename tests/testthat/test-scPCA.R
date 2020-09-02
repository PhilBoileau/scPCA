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


test_that("n_centers need not be passed when contrasts is a single number", {
  expect_silent(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = 1,
    penalties = 0,
  ))
  expect_silent(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = 1,
    penalties = 0.1
  ))
})

test_that("Users can provide cluster labels, bypassing clustering step", {
  expect_silent(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(0.5, 1, 2, 3, 5),
    penalties = 0,
    clusters = toy_df[, 31]
  ))
  expect_silent(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 5)),
    penalties = seq(0.1, 1, length.out = 3),
    alg = "rand_var_proj",
    clusters = toy_df[, 31]
  ))
})

test_that("`rotation` and `x` are always a matrices", {
  
  # cPCA, only one contrast
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = 1,
    penalties = 0,
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = 1,
    penalties = 0,
  )$x)[1],
  "matrix")
  
  # scPCA, only one contrast and L1 penalty
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = 1,
    penalties = 1,
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = 1,
    penalties = 1,
  )$x)[1],
  "matrix")
  
  # cPCA, abid et al version
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 10)),
    penalties = 0,
    n_centers = 1,
    n_medoids = 4
  )$rotation)[1],
  "list")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 10)),
    penalties = 0,
    n_centers = 1,
    n_medoids = 4
  )$rotation[[1]])[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 10)),
    penalties = 0,
    n_centers = 1,
    n_medoids = 4
  )$x)[1],
  "list")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = exp(seq(log(0.1), log(100), length.out = 10)),
    penalties = 0,
    n_centers = 1,
    n_medoids = 4
  )$x[[1]])[1],
  "matrix")
  
  # cPCA, Boileau et al
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = 0,
    n_centers = 4
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = 0,
    n_centers = 4
  )$x)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = 0,
    clusters = toy_df[, 31]
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = 0,
    clusters = toy_df[, 31]
  )$x)[1],
  "matrix")
  
  # scPCA, Boileau et al
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = 0.1,
    n_centers = 4
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = 0.1,
    n_centers = 4
  )$x)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = c(0.1, 1),
    clusters = toy_df[, 31]
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = c(0.1, 1),
    clusters = toy_df[, 31]
  )$x)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = c(0.1, 1),
    n_centers = 4
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = c(0.1, 1),
    n_centers = 4
  )$x)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = c(0.1, 1),
    clusters = toy_df[, 31]
  )$rotation)[1],
  "matrix")
  expect_equal(class(scPCA(
    target = toy_df[, 1:30],
    background = background_df,
    contrasts = c(1, 2, 3),
    penalties = c(0.1, 1),
    clusters = toy_df[, 31]
  )$x)[1],
  "matrix")
})
