context("Test checking of input arguments")
library(tibble)
library(Matrix)
library(DelayedArray)
library(sparseMatrixStats)
library(DelayedMatrixStats)

# set up dummy input for testing
target_df <- toy_df
background_df <- background_df
center <- TRUE
scale <- FALSE
n_eigen <- 2
contrasts <- c(1, 2)
penalties <- c(1, 2)
clust_method <- "kmeans"
linkage_method <- "complete"
clusters <- NULL
eigdecomp_tol <- 1e-10
eigdecomp_iter <- 1000
n_centers <- 2
scaled_matrix <- FALSE

test_that(paste0("Only data.frames, tibbles, matrices, sparse matrices and",
                 "DelayedMatrices pass"), {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    as.matrix(toy_df[, 1:30]), as.matrix(background_df),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix = TRUE
  ))
  expect_silent(checkArgs(
    as_tibble(toy_df[, 1:30]), as_tibble(background_df),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    as(as.matrix(toy_df[, 1:30]), "dgCMatrix"),
    as(as.matrix(background_df), "dgCMatrix"),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix = TRUE
  ))
  expect_silent(checkArgs(
    DelayedArray(toy_df[, 1:30]),
    DelayedArray(background_df),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
})

test_that(paste(
  "Catches target and background data with differing number",
  "of columns"
), {
  expect_error(checkArgs(
    toy_df, background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ),
  "ncol(target) not equal to ncol(background)",
  fixed = TRUE
  )
})

test_that("Center and scale arguments only handle Logical options", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method,clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, FALSE,
    TRUE, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, TRUE,
    FALSE, n_eigen, contrasts, penalties,
    clust_method, linkage_metho, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix = TRUE
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, FALSE,
    FALSE, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix = TRUE
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, 1,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, "a",
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    "scale", n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix = TRUE
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    12342, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scale_matrix = "T"
  ))
})

test_that("Argument n_eigen is set to an integer between 1 and ncol(target)", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 1, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 30, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 31, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 0, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, "n_eigen", contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 1.5, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix = TRUE
  ))
})

test_that("Contrasts is a non-negative, non-zero vector of length >= 1", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, 1, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, seq(1, 10, by = 0.1), penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, 0, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, c(seq(1, 10, by = 0.1), -1),
    penalties, clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, "contrasts", penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
})

test_that("Penalties is a non-negative vector of length at least 1", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, seq(1, 10, by = 0.1),
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, 0,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts,
    c(seq(1, 10, by = 0.1), -1),
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, "penalties",
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
})

test_that("ward.D linkage method cannot be used when clust_method is hclust", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "pam", linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "hclust", linkage_method = "single", clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "hclust", linkage_method = "ward.D2", clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "hclust", linkage_method = "ward.D", clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
})

test_that("Catches cluster assignments that don't match requirements", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method , clusters,
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method , clusters = toy_df[, 31],
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method , clusters = toy_df[1:25, 31],
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method , clusters = as.character(toy_df[, 31]),
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method , clusters = as.factor(toy_df[, 31]),
    eigdecomp_tol, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
})

test_that("Checks that RSpectra options are reasonable.", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol = 1e-3, eigdecomp_iter = 1000, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol = -1e-10, eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter = -1000, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol = -1e-10, eigdecomp_iter = -100, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol = "hello", eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol = c(1e-10, 1e-8), eigdecomp_iter, n_centers,
    scaled_matrix
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters,
    eigdecomp_tol, eigdecomp_iter = c(10, 100, 1000), n_centers,
    scaled_matrix
  ))
})

test_that("Check args picks up on missing n_centers", {
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters = NULL,
    eigdecomp_tol, eigdecomp_iter, n_centers = NULL,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters = toy_df[, 31],
    eigdecomp_tol, eigdecomp_iter, n_centers = NULL,
    scaled_matrix
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method, clusters = NULL,
    eigdecomp_tol, eigdecomp_iter, n_centers = 2,
    scaled_matrix
  ))
})
