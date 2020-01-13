context("Test checking of input arguments")
library(Matrix)
library(tibble)

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

test_that("Only data.frames, tibbles, matrices, and sparse matrices pass", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    as.matrix(toy_df[, 1:30]), as.matrix(background_df),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    as_tibble(toy_df[, 1:30]), as_tibble(background_df),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    as(as.matrix(toy_df[, 1:30]), "dgCMatrix"),
    as(as.matrix(background_df), "dgCMatrix"),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    as(as.matrix(toy_df[, 1:30]), "dgeMatrix"),
    as(as.matrix(background_df), "dgeMatrix"),
    center, scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
})

test_that(paste(
  "Catches target and background data with differing number",
  "of columns"
), {
  expect_error(checkArgs(
    toy_df, background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ),
  "ncol(target) not equal to ncol(background)",
  fixed = TRUE
  )
})

test_that("Center and scale arguments only handle Logical options", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, FALSE,
    TRUE, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, TRUE,
    FALSE, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, FALSE,
    FALSE, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, 1,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, "a",
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    "scale", n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    12342, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
})

test_that("Argument n_eigen is set to an integer between 1 and ncol(target)", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 1, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 30, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 31, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 0, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, "n_eigen", contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, 1.5, contrasts, penalties,
    clust_method, linkage_method
  ))
})

test_that("Contrasts is a non-negative, non-zero vector of length >= 1", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, 1, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, seq(1, 10, by = 0.1), penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, 0, penalties,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, c(seq(1, 10, by = 0.1), -1),
    penalties, clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, "contrasts", penalties,
    clust_method, linkage_method
  ))
})

test_that("Penalties is a non-negative vector of length at least 1", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, seq(1, 10, by = 0.1),
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, 0,
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts,
    c(seq(1, 10, by = 0.1), -1),
    clust_method, linkage_method
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, "penalties",
    clust_method, linkage_method
  ))
})

test_that("ward.D linkage method cannot be used when clust_method is hclust", {
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method, linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "pam", linkage_method
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "hclust", linkage_method = "single"
  ))
  expect_silent(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "hclust", linkage_method = "ward.D2"
  ))
  expect_error(checkArgs(
    toy_df[, 1:30], background_df, center,
    scale, n_eigen, contrasts, penalties,
    clust_method = "hclust", linkage_method = "ward.D"
  ))
})
