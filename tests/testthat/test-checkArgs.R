context("Test checkArgs functionality")
library(scPCA)
library(tidyverse)
library(Matrix)

# set up dummy input for testing
target_df <- toy_df
background_df <- background_df
center <- TRUE
scale <- FALSE
n_eigen <- 2
contrasts <- c(1, 2)
penalties <- c(1, 2)

test_that("only dataframes, tibbles, matrices, and sparse matrices pass", {
  expect_equal(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, n_eigen, contrasts, penalties), TRUE)
  expect_equal(checkArgs(as.matrix(toy_df[, 1:30]), as.matrix(background_df),
                         center, scale, n_eigen, contrasts, penalties), TRUE)
  expect_equal(checkArgs(as_tibble(toy_df[, 1:30]), as_tibble(background_df),
                         center, scale, n_eigen, contrasts, penalties), TRUE)
  expect_equal(checkArgs(as(as.matrix(toy_df[, 1:30]), "dgCMatrix"),
                         as(as.matrix(background_df), "dgCMatrix"),
                         center, scale, n_eigen, contrasts, penalties), TRUE)
  expect_equal(checkArgs(as(as.matrix(toy_df[, 1:30]), "dgeMatrix"),
                         as(as.matrix(background_df), "dgeMatrix"),
                         center, scale, n_eigen, contrasts, penalties), TRUE)
})

test_that("catches target and background data with differing num of cols", {
  expect_error(checkArgs(toy_df, background_df, center,
                         scale, n_eigen, contrasts, penalties),
               "ncol(target) not equal to ncol(background)", fixed = TRUE)
})

test_that("center and scale only handle bools", {
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, n_eigen, contrasts, penalties))
  expect_silent(checkArgs(toy_df[, 1:30], background_df, FALSE,
                          TRUE, n_eigen, contrasts, penalties))
  expect_silent(checkArgs(toy_df[, 1:30], background_df, TRUE,
                          FALSE, n_eigen, contrasts, penalties))
  expect_silent(checkArgs(toy_df[, 1:30], background_df, FALSE,
                          FALSE, n_eigen, contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, 1,
                          scale, n_eigen, contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, "a",
                          scale, n_eigen, contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         "scale", n_eigen, contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         12342, n_eigen, contrasts, penalties))
})

test_that("n_eigen is an integer between 1 and ncol(target)", {
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                          scale, 1, contrasts, penalties))
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                          scale, 30, contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, 31, contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, 0, contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, "n_eigen", contrasts, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, 1.5, contrasts, penalties))

})

test_that("contrasts is a non-negative, non-zero vector of length at least 1", {
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                          scale, n_eigen, 1, penalties))
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                          scale, n_eigen, seq(1, 10, by = 0.1), penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, n_eigen, 0, penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                          scale, n_eigen, c(seq(1, 10, by = 0.1), -1),
                         penalties))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, n_eigen, "contrasts", penalties))
})

test_that("pernalties is a non-negative vector of length at least 1", {
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                          scale, n_eigen, contrasts, penalties))
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                          scale, n_eigen, contrasts, seq(1, 10, by = 0.1)))
  expect_silent(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, n_eigen, contrasts, 0))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, n_eigen, contrasts,
                         c(seq(1, 10, by = 0.1), -1)))
  expect_error(checkArgs(toy_df[, 1:30], background_df, center,
                         scale, n_eigen, contrasts, "penalties"))
})
