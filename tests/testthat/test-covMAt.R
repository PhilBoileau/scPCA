context("Test covMat functionality")
library(scPCA)
library(tidyverse)

test_that("can be applied to dataframes, tibbles and matrices", {
  cov_test <- cov(scale(background_df))

  expect_equal(covMat(background_df), cov_test)
  expect_equal(covMat(as_tibble(background_df)), cov_test)
  expect_equal(covMat(as.matrix(background_df)), cov_test)
})

test_that("can center, scale, center and scale, or neither center and scale", {
  expect_silent(covMat(background_df, center = TRUE, scale = TRUE))
  expect_silent(covMat(background_df, center = TRUE, scale = FALSE))
  expect_silent(covMat(background_df, center = FALSE, scale = TRUE))
  expect_silent(covMat(background_df, center = FALSE, scale = FALSE))

  expect_silent(covMat(as_tibble(background_df), center = TRUE, scale = TRUE))
  expect_silent(covMat(as_tibble(background_df), center = TRUE, scale = FALSE))
  expect_silent(covMat(as_tibble(background_df), center = FALSE, scale = TRUE))
  expect_silent(covMat(as_tibble(background_df), center = FALSE, scale = FALSE))

  expect_silent(covMat(as.matrix(background_df), center = TRUE, scale = TRUE))
  expect_silent(covMat(as.matrix(background_df), center = TRUE, scale = FALSE))
  expect_silent(covMat(as.matrix(background_df), center = FALSE, scale = TRUE))
  expect_silent(covMat(as.matrix(background_df), center = FALSE, scale = FALSE))
})

test_that("variables with zero variance have zeros rows/cols in cov matrix", {
  test_df <- background_df %>%
    mutate(
      zero_var = 1
    )

  expect_equal(sum(covMat(test_df)[, 31]), 0)
  expect_equal(sum(covMat(test_df)[31, ]), 0)
})

