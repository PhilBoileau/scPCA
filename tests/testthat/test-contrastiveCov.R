context("Test contrastiveCov")
library(scPCA)

test_that("performs without issue when matrices are centered or scaled", {
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
                               center = TRUE, scale = TRUE))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
                               center = FALSE, scale = TRUE))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
                               center = TRUE, scale = FALSE))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
                               center = FALSE, scale = FALSE))
})

test_that(
"the num of contrastive cov matrices output is the num of contrastive
parameters",{
  expect_equal(length(contrastiveCov(toy_df[, -31], background_df, 100,
                                     center = TRUE, scale = TRUE)), 1)
  expect_equal(length(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
                              center = TRUE, scale = TRUE)), 3)
})

test_that(
"the contrastive cov matrices have the same num of cols as target and
background", {

  dim_cov <- ncol(contrastiveCov(toy_df[, -31], background_df, 4,
                                 center = TRUE, scale = TRUE)[[1]])
  expect_equal(dim_cov, ncol(toy_df[, -31]))
  expect_equal(dim_cov, ncol(background_df))
})

