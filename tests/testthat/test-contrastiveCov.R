context("Test contrastive covariance matrices")
library(BiocParallel)

test_that("performs without issue when matrices are centered or scaled", {
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = TRUE, scale = TRUE
  ))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = FALSE, scale = TRUE
  ))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = TRUE, scale = FALSE
  ))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = FALSE, scale = FALSE
  ))
})

test_that(paste("Number of contrastive covariance matrices matches number of",
                "contrastive parameters"), {
    expect_equal(length(contrastiveCov(toy_df[, -31], background_df, 100,
      center = TRUE, scale = TRUE
    )), 1)
    expect_equal(length(contrastiveCov(toy_df[, -31],
                                       background_df, c(0, 1, 2),
                                       center = TRUE, scale = TRUE)), 3)
})

test_that(paste("Number of columns of contrastive covariance matrix matches",
                "the number of columns in the target and background data"), {
    dim_cov <- ncol(contrastiveCov(toy_df[, -31], background_df, 4,
                                   center = TRUE, scale = TRUE)[[1]])
    expect_equal(dim_cov, ncol(toy_df[, -31]))
    expect_equal(dim_cov, ncol(background_df))
})

test_that(paste("Parallelized contrastive covariance routine matches",
                "sequential analog"), {
    register(SerialParam())
    expect_equal(contrastiveCov(toy_df[, -31],
                                background_df, c(0, 1, 2),
                                center = TRUE, scale = TRUE),
                 bpContrastiveCov(toy_df[, -31],
                                  background_df, c(0, 1, 2),
                                  center = TRUE, scale = TRUE))
})
