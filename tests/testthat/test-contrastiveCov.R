context("Test contrastive covariance matrices")
library(BiocParallel)
library(Matrix)
library(DelayedArray)
library(sparseMatrixStats)
library(DelayedMatrixStats)

dgC_back <- as(as.matrix(background_df), "dgCMatrix")
dgC_toy <- as(as.matrix(toy_df[, -31]), "dgCMatrix")
dm_back <- DelayedArray(background_df)
dm_toy <- DelayedArray(toy_df[, -31])

test_that(paste("performs without issue when matrices are centered or scaled,",
                "and when using ScaledMatrix objects"), {

  # normal dataframes
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = TRUE, scale = TRUE, scaled_matrix = TRUE
  ))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = FALSE, scale = TRUE
  ))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = TRUE, scale = FALSE, scaled_matrix = TRUE
  ))
  expect_silent(contrastiveCov(toy_df[, -31], background_df, c(0, 1, 2),
    center = FALSE, scale = FALSE
  ))

  # sparse matrices
  expect_silent(contrastiveCov(dgC_toy, dgC_back, c(0, 1, 2),
                               center = TRUE, scale = TRUE
  ))
  expect_silent(contrastiveCov(dgC_toy, dgC_back, c(0, 1, 2),
                               center = FALSE, scale = TRUE,
                               scaled_matrix = TRUE
  ))
  expect_silent(contrastiveCov(dgC_toy, dgC_back, c(0, 1, 2),
                               center = TRUE, scale = FALSE
  ))
  expect_silent(contrastiveCov(dgC_toy, dgC_back, c(0, 1, 2),
                               center = FALSE, scale = FALSE
  ))

  # DelayedMatrices
  expect_silent(contrastiveCov(dm_toy, dm_back, c(0, 1, 2),
                               center = TRUE, scale = TRUE,
                               scaled_matrix = TRUE
  ))
  expect_silent(contrastiveCov(dm_toy, dm_back, c(0, 1, 2),
                               center = FALSE, scale = TRUE
  ))
  expect_silent(contrastiveCov(dm_toy, dm_back, c(0, 1, 2),
                               center = TRUE, scale = FALSE,
                               scaled_matrix = TRUE
  ))
  expect_silent(contrastiveCov(dm_toy, dm_back, c(0, 1, 2),
                               center = FALSE, scale = FALSE
  ))
})

test_that(paste(
  "Number of contrastive covariance matrices matches number of",
  "contrastive parameters"
), {
  # normal df
  expect_equal(length(contrastiveCov(toy_df[, -31], background_df, 100,
    center = TRUE, scale = TRUE
  )), 1)
  expect_equal(length(contrastiveCov(toy_df[, -31],
    background_df, c(0, 1, 2),
    center = TRUE, scale = TRUE, scaled_matrix = TRUE
  )), 3)

  # sparse matrices
  expect_equal(length(contrastiveCov(dgC_toy, dgC_back, 100,
                                     center = TRUE, scale = TRUE,
                                     scaled_matrix = TRUE
  )), 1)
  expect_equal(length(contrastiveCov(dgC_toy, dgC_back, c(0, 1, 2),
                                     center = TRUE, scale = TRUE
  )), 3)

  # DelayedMatrices
  expect_equal(length(contrastiveCov(dm_toy, dm_back, 100,
                                     center = TRUE, scale = TRUE
  )), 1)
  expect_equal(length(contrastiveCov(dm_toy, dm_back, c(0, 1, 2),
                                     center = TRUE, scale = TRUE,
                                     scaled_matrix = TRUE
  )), 3)
})

test_that(paste(
  "Number of columns of contrastive covariance matrix matches",
  "the number of columns in the target and background data"
), {

  # data frames
  dim_cov <- ncol(contrastiveCov(toy_df[, -31], background_df, 4,
    center = TRUE, scale = TRUE
  )[[1]])
  expect_equal(dim_cov, ncol(toy_df[, -31]))
  expect_equal(dim_cov, ncol(background_df))

  # sparse matrices
  dim_cov <- ncol(contrastiveCov(dgC_toy, dgC_back, 4,
                                 center = TRUE, scale = TRUE,
                                 scaled_matrix = TRUE
  )[[1]])
  expect_equal(dim_cov, ncol(toy_df[, -31]))
  expect_equal(dim_cov, ncol(background_df))

  # DelayedMatrices
  dim_cov <- ncol(contrastiveCov(dm_toy, dm_back, 4,
                                 center = TRUE, scale = TRUE
  )[[1]])
  expect_equal(dim_cov, ncol(toy_df[, -31]))
  expect_equal(dim_cov, ncol(background_df))
})

test_that(paste(
  "Parallelized contrastive covariance routine matches",
  "sequential analog"
), {
  register(SerialParam())

  # regular data frame
  expect_equal(
    contrastiveCov(toy_df[, -31],
      background_df, c(0, 1, 2),
      center = TRUE, scale = TRUE
    ),
    bpContrastiveCov(toy_df[, -31],
      background_df, c(0, 1, 2),
      center = TRUE, scale = TRUE
    )
  )

  # sparse matrix
  expect_equal(
    contrastiveCov(dgC_toy, dgC_back, c(0, 1, 2),
                   center = TRUE, scale = TRUE, scaled_matrix = TRUE
    ),
    bpContrastiveCov(dgC_toy, dgC_back, c(0, 1, 2),
                     center = TRUE, scale = TRUE, scaled_matrix = TRUE
    )
  )

  # DelayedMatrix
  expect_equal(
    contrastiveCov(dm_toy, dm_back, c(0, 1, 2),
                   center = TRUE, scale = TRUE
    ),
    bpContrastiveCov(dm_toy, dm_back, c(0, 1, 2),
                     center = TRUE, scale = TRUE
    )
  )
})
