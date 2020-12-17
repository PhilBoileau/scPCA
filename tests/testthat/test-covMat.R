context("Test routines for generating covariance matrices")
library(tibble)
library(dplyr)
library(Matrix)
library(DelayedArray)

test_that(paste0("Routine can be applied to data.frame, tibble, matrix,",
                 " sparse matrices, and DelayedMatrix inputs"), {
  # make sure that all outputs agree
  cov_test <- cov(scale(background_df))
  expect_equal(covMat(background_df), cov_test)
  expect_equal(covMat(as_tibble(background_df)), cov_test)
  expect_equal(covMat(as.matrix(background_df)), cov_test)
  expect_equal(covMat(as(as.matrix(background_df), "dgCMatrix")), cov_test)
  expect_equivalent(as.matrix(covMat(DelayedArray(background_df))), cov_test)

  # check the ScaledMatrix variant
  expect_equal(covMat(background_df, scaled_matrix = TRUE), cov_test)
  expect_equal(covMat(as_tibble(background_df), scaled_matrix = TRUE), cov_test)
  expect_equal(covMat(as.matrix(background_df), scaled_matrix = TRUE), cov_test)
  expect_equal(
    covMat(as(as.matrix(background_df), "dgCMatrix"), scaled_matrix = TRUE),
    cov_test
  )
  expect_equivalent(
    as.matrix(covMat(DelayedArray(background_df), scaled_matrix = TRUE)),
    cov_test
  )
})

test_that("Routine can center, scale, center and scale, or neither", {

  # regular variant
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

  expect_silent(covMat(as(as.matrix(background_df), "dgCMatrix"),
                center = TRUE, scale = TRUE))
  expect_silent(covMat(as(as.matrix(background_df), "dgCMatrix"),
                       center = TRUE, scale = FALSE))
  expect_silent(covMat(as(as.matrix(background_df), "dgCMatrix"),
                       center = FALSE, scale = TRUE))
  expect_silent(covMat(as(as.matrix(background_df), "dgCMatrix"),
                       center = FALSE, scale = FALSE))

  expect_silent(covMat(DelayedArray(background_df),
                       center = TRUE, scale = TRUE))
  expect_silent(covMat(DelayedArray(background_df),
                       center = TRUE, scale = FALSE))
  expect_silent(covMat(DelayedArray(background_df),
                       center = FALSE, scale = TRUE))
  expect_silent(covMat(DelayedArray(background_df),
                       center = FALSE, scale = FALSE))

  # ScaledMatrix variant
  expect_silent(covMat(
      background_df, center = TRUE, scale = TRUE, scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(background_df, center = TRUE, scale = FALSE, scaled_matrix = TRUE)
  )
  expect_silent(
    covMat(background_df, center = FALSE, scale = TRUE, scaled_matrix = TRUE)
  )
  expect_silent(
    covMat(background_df, center = FALSE, scale = FALSE, scaled_matrix = TRUE)
  )

  expect_silent(
    covMat(
      as_tibble(background_df), center = TRUE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as_tibble(background_df), center = FALSE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as_tibble(background_df), center = TRUE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as_tibble(background_df), center = FALSE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )

  expect_silent(
    covMat(
      as.matrix(background_df), center = TRUE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as.matrix(background_df), center = TRUE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as.matrix(background_df), center = FALSE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as.matrix(background_df), center = FALSE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )

  expect_silent(
    covMat(
      as(as.matrix(background_df), "dgCMatrix"), center = TRUE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as(as.matrix(background_df), "dgCMatrix"), center = FALSE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as(as.matrix(background_df), "dgCMatrix"), center = TRUE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      as(as.matrix(background_df), "dgCMatrix"), center = FALSE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )

  expect_silent(
    covMat(
      DelayedArray(background_df), center = TRUE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      DelayedArray(background_df), center = TRUE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      DelayedArray(background_df), center = FALSE, scale = TRUE,
      scaled_matrix = TRUE
    )
  )
  expect_silent(
    covMat(
      DelayedArray(background_df), center = FALSE, scale = FALSE,
      scaled_matrix = TRUE
    )
  )

})

test_that(paste(
  "Variables with zero variance have zeros rows/columns in",
  "covariance matrix"), {
  # add a constant column to the background data
  test_df <- background_df %>%
    mutate(
      zero_var = 1
    )
  expect_equal(sum(covMat(test_df)[, 31]), 0)
  expect_equal(sum(covMat(test_df)[31, ]), 0)
})

test_that("Always returns a matrix object", {
 
  # regular variant
  expect_equal(class(covMat(background_df, center = TRUE, scale = TRUE)),
               c("matrix", "array"))
  expect_equal(class(covMat(as_tibble(background_df),
                            center = TRUE, scale = TRUE)),
               c("matrix", "array"))
  expect_equal(class(covMat(as.matrix(background_df),
                            center = TRUE, scale = TRUE)),
               c("matrix", "array"))
  expect_equal(class(covMat(as(as.matrix(background_df), "dgCMatrix"),
                       center = TRUE, scale = TRUE)),
               c("matrix", "array"))
  expect_equal(class(covMat(DelayedArray(background_df),
                       center = TRUE, scale = TRUE)),
               c("matrix", "array"))

  # ScaledMatrix Variant
  expect_equal(
    class(
      covMat(background_df, center = TRUE, scale = TRUE, scaled_matrix = TRUE)
    ),
    c("matrix", "array")
  )
  expect_equal(
    class(
      covMat(
        as_tibble(background_df), center = TRUE, scale = TRUE,
        scaled_matrix = TRUE
      )
    ),
    c("matrix", "array")
  )
  expect_equal(
    class(
      covMat(
        as.matrix(background_df), center = TRUE, scale = TRUE,
        scaled_matrix = TRUE
      )
    ),
    c("matrix", "array")
  )
  expect_equal(
    class(
      covMat(
        as(as.matrix(background_df), "dgCMatrix"), center = TRUE, scale = TRUE,
        scaled_matrix = TRUE
      )
    ),
    c("matrix", "array")
  )
  expect_equal(
    class(
      covMat(
        DelayedArray(background_df), center = TRUE, scale = TRUE,
        scaled_matrix = TRUE
      )
    ),
    c("matrix", "array")
  )
})
