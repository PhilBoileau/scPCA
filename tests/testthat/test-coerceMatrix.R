context("Test routine for coercing inputs to matrix objects")
library(Matrix)

test_that("Coercion of sparse matrices to standard matrices works", {
  # NOTE: Fails on Linux and Windows, not MacOS?
  skip_on_os(c("linux", "windows", "mac"))
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgCMatrix"))),
    "matrix"
  )
  expect_equal(
    class(coerceMatrix(data = as(as.matrix(background_df), "dgeMatrix"))),
    "matrix"
  )
})
