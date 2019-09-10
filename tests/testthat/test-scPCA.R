context("Test scPCA")
library(scPCA)

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
    n_centers = 4
  )

  expect_equal(length(scPCA_res), 6)
})

test_that(
  "scPCA outputs a list of length 6, where each element is a sublist of length
n_medoids when n_centers is set to 1", {
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
  }
)
