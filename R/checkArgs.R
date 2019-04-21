#' Check Arguments in scPCA Function
#'
#' @description Checks whether or not the all arguments in the scPCA
#' functions are input properly.
#'
#' @param target Target data
#' @param background Background data
#' @param center Indicates whether target and background should be centered.
#' @param scale Indicates whether target and background should be scaled.
#' @param num_eigen The number of eigen vectors to compute.
#' @param contrasts The vector of contrasts parameters
#' @param penalties The vector of penalty terms
#' @param num_medoids The number of medoids to find during spectral clustering.
#'
#' @import assertthat
#'
#' @author Philippe Boileau, \email{philippe_Boileau@@berkeley.edu}
#'
#' @return Whether all argument conditions are satisfied
#'
checkArgs <- function(target, background, center, scale,
                      num_eigen, contrasts, penalties, num_medoids) {

  # assert that the target and background data frames are of the right class
  assert_that(class(target) == "tbl_df" ||
    class(target) == "tbl" ||
    class(target) == "spec_tbl_df" ||
    class(target) == "data.frame" ||
    class(target) == "matrix")
  assert_that(class(background) == "tbl_df" ||
    class(background) == "tbl" ||
    class(background) == "spec_tbl_df" ||
    class(background) == "data.frame" ||
    class(background) == "matrix")

  # assert that target and background have the same number of variables
  assert_that(ncol(target) == ncol(background))

  # check the centering and scaling arguments
  assert_that(is.flag(center))
  assert_that(is.flag(scale))

  # check the number of eigenvectors to compute
  assert_that(is.count(num_eigen))
  assert_that(num_eigen <= ncol(target))

  # check the contrastive parameters
  if (see_if(!missing(contrasts))) {
    assert_that(length(contrasts) > 0)
    assert_that(is.numeric(contrasts))
    assert_that(all(contrasts > 0))
  }

  # check the num_medoids argument
  if (see_if(!missing(num_medoids))) {
    assert_that(is.count(num_medoids))
    assert_that(num_medoids < length(contrasts) * length(penalties))
  }

  # check penalty terms
  if (see_if(!is.null(penalties))) {
    assert_that(length(penalties) > 0)
    assert_that(is.numeric(penalties))
    assert_that(all(penalties >= 0))
  }
}
