.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "scPCA v",
    utils::packageDescription("scPCA")$Version, ": ",
    utils::packageDescription("scPCA")$Title
  ))
}
