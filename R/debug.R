#' Run gsPEN on a small sample of the provided data set (Only 100 samples)
#' @param ... Additional arguments to pass to gsPEN_R
#' @return The output of gsPEN_R
#' @importFrom utils data
#' @export
test_gsPEN <- function(...) {
  summaryZ <- Nvec <- plinkLD <- NULL
  data("summaryZ", envir = environment())
  data("Nvec", envir = environment())
  data("plinkLD", envir = environment())

  subset <- sample(nrow(summaryZ), 100)
  subset_summary_z <- summaryZ[subset, ]

  output <- gsPEN_R(
    summary_z = subset_summary_z, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    sub_tuning = 1,
    ...
  )

  return(output)
}

#' Run gsfPEN on a small sample of the provided data set (Only 100 samples)
#' @param ... Additional arguments to pass to gsfPEN_R
#' @return The output of gsfPEN_R
#' @importFrom utils data
#' @export
test_gsfPEN <- function(...) {
  summaryZ <- Nvec <- plinkLD <- funcIndex <- NULL
  data("summaryZ", envir = environment())
  data("Nvec", envir = environment())
  data("plinkLD", envir = environment())
  data("funcIndex", envir = environment())

  subset <- sample(nrow(summaryZ), 100)
  subset_summary_z <- summaryZ[subset, ]
  subset_func_index <- funcIndex[subset, ]

  output <- gsfPEN_R(
    summary_z = subset_summary_z, # nolint: object_usage_linter.
    n_vec = Nvec, # nolint: object_usage_linter.
    plinkLD = plinkLD, # nolint: object_usage_linter.
    func_index = subset_func_index, # nolint: object_usage_linter.
    sub_tuning = 1,
    ...
  )

  return(output)
}
