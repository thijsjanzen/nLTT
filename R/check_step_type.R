#' Check if the step type is valid
#'
#' Will \link{stop} if not
#' @inheritParams default_params_doc
#' @export
check_step_type <- function(step_type) {
  if (length(step_type) != 1 || (step_type != "lower" && step_type != "upper")) {
    stop(
      "stretch_nltt_matrix: step_type must be either 'lower' or 'upper', ",
      "step_type supplied was '", step_type, "' instead"
      )
  }
}
