#' Check if the time unit is valid
#'
#' Will \link{stop} if not
#' @inheritParams default_params_doc
#' @author Richèl J.C. Bilderbeek
#' @export
check_time_unit <- function(time_unit) {
  time_units <- c("since", "ago")
  if (!time_unit %in% time_units) {
    stop(
      "time_unit must be either 'since' or 'ago'. ",
      "Actual value: ", time_unit
    )
  }

}
