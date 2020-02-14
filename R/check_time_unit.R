#' Check if the time unit is valid
#'
#' Will \link{stop} if not
#' @author Richèl J.C. Bilderbeek
#' @examples
#' library(testthat)
#'
#' expect_silent(check_time_unit("since"))
#' expect_silent(check_time_unit("ago"))
#' expect_error(check_time_unit("nonsense"))
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

