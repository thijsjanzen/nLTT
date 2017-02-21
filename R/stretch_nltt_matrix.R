#' Stretch matrix 'm' with a timestep resolution of 'dt'.
#'
#' @param m A matrix of 2 columns and at least 2 rows
#' @param dt The resultion, a value e [0.0001, 1]. If 'dt' is set to a very small value, this function will stop
#' @param step_type can be 'lower' or 'upper'
#' @return The stretched matrix
#' @examples
#'   m <- matrix( c(c(0.0, 1.0), c(0.5, 1.0)), ncol = 2, nrow = 2)
#'   expected <- matrix(
#'     c(
#'       c(0.0, 0.5, 1.0),  # Timepoints
#'       c(0.5, 0.5, 1.0)   # Values
#'     ),
#'     ncol = 2, nrow = 3
#'   )
#'   result <- stretch_nltt_matrix(m = m, dt = 0.5, step_type = "lower")
#'   testit::assert(identical(result, expected))
#'
#' @author Richel Bilderbeek
#' @export
stretch_nltt_matrix <- function(
  m,
  dt,
  step_type
) {
  if (!is.matrix(m)) {
    stop("stretch_nltt_matrix: m must be a matrix, ",
      "m is of class '", class(m), "' instead")
  }
  if (ncol(m) != 2) {
    stop("stretch_nltt_matrix: m must have two columns, ",
      "m has ", ncol(m), " columns instead")
  }
  if (step_type != "lower" && step_type != "upper") {
    stop("stretch_nltt_matrix: step_type must be either 'lower' or 'upper', ",
      "step_type supplied was '", step_type, "' instead")
  }

  # Remove rows with same t's, take the first
  rows_to_delete <- NULL
  for (i in seq(1, nrow(m) - 1)) { # -1 because in the body i+1 will be used # nolint
    if (m[i, 1] == m[i + 1, 1]) {
      rows_to_delete <- c(rows_to_delete, i + 1)
    }
  }
  if (!is.null(rows_to_delete)) {
    m <- m[ -rows_to_delete, ]
  }

  # Prepare a new matrix called n
  n_nrow <- 1 + (1.0 / dt)
  testit::assert(all.equal(1.0 / (n_nrow - 1), dt))
  n_ts <- seq(0.0, 1.0, length.out = n_nrow)

  # I am unsure why seq cannot fulfill its promise to create
  # an output of length.out.
  # Stop when this happens
  if (length(n_ts) != n_nrow) {
    stop("dt too small")
  }
  testit::assert(length(n_ts) == n_nrow)

  n_ns <- rep(NA, times = n_nrow)
  testit::assert(length(n_ns) == n_nrow)
  n <- matrix( c(n_ts, n_ns), ncol = 2, nrow = n_nrow)

  testit::assert(nrow(n) == n_nrow)

  # Add endtime at the bottom
  m <- rbind(m, c(1e99, 1.0))

  # Fill in the nLTT values (the second column) in n
  # (n already has the time values in its first column):
  # - go through the values of m until desired timepoint is exceded
  # - copy that nLTT value to n
  #
  #       m                 n
  #
  # |
  # |        +-*     |
  # |      +-*       |    +----*
  # |    +-*      -> |    |
  # |  +-*           +----*
  # |+-*             |
  # |*               *
  # +-----------     +----------------
  # ^    ^    ^      ^    ^    ^
  #

  m_row_index <- 1
  for (n_row_index in 1:n_nrow) {
    # Find the value in m that has a time later than the time value in n
    while (n[n_row_index, 1] >= m[m_row_index + 1, 1]) {
      # Work further through m
      m_row_index <- m_row_index + 1
    }
    # Copy the nLTT value from m
    n[n_row_index, 2] <- m[m_row_index + ifelse(step_type == "lower", 0, 1), 2]  # nolint
  }
  n
}
