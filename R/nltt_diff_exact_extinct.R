#' Checks that event times are correct
#'
#' @description Checks \code{event_times} and \code{event_times2} are of the
#'   appropriate class and have expected characteristics for correct calculation
#'   of NLTT in \code{\link{nltt_diff_exact_extinct}}.
#' @author Pedro Neves and Richel Bilderbeek and Thijs Janzen
#'
#' @param event_times event times of the first phylogeny
#' @param event_times2 event times of the second phylogeny
#' @param time_unit the time unit of the branching times
#' \itemize{
#'  \item{"ago: "}{the branching times are postive,
#'    as these are in time units ago}
#'  \item{"since: "}{the branching times are negative,
#'    as these are in time units since present}
#' }
#'
#' @return Nothing. Throws error with helpful error message if
#' \code{event_times} and \code{event_times2} are not correct.
#' @export
check_input_event_times <- function(event_times, event_times2, time_unit) {

  # This function doesn't handle phylo objects
  if (!is.numeric(event_times) || !is.numeric(event_times2)) {
    stop("event times should be a numeric vector of event times, not a
       phylogeny. Use nltt_diff_exact_norm_brts for phylo objects instead.")
  }

  if (time_unit == "since") {
    if (!all(event_times <= 0.0) || !all(event_times2 <= 0.0)) {
      stop("event times must be negative, ",
           "for example -3 time units since the present")
    }
  } else if (!all(event_times >= 0.0) || !all(event_times2 >= 0.0)) {
    stop("event times must be positive, ",
         "for example 3 time units ago")
  }
}

#' Calculates the exact difference between the nLTT
#' curves of the event times. This includes extinction events.
#' @description Takes branching times such as (for example) as returned by the
#'   DDD package.
#' @author Pedro Neves and Richel Bilderbeek and Thijs Janzen
#' @inheritParams default_params_doc
#' @param event_times event times of the first phylogeny
#' @param species_number the number of species at each event time of the first
#' phylogeny
#' @param event_times2 event times of the second phylogeny
#' @param species_number2 the number of species at each event time of the second
#' phylogeny
#' @param distance_method how the difference between the two nLTTs is summed
#' \itemize{
#'  \item{"abs: "}{the absolute distance between the two nLTTs is summed}
#'  \item{"squ: "}{the squared distance between the two nLTTs is summed}
#' }
#' @param normalize should the output be normalized? Default is TRUE.
#'
#' @examples
#'
#' # Generate data
#' n <- 10
#' b_times_n <- (seq(1, n) / n)
#' lineages_n <- b_times_n
#' b_times2_n <- b_times_n * b_times_n
#' lineages2_n <- b_times2_n
#'
#' # Calculate nLTT
#' out <- nLTT::nltt_diff_exact_extinct(
#'   event_times = b_times_n,
#'   species_number  = lineages_n,
#'   event_times2 = b_times2_n,
#'   species_number2 = lineages2_n,
#'   time_unit = "ago",
#'   distance_method = "abs"
#' )
#'
#' @export
nltt_diff_exact_extinct <- function(
  event_times,
  species_number,
  event_times2,
  species_number2,
  distance_method = "abs",
  time_unit = "since",
  normalize = TRUE) {
  nLTT::check_time_unit(time_unit)
  if (!is.numeric(event_times)) {
    stop(
    "event times should be a numeric vector of event times, not a
       phylogeny. Use nltt_diff_exact_norm_brts for phylo objects instead."
      )

  }
  testthat::expect_true(is.numeric(event_times2))
  nLTT::check_input_event_times(
    event_times = event_times,
    event_times2 = event_times2,
    time_unit = time_unit
  )

  # Each branching time must have a number of lineages to accompany it
  testit::assert(length(event_times) == length(species_number))
  testit::assert(length(event_times2) == length(species_number2))

  if (time_unit == "ago") {
    # Conformize te b_times for the classic calculation
    event_times <- c(-1.0 * rev(sort(event_times)), 0)
    event_times2 <- c(-1.0 * rev(sort(event_times2)), 0)
    species_number  <- c(species_number, utils::tail(species_number, n = 1))
    species_number2 <- c(species_number2, utils::tail(species_number2, n = 1))
  }

  # Each branching time must have a number of species to accompany it
  testit::assert(length(event_times) == length(species_number))
  testit::assert(length(event_times2) == length(species_number2))
  # We calculate with 'since the present' as a time unit
  testit::assert(all(event_times <= 0.0))
  testit::assert(all(event_times2 <= 0.0))
  testit::assert(all(species_number >= 0.0))
  testit::assert(all(species_number2 >= 0.0))

  if (normalize) {
    #normalize event times
    event_times_n <- 1.0 - event_times / min(event_times)

    #normalize lineages
    species_number_n <- species_number / max(species_number)
    # Normalizations must have worked
    testit::assert(all(event_times_n >= 0.0 & event_times_n <= 1.0))
    testit::assert(all(species_number_n >= 0.0 & species_number_n <= 1.0))

    #normalize branching times
    event_times2_n <- 1 - event_times2 / min(event_times2)
    #normalize lineages
    species_number2_n <- species_number2 / max(species_number2)

    nltt_diff_exact_calc_extinct(
      event_times = event_times_n,
      species_number = species_number_n,
      event_times2 = event_times2_n,
      species_number2 = species_number2_n,
      distance_method = distance_method
    )
  } else {
    nltt_diff_exact_calc_extinct(
      event_times = event_times,
      species_number = species_number,
      event_times2 = event_times2,
      species_number2 = species_number2,
      distance_method = distance_method
    )
  }
}

#' Calculates the exact difference between the nLTT
#' curves of the event times. This includes extinction events.
#' @author Thijs Janzen and Richel Bilderbeek and Pedro Neves
#' @param event_times event times of the first phylogeny
#' @param species_number the number of species at each event time of the first
#' phylogeny
#' @param event_times2 event times of the second phylogeny
#' @param species_number2 the number of species at each event time of the second
#' phylogeny
#' @param distance_method (string) absolute, or squared distance?
#'
#' @examples
#'
#' # Generate data
#' n <- 10
#' b_times_n <- (seq(1, n) / n)
#' lineages_n <- b_times_n
#' b_times2_n <- b_times_n * b_times_n
#' lineages2_n <- b_times2_n
#'
#' # Calculate nLTT
#' out <- nLTT::nltt_diff_exact_calc_extinct(
#'   event_times = b_times_n,
#'   species_number  = lineages_n,
#'   event_times2 = b_times2_n,
#'   species_number2 = lineages2_n,
#'   distance_method = "abs"
#' )
#' #'
#' @export
nltt_diff_exact_calc_extinct <- function(
  event_times,
  species_number,
  event_times2,
  species_number2,
  distance_method) {

  # Are equally sized?
  testit::assert(length(event_times) == length(species_number))
  testit::assert(length(event_times2) == length(species_number2))
  # Are non-trivial?
  testit::assert(length(event_times) > 0)
  testit::assert(length(species_number) > 0)
  testit::assert(length(event_times2) > 0)
  testit::assert(length(species_number) > 0)

  # Are sorted?
  testit::assert(all(sort(event_times) == event_times))
  testit::assert(all(sort(event_times2) == event_times2))

  # Make a list of all branching times, and remove duplicates
  # TODO: Use more effecient merge, as both ranges are sorted
  all_event_times <- unique(sort(c(event_times, event_times2)))

  diff <- 0.0
  #iterate through all branching times
  for (k in 2:length(all_event_times)) {
    tim <- all_event_times[k]
    #find the index of the first branching time
    #that is up to the focal branching time
    suppressWarnings(
      index1 <- max(which(event_times < tim))
    )
    suppressWarnings(
      index2 <- max(which(event_times2 < tim))  #same for the other tree
    )
    #find the number of lineages at time "tim" for tree 1
    lins1 <- species_number[max(index1, 1)]
    #find the number of lineages at time "tim" for tree 2
    lins2 <- species_number2[max(index2, 1)]

    #the amount of time that this difference persisted
    dt <- all_event_times[k] - all_event_times[k - 1]
    if (distance_method == "abs") {
      diff <- diff + dt * abs(lins1 - lins2)     #update the difference
    }
    if (distance_method == "squ")  {
      diff <- diff + dt * (lins1 - lins2) * (lins1 - lins2)
    }
  }
  diff
}
