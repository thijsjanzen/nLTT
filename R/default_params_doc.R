#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param dt The timestep resolution,
#'   a value bigger than zero and less or equal to one.
#'   1/dt is the number of points that will be evaluated
#' @param phylogenies a collection of one or more phylogenies,
#'   where the phylogenies are of type \link[ape]{phylo}.
#'   This collection can both be a list of \link[ape]{phylo} or
#'   a \link[ape]{multiphylo}.
#' @param step_type when between two points, where the second point
#'   has both a higher x and y coordinat, which y coordinat to follow.
#'  'step_type' can be:
#'  \itemize{
#'     \item lower maintain the y-coordinat of the leftmost point
#'     \item upper already use the y-coordinat of the rightmost point
#'  }
#' @param time_unit the time unit of the branching times
#' \itemize{
#'  \item{"ago: "}{the branching times are postive,
#'    as these are in time units ago}
#'  \item{"since: "}{the branching times are negative,
#'    as these are in time units since present}
#' }
#' @author Richèl J.C. Bilderbeek
#' @note This is an internal function, so it should be marked with
#'   \code{@noRd}. This is not done, as this will disallow all
#'   functions to find the documentation parameters
default_params_doc <- function(
  dt,
  phylogenies,
  time_unit
) {
  # Nothing
}
