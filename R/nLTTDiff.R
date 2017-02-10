#' Calculates the exact, difference between the lineage through time
#' curves of tree1 & tree2 (normalized in time and for the number of lineages)
#' @author Thijs Janzen
#' @param tree1 (phylo) First phylogenetic tree
#' @param tree2 (phylo) Second phylogenetic tree
#' @param distance_method (string) absolute, or squared distance?
#' @param ignore_stem (logical) Should the phylogeny its stem be ignored?
#' @return (scalar) normalized Lineage-Through-Time difference between tree1 & tree2
#' @export
nltt_diff_exact <- function(
  tree1,
  tree2,
  distance_method = "abs",
  ignore_stem = TRUE
) {
  return (nLTT::nltt_diff( tree1, tree2, distance_method, ignore_stem))
}

#' Calculates the exact difference between the nLTT
#' curves of the branching times
#' @author Thijs Janzen and Richel Bilderbeek
#' @param b_times branching times of the first phylogeny,
#' @param lineages the number of lineages, usually one to the number of lineages
#' @param b_times2 branching times of the first phylogeny
#' @param lineages2 the number of lineages, usually one to the number of lineages
#' @param distance_method how the difference between the two nLTTs is summed
#' \itemize{
#'  \item{"abs: "}{the absolute distance between the two nLTTs is summed}
#'  \item{"squ: "}{the squared distance between the two nLTTs is summed}
#' }
#' @param time_unit the time unit of the branching times
#' \itemize{
#'  \item{"ago: "}{the branching times are postive, as these are in time units ago}
#'  \item{"since: "}{the branching times are negative, as these are in time units since present}
#' }
#' @export
nltt_diff_exact_brts <- function(
  b_times,
  lineages,
  b_times2,
  lineages2,
  distance_method = "abs",
  time_unit = "since") {
  if (time_unit != "since" && time_unit != "ago") {
    stop("time_unit must be either 'since' or 'ago'")
  }
  # Each branching time must have a number of lineages to accompany it
  testit::assert(length(b_times ) == length(lineages ))
  testit::assert(length(b_times2) == length(lineages2))

  if (time_unit == "since"){
    if (!all(b_times <= 0.0) || !all(b_times2 <= 0.0)) {
      stop("branching times must be negative, ",
        "for example -3 time units since the present")
    }

  }
  if (time_unit == "ago"){
    if (!all(b_times >= 0.0) || !all(b_times2 >= 0.0)) {
      stop("branching times must be positive, ",
        "for example 3 time units ago")
    }
   # Conformize te b_times for the classic calculation
    b_times <- c(-1.0 * rev(sort(b_times)), 0)
    b_times2 <- c(-1.0 * rev(sort(b_times2)), 0)
    lineages  <- c(lineages, utils::tail(lineages, n = 1))
    lineages2 <- c(lineages2, utils::tail(lineages2, n = 1))
  }

  # Each branching time must have a number of lineages to accompany it
  testit::assert(length(b_times ) == length(lineages ))
  testit::assert(length(b_times2) == length(lineages2))
  # We calculate with 'since the present' as a time unit
  testit::assert(all(b_times <= 0.0))
  testit::assert(all(b_times2 <= 0.0))
  testit::assert(all(lineages >= 0.0))
  testit::assert(all(lineages2 >= 0.0))


  b_times_N <- 1.0 - b_times / min(b_times) #normalize branching times
  lineages_N <- lineages / max(lineages)  #normalize lineages
  # Normalizations must have worked
  testit::assert(all(b_times_N >= 0.0 & b_times_N <= 1.0))
  testit::assert(all(lineages_N >= 0.0 & lineages_N <= 1.0))

  b_times2_N <- 1 - b_times2 / min(b_times2) #normalize branching times
  lineages2_N <- lineages2 / max(lineages2)  #normalize lineages
  # Normalizations must have worked

  return (nltt_diff_exact_norm_brts(
    b_times_N = b_times_N,
    lineages_N = lineages_N,
    b_times2_N = b_times2_N,
    lineages2_N = lineages2_N,
    distance_method = distance_method)
  )
}



#' Calculates the exact difference between the nLTT
#' curves of the branching times
#' @author Thijs Janzen and Richel Bilderbeek
#' @param b_times_N branching times of the first phylogeny
#' @param lineages_N the number of lineages, usually one to the number of lineages
#' @param b_times2_N branching times of the first phylogeny
#' @param lineages2_N the number of lineages, usually one to the number of lineages
#' @param distance_method (string) absolute, or squared distance?
#' @export
nltt_diff_exact_norm_brts <- function(
  b_times_N,
  lineages_N,
  b_times2_N,
  lineages2_N,
  distance_method) {

  # Are equally sized?
  testit::assert(length(b_times_N) == length(lineages_N))
  testit::assert(length(b_times2_N) == length(lineages2_N))
  # Are non-trivial?
  testit::assert(length(b_times_N) > 0)
  testit::assert(length(lineages_N) > 0)
  testit::assert(length(b_times2_N) > 0)
  testit::assert(length(lineages2_N) > 0)
  # Start with zeroes? Not necessarily
  testit::assert(b_times_N[1] >= 0.0)
  testit::assert(lineages_N[1] >= 0.0)
  testit::assert(b_times2_N[1] >= 0.0)
  testit::assert(lineages2_N[1] >= 0.0)
  # Ends with ones
  testit::assert(tail(b_times_N, n = 1) == 1.0)
  testit::assert(tail(lineages_N, n = 1) == 1.0)
  testit::assert(tail(b_times2_N, n = 1) == 1.0)
  testit::assert(tail(lineages2_N, n = 1) == 1.0)
  # Are in range?
  testit::assert(all(b_times_N >= 0.0 & b_times_N <= 1.0))
  testit::assert(all(lineages_N >= 0.0 & lineages_N <= 1.0))
  testit::assert(all(b_times2_N >= 0.0 & b_times2_N <= 1.0))
  testit::assert(all(lineages2_N >= 0.0 & lineages2_N <= 1.0))
  # Are sorted?
  testit::assert(all(sort(b_times_N) == b_times_N))
  testit::assert(all(sort(lineages_N) == lineages_N))
  testit::assert(all(sort(b_times2_N) == b_times2_N))
  testit::assert(all(sort(lineages2_N) == lineages2_N))

  # Make a list of all branching times, and remove duplicates
  # TODO: Use more effecient merge, as both ranges are sorted
  all_b_times <- unique(sort(c(b_times_N, b_times2_N)))

  diff <- 0.0
  #iterate through all branching times
  for (k in 2:length(all_b_times)) {
      tim <- all_b_times[k]
      #find the index of the first branching time
      #that is up to the focal branching time
      index1 <- max(which(b_times_N < tim))
      index2 <- max(which(b_times2_N < tim))  #same for the other tree

      #find the number of lineages at time "tim" for tree 1
      lins1 <- lineages_N[max(index1, 1)]
      #find the number of lineages at time "tim" for tree 2
      lins2 <- lineages2_N[max(index2, 1)]

      #the amount of time that this difference persisted
      dt <- all_b_times[k] - all_b_times[k - 1]
      if (distance_method == "abs") {
        diff <- diff + dt * abs( lins1 - lins2)     #update the difference
      }
      if (distance_method == "squ")  {
        diff <- diff + dt * ( lins1 - lins2) * ( lins1 - lins2)
      }
  }
  return ( diff)
}

#' Calculates the exact difference between the lineage through time curves of tree1 & tree2 (normalized in time and for the number of lineages)
#' @author Thijs Janzen
#' @param tree1 (phylo) First phylogenetic tree
#' @param tree2 (phylo) Second phylogenetic tree
#' @param distance_method (string) absolute, or squared distance?
#' @param ignore_stem logical    Should the phylogeny its stem be ignored?
#' @return (scalar) normalized Lineage-Through-Time difference between tree1 & tree2
#' @export
nltt_diff <- function(
 tree1, tree2, distance_method = "abs", ignore_stem = TRUE) {

  if (!ape::is.binary(tree1) || !ape::is.binary(tree2)) {
    stop("phylogenies must both be binary")
  }

  if (any(ape::branching.times(tree1) < 0.0)) {
    stop("tree1 cannot have negative branching times")
  }
  if (any(ape::branching.times(tree2) < 0.0)) {
    stop("tree2 cannot have negative branching times")
  }

  #branching times of tree1, including the present time (0)
  b_times <- c(-1.0 * rev(sort(ape::branching.times(tree1))), 0)
  if (!ignore_stem) {
    stem_length1 <- ifelse(is.null(tree1$root.edge), 0.0, tree1$root.edge)
    b_times <- c(b_times[1] - stem_length1, b_times)
    testit::assert(all(b_times <= 0.0))
  }

  # Same for other tree
  b_times2 <- c(-1 * rev(sort(ape::branching.times(tree2))), 0)
  options(scipen = 6)
  if (!ignore_stem) {
    stem_length2 <- ifelse(is.null(tree2$root.edge), 0.0, tree2$root.edge)
    b_times2 <- c(b_times2[1] - stem_length2, b_times2)
    testit::assert(all(b_times2 <= 0.0))
  }
  # the number of lineages per branching time
  first_n_lineages1 <- ifelse(ignore_stem, 2, 1)
  n_taxa1 <- length(tree1$tip.label)
  lineages <- c(first_n_lineages1:n_taxa1, n_taxa1)
  # Each branching time must have a number of lineages to accompany it
  testit::assert(length(b_times) == length(lineages))

  # the number of lineages per branching time
  first_n_lineages2 <- ifelse(ignore_stem, 2, 1)
  n_taxa2 <- length(tree2$tip.label)
  lineages2 <- c(first_n_lineages2:n_taxa2, n_taxa2)
  # Each branching time must have a number of lineages to accompany it
  testit::assert(length(b_times2) == length(lineages2))

  return (
    nltt_diff_exact_brts(
      b_times = b_times,
      lineages = lineages,
      b_times2 = b_times2,
      lineages2 = lineages2,
      distance_method = distance_method
    )
  )
}

################################################################################
#
# @brief Wrapper to calculate the nLTT statistic
#
# @date Last modified: 2015-21-04
# @author Thijs Janzen
# @since 2015-21-04, version 1.1
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @param    distance_method        string     Method to calculate the difference, either absolute, or squared
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################
#' This function takes two ultrametric phylogenetic trees, calculates the normalized Lineage-Through-Time statistic for both trees and then calculates the difference between the two statistics.
#' @title Calculate the difference between two normalized Lineage-Through-Time curves, given two phylogenetic trees.
#' @usage nLTTstat(tree1, tree2, distance_method = "abs", ignore_stem = TRUE)
#' @param tree1 an object of class \code{"phylo"}
#' @param tree2 an object of class \code{"phylo"}
#' @param distance_method Chosen measurement of distance between the two nLTT curves, options are (case sensitive):\cr
#'   - "abs": use the absolute distance\cr
#'   - "squ": use the squared distance;\cr
#' @param ignore_stem a boolean whether to ignore the stem length
#' @return The difference between the two nLTT statistics
#' @author Thijs Janzen
#' @examples
#'   data(exampleTrees)
#'   nltt_plot(exampleTrees[[1]])
#'   nltt_lines(exampleTrees[[2]], lty=2)
#'   nLTTstat(
#'     exampleTrees[[1]], exampleTrees[[2]],
#'     distance_method = "abs", ignore_stem = TRUE)
#' @export
nLTTstat <- function(  # nolint keep function name non-all-lowercase, due to backwards compatibility
  tree1,
  tree2,
  distance_method = "abs",
  ignore_stem = TRUE
) {
  if (!inherits(tree1, "phylo")) {
    # Just checking
    stop("nLTTstat: ",
      "tree1 must be of class 'phylo', ",
      "but was of type '", class(tree1), "' instead")
  }
  if (!inherits(tree2, "phylo")) {
    # Just checking
    stop("nLTTstat: ",
      "tree2 must be of class 'phylo', ",
      "but was of type '", class(tree2), "' instead")
  }

  if ( distance_method != "abs" && distance_method != "squ") {
    stop("nLTTstat: distance method unknown")
  }

  if (!is.logical(ignore_stem)) {
    stop("nLTTstat: ignore_stem must be logical")
  }

  diff <- nLTT::nltt_diff( tree1, tree2, distance_method, ignore_stem)
  return (diff)
}

################################################################################
#
# @brief Wrapper to calculate the nLTT statistic - exact version
#
# @date Last modified: 2016-26-04
# @author Thijs Janzen
# @since 2016-26-04, version 1.2
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @param    distance_method        string     Method to calculate the difference, either absolute, or squared
# @param    ignore_stem            logical    Should the phylogeny its stem be ignored?
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################
#' Calculate the exact difference between two normalized Lineage-Through-Time curves, given two phylogenetic trees.
#' @description This function takes two ultrametric phylogenetic trees, calculates the normalized Lineage-Through-Time statistic for both trees and then calculates the exact difference between the two statistics. Whereas the function \code{nLTTstat} uses an approximation to calculate the difference (which is faster for large trees), the function \code{nLTTstat_exact} calculates the exact difference, and should generally be preferred. Although the estimates are highly similar, \code{nLTTstat_exact} tends to return slightly higher values.
#' @usage
#'   nLTTstat_exact(tree1, tree2, distance_method = "abs", ignore_stem = TRUE)
#' @param tree1 an object of class \code{"phylo"}
#' @param tree2 an object of class \code{"phylo"}
#' @param distance_method
#'   Chosen measurement of distance between the two nLTT curves, options are (case sensitive):\cr
#'   - "abs": use the absolute distance.\cr
#'   - "squ": use the squared distance
#' @param ignore_stem a boolean whether to ignore the stem length
#' @return The exact difference between the two nLTT statistics
#' @author Thijs Janzen
#' @examples
#'   data(exampleTrees)
#'   nltt_plot(exampleTrees[[1]])
#'   nltt_lines(exampleTrees[[2]], lty = 2)
#'   nLTTstat_exact(
#'     exampleTrees[[1]],
#'     exampleTrees[[2]],
#'     distance_method = "abs",
#'     ignore_stem = TRUE
#'   )
#' @export
nLTTstat_exact <- function( # nolint keep function name non-all-lowercase, due to backwards compatibility
  tree1,
  tree2,
  distance_method = "abs",
  ignore_stem = TRUE
) {
  if (!inherits(tree1, "phylo")) {
    # Just checking
    stop("nLTTstat_exact: tree1 must be of class 'phylo', ",
      "but was of type '", class(tree1), "' instead")
  }
  if (!inherits(tree2, "phylo")) {
    # Just checking
    stop("nLTTstat_exact: tree2 must be of class 'phylo', ",
      "but was of type '", class(tree2), "' instead")
  }
  if (distance_method != "abs" && distance_method != "squ") {
    stop("nLTTstat_exact: distance method unknown")
  }
  if (!is.logical(ignore_stem)) {
    stop("nLTTstat_exact: ignore_stem must be logical")
  }
  diff <- nLTT::nltt_diff_exact( tree1, tree2, distance_method, ignore_stem)
  return (diff)
}
