#' Collect the normalized number of lineages from the stem age
#' @param phylogeny a phylogeny of class 'phylo'
#' @return branching times, in time units before the present
#' @examples
#'   phylogeny <- ape::read.tree(text = "((a:2,b:2):1,c:3);")
#'   phylogeny$root.edge <- 2 # nolint ape variable name
#'   testthat::expect_true(
#'     all(nLTT::get_branching_times(phylogeny) == c(5, 3, 2)))
#' @author Richel Bilderbeek
#' @export
get_norm_n <- function(phylogeny) {
  ns <- nLTT::get_n_lineages(phylogeny)
  # Repeat the last value to have te ns at the present
  ns <- c(ns, utils::tail(ns, n = 1))
  # Reverse to have 'times ago' (e.g. 4,3,0) -> (0, -1, -4) -> (0,1,4)
  ns <- ns - ns[1]
  ns <- -ns
  # Normalize
  ns <- ns / utils::tail(ns, n = 1)
  return(ns)
}
