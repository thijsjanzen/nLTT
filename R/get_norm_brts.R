#' Collect the normalized branching times from the stem age
#' @param phylogeny a phylogeny of class 'phylo'
#' @return branching times, in time units before the present
#' @examples
#'   phylogeny <- ape::read.tree(text = "((a:2,b:2):1,c:3);")
#'   phylogeny$root.edge <- 2
#'   testthat::expect_true(
#'     all(nLTT::get_branching_times(phylogeny) == c(5, 3, 2)))
#' @author Richel Bilderbeek
#' @export
get_norm_brts <- function(phylogeny) {
  brts <- nLTT::get_branching_times(phylogeny)
  # Repeat the last value to have te brts at the present
  brts <- c(brts, utils::tail(brts, n = 1))
  # Reverse to have 'times ago' (e.g. 4,3,0) -> (0, -1, -4) -> (0,1,4)
  brts <- brts - brts[1]
  brts <- -brts
  # Normalize
  brts <- brts / utils::tail(brts, n = 1)
 return (brts)
}
