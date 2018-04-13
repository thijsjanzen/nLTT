#' Collect the branching times from the stem age
#' @param phylogeny a phylogeny of class 'phylo'
#' @return branching times, in time units before the present
#' @examples
#'   phylogeny <- ape::read.tree(text = "((a:2,b:2):1,c:3);")
#'   phylogeny$root.edge <- 2 # nolint ape variable name
#'   testthat::expect_true(
#'     all(nLTT::get_branching_times(phylogeny) == c(5, 3, 2)))
#' @author Richel Bilderbeek
#' @export
get_branching_times <- function(phylogeny) {
  brts <- ape::branching.times(phylogeny)
  if (!is.null(phylogeny$root.edge)) { # nolint ape variable name
    # Put the stem age first
    return(c(brts[1] + phylogeny$root.edge, brts))
  }
  brts
}
