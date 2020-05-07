#' Collect the number of lineages from the stem age
#' @param phylogeny a phylogeny of class 'phylo'
#' @return number of lineages, will go from 1 to the number of tips,
#'   if there is a stem, will go from 2 to the number of tips
#'   if there is no stem
#' @examples
#'   phylogeny <- ape::read.tree(text = "((a:2,b:2):1,c:3);")
#'   testthat::expect_true(
#'     all(nLTT::get_n_lineages(phylogeny) == c(2, 3)))
#'   phylogeny$root.edge <- 2 # nolint ape variable name
#'   testthat::expect_true(
#'     all(nLTT::get_n_lineages(phylogeny) == c(1, 2, 3)))
#' @author Richèl Bilderbeek
#' @export
get_n_lineages <- function(phylogeny) {
  if (!is.null(phylogeny$root.edge)) { # nolint ape variable name
    return(seq_along(phylogeny$tip.label))
  }
  2:length(phylogeny$tip.label)
}
