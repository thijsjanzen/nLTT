#' Get the average nLTT from a collection of phylogenies
#'
#' @param phylogenies the phylogenies, supplied as either a list
#'   or a multiPhylo object, where the phylogenies are of type 'phylo'
#' @param dt The timestep resolution, where 1/dt is
#'   the number of points evaluated
#' @return A matrix of timepoints with the average number of (normalized)
#'   lineages through (normalized) time
#' @examples
#'   get_average_nltt_matrix(c(ape::rcoal(10), ape::rcoal(20)))
#'
#' @author Richèl J.C. Bilderbeek
#' @export
get_average_nltt_matrix <- function(
  phylogenies,
  dt = 0.001
) {

  if (length(phylogenies) < 1) {
    stop("there must be at least one phylogeny supplied")
  }
  if (class(phylogenies) != "multiPhylo" && class(phylogenies) != "list") {
    stop("phylogenies must be of class 'multiPhylo' or 'list'")
  }
  if (!inherits(phylogenies[[1]], "phylo")) {
    # Stop imposed by ape::ltt.plot.coords
    stop("phylogenies must be of type phylo")
  }
  if (dt <= 0.0 || dt >= 1.0) {
    stop("dt must be between (not including) zero and one")
  }

  sz <- length(phylogenies)

  nltts <- NULL
  for (phylogeny in phylogenies) {
    testit::assert(length(phylogeny$tip.label) > 0)
    nltts <- c(nltts, list(nLTT::get_phylogeny_nltt_matrix(phylogeny)))
  }
  testit::assert(length(nltts) == length(phylogenies))

  stretch_matrices <- NULL
  for (nltt in nltts) {
    stretch_matrix <- nLTT::stretch_nltt_matrix(
      m = nltt, dt = dt, step_type = "upper"
    )
    stretch_matrices <- c(stretch_matrices, list(stretch_matrix))
  }
  testit::assert(length(stretch_matrices) == length(nltts))

  xy <- stretch_matrices[[1]]
  if (sz > 1) {
    for (i in seq(2, sz)) {
      xy <- (xy + stretch_matrices[[i]])
    }
  }
  xy <- (xy / sz)

  xy
}
