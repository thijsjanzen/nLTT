#' Extract the nLTT matrix from a phylogeny
#' @param phylogeny A phylogeny of type phylo
#' @return a matrix
#' @author Richèl Bilderbeek
#' @export
get_phylogeny_nltt_matrix <- function(phylogeny) {
  if (!inherits(phylogeny, "phylo")) {
    # Stop imposed by ape::ltt.plot.coords
    stop("get_phylogeny_nltt_matrix: ",
      "phylogeny must be of class 'phylo', ",
      "but was of type '", class(phylogeny), "' instead")
  }
  xy <- ape::ltt.plot.coords(phylogeny, backward = TRUE, tol = 1e-06)
  xy[, 2] <- xy[, 2] / max(xy[, 2])
  xy[, 1] <- xy[, 1] + abs(min(xy[, 1]))
  xy[, 1] <- xy[, 1] / max(xy[, 1])
  xy
}
