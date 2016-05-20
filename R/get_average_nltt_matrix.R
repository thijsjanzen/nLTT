#' Get the average nLTT from a collection of phylogenies
#'
#' @param phylogenies the phylogenies, supplied as either a list or a multiPhylo object, where the phylogenies are of type 'phylo'
#' @param dt The timestep resolution, where 1/dt is the number of points evaluated
#' @return A matrix of timepoints with the average number of (normalized) lineages through (normalized) time
#' @examples
#'   get_average_nltt(c(ape::rcoal(10), ape::rcoal(20)))
#'
#' @author Richel Bilderbeek
#' @export
get_average_nltt_matrix <- function(
  phylogenies,
  dt = 0.001) {
  ribir::get_average_nltt_matrix_impl_1(phylogenies = phylogenies, dt = dt)
}



#' Internal function. First implementation of get_average_nltt_matrix.
#'
#' @param phylogenies the phylogenies, supplied as either a list or a multiPhylo object, where the phylogenies are of type 'phylo'
#' @param dt The timestep resolution, where 1/dt is the number of points evaluated
#' @return A matrix of timepoints with the average number of (normalized) lineages through (normalized) time
#' @author Richel Bilderbeek
#' @export
get_average_nltt_matrix_impl_1 <- function(phylogenies, dt) {
  if (length(phylogenies) < 1) {
    stop("get_average_nltt_matrix: ",
         "there must be at least one phylogeny supplied")
  }
  if (class(phylogenies) != "multiPhylo" && class(phylogenies) != "list") {
    stop("get_average_nltt_matrix: ",
      "phylogenies must be of class 'multiPhylo' or 'list', ",
      "used '", class(phylogenies), "' instead")
  }
  if (!inherits(phylogenies[[1]], "phylo")) {
    # Stop imposed by ape::ltt.plot.coords
    stop("get_average_nltt_matrix: ",
      "phylogenies must be of type phylo, ",
      "instead of '", class(phylogenies[[1]]), "'")
  }
  if (dt <= 0.0 || dt >= 1.0) {
    stop("get_average_nltt_matrix: ",
      "dt must be between (not including) zero and one, ",
      "dt was ", dt, " instead")
  }

  sz <- length(phylogenies)

  nltts <- NULL
  for (phylogeny in phylogenies) {
    nltts <- c(nltts, list(ribir::get_phylogeny_nltt_matrix(phylogeny)))
  }
  testit::assert(length(nltts) == length(phylogenies))

  stretch_matrices <- NULL
  for (nltt in nltts) {
    stretch_matrix <- ribir::stretch_nltt_matrix(
      nltt, dt = dt, step_type = "upper"
    )
    stretch_matrices <- c(stretch_matrices, list(stretch_matrix))
  }
  testit::assert(length(stretch_matrices) == length(nltts))

  xy <- stretch_matrices[[1]]
  for (i in seq(2, sz)) {
    xy <- (xy + stretch_matrices[[i]])
  }
  xy <- (xy / sz)

  xy
}

#' Internal function. Second implementation of get_average_nltt_matrix. A stub now.
#'
#' @param phylogenies the phylogenies, supplied as either a list or a multiPhylo object, where the phylogenies are of type 'phylo'
#' @param dt The timestep resolution, where 1/dt is the number of points evaluated
#' @return A matrix of timepoints with the average number of (normalized) lineages through (normalized) time
#' @author Richel Bilderbeek
#' @export
get_average_nltt_matrix_impl_2 <- function(phylogenies, dt) {
  xy <- ribir::get_average_nltt_matrix_impl_1(
    phylogenies = phylogenies, dt = dt
  )
  xy
}
