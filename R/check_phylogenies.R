#' Check if the input is a valid collection of one or more phylogenies
#'
#' Will \link{stop} if not
#' @inheritParams default_params_doc
#' @export
check_phylogenies <- function(phylogenies) {
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
}
