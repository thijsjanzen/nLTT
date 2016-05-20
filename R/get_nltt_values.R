#' Collects the nLTT values of all phylogenies in the melted/uncast/long form
#'
#' @param phylogenies the phylogenies, supplied as either a list or a multiPhylo object, where the phylogenies are of type 'phylo'
#' @param dt The timestep resolution, where 1/dt is the number of points evaluated
#' @return A dataframe of timepoints with the nLTT value of each phylogeny in time
#' @examples
#'   library(ape)
#'   library(ggplot2)
#'   library(ribir)
#'
#'   # Create some random phylogenies
#'   phylogeny1 <- rcoal(10)
#'   phylogeny2 <- rcoal(20)
#'   phylogeny3 <- rcoal(30)
#'   phylogeny4 <- rcoal(40)
#'   phylogeny5 <- rcoal(50)
#'   phylogeny6 <- rcoal(60)
#'   phylogeny7 <- rcoal(70)
#'   phylogenies <- c(phylogeny1, phylogeny2, phylogeny3,
#'     phylogeny4, phylogeny5, phylogeny6, phylogeny7
#'   )
#'
#'   # Obtain the nLTT values
#'   dt <- 0.2
#'   nltt_values <- get_nltt_values(phylogenies, dt = dt)
#'
#'   # Check properties of nltt_values
#'   testit::assert(names(nltt_values) == c("id", "t", "nltt"))
#'   nltt_values_per_phylogeny <- (1 + (1 / dt))
#'   n_phylogenies <- length(phylogenies)
#'   testit::assert(nrow(nltt_values)
#'     == nltt_values_per_phylogeny * n_phylogenies
#'   )
#'
#'   # Plot the phylognies, where the individual nLTT values are visible
#'   qplot(t, nltt, data = nltt_values, geom = "point",
#'     ylim = c(0,1),
#'     main = "Average nLTT plot of phylogenies", color = id, size = I(0.1)
#'   ) + stat_summary(
#'     fun.data = "mean_cl_boot", color = "red", geom = "smooth"
#'   )
#'
#'   # Plot the phylognies, where the individual nLTT values are omitted
#'   qplot(
#'     t, nltt, data = nltt_values, geom = "blank", ylim = c(0,1),
#'     main = "Average nLTT plot of phylogenies"
#'   ) + stat_summary(
#'     fun.data = "mean_cl_boot", color = "red", geom = "smooth"
#'   )
#' @author Richel Bilderbeek
#' @export
get_nltt_values <- function(phylogenies, dt) {
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

  n_cols <- 3 # ID, t, nLTT(t)
  n_rows_per_phylogeny <- (1.0 / dt) + 1
  n_phylogenies <- length(phylogenies)
  n_rows <- n_rows_per_phylogeny * n_phylogenies
  m <- matrix(nrow = n_rows, ncol = n_cols)
  for (i in seq(1, n_phylogenies)) {
    testit::assert(i >= 1 && i <= n_phylogenies)
    new_col <- ribir::stretch_nltt_matrix(
      m = ribir::get_phylogeny_nltt_matrix(phylogenies[[i]]),
      dt = dt,
      step_type = "upper"
    )
    testit::assert(length(new_col[, 2]) == n_rows_per_phylogeny)
    row_from_index <- (i - 1) * n_rows_per_phylogeny + 1
    row_to_index <- (i - 0) * n_rows_per_phylogeny
    testit::assert(row_from_index >= 1)
    testit::assert(row_from_index <= n_rows)
    testit::assert(row_to_index >= 1)
    testit::assert(row_to_index <= n_rows)
    m[row_from_index:row_to_index, 1] <- i # Recycling
    m[row_from_index:row_to_index, 2] <- seq(0, 1, dt)
    m[row_from_index:row_to_index, 3] <- new_col[, 2, drop = FALSE]
  }
  z <- as.data.frame(x = m)
  colnames(z) <- c("id", "t", "nltt")
  z[, 1] <- sapply(z[, 1], as.integer)
  z
}
