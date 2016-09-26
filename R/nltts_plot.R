#' Get the average nLTT from a collection of phylogenies
#'
#' @param phylogenies the phylogenies, where the phylogenies are of type 'phylo'
#' @param dt The timestep resolution, where 1/dt is the number of points evaluated
#' @param plot_nltts Also plot each nLLT line
#' @param xlab Label on the x axis
#' @param ylab Label on the y axis
#' @param replot If false, start a clean plot. If true, plot the new data over the current
#' @param ... Plotting options
#' @return Nothing
#' @examples
#'   nltts_plot(c(ape::rcoal(10), ape::rcoal(10)))
#'   nltts_plot(c(ape::rcoal(10), ape::rcoal(20)), dt = 0.1)
#'
#' @author Richel Bilderbeek
#' @export
nltts_plot <- function(
  phylogenies,
  dt = 0.001,
  plot_nltts = FALSE,
  xlab = "Normalized Time",
  ylab = "Normalized Lineages",
  replot = FALSE,
  ...
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


  xy <- nLTT::get_average_nltt_matrix(
    phylogenies = phylogenies,
    dt = dt
  )

  # Set the shape of the plot
  if (replot == FALSE) {
    graphics::plot.default(
      xy,
      xlab = "Normalized Time",
      ylab = "Normalized Lineages",
      xaxs = "r",
      yaxs = "r",
      type = "S",
      xlim = c(0, 1),
      ylim = c(0, 1),
      ...
    )
  }

  # Draw the nLTTS plots used
  if (plot_nltts == TRUE) {
    # Copied

    nltts <- NULL
    for (phylogeny in phylogenies) {
      nltts <- c(nltts, list(nLTT::get_phylogeny_nltt_matrix(phylogeny)))
    }
    testit::assert(length(nltts) == length(phylogenies))

    stretch_matrices <- NULL
    for (nltt in nltts) {
      stretch_matrix <- nLTT::stretch_nltt_matrix(
        nltt, dt = dt, step_type = "upper"
      )
      stretch_matrices <- c(stretch_matrices, list(stretch_matrix))
    }
    testit::assert(length(stretch_matrices) == length(nltts))
    # End of copy

    for (stretch_matrix in stretch_matrices) {
      graphics::lines.default(
        stretch_matrix,
        xaxs = "r",
        yaxs = "r",
        type = "S",
        col = "grey",
        xlim = c(0, 1),
        ylim = c(0, 1)
      )
    }
  }

  # Redraw the average nLTT plot
  graphics::lines.default(
    xy,
    xaxs = "r",
    yaxs = "r",
    type = "S",
    xlim = c(0, 1),
    ylim = c(0, 1),
    ...
  )
}
