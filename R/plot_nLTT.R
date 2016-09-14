################################################################################
#
# @brief plot the normalized nLTT plot of a tree
#
# @date Last modified: 2014-20-09
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#
# @param    tree1                  phylo      Phylogenetic tree
# @param    xlab                   string     Label on the x-axis
# @param    ylab                   string     Label on the y-axis
#
#
################################################################################
#' This function uses a modified version of the ltt.plot function from \code{"ape"} to plot the normalized number of lineages through normalized time, where the number of lineages is normalized by dividing by the number of tips of the tree, and the time is normalized by the total time between the most common recent ancestor and the present, such that t(MRCA) = 0 & t(present) = 1.
#' @title Normalized version of the ape function ltt.plot
#' @usage nltt_plot(phy, xlab = "Normalized Time", ylab = "Normalized Lineages", ...)
#' @param phy an object of class \code{"phylo"};
#' @param xlab a character string (or a variable of mode character)
#'   giving the label for the \eqn{x}-axis (default is "Normalized Time").
#' @param ylab a character string (or a variable of mode character)
#'   giving the label for the \eqn{y}-axis (default is "Normalized Lineages").
#' @param \dots further graphical arguments that can be passed to \code{plot()}
#' @author Thijs Janzen
#' @examples
#'   data(exampleTrees)
#'   nltt_plot(exampleTrees[[1]])
#' @export
nltt_plot <- function( phy, xlab = "Normalized Time",
  ylab = "Normalized Lineages", ...) {

  if (!inherits(phy, "phylo")) {
    # Stop imposed by ape::ltt.plot.coords
    stop("nltt_plot: ",
      "phylogeny must be of class 'phylo', ",
      "but was of type '", class(phy), "' instead")
  }

  #we use the ltt.plot.coords function from the package ape
  xy <- ape::ltt.plot.coords( phy, backward = TRUE, tol = 1e-6)
  xy[, 2] <- xy[, 2] / max(xy[, 2]) #normalize number lineages

  xy[, 1] <- xy[, 1] + abs( min( xy[, 1])) #make sure time runs from 0..T
  xy[, 1] <- xy[, 1] / max( xy[, 1])      #normalize time

  graphics::plot.default(xy, xlab = xlab, ylab = ylab, xaxs = "r", yaxs = "r",
    type = "S", ...) #type = "S" ensures a stepwise function
}

################################################################################
#
# @brief add a nLTT line to an existing nLTT plot
# @date Last modified: 2014-20-09
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#
# @param    tree1                  phylo      Phylogenetic tree
#
#
################################################################################
#' This is a modified version of the \code{ape} function ltt.lines: add the normalized Lineage-Through-Time statistic of a phylogenetic tree to an already existing plot
#' @title Normalized version of the ape function ltt.lines.
#' @usage nltt_lines(phy, ...)
#' @param phy an object of class \code{"phylo"};
#' @param \dots further graphical arguments that can be passed to \code{lines()}
#' @author Thijs Janzen
#' @examples
#'   data(exampleTrees)
#'   nltt_plot(exampleTrees[[1]])
#'   nltt_lines(exampleTrees[[2]], lty=2)
#' @export
nltt_lines <- function(phy, ...) {

  if (!inherits(phy, "phylo")) {
    # Stop imposed by ape::ltt.plot.coords
    stop("nltt_lines: ",
      "phylogeny must be of class 'phylo', ",
      "but was of type '", class(phy), "' instead")
  }
  xy <- ape::ltt.plot.coords( phy, backward = TRUE, tol = 1e-6)
  xy[, 2] <- xy[, 2] / max( xy[, 2]) #normalize number lineages

  xy[, 1] <- xy[, 1] + abs( min( xy[, 1])) #make sure time runs from 0..T
  xy[, 1] <- xy[, 1] / max( xy[, 1])      #normalize time
  graphics::lines(xy, type = "S", ...)
}
