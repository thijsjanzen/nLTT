
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


nLTT.plot <- function(phy,xlab="Normalized Time",ylab="Normalized Lineages", ...) {

 xy <- ltt.plot.coords(phy,backward=TRUE,tol=1e-6);  #we use the ltt.plot.coords function from the package ape
 xy[,2] <- xy[,2]/max(xy[,2]); #normalize number lineages
 
 xy[,1] <- xy[,1] + abs(min(xy[,1])); #make sure time runs from 0..T
 xy[,1] <- xy[,1] / max(xy[,1]);      #normalize time
 
 plot.default(xy, xlab = xlab, ylab = ylab, xaxs = "r", yaxs = "r", 
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



nLTT.lines <- function(phy, ...) {

  xy <- ltt.plot.coords(phy, backward=TRUE, tol=1e-6)
  xy[,2] <- xy[,2]/max(xy[,2]); #normalize number lineages
 
  xy[,1] <- xy[,1] + abs(min(xy[,1])); #make sure time runs from 0..T
  xy[,1] <- xy[,1] / max(xy[,1]);      #normalize time
  lines(xy, type = "S", ...)
}
		