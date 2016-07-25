################################################################################
#
# @brief Calculates the exact, difference between the lineage through time curves of tree1 & tree2 (normalized in time and for the number of lineages)
#
# @date Last modified: 2016-05-20
# @author Thijs Janzen
# @since 2016-05-20, version 1.2.1
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @param    distance_method        string     absolute, or squared distance?
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################

nltt_diff_exact <- function(tree1, tree2, distance_method = "abs") { # nolint
  #branching times of tree1, including the present time (0)
  b_times <- c(-1 * rev(sort(ape::branching.times(tree1))), 0);

  #the number of lineages, we assume that the
  # first branching time indicates the crown age.
  lineages <- c(2:length(b_times), length(b_times));

  b_times_N <- 1 - b_times / min(b_times); #normalize branching times
  lineages_N <- lineages / max(lineages);  #normalize lineages

  b_times2 <- c(-1 * rev(sort(ape::branching.times(tree2))), 0);
  lineages2 <- c(2:length(b_times2), length(b_times2));
  b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times
  lineages2_N <- lineages2 / max(lineages2);  #normalize lineages

  #make a list of all branching times, and remove duplicates
  all_b_times <- unique(sort(c(b_times_N, b_times2_N)));
  diff <- 0;
  #iterate through all branching times
  for (k in 2:length( all_b_times)) {
      tim <- all_b_times[k];
      #find the index of the first branching time
      #that is up to the focal branching time
      index1 <- max(which(b_times_N < tim));
      index2 <- max(which(b_times2_N < tim));  #same for the other tree

      #find the number of lineages at time "tim" for tree 1
      lins1 <- lineages_N[max(index1, 1)];
      #find the number of lineages at time "tim" for tree 2
      lins2 <- lineages2_N[max(index2, 1)];

      #the amount of time that this difference persisted
      dt <- all_b_times[k] - all_b_times[k - 1]
      if (distance_method == "abs") {
        diff <- diff + dt * abs( lins1 - lins2);     #update the difference
      }
      if (distance_method == "squ")  {
        diff <- diff + dt * ( lins1 - lins2) * ( lins1 - lins2);
      }
  }
  return ( diff);
}
################################################################################
#
# @brief Calculates the exact difference between the lineage through time curves of tree1 & tree2 (normalized in time and for the number of lineages)
#
# @date Last modified: 2016-05-20
# @author Thijs Janzen
# @since 2016-05-20, version 1.2.1
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @param    distance_method        string     absolute, or squared distance?
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################

nltt_diff <- function(tree1, tree2, distance_method = "abs")  {
  #branching times of tree1, including the present time (0)
  b_times    <- c(-1 * rev( sort( ape::branching.times( tree1))), 0);

  #the number of lineages
  #we assume that the first branching time indicates the crown age.
  lineages   <- c( 2:length( b_times), length( b_times));
  b_times_N  <- 1 - b_times / min( b_times); #normalize branching times
  lineages_N <- lineages / max( lineages); #normalize lineages

  #method = constant ensures a step function
  ltt1       <- approxfun( b_times_N, lineages_N, method = "constant");

  #branching times of tree2, including the present time (0)
  b_times2    <- c(-1 * rev( sort( ape::branching.times( tree2))), 0);

  #the number of lineages
  #we assume that the first branching time indicates the crown age.
  lineages2   <- c( 2:length( b_times2), length( b_times2));

  b_times2_N  <- 1 - b_times2 / min( b_times2); #normalize branching times
  lineages2_N <- lineages2 / max( lineages2);  #normalize lineages
  #method = constant ensures a step function
  ltt2        <- approxfun( b_times2_N, lineages2_N, method = "constant");

  #function f is the absolute difference in time t: 0 <= t < 1
  f <- function( t, x, p) {
       if ( distance_method == "abs" ) {
         output <- abs( ltt1(t) - ltt2(t) );
       }
       if ( distance_method == "squ" ) {
         output <- (ltt1(t) - ltt2(t)) * (ltt1(t) - ltt2(t))
       }
       return( list( output) );
  }

  #evaluation points of the integration function below
  #more points ensures higher precision.
  times <- ( 0:100 ) / 100;
  #integrate over t: 0 < t < 1, notice tcrit=0 indicating t
  #should never be larger than 0:
  int_1 <- deSolve::lsoda( 0, times, func = f, tcrit = c( 1 ) );
  total_area <- int_1[length(times), 2]
  return( total_area[[1]])
}

################################################################################
#
# @brief Wrapper to calculate the nLTT statistic
#
# @date Last modified: 2015-21-04
# @author Thijs Janzen
# @since 2015-21-04, version 1.1
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @param    distance_method        string     Method to calculate the difference, either absolute, or squared
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################

nLTTstat <- function( tree1, tree2, distance_method = "abs") { # nolint keep function name non-all-lowercase, due to backwards compatibility
  if (!inherits(tree1, "phylo")) {
    # Just checking
    stop("nLTTstat: ",
      "tree1 must be of class 'phylo', ",
      "but was of type '", class(tree1), "' instead")
  }
  if (!inherits(tree2, "phylo")) {
    # Just checking
    stop("nLTTstat: ",
      "tree1 must be of class 'phylo', ",
      "but was of type '", class(tree2), "' instead")
  }

  if ( distance_method != "abs" && distance_method != "squ") {
    cat( "chosen unknown distance method!\n" );
    flush.console();
  }
  diff <- nLTT::nltt_diff( tree1, tree2, distance_method);
  return (diff);
}

################################################################################
#
# @brief Wrapper to calculate the nLTT statistic - exact version
#
# @date Last modified: 2016-26-04
# @author Thijs Janzen
# @since 2016-26-04, version 1.2
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @param    distance_method        string     Method to calculate the difference, either absolute, or squared
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################

nLTTstat_exact <- function(tree1, tree2, distance_method = "abs") { # nolint keep function name non-all-lowercase, due to backwards compatibility
  if (!inherits(tree1, "phylo")) {
    # Just checking
    stop("nLTTstat_exact: ",
      "tree1 must be of class 'phylo', ",
      "but was of type '", class(tree1), "' instead")
  }
  if (!inherits(tree2, "phylo")) {
    # Just checking
    stop("nLTTstat_exact: ",
      "tree2 must be of class 'phylo', ",
      "but was of type '", class(tree2), "' instead")
  }
  if (distance_method != "abs" && distance_method != "squ") {
    cat("chosen unknown distance method!\n");
    flush.console()
  }
  diff <- nLTT::nltt_diff_exact( tree1, tree2, distance_method)
  return (diff)
}
