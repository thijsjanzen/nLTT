
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
# @param    distanceMethod         string     Method to calculate the difference, either absolute, or squared
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################

nLTTstat <- function(tree1, tree2, distanceMethod = "abs")
{
 	diff <- -10;
 	if(distanceMethod == "abs") diff = normLTTdiffABS(tree1,tree2); #nLTT statistic using the absolute difference
 	if(distanceMethod == "squ") diff = normLTTdiffSQ(tree1,tree2);  #nLTT statistic using the squared difference
 	 	
 	if(diff < 0) {cat("chosen unknown distance method!\n"); flush.console();}
 	return(diff);	
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
# @param    distanceMethod         string     Method to calculate the difference, either absolute, or squared
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################


nLTTstat_exact <- function(tree1, tree2, distanceMethod = "abs")
{
 	diff <- -10;
 	if(distanceMethod == "abs") diff = normLTTdiffexactABS(tree1,tree2); #nLTT statistic using the absolute difference
 	if(distanceMethod == "squ") diff = normLTTdiffexactSQ(tree1,tree2);  #nLTT statistic using the squared difference
 	 	
 	if(diff < 0) {cat("chosen unknown distance method!\n"); flush.console();}
 	return(diff);	
}

################################################################################
# 
# @brief Calculates the absolute, exact, difference between the lineage through time curves of tree1 & tree2 (normalized in time and for the number of lineages)
#
# @date Last modified: 2016-26-04
# @author Thijs Janzen
# @since 2016-26-04, version 1.2
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################


normLTTdiffexactABS <- function(tree1,tree2)
{
	b_times <- c(-1 * rev(sort(branching.times(tree1))),0); #branching times of tree1, including the present time (0)
  	lineages <- c(2:length(b_times),length(b_times));     #the number of lineages, we assume that the first branching time indicates the crown age.
  	b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  	lineages_N <- lineages / max(lineages);  #normalize lineages  	
	
	b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);
    lineages2 <- c(2:length(b_times2),length(b_times2));  
    b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
    lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
	
	allBtimes <- unique(sort(c(b_times_N,b_times2_N))); #make a list of all branching times, and remove duplicates
	diff <- 0;
	for(k in 2:length(allBtimes)) #iterate through all branching times
	{
			tim <- allBtimes[k];
			index1 <- max(which(b_times_N <= tim));   #find the index of the first branching time that is up to the focal branching time
			index2 <- max(which(b_times2_N <= tim));  #same for the other tree
			lins1 <- lineages_N[index1];              #find the number of lineages at time "tim" for tree 1
			lins2 <- lineages2_N[index2];             #find the number of lineages at time "tim" for tree 2
			dt <- allBtimes[k] - allBtimes[k-1]       #the amount of time that this difference persisted
			diff <- diff + dt * abs(lins1-lins2);     #update the difference
	}                                             
	return(diff);
}

################################################################################
# 
# @brief Calculates the squared, exact, difference between the lineage through time curves of tree1 & tree2 (normalized in time and for the number of lineages)
#
# @date Last modified: 2016-26-04
# @author Thijs Janzen
# @since 2016-26-04, version 1.2
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################


normLTTdiffexactSQ <- function(tree1,tree2)
{
	b_times <- c(-1 * rev(sort(branching.times(tree1))),0); #branching times of tree1, including the present time (0)
  	lineages <- c(2:length(b_times),length(b_times));     #the number of lineages, we assume that the first branching time indicates the crown age.
  	b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  	lineages_N <- lineages / max(lineages);  #normalize lineages  	
	
	b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);
    lineages2 <- c(2:length(b_times2),length(b_times2));  
    b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
    lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
	
	allBtimes <- unique(sort(c(b_times_N,b_times2_N))); #make a list of all branching times, and remove duplicates
  diff <- 0;
  for(k in 2:length(allBtimes)) #iterate through all branching times
  {
      tim <- allBtimes[k];
      index1 <- max(which(b_times_N <= tim));   #find the index of the first branching time that is up to the focal branching time
      index2 <- max(which(b_times2_N <= tim));  #same for the other tree
      lins1 <- lineages_N[index1];              #find the number of lineages at time "tim" for tree 1
      lins2 <- lineages2_N[index2];             #find the number of lineages at time "tim" for tree 2
      dt <- allBtimes[k] - allBtimes[k-1]       #the amount of time that this difference persisted
      diff <- diff + dt * (lins1-lins2)*(lins1-lins2);    #update the difference
  }       
	return(diff);
}


################################################################################
# 
# @brief Calculates the absolute, difference between the lineage through time curves of tree1 & tree2 (normalized in time and for the number of lineages). Uses an approximation, which is faster for large trees.
#
# @date Last modified: 2014-20-09
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################




normLTTdiffABS <- function(tree1, tree2) {

  b_times <- c(-1 * rev(sort(branching.times(tree1))),0);  #branching times of tree1, including the present time (0)
  lineages <- c(2:length(b_times),length(b_times));  #the number of lineages, we assume that the first branching time indicates the crown age.
  b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  lineages_N <- lineages / max(lineages);  #normalize lineages  
  ltt1 <- approxfun(b_times_N,lineages_N,method="constant"); #method = constant ensures a step function

  b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);  #branching times of tree2, including the present time (0)
  lineages2 <- c(2:length(b_times2),length(b_times2));      #the number of lineages, we assume that the first branching time indicates the crown age.
  b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
  lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
  ltt2 <- approxfun(b_times2_N,lineages2_N,method="constant"); #method = constant ensures a step function

  f <- function(t,x,p) { #function f is the absolute difference in time t: 0 <= t < 1
       output <- abs( ltt1(t) - ltt2(t));
       return(list(output));
  }
  
  times <- (0:100)/100; #evaluation points of the integration function below - more points ensures higher precision.
  int_1 <- lsoda(0,times,func=f,tcrit=c(1)); #integrate over t: 0 < t < 1, notice tcrit=0 indicating t should never be larger than 0.
  total_area <- int_1[length(times),2]
  
  return(total_area);
}


################################################################################
# 
# @brief Calculates the squared, difference between the lineage through time curves of tree1 & tree2 (normalized in time and for the number of lineages). Uses an approximation, which is faster for large trees.
#
# @date Last modified: 2014-20-09
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#
# @param    tree1                  phylo      First phylogenetic tree
# @param    tree2                  phylo      Second phylogenetic tree
# @return                          scalar     normalized Lineage-Through-Time difference between tree1 & tree2
#
################################################################################


normLTTdiffSQ <- function(tree1, tree2) {

  b_times <- c(-1 * rev(sort(branching.times(tree1))),0);  #branching times of tree1, including the present time (0)
  lineages <- c(2:length(b_times),length(b_times)); #the number of lineages, we assume that the first branching time indicates the crown age.
  b_times_N <- 1 - b_times / min(b_times); #normalize branching times  
  lineages_N <- lineages / max(lineages);  #normalize lineages  
  ltt1 <- approxfun(b_times_N,lineages_N,method="constant"); #method = constant ensures a step function

  b_times2 <- c(-1 * rev(sort(branching.times(tree2))),0);  #branching times of tree2, including the present time (0)
  lineages2 <- c(2:length(b_times2),length(b_times2));   #the number of lineages, we assume that the first branching time indicates the crown age.
  b_times2_N <- 1 - b_times2 / min(b_times2); #normalize branching times  
  lineages2_N <- lineages2 / max(lineages2);  #normalize lineages  
  ltt2 <- approxfun(b_times2_N,lineages2_N,method="constant"); #method = constant ensures a step function

  f <- function(t,x,p) { #function f is the absolute difference in time t: 0 <= t < 1
       output <- ltt1(t) - ltt2(t);
	   output <- output * output;
       return(list(output));
  }
  
  times <- (0:100)/100; #evaluation points of the integration function below - more points ensures higher precision.
  int_1 <- lsoda(0,times,func=f,tcrit=c(1)); #integrate over t: 0 < t < 1, notice tcrit=0 indicating t should never be larger than 0.
  total_area <- int_1[length(times),2]
  
  return(total_area);
}


