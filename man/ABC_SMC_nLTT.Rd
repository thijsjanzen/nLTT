\name{ABC_SMC_nLTT}
\alias{ABC_SMC_nLTT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to perform Approximate Bayesian Computation within an Sequential Markov Chain (ABC-SMC), for diversification analysis of phylogenetic trees.

}
\description{
This function performs ABC-SMC as described in Toni 2009 for given diversification model, provided a phylogenetic tree. ABC-SMC is not limited to only using the normalized LTT as statistic.
}
\usage{
ABC_SMC_nLTT(tree, statistics, simFunc, initEpsilon, 
		     PRIOR_GEN, PRIOR_DENS, numParticles, sigma = 0.05, stopRate = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
	an object of class \code{"phylo"}; the tree upon which we want to fit our diversification model
}
  \item{statistics}{
	A vector containing functions that take a tree as an argument and return a single scalar value (the statistic).
}
  \item{simFunc}{
	A function that implements the diversification model and returns an object of class \code{"phylo"}.
}
  \item{initEpsilon}{
	A vector containing the initial threshold values for the summary statistics from the vector \code{statistics}.
}
  \item{PRIOR_GEN}{
	Function to generate parameters from the prior distribution of these parameters (e.g. a function returning lambda and mu in case of the birth-death model)
}
  \item{PRIOR_DENS}{
	Function to calculate the prior probability of a set of parameters.
}
  \item{numParticles}{
	Number of particles to be used per iteration of the ABC-SMC algorithm.
}
  \item{sigma}{
	Standard deviation of the perturbance distribution (perturbance distribution is a gaussian with mean 0).
}
  \item{stopRate}{
	If the acceptance rate drops below \code{stopRate}, stop the ABC-SMC algorithm  and assume convergence.
}
}
\value{
   A matrix with \code{n} columns, where \code{n} is the number of parameters you are trying to estimate.
}
\references{
 Toni, T., Welch, D., Strelkowa, N., Ipsen, A., & Stumpf, M.P.H. (2009). Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems. Journal of the Royal Society Interface, 6(31), 187-202.
}
\author{
	Thijs Janzen
}
\examples{

############# Code to run the ABC-SMC algorithm would involve a lot 
############# of dependencies on other packages, and be computationally demanding
############# here we sketch out an example, but do not execute it:

\dontrun{


 prior_gen <- function() {
 	return( rexp(n=2,rate=0.1) );
 }

 prior_dens <- function(val) {
 	return( dexp(val[1],rate=0.1) * dexp(val[2],rate=0.1));
 }

 require(TESS);
 treeSim <- function(params) {
 	t <- sim.globalBiDe.age(n=1,lambda=params[1],mu=params[2],age=10)[[1]];
	return(t);
 }

 obs <- treeSim(c(0.5,0.1));
 
 statWrapper <- function(tree1) {
	return(nLTT(tree1,obs,"abs"));
 }
	
 stats <- c(statWrapper);
 
 results <- ABC_SMC_nLTT(obs,stats,treeSim,initEpsilon = 0.2,
		     PRIOR_GEN = prior_gen,PRIOR_DENS <- prior_dens,
  		     numParticles = 1000, sigma = 0.05, stopRate = 1e-5);

} 


}
