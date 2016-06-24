\name{abc_smc_nltt}
\alias{abc_smc_nltt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to perform Approximate Bayesian Computation within an Sequential Markov Chain (ABC-SMC), for diversification analysis of phylogenetic trees.

}
\description{
This function performs ABC-SMC as described in Toni 2009 for given diversification model, provided a phylogenetic tree. ABC-SMC is not limited to only using the normalized LTT as statistic.
}
\usage{
abc_smc_nltt(tree, statistics, simulation_function, init_epsilon_values,
		     prior_generating_function, prior_density_function,
		     number_of_particles = 1000, sigma = 0.05, stop_rate = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
	an object of class \code{"phylo"}; the tree upon which we want to fit our diversification model
}
  \item{statistics}{
	A vector containing functions that take a tree as an argument and return a single scalar value (the statistic).
}
  \item{simulation_function}{
	A function that implements the diversification model and returns an object of class \code{"phylo"}.
}
  \item{init_epsilon_values}{
	A vector containing the initial threshold values for the summary statistics from the vector \code{statistics}.
}
  \item{prior_generating_function}{
	Function to generate parameters from the prior distribution of these parameters (e.g. a function returning lambda and mu in case of the birth-death model)
}
  \item{prior_density_function}{
	Function to calculate the prior probability of a set of parameters.
}
  \item{number_of_particles}{
	Number of particles to be used per iteration of the ABC-SMC algorithm.
}
  \item{sigma}{
	Standard deviation of the perturbance distribution (perturbance distribution is a gaussian with mean 0).
}
  \item{stop_rate}{
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
 	return ( rexp(n=2, rate=0.1) );
 }

 prior_dens <- function(val) {
 	return ( dexp( val[1], rate = 0.1) * dexp( val[2], rate = 0.1) );
 }

 require(TESS);
 treeSim <- function(params) {
 	t <- TESS.sim.age(n=1, lambda = params[1], mu = params[2], age = 10)[[1]];
	return (t);
 }

 obs <- treeSim(c(0.5,0.1));

 statWrapper <- function(tree1) {
	return ( nLTTstat_exact(tree1, obs, "abs"));
 }

 stats <- c(statWrapper);

 results <-  abc.smc.nltt(obs, stats, treeSim, init_epsilon_values = 0.2,
		     prior_generating_function = prior_gen,
		     prior_density_function = prior_dens,
  		     number_of_particles = 1000, sigma = 0.05, stop_rate = 1e-5);

}


}
