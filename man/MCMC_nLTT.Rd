\name{mcmc_nltt}
\alias{mcmc_nltt}
\title{
  Code to perform Metropolis-Hastings MCMC for a diversification model, given a phylogenetic tree. This function was used in the MEE paper to calculate the likelihood reference estimates.
}
\description{
  This function performs Metropolis-Hastings MCMC, where the user provides a likelihood function and a phylogenetic tree.
}
\usage{
mcmc_nltt(phy, likelihood_function, parameters, logtransforms,
		  iterations, burnin = round(iterations / 3), thinning = 1, sigma=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
	an object of class \code{"phylo"}; the tree upon which we want to fit our diversification model
  }
  \item{likelihood_function}{
	Function that calculates the likelihood of our diversification model, given the tree. Function should me of the format function(parameters,phy).
}
  \item{parameters}{
	Initial parameters to start the chain.
}
  \item{logtransforms}{
	Whether to perform jumps on logtransformed parameters (TRUE) or not (FALSE)
}
  \item{iterations}{
	Length of the chain
}
  \item{burnin}{
	Length of the burnin, default is 30% of \item{iterations}
}
  \item{thinning}{
	Size of thinning, default = 1
}
  \item{sigma}{
	Standard deviation of the jumping distribution, which is N(0, sigma).
}
}
\value{
  An \code{MCMC} object, as used by the package \code{"coda"}.
}

\author{
	Sebastian Hoehna & Thijs Janzen
}
\examples{

## MCMC examples are typically very slow ####
\dontrun{

require(TESS);

obs <- TESS.sim.age(n = 1, lambda = 0.5, mu = 0.1, age = 10)[[1]];

LL_BD <- function(params, phy) {
 lnl <- tess.likelihood(phy, lambda = params[1], mu = params[2],
 							  samplingProbability = 1, log = TRUE);
 prior1 <- dunif( params[1], 0, 100, log = TRUE)
 prior2 <- dunif( params[2], 0, 100, log = TRUE);
 return(lnl + prior1 + prior2);
}

require(coda);

mcmc_out <- mcmc_nltt(obs, LL_BD, c(0.5, 0.1), c(TRUE, TRUE),
			iterations = 1000, burnin = 100, thinning = 10, sigma = 1)
plot(mcmc_out);
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MCMC}
\keyword{ Likelihood}% __ONLY ONE__ keyword per line
