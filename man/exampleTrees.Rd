\name{exampleTrees}
\alias{exampleTrees}
\docType{data}
\title{ example trees to test the functionality of the package}
\description{100 phylogenetic trees of class \code{phylo}, generated using the sim.globalBiDe.age function from the \code{TESS} package, with lambda = 0.3, mu = 0.1, age = 10.}
\usage{data(exampleTrees)}
\format{  
	A list containing objects of class \code{phylo}.
}

\examples{
 data(exampleTrees);
 obs <- exampleTrees[[1]];
 nLTT.plot(obs);
}

\keyword{datasets}
