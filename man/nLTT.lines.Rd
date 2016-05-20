\name{nLTT.lines}
\alias{nLTT.lines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Normalized version of the ape function ltt.lines.
}
\description{
This is a modified version of the \code{ape} function ltt.lines: add the normalized Lineage-Through-Time statistic of a phylogenetic tree to an already existing plot
}
\usage{
nLTT.lines(phy, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
	an object of class \code{"phylo"};

}
  \item{\dots}{
further graphical arguments that can be passed to \code{lines()}
}
}
\author{
  Thijs Janzen
}

\examples{
	data(exampleTrees);
	nLTT.plot(exampleTrees[[1]]);
	nLTT.lines(exampleTrees[[2]],lty=2);
}