\name{nLTTstat}
\alias{nLTTstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the difference between two normalized Lineage-Through-Time curves, given two phylogenetic trees.
}
\description{
  This function takes two ultrametric phylogenetic trees, calculates the normalized Lineage-Through-Time statistic for both trees and then calculates the difference between the two statistics.
}
\usage{
nLTTstat(tree1, tree2,distance_method = "abs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{tree1}{
an object of class \code{"phylo"};
}

  \item{tree2}{
an object of class \code{"phylo"};
}

  \item{distance_method}{
Chosen measurement of distance between the two nLTT curves, options are (case sensitive):\cr
	- "abs": use the absolute distance\cr
	- "squ": use the squared distance;\cr
  }
}
\value{
The difference between the two nLTT statistics
}
\author{
	Thijs Janzen
}

\examples{
	data(exampleTrees);
	nltt_plot(exampleTrees[[1]]);
	nnltt_lines(exampleTrees[[2]], lty=2);
	nLTTstat(exampleTrees[[1]], exampleTrees[[2]], distance_method = "abs");
}
