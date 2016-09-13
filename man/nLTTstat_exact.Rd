\name{nLTTstat_exact}
\alias{nLTTstat_exact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the exact difference between two normalized Lineage-Through-Time curves, given two phylogenetic trees.
}
\description{
  This function takes two ultrametric phylogenetic trees, calculates the normalized Lineage-Through-Time statistic for both trees and then calculates the exact difference between the two statistics. Whereas the function \code{nLTTstat} uses an approximation to calculate the difference (which is faster for large trees), the function \code{nLTTstat_exact} calculates the exact difference, and should generally be preferred. Although the estimates are highly similar, \code{nLTTstat_exact} tends to return slightly higher values.
}
\usage{
nLTTstat_exact(tree1, tree2, distance_method = "abs", ignore_stem = TRUE)
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
  \item{ignore_stem}{
    a boolean whether to ignore the stem length
  }

}
\value{
The exact difference between the two nLTT statistics
}
\author{
	Thijs Janzen
}

\examples{
	data(exampleTrees);
	nltt_plot(exampleTrees[[1]]);
	nltt_lines(exampleTrees[[2]], lty=2);
	nLTTstat_exact(exampleTrees[[1]], exampleTrees[[2]], distance_method = "abs");
}
