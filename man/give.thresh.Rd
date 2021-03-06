% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/give.thresh.R
\name{give.thresh}
\alias{give.thresh}
\title{Compute the maximum distance of the minimum spanning tree based on a distance
matrix}
\usage{
give.thresh(matdist)
}
\arguments{
\item{matdist}{A distance matrix (class \code{dist} or \code{matrix})}
}
\value{
The maximum distance in the minimum spanning tree.
}
\description{
It is used to select a truncation value for the dbMEM approach. It
returns the minimum value that keep all samples connected.
}
\examples{
xy <- matrix(rnorm(60),30,2)
dxy <- dist(xy)
th <- give.thresh(dxy)
}
\author{
Stéphane Dray \email{stephane.dray@univ-lyon1.fr}
}
