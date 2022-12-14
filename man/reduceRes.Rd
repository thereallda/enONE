% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{reduceRes}
\alias{reduceRes}
\title{Combine list of DE results}
\usage{
reduceRes(res.ls, logfc.col, levels = names(res.ls))
}
\arguments{
\item{res.ls}{Named list of differential analysis results tables.
Each elements in the list correspond to a table of differential analysis
results between two groups of samples.}

\item{logfc.col}{Column name of the log fold-change.}

\item{levels}{Factor levels of the groups, default order by the element order of \code{res.ls}.}
}
\value{
data.frame
}
\description{
Combine list of DE results
}
