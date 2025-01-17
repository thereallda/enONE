% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{Enone-class}
\alias{Enone-class}
\alias{createEnone}
\title{Enone object and constructor}
\usage{
createEnone(
  data,
  col.data,
  spike.in.prefix = NULL,
  input.id = "Input",
  enrich.id = "Enrich",
  synthetic.id = NULL
)
}
\arguments{
\item{data}{A un-normalized count data matrix of shape n x p, where n is the
number of samples and p is the number of features.}

\item{col.data}{\code{data.frame} with at least two columns (indicate condition and enrich groups).
Rows of \code{col.data} correspond to columns of \code{data}.}

\item{spike.in.prefix}{A character specify the prefix of spike-in id, e.g., "^FB" stands for fly spike-in id, default: NULL.}

\item{input.id}{Input library id, must be consistent with the enrich column of \code{col.data}, e.g., "Input".}

\item{enrich.id}{Enrich library id, must be consistent with the enrich column of \code{col.data}, e.g., "Enrich".}

\item{synthetic.id}{Vector of synthetic RNA id, e.g. c("Syn1","Syn2"), default: NULL.}
}
\value{
Enone object
}
\description{
\code{Enone} object extends the \code{SummarizedExperiment} class.
The \code{createEnone} is a easy constructor of \code{Enone} object
}
\details{
Description of each slot: \cr
\code{assay} \code{SummarizedExperiment::Assays} object, contains all counts. \cr
\code{counts} list for holding raw/normalized counts of sample and spike_in. \cr
\code{enone_factor} list for holding normalization factors of sample and spike_in. \cr
\code{enone_metrics} data frame with normalization methods in row and metrics in columns. \cr
\code{enone_score} data frame with normalization methods in row and scores in columns. \cr
\code{enrichment} list for holding all differential analysis results. \cr
\code{enrichment_filtered} lists for holding filtered differential analysis results. \cr
\code{parameter} list of parameters.\cr
}
