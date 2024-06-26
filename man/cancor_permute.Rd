% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cancor_permute.R
\name{cancor_permute}
\alias{cancor_permute}
\title{Calculate summary statistic and perform permutation test}
\usage{
cancor_permute(x, y, method, num_permutations, alternative, component_id = 1)
}
\arguments{
\item{x}{A vector or data frame.}

\item{y}{A vector or data frame.}

\item{method}{Specify analysis_type: cancor, pearson, spearman, regression}

\item{num_permutations}{Number of permutations (by default: 100).}

\item{alternative}{A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}

\item{componen_id}{An index to extract canonical correlation component (by default: 1).}
}
\value{
A data frame of summary statistic and p-value
}
\description{
Compute summary statistic (e.g. canonical correlation, F-statistic, likelihood ratio test statistic) and obtain p-value based on permutation test.
}
\author{
Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
}
