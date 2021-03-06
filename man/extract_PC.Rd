% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_PC.R
\name{extract_PC}
\alias{extract_PC}
\title{Extract top n PCs}
\usage{
extract_PC(data, num_PCs = 2)
}
\arguments{
\item{data}{A data frame to conduct PCA}

\item{num_PCs}{Number of PCs to be extracted}
}
\value{
A PCA matrix with \code{num_PCs} columns
}
\description{
This function performs PCA on a given matrix and returns a PCA matrix.
}
\author{
Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
}
