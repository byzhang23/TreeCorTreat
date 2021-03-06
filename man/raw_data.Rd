% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raw_data.R
\docType{data}
\name{raw_data}
\alias{raw_data}
\title{Raw data for TreeCorTreat package}
\format{
A list that contains three elements: sample-level metadata, cell-level metadata and gene expression matrix (raw count)
}
\source{
\href{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9357/}{ArrayExpress}
}
\usage{
data(raw_data)
}
\description{
Processed single-cell RNA-seq data downloaded from ArrayExpress: E-MTAB-9357. It includes 8 PMBC samples, 4 healthy donors and 4 severe COVID-19 patients.
}
\examples{
data(raw_data)
}
\references{
Su Y, Chen D, Yuan D, et al. Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19. Cell. 2020;183(6):1479-1495.e20. doi:10.1016/j.cell.2020.10.037.
(\href{https://pubmed.ncbi.nlm.nih.gov/33171100/}{PubMed})
}
\keyword{datasets}
