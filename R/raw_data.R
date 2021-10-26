#' Raw data for TreeCorTreat package
#'
#' Processed single-cell RNA-seq data downloaded from ArrayExpress: E-MTAB-9357. It includes 8 PMBC samples, 4 healthy donors and 4 severe COVID-19 patients.
#'
#' @docType data
#' @usage data(raw_data)
#' @keywords datasets
#' @format A list that contains three elements: sample-level metadata, cell-level metadata and gene expression matrix (raw count)
#' @references Su Y, Chen D, Yuan D, et al. Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19. Cell. 2020;183(6):1479-1495.e20. doi:10.1016/j.cell.2020.10.037.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/33171100/}{PubMed})
#'
#' @source \href{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9357/}{ArrayExpress}
#'
#' @examples
#' data(raw_data)
"raw_data"