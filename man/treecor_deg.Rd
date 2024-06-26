% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treecor_deg.R
\name{treecor_deg}
\alias{treecor_deg}
\title{Differential gene expression (DGE) analysis}
\usage{
treecor_deg(
  expr,
  hierarchy_list,
  cell_meta,
  sample_meta,
  response_variable,
  separate = T,
  weight = NULL,
  formula = NULL,
  coef = 2,
  fdr_cutoff = 0.05,
  filter_prop = 0.1,
  pseudobulk_list = NULL,
  ncores = parallel::detectCores(),
  save_as_csv = T,
  verbose = T
)
}
\arguments{
\item{expr}{A raw count gene expression matrix with genes on rows and cells on columns. Note that cell barcode shall use ':' to separate sample name and barcode (i.e. "sample:barcode")}

\item{hierarchy_list}{A hierarchy list by running \code{'extract_hrchy_text()'} or \code{'extract_hrchy_seurat()'} function. Contains 'edges', 'layout', 'immediate_children' and 'leaves_info' as 4 elements.}

\item{cell_meta}{Cell-level metadata, where each row is a cell. Must contain these columns: barcode, celltype and sample.}

\item{sample_meta}{Sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables to be used in the analysis, such as covariates or outcomes of interest.}

\item{response_variable}{A vector of response variables.}

\item{separate}{A TRUE (default) or FALSE indicator, specifying how to evaluate multivariate outcomes.
\itemize{
\item TRUE: evaluate multivariate phenotype separately (it is equivalent to run this pipeline for each univariate phenotype).
\item FALSE: evaluate multivariate phenotype jointly.
}}

\item{weight}{A weight matrix to combine multivariate phenotype. The dimension should be number_phenotype * 1 If none is provided, then PC1 will be used as a joint univariate phenotype.}

\item{formula}{An object of class 'formula': a symbolic description of adjustment formula (i.e. only includes covariates other than response variable(s))}

\item{coef}{A column number or column name specifying which coefficient to be extracted (by default: 2).}

\item{fdr_cutoff}{Cutoff value for FDR. Only genes with lower FDR are listed. Default is 0.05.}

\item{filter_prop}{A number ranges from 0 to 1, to filter low expressed genes across samples (by default: 0.1). Genes with at least this proportion of samples with log2-normalized count greater than 0.01 are retained.}

\item{pseudobulk_list}{A list of sample-level (adjusted) pseudobulk for each node. Default is NULL. Users can provide their processed pseudobulk list (e.g. after covariate adjustment) via this parameter. Note that the names of list shall be matched with \code{`id`} extracted from \code{`hierarchy_list`}.}

\item{ncores}{Number of cores to be used. If ncores > 1, it will be implemented in a parallel mode.}

\item{save_as_csv}{An indicator to save identified DEGs in csv files. DEGs for each cell cluster is saved as 'responsevariable_celltype_DEG.csv' and a summary file of combining DEGs from all cell clusters is saved as 'responsevariable_combinedDEG.csv'.}

\item{verbose}{Show progress}
}
\value{
A list of three elements:
\itemize{
\item dge.summary: A summary table of number of DEGs for each tree node.
\item dge.ls: A comprehensive list of outcome(s)-associated DEGs for each tree node. Use \code{`result$dge.ls$response_variable[[celltype]]`} to extract DEGs for a specific cell type
\item pseudobulk.ls: A list of sample-level pseudobulk gene expression matrix for each cell cluster. Use \code{`result$pseudobulk.ls[[celltype]]`} to extract.
}
}
\description{
Constructs sample-level normalized pseudobulk gene expression matrix and uses LIMMA to identify differentially expressed genes (DEGs) for each cell cluster.
}
\examples{
# default setting
result <- treecor_deg(expr,hierarchy_list, cell_meta, sample_meta, response_variable = 'severity')
# obtain summary table
result$dge.summary
# extract DEGs of severity in all cell types
result$dge.ls$severity
# extract DEGs for celltype 'T'
result$dge.ls$severity[['T']]
# extract sample-level pseudobulk for all cell clusters
result$pseudobulk.ls
# extract sample-level pseudobulk for celltype 'T'
result$pseudobulk.ls[['T']]
}
\author{
Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
}
