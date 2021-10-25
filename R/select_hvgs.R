#' Select highly variable genes (HVGs)
#'
#' This function takes in normalized sample-level pseudobulk as input and identifies highly variable genes.
#' @param expr      Sample-level pseudobulk gene expression matrix
#' @param num_hvgs  Number of highly variable genes extracted from sample-level pseudobulk. Highly variable genes (HVGs) are defined as genes with positive residuals by fitting a gene-specific lowess (\code{loessFit}) function of standard deviation against its mean. The HVGs are ordered by residuals in a descending order. Default is NULL, which includes all genes with positive residuals. Can specify a reasonable positive number (e.g. 1000).
#' @return A vector of highly variable gene names
#' @importFrom limma loessFit
#' @export
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji

select_hvgs <- function(expr,
                        num_hvgs = NULL){

    cm <- rowMeans(expr)
    csd <- sqrt((rowMeans(expr*expr) - cm^2) / (ncol(expr) - 1) * ncol(expr))
    mod <- limma::loessFit(csd,cm)
    residual <- resid(mod)
    names(residual) <- rownames(expr)
    residual <- sort(residual,decreasing = T)

    if(is.null(num_hvgs)){
        hvgs <- names(residual)[which(residual > 0)]
    }else{
        hvgs <- names(residual)[1:num_hvgs]
    }

    return(hvgs)
}
