#' Extract top n PCs
#'
#' This function performs PCA on a given matrix and returns a PCA matrix.
#' @param data A data frame to conduct PCA
#' @param num_PCs Number of PCs to be extracted
#' @return A PCA matrix with \code{num_PCs} columns
#' @export
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji

extract_PC <- function(data,num_PCs = 2){

    if(!is.null(data)){
        dat_pca <- prcomp(t(data),scale. = T)$x
        res <- dat_pca[,1:min(num_PCs,ncol(dat_pca))]
    }else{
        res <- NULL
    }
    res
}
