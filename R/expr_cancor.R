#' Compute summary statistic based on gene expression features
#'
#' This function takes in a list of cluster-level PCA matrices as input and generates canonical correlation and p-value based on permutation test.
#' @param pca.ls            A list of cell cluster-specific PCA matrices
#' @param sample_meta       Sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables such as covariates or outcome of interest.
#' @param response_variable A vector of response variables.
#' @param analysis_type     Either 'cancor' (canonical correlation) or 'regression' (F-statistic).
#' @param num_permutations  Number of permutations (by default: 1000).
#' @param alternative       A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param componen_id       An index to extract canonical correlation component (by default: 1).
#' @param ncores            Number of cores
#' @return A data frame for summary statistic (canonical correlation or F-stat), p-value and adjusted p-value for each node.
#' @export
#' @import dplyr
#' @importFrom parallel mclapply
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji

expr_cancor <- function(pca.ls, sample_meta, response_variable, analysis_type, num_permutations, alternative, component_id = 1, ncores){

    summary <- do.call(rbind,mclapply(pca.ls,function(pc){
        if(is.null(pc)){
            cancor_p <- data.frame(cancor = NA,p = NA)
        }else{
            tmp <- suppressMessages(inner_join(data.frame(sample = rownames(pc)),sample_meta))
            response <- tmp[,response_variable,drop = F]
            cancor_p <- cancor_permute(pc,
                           response,
                           method = analysis_type,
                           num_permutations,
                           alternative,
                           component_id)
        }
        cancor_p
    },mc.cores = ncores))

    summary$adjp <- p.adjust(summary$p,method = 'BY') # BY procedure - dependent test
    summary$id <- names(pca.ls) %>% as.numeric
    res <- summary %>% select(id,cancor,p,adjp) %>% mutate(direction = ifelse(cancor>0,'+','-'),
                                                          p.sign = ifelse(p<0.05,'sig','ns'),
                                                          adjp.sign = ifelse(adjp<0.05,'sig','ns'))

    if(length(setdiff(unique(res$direction),NA))>1){
        res <- res %>% mutate(cancor = abs(cancor))
    }
    if(analysis_type == 'regression_concat'){
        colnames(res) <- gsub('cancor','Fstat',colnames(res))
    }
    return(res)
}
