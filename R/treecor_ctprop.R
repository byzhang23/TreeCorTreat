#' Pipeline for TreeCorTreat using celltype proportion as features
#'
#' @param hierarchy_list        A hierarchy list by running \code{'extract_hrchy_string()'} or \code{'extract_hrchy_seurat()'} function. Contains four elements: 'edges', 'layout', 'immediate_children' and 'leaves_info'.
#' @param cell_meta             A data frame for cell-level metadata, with cells on the row. Must contain these columns: 'barcode', 'celltype' and 'sample'.
#' @param sample_meta           A data frame for sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables to be used in the analysis, such as covariates or outcomes of interest.
#' @param response_variable     A vector of response variables. Extract these columns from 'sample_meta'.
#' @param method                A character string indicating which approach is used to summarize features. One of 'aggregate' (default) or 'concat_leaf' or 'concat_immediate_children'.
#' @param formula               An object of class 'formula': a symbolic description of the model to be fitted, adjusting for confounders.
#' @param separate              A TRUE (default) or FALSE indicator, specifying how to evaluate multivariate outcomes.
#'                              \itemize{
#'                              \item TRUE: evaluate multivariate phenotype separately (it is equivalent to run this pipeline for each univariate phenotype).
#'                              \item FALSE: evaluate multivariate phenotype jointly.
#'                              }
#' @param analysis_type         Specify a summary statistic for analysis, must be one of
#'                              \itemize{
#'                              \item pearson (default): Pearson correlation (only for \code{method = 'aggregate'})
#'                              \item spearman: Spearman correlation (only for \code{method = 'aggregate'})
#'                              \item cancor: Canonical correlation
#'                              \item regression: Fit a linear regression with cell type proportion feature vectors as dependent variable and phenotype(s) as explanatory variables (full model). Compare between full model (with phenotypes in explanatory variables) and reduced model (without phenotypes).
#'                              }
#' @param num_cancor_components Number of canonical components to be extracted. Only works for multivariate phenotypes when specifying \code{'separate = F'} and \code{'analysis_type = cancor'}.
#' @param num_permutations      Number of permutations (by default: 1000).
#' @param alternative           A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param num_PCs               Number of PCs
#' @param verbose               Show progress
#' @return A list of two elements: the first element is a table of summary statistic for each tree node and second element is a list of PC matrices for each tree node.
#' @export
#' @import dplyr
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#'
#' @examples
#' # default setting
#' result <- treecor_ctprop(hierarchy_list, cell_meta, sample_meta, response_variable = 'severity')
#' # extract a table of summary statistic for each cell cluster
#' result[[1]] # or result$canonical_corr
#' # extract a list of PC matrices for each cell cluster
#' result[[2]] # or result$pc.ls
#' # extract PC matrix for celltype 'T'
#' result$pc.ls[['T']]

treecor_ctprop <- function(hierarchy_list, cell_meta, sample_meta, response_variable, method = 'aggregate', formula = NULL, separate = T, analysis_type = 'pearson',num_cancor_components = 1, num_permutations = 1000, alternative = 'two.sided', num_PCs = 2,verbose = T){

    ## Extract info
    label_info <- hierarchy_list$layout
    leaves_info <- hierarchy_list$leaves_info
    if(method == 'concat_immediate_children'){
        immediate_children <- hierarchy_list$immediate_children
    }else{
        immediate_children <- NULL
    }

    ## Check
    if(!separate & length(response_variable)>nrow(sample_meta)){
        # multivariate outcome when p > n; combine
        pc_phenotype <- prcomp(sample_meta[,response_variable],scale. = T)$x[,1:num_PCs]
        warning(paste0('Number of phenotypes > number of samples.\nProceed with PCA of response variables.\nTop ',num_PCs,' PCs are extracted and phenotypes are renamed as combined_1,...,combined_',num_PCs))
        for(i in 1:num_PCs){
            sample_meta[,paste0('combined_',i)] <- pc_phenotype[,i]
        }
        response_variable <- paste0('combined_',1:num_PCs)
    }
    if(min(length(response_variable),num_PCs) < num_cancor_components & !separate){
        warning(paste0('Number of canonical components exceeds ', min(length(response_variable),num_PCs),'. Set num_cancor_components = ',min(length(response_variable),num_PCs)))
        num_cancor_components <- min(length(response_variable),num_PCs)
    }
    if(method == 'aggregate' & num_cancor_components>1){
        warning('Number of canonical components > 1 when method = `aggregate`. . Set num_cancor_components = 1')
        num_cancor_components <- 1
    }
    if(length(response_variable)==1 & separate){
        separate <- F
    }
    if(sum(analysis_type %in% c('regression','cancor','pearson','spearman'))==0){
        stop('Invalid analysis type. Must be one of `regression`,`cancor`,`pearson`,`spearman`.')
    }



    ## =============================== ##
    ## 1. Summarize features
    ## =============================== ##
    if(verbose){
        message(paste0('=== Use ',method,' method to summarize celltype proportion features and Perform PCA ==='))
        if(!is.null(formula)){
            message(paste0('Adjusted for confounders: ',formula))
        }

        message('Response variable(s):\n',paste(response_variable,collapse = ' '))
        message(paste0('Evaluate response variable(s) ',ifelse(separate,'separately','jointly')))
    }
    if(!separate){
            if(num_cancor_components == 1){
                cancor.dat <- ctprop_cancor(cell_meta,
                                            sample_meta,
                                            response_variable,
                                            method,
                                            formula,
                                            analysis_type,
                                            leaves_info,
                                            immediate_children,
                                            num_permutations,
                                            alternative,
                                            num_PCs,
                                            component_id = 1,
                                            verbose)$canonical_corr
                colnames(cancor.dat)[-1] = paste0(paste(response_variable,collapse = '_'),'.',colnames(cancor.dat)[-1])
            }else{

                cancor.dat <- (Reduce(inner_join,lapply(1:num_cancor_components,function(c_id){
                    tmp <- ctprop_cancor(cell_meta,
                                         sample_meta,
                                         response_variable,
                                         method,
                                         formula,
                                         analysis_type,
                                         leaves_info,
                                         immediate_children,
                                         num_permutations,
                                         alternative,
                                         num_PCs,
                                         component_id = c_id,
                                         verbose)$canonical_corr
                    colnames(tmp)[-1] = paste0(paste(response_variable,collapse = '_'),'_Component',c_id,'.',colnames(tmp)[-1])
                    tmp
                })))
            }
        }else{
            ## evaluate response variables separately
            cancor.dat <- (Reduce(inner_join,lapply(response_variable,function(r){
                tmp <- ctprop_cancor(cell_meta,
                                     sample_meta,
                                     r,
                                     method,
                                     formula,
                                     analysis_type,
                                     leaves_info,
                                     immediate_children,
                                     num_permutations,
                                     alternative,
                                     num_PCs,
                                     component_id = 1,
                                     verbose)$canonical_corr
                colnames(tmp)[-1] = paste0(r,'.',colnames(tmp)[-1])
                tmp
            })))
        }

    pc.ls <- ctprop_cancor(cell_meta,
                           sample_meta,
                           response_variable[1],
                           method,
                           formula,
                           analysis_type,
                           leaves_info,
                           immediate_children,
                           num_permutations,
                           alternative,
                           num_PCs,
                           component_id = 1,
                           verbose)$pc_ls

    res <- suppressMessages(inner_join(cancor.dat,label_info))
    if(!is.null(pc.ls)) names(pc.ls) <- label_info$label[as.integer(names(pc.ls))]

    res.ls <- list(canonical_corr = res,
                   pc_ls = pc.ls)
}
