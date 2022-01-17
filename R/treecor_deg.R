#' Differential gene expression (DGE) analysis
#'
#' Constructs sample-level normalized pseudobulk gene expression matrix and uses LIMMA to identify differentially expressed genes (DEGs) for each cell cluster.
#' @param expr              A raw count gene expression matrix with genes on rows and cells on columns. Note that cell barcode shall use ':' to separate sample name and barcode (i.e. "sample:barcode")
#' @param hierarchy_list    A hierarchy list by running \code{'extract_hrchy_text()'} or \code{'extract_hrchy_seurat()'} function. Contains 'edges', 'layout', 'immediate_children' and 'leaves_info' as 4 elements.
#' @param cell_meta         Cell-level metadata, where each row is a cell. Must contain these columns: barcode, celltype and sample.
#' @param sample_meta       Sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables to be used in the analysis, such as covariates or outcomes of interest.
#' @param response_variable A vector of response variables.
#' @param separate          A TRUE (default) or FALSE indicator, specifying how to evaluate multivariate outcomes.
#'                              \itemize{
#'                              \item TRUE: evaluate multivariate phenotype separately (it is equivalent to run this pipeline for each univariate phenotype).
#'                              \item FALSE: evaluate multivariate phenotype jointly.
#'                              }
#' @param weight            A weight matrix to combine multivariate phenotype. The dimension should be number_phenotype * 1 If none is provided, then PC1 will be used as a joint univariate phenotype.
#' @param formula           An object of class 'formula': a symbolic description of adjustment formula (i.e. only includes covariates other than response variable(s))
#' @param coef              A column number or column name specifying which coefficient to be extracted (by default: 2).
#' @param fdr_cutoff        Cutoff value for FDR. Only genes with lower FDR are listed. Default is 0.05.
#' @param filter_prop       A number ranges from 0 to 1, to filter low expressed genes across samples (by default: 0.1). Genes with at least this proportion of samples with log2-normalized count greater than 0.01 are retained.
#' @param pseudobulk_list   A list of sample-level (adjusted) pseudobulk for each node. Default is NULL. Users can provide their processed pseudobulk list (e.g. after covariate adjustment) via this parameter. Note that the names of list shall be matched with \code{`id`} extracted from \code{`hierarchy_list`}.
#' @param ncores            Number of cores to be used. If ncores > 1, it will be implemented in a parallel mode.
#' @param save_as_csv       An indicator to save identified DEGs in csv files. DEGs for each cell cluster is saved as 'responsevariable_celltype_DEG.csv' and a summary file of combining DEGs from all cell clusters is saved as 'responsevariable_combinedDEG.csv'.
#' @param verbose           Show progress
#' @return A list of three elements:
#' \itemize{
#' \item dge.summary: A summary table of number of DEGs for each tree node.
#' \item dge.ls: A comprehensive list of outcome(s)-associated DEGs for each tree node. Use \code{`result$dge.ls$response_variable[[celltype]]`} to extract DEGs for a specific cell type
#' \item pseudobulk.ls: A list of sample-level pseudobulk gene expression matrix for each cell cluster. Use \code{`result$pseudobulk.ls[[celltype]]`} to extract.
#' }
#' @export
#' @import dplyr parallel limma
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#'
#' @examples
#' # default setting
#' result <- treecor_deg(expr,hierarchy_list, cell_meta, sample_meta, response_variable = 'severity')
#' # obtain summary table
#' result$dge.summary
#' # extract DEGs of severity in all cell types
#' result$dge.ls$severity
#' # extract DEGs for celltype 'T'
#' result$dge.ls$severity[['T']]
#' # extract sample-level pseudobulk for all cell clusters
#' result$pseudobulk.ls
#' # extract sample-level pseudobulk for celltype 'T'
#' result$pseudobulk.ls[['T']]

treecor_deg <- function(expr,hierarchy_list, cell_meta, sample_meta, response_variable, separate = T, weight = NULL, formula = NULL, coef = 2, fdr_cutoff = 0.05,
                       filter_prop = 0.1, pseudobulk_list = NULL,ncores = parallel::detectCores(), save_as_csv = T, verbose = T){

    if (Sys.info()[['sysname']]=='Windows') {
        message('Parallel is disabled for Windows. Running with one core')
        ncores <- 1
    }

    ## Extract info
    label_info <- hierarchy_list$layout
    if(sum(duplicated(label_info$label))) stop('Duplicated cell type name.')
    leaves_info <- hierarchy_list$leaves_info

    if(!is.null(formula)){
        formula = gsub('\\~|\\ ','',formula)
    }

    ## =============================== ##
    ## 1. Pseudobulk
    ## =============================== ##
    if(verbose){
        message('=== Construct sample-level pseudobulk ===')
    }
    unq_id <- unique(leaves_info$id)

    if(is.null(pseudobulk_list)){
        pb.ls <- mclapply(unq_id,function(tid){
            if(verbose) message(paste0('node ',tid))
            node_info <- leaves_info %>% filter(id==tid)
            sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
            pb <- construct_pseudobulk(expr,
                                       sub_meta,
                                       filter_prop = filter_prop)
        },mc.cores = ncores)
        names(pb.ls) <- unq_id
    }else{
        # input pseudobulk list (user-specify)
        message('Input pseudobulk list:')
        pb_names <- intersect(as.character(unq_id),names(pseudobulk_list))
        if(length(pb_names) != length(unq_id)){
            stop('Invalid pseudobulk list, please check if the following nodes exist: \n', paste(setdiff(as.character(unq_id),pb_names),collapse = ','))
        }
        pb.ls <- pseudobulk_list
    }

    if(!separate & length(response_variable)>1){
        # multivariate outcome; combine
        if(!is.null(weight)){
            sample_meta$combined_phenotype <- as.matrix(sample_meta[,response_variable]) %*% weight
        }else{
            # default: PC1
            sample_meta$combined_phenotype <- prcomp(sample_meta[,response_variable],scale. = T)$x[,1]
            warning('Multiple response variables are provided. Proceed with 1st PC of response variables and rename as `combined_phenotype`')
        }
        response_variable <- 'combined_phenotype'
    }

    ## =============================== ##
    ## 2. LIMMA
    ## =============================== ##
    if(verbose){
        message('=== Running LIMMA ===')
    }
    resp.ls <- lapply(response_variable,function(r){

        if(verbose) message(paste0('## Response variable: ',r))

        dge.ls <- mclapply(unq_id,function(tid){

            if(verbose) message(paste0('node ',tid))
            pb <- pb.ls[[as.character(tid)]]
            if(ncol(pb)>0){
                tmp_meta <- suppressMessages(inner_join(data.frame(sample = colnames(pb)),sample_meta))

                if(nrow(tmp_meta)>1){
                    # drop factor level
                    factor.id <- which(sapply(1:ncol(tmp_meta),function(x) is.factor(tmp_meta[,x])))
                    for(k in factor.id){
                        tmp_meta[,k] <- droplevels(tmp_meta[,k])
                    }

                    # design matrix
                    design <- try(model.matrix(as.formula(paste(c(paste0('~',r),formula),collapse = '+')),data = tmp_meta))
                    if(verbose) message(paste0('Checking model: ',paste(c(paste0('~',r),formula),collapse = '+')))

                    # check collinearity/error in building design matrix
                    if(class(design)[1] == 'try-error'|qr(design)$rank<ncol(design)){
                        warning(paste0('Error occurs, proceed with dropping the covariates: `formula = ~',r,'`'))
                        design <- model.matrix(as.formula(paste0('~',r)),data = tmp_meta)
                        if(verbose & !(qr(design)$rank<ncol(design))) message(paste0('Final linear model fit: ~',r))
                    }
                    if(verbose){
                        message('Column names of design matrix:\n', paste(colnames(design),collapse  = ' '))
                        if(is.numeric(coef)){
                            if(coef <= ncol(design)){
                                message('Extract coefficient by column name: ', colnames(design)[coef])
                            }else{
                                stop(paste0('Incorrect column id specifying which coefficient to be extracted. Must be in the range of 1-',ncol(design)))
                            }
                        }else{
                            if(sum(colnames(design)==coef)>0){
                                message('Extract coefficient by column name: ', coef)
                            }else{
                                stop(paste0('Incorrect column name specifying which coefficient to be extracted. Must be one of :',paste(colnames(design),collapse  = ' ')))
                            }
                        }
                    }
                    # limma
                    if(qr(design)$rank<ncol(design)){
                        if(verbose) warning('Check model again. Proceed with num_DEG = 0.')
                        diff <- NULL
                    }else{
                        fit <- lmFit(pb,design)
                        eb <- eBayes(fit)
                        diff <- topTable(eb,coef = coef,n=nrow(pb),sort.by = "P", p.value = fdr_cutoff)
                        if(nrow(diff)>0){
                            diff <- diff[,c('logFC','t','P.Value','adj.P.Val')]
                            if(is.numeric(coef)){
                                colnames(diff) <- paste0(colnames(design)[coef],'.',c('logFC','t','p','fdr'))
                            }else{
                                colnames(diff) <- paste0(coef,'.',c('logFC','t','p','fdr'))
                            }
                        }else{
                            diff <- NULL
                        }
                    }

                }else{
                    diff <- NULL
                }
            }else{
                diff <- NULL
            }
            diff
        },mc.cores = ncores)
        names(dge.ls) <- label_info$label[unq_id]
        dge.ls

    })
    names(resp.ls) <- response_variable
    names(pb.ls) <- label_info$label[unq_id]

    ## =============================== ##
    ## 3. Final result
    ## =============================== ##
    if(save_as_csv){
        for(resp in response_variable){
            rdt <- resp.ls[[resp]]
            for(j in names(rdt)){
                rjdt <- rdt[[j]]
                if(!is.null(rjdt)){
                    write.csv(rjdt,paste0(resp,'_',j,'_DEG.csv'),row.names = T)
                }
            }
            # all in 1 csv
            rdt.comb <- do.call(rbind,lapply(1:length(rdt),function(j){
                dt <- rdt[[j]]
                if(!is.null(dt)){
                    dt$Gene <- rownames(dt)
                    dt$celltype <- names(rdt)[j]
                }
                dt
            }))
            if(nrow(rdt.comb)>0) write.csv(rdt.comb,paste0(resp,'_combinedDEG.csv'),row.names = T)
        }
    }

    if(verbose){
        message('=== Summarize final result ===')
    }
    dge.dat <- suppressMessages(Reduce(inner_join,lapply(1:length(resp.ls),function(i){

        dt <- resp.ls[[i]]
        df <- sapply(dt,function(d) ifelse(is.null(d),0,nrow(d))) %>% unlist
        setNames(data.frame(label = names(df),df),
                 c('label',paste0(names(resp.ls)[i],'.num_deg')))


    })))
    res <- suppressMessages(inner_join(dge.dat,label_info))
    final.ls <- list(dge.summary = res,
                     dge.ls = resp.ls,
                     pseudobulk.ls = pb.ls)
}
