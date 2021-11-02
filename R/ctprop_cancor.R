#' Core function for running treecor_ctprop() pipeline.
#'
#' This function calculates celltype proportion for each sample, performs Principal Component Analysis (PCA), calculates correlation/summary statistic based on top K chosen PCs and obtains permutation based p-value.
#' @param cell_meta             A data frame for cell-level metadata, where each row is a cell. Must contain these columns: barcode, celltype and sample.
#' @param sample_meta           A data frame for sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables such as covariates or outcome of interest.
#' @param response_variable     A vector of response variables. Will be used to calculate canonical correlation.
#' @param method                A character string indicating which approach is used to summarize features. One of 'concat_leaf (default)' or 'concat_immediate_children' or 'aggregate'.
#' @param formula               An object of class 'formula': a symbolic description of the model to be fitted, adjusting for confounders.
#' @param analysis_type         \itemize{
#' \item pearson: Pearson correlation (default; only works when method = 'aggregate')
#' \item spearman: Spearman correlation (only works when method = 'aggregate')
#' \item cancor: Canonical correlation
#' \item regression: Regression framework
#' }
#' @param leaves_info           A data frame that encodes leaf children for each node. Can be extracted from 'hierarchy_list' by running \code{'extract_hrchy_string()'} or \code{'extract_hrchy_seurat()'} function
#' @param immediate_children    A data frame that contains immediate children for each node. Can be extracted from 'hierarchy_list' by running \code{'extract_hrchy_string()'} or \code{'extract_hrchy_seurat()'} function
#' @param component_id          An index to extract canonical correlation component (by default: 1).
#' @param num_permutations      Number of permutations (by default: 1000).
#' @param alternative           Specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param num_PCs               Number of PCs
#' @param verbose               Show progress
#' @return A list of two elements: the first element is a table of summary statistic for each tree node and second element is a list of PC matrices for each tree node.
#' @export
#' @import dplyr
#' @importFrom lmtest lrtest
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji


ctprop_cancor <- function(cell_meta, sample_meta, response_variable, method,formula = NULL, analysis_type = 'pearson',leaves_info, immediate_children, component_id = 1, num_permutations = 100, alternative = 'two.sided', num_PCs = 2,verbose = T){

    count_ct_sample <- as.matrix(table(cell_meta$celltype,cell_meta$sample))
    ctprop <- t(t(count_ct_sample)/colSums(count_ct_sample))

    unq_leaves_info <- leaves_info %>% filter(leaf) %>% select(id,label) %>% unique
    miss_label <- setdiff(unq_leaves_info$label,rownames(ctprop))
    unq_id <- unique(leaves_info$id)

    if(length(miss_label)>0){
        warning('Does not contain cells from the following celltype(s): ', paste(miss_label,collapse = ' '))
    }
    
    # check response_variable numeric
    y <- sample_meta[,response_variable,drop = F]
    convert_column <- which(!sapply(y,is.numeric))
    response_unique <- unique(y)
    rank_zero <- ifelse(is.null(ncol(y)), length(response_unique)==1,nrow(response_unique)==1)
    
    ## Check if all values are numeric
    if(length(convert_column)>0 & analysis_type %in% c('pearson','spearman','cancor')){
        warning(paste0(paste(colnames(y)[convert_column],collapse = ','),' is not numeric. Proceed with factorized numeric values.'))
    }

    ## m1: aggregate (analysis_type)
    if(method == 'aggregate'){
        summary <- c()

        for(i in 1:length(unq_id)){
            tid <- unq_id[i]
            if(verbose) message(paste0('node ',tid))
            node_info <- leaves_info %>% dplyr::filter(id == tid) %>% filter(children %in% rownames(ctprop))

            if(nrow(node_info)>0){
                if(node_info$leaf[1]){
                    # leaf
                    pb <- ctprop[node_info$label,]
                }else{
                    # non-leaf
                    children <- node_info$children
                    pb <- ctprop[children,,drop=F] %>% colSums
                }

                # covariate adjustment
                tmp_meta <- suppressMessages(inner_join(data.frame(sample = names(pb)),sample_meta))
                response <- tmp_meta[,response_variable,drop = F]
                if(length(unique(pb))>1){
                    if(!is.null(formula) & nrow(tmp_meta)>1){
                        # drop factor level
                        factor.id <- which(sapply(1:ncol(tmp_meta),function(x) is.factor(tmp_meta[,x])))
                        for(k in factor.id){
                            tmp_meta[,k] <- droplevels(tmp_meta[,k])
                        }
                        fit <- lm(as.formula(paste0('pb',formula)),data = tmp_meta)
                        pb <- residuals(fit)
                        names(pb) <- tmp_meta$sample
                    }
                }

                if(length(response_variable) == 1){
                    # Switch to pearson correlation (absolute value = cancor)
                    if(length(unique(pb))>1){
                        summary <- rbind(summary,
                                         cancor_permute(pb,
                                                        response,
                                                        method = analysis_type,
                                                        num_permutations,
                                                        alternative))
                    }else{
                        summary <- rbind(summary,
                                         data.frame(cancor = NA,p = NA))
                    }
                    corr_type <- analysis_type
                }else{
                    # Canonical correlation
                    if(length(unique(pb))>1){
                        if(!(analysis_type %in% c('regression','cancor'))){
                            stop('Cannot compute Pearson/Spearman correlation of multiple phenotypes jointly.')
                        }

                        summary <- rbind(summary,
                                         cancor_permute(pb,
                                                        response,
                                                        method = analysis_type,
                                                        num_permutations,
                                                        alternative))
                    }else{
                        summary <- rbind(summary,
                                         data.frame(cancor = NA,p = NA))
                    }
                    corr_type <- ifelse(analysis_type == 'regression','regression','cancor')
                }
            }else{
                summary <- rbind(summary,
                                 data.frame(cancor = NA,p = NA))
                corr_type <- ifelse(analysis_type %in% c('pearson','spearman','regression'),analysis_type,'cancor')
            }
        }
        summary$method <- 'aggregate'
        summary$analysis_type <- corr_type
        pc.ls <- NULL

    }else if(method=='concat_leaf'){
        # m2: concatenate leaf
        if(analysis_type %in% c('pearson','spearman')){
            stop('In concatenating leaves, cannot calculate pearson/spearman correlation.')
        }
        summary <- c()
        pc.ls <- list()
        for(i in 1:length(unq_id)){
            tid <- unq_id[i]
            if(verbose) message(paste0('node ',tid))
            node_info <- leaves_info %>% dplyr::filter(id==tid) %>% filter(children %in% rownames(ctprop))

            if(nrow(node_info)>0){
                pb <- ctprop[node_info$children,,drop=F]

                # covariate adjustment
                tmp_meta <- suppressMessages(inner_join(data.frame(sample = colnames(pb)),sample_meta))
                response <- tmp_meta[,response_variable,drop = F]

                if(!is.null(formula) & nrow(tmp_meta)>1){
                    factor.id <- which(sapply(1:ncol(tmp_meta),function(x) is.factor(tmp_meta[,x])))
                    for(k in factor.id){
                        tmp_meta[,k] <- droplevels(tmp_meta[,k])
                    }
                    fit <- lm(as.formula(paste0('t(pb)',formula)),data = tmp_meta)
                    pb <- t(residuals(fit))
                }

                if(node_info$leaf[1]){
                    # leaf
                    pc <- matrix(pb,ncol = 1,dimnames = list(tmp_meta$sample))
                }else{
                    # PCA
                    pc <- extract_PC(pb,num_PCs)
                }
                pc.ls[[i]] <- pc
                summary <- rbind(summary,
                                 cancor_permute(pc,
                                                response,
                                                method = ifelse(analysis_type == 'regression','regression_concat','cancor'),
                                                num_permutations,
                                                alternative,
                                                component_id))
            }else{
                summary <- rbind(summary,
                                 data.frame(cancor = NA, p = NA))
                pc.ls[[i]] <- NULL
            }
        }
        summary$method <- 'concat_leaf'
        summary$analysis_type <- analysis_type
        names(pc.ls) <- unq_id

    }else if(method=='concat_immediate_children'){
        # m3: concatenate immediate children
        if(analysis_type %in% c('pearson','spearman')){
            stop('In concatenating immediate children, cannot calculate pearson/spearman correlation.')
        }
        summary <- c()
        pc.ls <- list()
        annot_data <- leaves_info %>% select(label,id) %>% unique

        for(i in 1:length(unq_id)){
            mid <- unq_id[i]
            if(verbose) message(paste0('node ',mid))
            mid_children <- immediate_children$children[immediate_children$id==mid]
            tid <- annot_data %>% filter(label %in% mid_children) %>% select(id) %>% unlist

            node_info <- leaves_info %>% dplyr::filter(id %in% tid) %>% filter(children %in% rownames(ctprop))
            if(nrow(node_info)>0){
                node_info_ls <- split(node_info,factor(node_info$id))
                pb <- do.call(rbind,lapply(node_info_ls,function(r){
                    if(nrow(r)>1){
                        ctprop[r$children,,drop=F] %>% colSums()
                    }else{
                        ctprop[r$children,,drop=F]
                    }
                }))

                # covariate adjustment
                tmp_meta <- suppressMessages(inner_join(data.frame(sample = colnames(pb)),sample_meta))
                response <- tmp_meta[,response_variable,drop = F]
                if(!is.null(formula) & nrow(tmp_meta)>1){
                    # drop factor level
                    factor.id <- which(sapply(1:ncol(tmp_meta),function(x) is.factor(tmp_meta[,x])))
                    for(k in factor.id){
                        tmp_meta[,k] <- droplevels(tmp_meta[,k])
                    }
                    fit <- lm(as.formula(paste0('t(pb)',formula)),data = tmp_meta)
                    pb <- t(residuals(fit))
                }

                # PCA
                if(nrow(pb)>1){
                    pc <- extract_PC(pb,num_PCs)
                }else{
                    pc <- matrix(pb,ncol = 1,dimnames = list(tmp_meta$sample))
                }
                pc.ls[[i]] <- pc
                summary <- rbind(summary,
                                 cancor_permute(pc,
                                                response,
                                                method = ifelse(analysis_type == 'regression','regression_concat','cancor'),
                                                num_permutations,
                                                alternative,
                                                component_id))
            }else{
                summary <- rbind(summary,
                                 data.frame(cancor = NA,p = NA))
                pc.ls[[i]] <- NULL
            }
        }
        summary$method <- 'concat_immediate_children'
        summary$analysis_type <- analysis_type
        names(pc.ls) <- unq_id
    }else{
        stop('Incorrect `method` input.')
    }
    summary$adjp <- p.adjust(summary$p,method = 'BY') # BY procedure - dependent test
    summary$id <- unq_id
    res <- summary %>%
        select(id,method,analysis_type,cancor,p,adjp) %>%
        mutate(direction = ifelse(cancor>0,'+','-'),
               p.sign = ifelse(p<0.05,'sig','ns'),
               adjp.sign = ifelse(adjp<0.05,'sig','ns'))

    if(analysis_type=='regression'){
        if(method == 'aggregate'){
            colnames(res) <- gsub('cancor','chisq',colnames(res))
            res$chisq <- abs(res$chisq)
        }else{
            colnames(res) <- gsub('cancor','Fstat',colnames(res))
        }
    }else{
        colnames(res) <- gsub('cancor',res$analysis_type[1],colnames(res))
    }

    res.ls <- list(canonical_corr = res,
                   pc_ls = pc.ls)

}

