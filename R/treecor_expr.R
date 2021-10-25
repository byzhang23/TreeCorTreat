#' Pipeline for TreeCorTreat using gene expression as features
#'
#' @param expr                  A count gene expression matrix. Can directly extract from Seurat object using \code{seurat_object@assays$RNA@counts}.
#' @param hierarchy_list        A hierarchy list by running \code{extract_hrchy_string()} or \code{'extract_hrchy_seurat()'} functions, which contains four elements: 'edges', 'layout', 'immediate_children' and 'leaves_info'.
#' @param cell_meta             Cell-level metadata, where each row is a cell. Must contain these columns: barcode, celltype and sample.
#' @param sample_meta           Sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables such as covariates or outcome of interest.
#' @param response_variable     A vector of response variables.
#' @param method                A character string indicating which approach is used to summarize features. One of 'concat_leaf' or 'concat_immediate_children' or 'aggregate'(default).
#' @param formula               An object of class 'formula': a symbolic description of the model to be fitted, adjusting for confounders.
#' @param separate              A TRUE (default) or FALSE indicator, specifying how to evaluate multivariate outcomes.
#'                              \itemize{
#'                              \item TRUE: evaluate multivariate phenotype separately (it is equivalent to run this pipeline for each univariate phenotype).
#'                              \item FALSE: evaluate multivariate phenotype jointly.
#'                              }
#' @param analysis_type         Either 'cancor' (canonical correlation, by default) or 'regression' (F statistics) to evaluate association between gene expression and samples' phenotype.
#' @param num_cancor_components Number of canonical components to be extracted. Only works when separate = F.
#' @param num_permutations      Number of permutations (by default: 1000).
#' @param alternative           A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param num_PCs               Number of PCs used in canonical correlation calculation
#' @param num_hvgs              Number of highly variable genes extracted from sample-level pseudobulk. Highly variable genes (HVGs) are defined as genes with positive residuals by fitting a gene-specific lowess function of standard deviation against its mean. The HVGs are ordered by residuals in a descending order. Default is NULL, which includes all genes with positive residuals. Can specify a reasonable positive number (e.g. 1000).
#' @param threshold_pct_samples Threshold for selecting cell clusters with at least threshold_pct_samples\% of total samples. Ranges from 0 to 100, where '100' indicates to include cell clusters with all samples; '0' indicates to include all cell clusters. For any non-leaf node, there will be a trade-off between number of cell clusters and number of samples to be included.
#' @param method_threshold      A character string specifying method to apply the threshold (i.e. 'threshold_pct_samples'), must be one of 'global' or 'local'.
#'                              \itemize{
#'                              \item global (default): indicates a minimum sample size is calculated at global level and only affects non-leaf node with sample size that exceeds global cutoff.
#'                              \item local: indicates that threshold_pct_samples is applied to every hierarchy and affect every non-leaf node.
#'                              }
#' @param filter_prop           A number ranges from 0 to 1, to filter low expressed genes across samples (by default: 0.1). Genes with at least this proportion of samples with log2-normalized count greater than 0.01 are retained.
#' @param pseudobulk_list       A list of sample-level pseudobulk for each node.  Default is NULL. Each element is a dataframe with rows representing genes and columns representing samples. Users can provide their processed pseudobulk list (e.g. use combat or other batch correction methods) via this parameter. Note that the names of list shall be matched with id extracted from \code{hierarchy_list$layout[,id]}.
#' @param ncores                Number of cores to be used. If ncores > 1, it will be implemented in a parallel mode.
#' @param verbose               Show progress
#' @return A list that contains:
#' \itemize{
#' \item canonical_corr: A data frame for summary statistic (e.g. canonical correlation or F-statistic), p-value, adjusted p-value and label information for each node.
#' \item pc_ls: A list for top n PC matrices for each node
#' }
#' @export
#' @import dplyr parallel
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#' @examples
#'# default setting
#' result <- treecor_expr(expr,hierarchy_list, cell_meta, sample_meta, response_variable = 'severity')
#' # obtain summary statistic for each cell type
#' result$canonical_corr # or result[[1]]
#' # extract PC matrix for celltype 'T'
#' result$pc.ls[['T']]

treecor_expr <- function(expr,hierarchy_list, cell_meta, sample_meta, response_variable, method = 'aggregate',formula = NULL, separate = T, analysis_type = 'cancor',num_cancor_components = 1, num_permutations = 1000, alternative = 'two.sided', num_PCs = 2,num_hvgs = NULL, threshold_pct_samples = 100, method_threshold = 'global', filter_prop = 0.1, pseudobulk_list = NULL, ncores = parallel::detectCores(), verbose = T){

    ## add progress bar
    if (Sys.info()[['sysname']]=='Windows') {
        message('Parallel is disabled for Windows. Running with one core')
        ncores <- 1
    }

    ## Extract info
    label_info <- hierarchy_list$layout
    leaves_info <- hierarchy_list$leaves_info
    if(method == 'concat_immediate_children'){
        immediate_children <- hierarchy_list$immediate_children
    }else{
        immediate_children <- NULL
    }

    ## Check (ADD: check NA in cell_meta, sample_meta; required columns in both meta)
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
        warning('Number of canonical components > 1 when method = `aggregate`. Set num_cancor_components = 1')
        num_cancor_components <- 1
    }
    if(length(response_variable)==1 & separate){
        separate <- F
    }
    # if(length(response_variable)>1 & separate==F & analysis_type == 'regression'){
    #     stop(paste0('Please provide a new column that combine ',length(response_variable), ' phenotypes jointly'))
    # }
    if(sum(analysis_type %in% c('regression','cancor'))==0){
        stop('Invalid analysis type. Must be one of regression or cancor.')
    }
    if(analysis_type == 'regression'){
        analysis_type <- 'regression_concat'
    }

    unq_leaves_info <- leaves_info %>% filter(leaf) %>% select(id,label) %>% unique
    all_celltypes <- unique(cell_meta$celltype)
    miss_label <- setdiff(unq_leaves_info$label,all_celltypes)
    if(length(miss_label)>0){
        warning('Does not contain cells from the following celltype(s): ', paste(miss_label,collapse = ' '))
    }
    ## =============================== ##
    ## 1. Pseudobulk, Regression, HVGs
    ## =============================== ##
    unq_id <- unique(leaves_info$id)
    leaf_id <- leaves_info$id[leaves_info$leaf]

    if(method == 'concat_leaf'){
        cal_id <- leaf_id
    }else{
        cal_id <- rev(unq_id) # first leaf then non-leaf node
    }

    ## pb.ls
    if(!is.null(pseudobulk_list)){
        # input pseudobulk list (user-specify)
        message('Input pseudobulk list:')
        pb_names <- intersect(as.character(cal_id),names(pseudobulk_list))
        if(length(pb_names) != length(cal_id)){
            stop('Invalid pseudobulk list, please check if the following nodes exist: \n', paste(setdiff(as.character(cal_id),pb_names),collapse = ','))
        }
        pb.ls <- pseudobulk_list
    }else{
        # default pipeline
        if(verbose){
            message('=== Calculate Pseudobulk and Extract Highly Variable Genes ===')
            message(paste0('Calculate pseudobulk on ',length(cal_id),' nodes'))
            if(!is.null(formula)){
                message(paste0('Adjusted for confounders: ',formula))
            }
        }
        pb.ls <- mclapply(cal_id,function(tid){
            if(verbose) message(paste0('node ',tid))
            node_info <- leaves_info %>% filter(id==tid) %>% filter(children %in% all_celltypes)
            if(nrow(node_info)>0){
                sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
                pb <- construct_pseudobulk(expr,
                                           sub_meta,
                                           filter_prop = filter_prop)
                if(ncol(pb) > 1){
                    if(!is.null(formula)){
                        tmp_meta <- suppressMessages(inner_join(data.frame(sample = colnames(pb)),sample_meta))
                        # drop factor level
                        factor.id <- which(sapply(1:ncol(tmp_meta),function(x) is.factor(tmp_meta[,x])))
                        for(k in factor.id){
                            tmp_meta[,k] <- droplevels(tmp_meta[,k])
                        }
                        fit <- lm(as.formula(paste0('t(pb)',formula)),data = tmp_meta)
                        pb <- t(residuals(fit))
                    }
                }else{
                    pb <- NULL
                }
            }else{
                pb <- NULL
            }
            pb
        },mc.cores = ncores)
        names(pb.ls) <- cal_id
    }

    ## HVG selection
    if(verbose) message('Selecting highly variable genes...')
    pb.ls <- mclapply(cal_id,function(tid){
        if(verbose) message(paste0('node ',tid))
        pb <- pb.ls[[as.character(tid)]]
        if(!is.null(pb)){
            if(ncol(pb)>1){
                hvgs <- select_hvgs(pb,num_hvgs)
                pb <- pb[hvgs,]
            }else{
                pb <- NULL
            }
        }else{
            pb <- NULL
        }
        pb
    },mc.cores = ncores)
    names(pb.ls) <- cal_id

    ## =============================== ##
    ## 2. Calculate PC (method)
    ## =============================== ##
    global_cutoff <- nrow(sample_meta) * threshold_pct_samples/100
    if(verbose){
        message(paste0('=== Use ',method,' method to summarize gene expression features and Perform PCA ==='))
        message(paste0('Global sample size cutoff: ', global_cutoff))
        message(paste0('Use ',method_threshold,' method to select cell clusters at each non-leaf node'))
    }

    if(method == 'aggregate'){

        pc.ls <- mclapply(unq_id,function(tid){
            if(verbose) message(paste0('node ',tid))
            pc <- extract_PC(pb.ls[[as.character(tid)]],num_PCs)
            pc
        },mc.cores = ncores)
        names(pc.ls) <- unq_id

    }else if(method == 'concat_leaf'){

        pc.ls <- mclapply(unq_id,function(tid){
            if(verbose) message(paste0('node ',tid))
            if(tid %in% leaf_id){
                pc <- extract_PC(pb.ls[[as.character(tid)]],num_PCs)
            }else{
                node_info <- leaves_info %>% filter(id==tid)
                children_id <- leaves_info %>% filter(leaf==T & label %in% node_info$children) %>% select(id) %>% unlist %>% as.character
                nsample <- sapply(pb.ls[children_id],ncol) %>% unlist
                local_max <- max(nsample)
                if(method_threshold == 'global'){
                    # Only retain cell types with threshold_pct_samples/100 * nrow(sample_meta) or largest sample size
                    final_id <- names(nsample)[nsample >= min(global_cutoff,local_max)]
                }else if(method_threshold == 'local'){
                    local_cutoff <- local_max * threshold_pct_samples/100
                    final_id <- names(nsample)[nsample >= local_cutoff]
                }else{
                    stop('Incorrect `method_threshold` input.')
                }

                sub.pb.ls <- pb.ls[final_id]
                # extract intersected samples
                intersect.cols <- Reduce(intersect,lapply(sub.pb.ls,colnames))
                pc <- extract_PC(do.call(rbind,lapply(sub.pb.ls,function(x) x[,intersect.cols,drop=F])),num_PCs)
            }
            pc
        },mc.cores = ncores)
        names(pc.ls) <- unq_id

    }else if(method == 'concat_immediate_children'){

        pc.ls <- mclapply(unq_id,function(tid){
            if(verbose) message(paste0('node ',tid))
            if(tid %in% leaf_id){
                pc <- extract_PC(pb.ls[[as.character(tid)]],num_PCs)
            }else{
                tid_child <- immediate_children$children[immediate_children$id==tid]
                children_id <- immediate_children %>% filter(label %in% tid_child) %>% select(id) %>% unique %>% unlist %>% as.character
                nsample <- sapply(pb.ls[children_id],ncol) %>% unlist
                local_max <- max(nsample)
                if(method_threshold == 'global'){
                    # Only retain cell types with threshold_pct_samples/100 * nrow(sample_meta) or largest sample size
                    final_id <- names(nsample)[nsample >= min(global_cutoff,local_max)]
                }else if(method_threshold == 'local'){
                    local_cutoff <- local_max * threshold_pct_samples/100
                    final_id <- names(nsample)[nsample >= local_cutoff]
                }else{
                    stop('Incorrect `method_threshold` input.')
                }

                sub.pb.ls <- pb.ls[final_id]
                # extract intersected samples
                intersect.cols <- Reduce(intersect,lapply(sub.pb.ls,colnames))
                pc <- extract_PC(do.call(rbind,lapply(sub.pb.ls,function(x) x[,intersect.cols,drop=F])),num_PCs)
            }
            pc
        },mc.cores = ncores)
        names(pc.ls) <- unq_id

    }else{
        stop('Incorrect `method` input.')
    }

    ## =============================== ##
    ## 3. Cancor and p-value
    ## =============================== ##
    if(verbose){
        message('=== Calculate statistics and Perform permutation test ===')
        message('Response variable(s):\n',paste(response_variable,collapse = ' '))
        message(paste0('Evaluate response variable(s) ',ifelse(separate,'separately','jointly')))
    }

   if(!separate){
        if(num_cancor_components == 1){
            cancor.dat <- expr_cancor(pc.ls,
                                      sample_meta,
                                      response_variable,
                                      analysis_type,
                                      num_permutations,
                                      alternative,
                                      component_id = 1,
                                      ncores)
            colnames(cancor.dat)[-1] <- paste0(paste(response_variable,collapse = '_'),'.',colnames(cancor.dat)[-1])
        }else{

            cancor.dat <- suppressMessages(Reduce(inner_join,lapply(1:num_cancor_components,function(c_id){
                  tmp <- expr_cancor(pc.ls,
                                     sample_meta,
                                     response_variable,
                                     analysis_type,
                                     num_permutations,
                                     alternative,
                                     component_id = c_id,
                                     ncores)
                  colnames(tmp)[-1] <- paste0(paste(response_variable,collapse = '_'),'_Component',c_id,'.',colnames(tmp)[-1])
                  tmp
            })))
        }
    }else{
        ## evaluate response variables separately
        cancor.dat <- suppressMessages(Reduce(inner_join,lapply(response_variable,function(r){
              tmp <- expr_cancor(pc.ls,
                                 sample_meta,
                                 r,
                                 analysis_type,
                                 num_permutations,
                                 alternative,
                                 component_id = 1,
                                 ncores)
              colnames(tmp)[-1] <- paste0(r,'.',colnames(tmp)[-1])
              tmp
        })))
    }

    res <- suppressMessages(inner_join(cancor.dat,label_info))
    names(pc.ls) <- label_info$label[unq_id]

    res.ls <- list(canonical_corr = res,
                   pc_ls = pc.ls)
}
