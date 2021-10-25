#' Plot heatmap for top n differentially expressed genes
#'
#' @param sample_meta       A data frame for sample-level metadata, where each row is a sample. Must contain 'sample' column and additional variables such as covariates or outcome of interest.
#' @param pseudobulk        A sample-level pseudobulk, which can be obtained by running `limma_expr` pipeline and can be extracted via `result$pseudobulk.ls[[cell_type_name]]`.
#' @param deg_result        Differential expression output obtained by running limma_expr` pipeline` and can be extracted via `result$dge.ls[[phenotype]][[cell_type_name]]`.
#' @param deg_logFC         Filter `deg_result` by logFC. One of the following:
#' #' \itemize{
#' \item positive: limit to genes with logFC>0.
#' \item negative: limit to genes with logFC<0.
#' \item both (default): keep genes with both positive and negative logFC.
#' }
#' @param top_n             A number to extract top n differentially expressed genes. If \code{deg_logFC = 'both'}, then top n DEGs with positive logFC will be combined with top n DEGs with negative logFC, resulting in a total of 2n genes.
#' @param scale             A character indicates if the values should be centered and scaled in either the `row` direction or the `column` direction, or `none`.
#' @param cluster_rows      A True of False indicator to determine if rows should be clustered.
#' @param cluster_cols      A True of False indicator to determine if columns should be clustered.
#' @param show_rows         A True of False indicator to determine if row names should be displayed.
#' @param cluster_cols      A True of False indicator to determine if column names should be displayed.
#' @param annotation_col    A vector specifying column names to be annotated in the sidebars of a heatmap.
#' @param annotation_colors List for manually specifying colors for annotation_col.
#' @param color             A vector of colors used in heatmap.
#' @param breaks            A sequence of numbers specifying range of values in the matrix. Must be an element longer than color vector. See more in `?pheatmap`.
#' @param fontsize          Font size.
#' @return Heatmap with DEGs on the row and samples on the column.
#' @export
#' @import dplyr pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#'
#' @examples
#' # Assume obtain pseudobulk of celltype 'T' from treecor_deg()
#' result <- treecor_deg(expr,hierarchy_list, cell_meta, sample_meta, response_variable = 'severity')
#' pseudobulk <- result$pseudobulk.ls[['T']]
#' deg_result <- result$dge.ls$severity[['T']]
#' # Top 10 DEGs with logFC>0
#' treecor_degheatmap(sample_meta,pseudobulk, deg_result, top_n = 10, deg_logFC = 'positive')
#' # Top 10 DEGs with logFC<0
#' treecor_degheatmap(sample_meta,pseudobulk, deg_result, top_n = 10, deg_logFC = 'negative')
#' # Combined the above two options: top 10 positively-associated DEGs and top 10 negatively-associated DEGs
#' treecor_degheatmap(sample_meta,pseudobulk, deg_result, top_n = 10, deg_logFC = 'both')


treecor_degheatmap <- function(sample_meta,pseudobulk, deg_result, top_n = 10, deg_logFC = 'both', scale = 'row',cluster_rows = F,cluster_cols = F, show_colnames = F, show_rownames = T, annotation_col = NA, annotation_colors = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),breaks = NA,fontsize = 12){

    if(deg_logFC=='both'){
        res <- deg_result
    }else if(deg_logFC=='positive'){
        res <- deg_result[deg_result[,1]>0,]
        if(nrow(res)==0) res <- NULL
    }else if(deg_logFC=='negative'){
        res <- deg_result[deg_result[,1]<0,]
        if(nrow(res)==0) res <- NULL
    }

    if(is.null(res)){
        stop('No differentiallly expressed genes.')
    }else{
        # extract top n
        if(deg_logFC=='both'){
            pos <- res[res[,1]>0,]
            neg <- res[res[,1]<0,]
            pos_topn <- pos %>% arrange(.[[3]],desc(abs(.[[1]]))) %>% head(top_n)
            neg_topn <- neg %>% arrange(.[[3]],desc(abs(.[[1]]))) %>% head(top_n)

            genes <- c(rownames(pos_topn),rownames(neg_topn))
        }else{
            res_topn <- res %>% arrange(.[[3]],desc(abs(.[[1]]))) %>% head(top_n)
            genes <- rownames(res_topn)
        }
        genes <- setdiff(genes,NA)

        # annotation
        if(!is.na(annotation_col[1])){
            annot <- suppressMessages(inner_join(data.frame(sample = colnames(pseudobulk)),sample_meta[,c('sample',annotation_col)]) %>% arrange(across(annotation_col)))
            rownames(annot) <- annot$sample
            mat <- pseudobulk[genes,annot$sample]
            annot <- annot %>% select(-sample)
        }else{
            annot <- NA
            mat <- pseudobulk[genes,]
        }

        if(!is.na(breaks[1])){
           if(scale=='none'){
               if(min(breaks) > min(mat)) breaks <- c(min(mat),breaks)
               if(max(breaks) < max(mat)) breaks <- c(breaks,max(mat))
           }else if(scale=='col'){
               mat_scale <- scale(mat)
               if(min(breaks) > min(mat_scale)) breaks <- c(min(mat_scale),breaks)
               if(max(breaks) < max(mat_scale)) breaks <- c(breaks,max(mat_scale))

           }else if(scale=='row'){
               mat_scale <- t(scale(t(mat)))
               if(min(breaks) > min(mat_scale)) breaks <- c(min(mat_scale),breaks)
               if(max(breaks) < max(mat_scale)) breaks <- c(breaks,max(mat_scale))

           }
        }

        nbreaks <- length(breaks)
        ncolor <- length(color)
        if(nbreaks != (ncolor + 1)){
            warning(paste(paste0('Length of `breaks` vector: ', nbreaks),
                          paste0('Length of `color` vector: ', ncolor),
                          paste0('length(breaks) ', ifelse(nbreaks<ncolor+1,'<','>'),' (length(color) + 1)'),
                          sep = '\n'))
        }
        # heatmap
        pheatmap(mat,
                 color = color,
                 breaks = breaks,
                 scale = scale,
                 show_colnames = show_colnames,
                 show_rownames = show_rownames,
                 cluster_rows = cluster_rows,
                 cluster_cols = cluster_cols,
                 annotation_col = annot,
                 annotation_colors = annotation_colors,
                 fontsize = fontsize)
    }

}
