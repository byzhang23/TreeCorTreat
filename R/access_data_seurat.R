#' Extract cell-level meta data and count matrix from Seurat object
#'
#' @param seurat_obj        A Seurat object
#' @param reduction         Extract x,y coordinates from a type of dimension reduction. If run default \code{harmony_integration()} function, can be one of 'umap' (default) or 'harmony'.
#' @param output_dir        An output directory
#' @return A list of filtered inputs:
#'         \itemize{
#'         \item count: a filtered count matrix
#'         \item cell_meta: a cell-level metadata
#'         \item phylo_tree: phylogenetic tree
#'         }
#'
#' @export
#' @import dplyr Seurat
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#' @examples
#' \dontrun{
#' access_data_seurat(seurat_object,reduction = 'umap')
#' }
#'

access_data_seurat <- function(seurat_obj, reduction = 'umap', output_dir){

    if(!dir.exists(output_dir)) dir.create(output_dir)

    ## coord
    reduction <- tolower(reduction)
    if(length(reduction) > 1){
        warning(paste0('There are multiple reductions provided. Will proceed with first one: ', reduction[1]))
        reduction <- reduction[1]
    }
    red_ls <- seurat_obj@reductions
    if(sum(names(red_ls) %in% reduction)>0){

        coord <- red_ls[[reduction]]@cell.embeddings %>% data.frame %>% mutate(barcode = rownames(.))

    }else{
        stop(paste0('Not found in `seurat_obj@reductions`: ',reduction))
    }

    ## cell meta
    cell_meta <- suppressMessages(seurat_obj[[]] %>% mutate(barcode = rownames(.)) %>% inner_join(coord) %>% mutate(celltype = as.character(seurat_clusters)))
    ## count
    count <- seurat_obj@assays$RNA@counts
    saveRDS(cell_meta, paste0(output_dir,'cell_meta.rds'))
    saveRDS(count, paste0(output_dir,'count.rds'))

    ## Extract Phylogenetic tree
    phylo_tree <- Seurat::Tool(object = seurat_obj, slot = 'BuildClusterTree')
    if(is.null(phylo_tree)){
        warning('`BuildClusterTree` slot is empty. Please run `BuildClusterTree`.')
        phylo_tree <- NULL

    }else{
        if(sum(phylo_tree$tip.label=='0')>0) phylo_tree$tip.label <- as.character(as.integer(phylo_tree$tip.label) + 1) # start from 1 (leaf node)
        saveRDS(phylo_tree,paste0(output_dir,'phylo_tree.rds'))
    }

    res.ls <- list(count = count,
                   cell_meta = cell_meta,
                   phylo_tree = phylo_tree)

}


