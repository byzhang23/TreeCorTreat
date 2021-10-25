#' Extract hierarchy from Seurat object (Data-driven approach)
#'
#' This function takes in either Seurat object or phylogenetic tree and generates hierarchy information.
#' @param input   An input data
#' @param type    Input data type, must be one of 'phylo_tree' (default) or 'seurat_object'.
#' @return A list of hierarchy information:
#'         \itemize{
#'         \item edges: a data frame that documents that parent-children relationship.
#'         \item layout: a data frame that contains id, xy coordiates, label and a leaf indicator for each node.
#'         \item immediate_children: a data frame that captures immediate children for each node.
#'         \item leaves_info: a data frame that summarizes leaf children for each node.
#'         }
#'
#' @export
#' @import dplyr
#' @importFrom Seurat Tool
#' @importFrom ggraph create_layout
#' @importFrom igraph graph_from_data_frame all_simple_paths
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#' @examples
#' # Assume seurat_object is obtained after running treecor_harmony()
#' extract_hrchy_seurat(seurat_object, type = 'seurat_object')


extract_hrchy_seurat <- function(input, type = 'phylo_tree'){

    if(type == 'seurat_object'){
        ## Extract Phylogenetic tree
        phylo_tree <- Seurat::Tool(object = input, slot = 'BuildClusterTree')
        if(is.null(phylo_tree)){
            stop('`BuildClusterTree` slot is empty. Please run `BuildClusterTree` and try again.')

        }
    }else{
        phylo_tree <- input
    }
    if(sum(phylo_tree$tip.label=='0')>0) phylo_tree$tip.label <- as.character(as.integer(phylo_tree$tip.label) + 1) # start from 1 (leaf node)

    ## Return 1: edges
    edges <- data.frame(parent = as.character(phylo_tree$edge[,1]),children = as.character(phylo_tree$edge[,2]))

    ## Return 2: layout
    graph <- graph_from_data_frame(edges %>% dplyr::rename(from = parent,to = children))
    layout <- create_layout(graph,layout = 'dendrogram') %>%
        dplyr::rename(label = name,id = .ggraph.index) %>%
        select(x,y,id,label,leaf)
    leaf_leaf <- data.frame(parent = layout$label[layout$leaf]) %>%
        mutate(children = parent)

    ## Return 3: immediate_children
    immediate_children <- suppressMessages(
        rbind(edges,leaf_leaf) %>%
            dplyr::rename(label=parent) %>%
            inner_join(layout))

    ## Return 4: leaves_info
    leaves <- layout$label[layout$leaf]
    leaves_info <- suppressMessages(do.call(rbind,lapply(unique(edges$parent),function(p){
        path_list <- all_simple_paths(graph,from = p)
        last <- sapply(path_list,function(x) x[length(x)] %>% names)
        data.frame(label = p,children = intersect(leaves,last))
    })) %>% rbind(leaf_leaf %>% dplyr::rename(label = parent)) %>% inner_join(layout))

    ## return
    res <- list(edges = edges,
                layout = layout,
                immediate_children = immediate_children,
                leaves_info = leaves_info)


}
