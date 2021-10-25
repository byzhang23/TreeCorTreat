#' Modify label for tree nodes
#'
#' Rename the node labels
#' @param new_label      A data frame specifying old and new names for each node. Must contain two columns named as 'old' (corresponds to old celltype name) and 'new' (corresponds to new celltype name).
#' @param hierarchy_list A hierarchy list by running \code{extract_hrchy_string()} or \code{extract_hrchy_seurat()} functions, which contains four elements: 'edges', 'layout', 'immediate_children' and 'leaves_info'.
#' @param cell_meta      Cell-level metadata, where each row is a cell. Must contain these columns: barcode, celltype and sample.
#' @param verbose        Show progress
#' @return A list that contains modified \code{hierarchy_list} and/or \code{cell_meta}
#' @export
#' @import dplyr
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji

modify_label <- function(new_label, hierarchy_list, cell_meta, verbose = T){


    new_label <- new_label %>% mutate(new = ifelse(is.na(new),old,new))
    find_new_name <- function(col){
        suppressMessages(inner_join(data.frame(old = col),new_label))[,2]
    }

    if(!is.null(hierarchy_list)){
        if(nrow(new_label) != nrow(hierarchy_list$layout)){
            stop('Incomplete label information. Please provide labels of all nodes (both old and new names).')
        }

        ## 1. Hierarchy_list
        if(verbose){
            message('Modify label in hierarchy_list')
        }
        # (1) edges
        hierarchy_list$edges <- hierarchy_list$edges %>% mutate(parent = find_new_name(parent),
                                                                children = find_new_name(children))
        # (2) layout
        hierarchy_list$layout <- hierarchy_list$layout %>% mutate(label = find_new_name(label))
        # (3) immediate_children
        hierarchy_list$immediate_children <- hierarchy_list$immediate_children %>% mutate(label = find_new_name(label),
                                                                                          children = find_new_name(children))
        # (4) leaves_info
        hierarchy_list$leaves_info <- hierarchy_list$leaves_info %>% mutate(label = find_new_name(label),
                                                                            children = find_new_name(children))
    }

    if(!is.null(cell_meta)){
        ## 2. cell_meta
        if(verbose){
            message('Modify label in cell_meta')
        }
        cell_meta <- cell_meta %>% mutate(celltype = find_new_name(celltype))
    }

    list(hierarchy_list = hierarchy_list,
         cell_meta = cell_meta)
}
