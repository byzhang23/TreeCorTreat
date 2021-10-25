#' Extract hierarchy from a string (based on prior knowledge)
#'
#' Parse an input string by extracting parent-children relationship.
#' @param input_string An input string with nested labeled parentheses that characterizes hierarchy of annotated cell types. Avoid using comma (','), parentheses ('()') and special character (e.g.'@') in celltype names.
#' @param special_character A character denotes non-leaf node (including root node and internal nodes).
#' @param plot  Plot tree structure or not.
#' @return A list of hierarchy information:
#'         \itemize{
#'         \item edges: a data frame that documents that parent-children relationship.
#'         \item layout: a data frame that contains id, xy coordiates, label and a leaf indicator for each node.
#'         \item immediate_children: a data frame that captures immediate children for each node.
#'         \item leaves_info: a data frame that summarizes leaf children for each node.
#'         }
#' @export
#' @import dplyr ggplot2
#' @importFrom ggraph create_layout
#' @importFrom igraph graph_from_data_frame all_simple_paths
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#' @examples
#' input_string = '@T(@CD4+(CD4_c1,CD4_c2),@CD8+(CD8_c1,CD8_c2))'
#' extract_hrchy_string(input_string, special_character = '@')

extract_hrchy_string <- function(input_string,special_character = '@', plot = T){

    if(special_character %in% c('(',')','()',',')){
        stop('Avoid using parentheses or comma as special character. Please choose another symbol instead.')
    }
    message(paste0('Make sure these symbols are not used in cell type names: parentheses (`(`,`)`), comma (`,`) and special character `',special_character,'`'))
    ## Find matched parentheses
    open <- gregexpr(paste0('\\('), input_string)[[1]]
    close <- gregexpr(paste0('\\)'), input_string)[[1]]
    if(length(open) != length(close) | sum(open<0)>0 | sum(close<0)>0){
        stop('Parenetheses are not matched. Please check the input_string again.')
    }

    ## Sort ')' positions
    match_parentheses <- function(input_string,open){

        close <- open
        counter <- 1
        while (counter > 0){
            close <- close + 1
            m <- substr(input_string,close,close)
            if(m == '('){counter <- counter + 1}else if(m == ')'){counter <- counter - 1}
        }
        return(close)
    }
    close <- sapply(open,function(o) match_parentheses(input_string,o))

    ## Check input_string
    nonleaf_idx <- gregexpr(paste0('\\',special_character), input_string)[[1]]
    if(length(nonleaf_idx) != length(open) | sum(nonleaf_idx<0)>0){
        stop(paste0('Special character ',special_character,' is not matched. Please check the input_string again.'))
    }

    ## Direct extract parent(children) format
    rm_most_inner_parantheses <- function(text){
        gsub("\\([^()]*\\)", "", text)
    }

    extract_format <- function(text){
        old <- text
        new <- NULL
        repeat {
            new <- rm_most_inner_parantheses(old)
            if (identical(new,old)|!grepl('\\(',new)) break
            old <- new
        }
        return(old) # parent(children node)
    }

    hierarchy_levels <- sapply(1:length(nonleaf_idx),function(i){
        sub_text <- substr(input_string,nonleaf_idx[i],close[i])
        extract_format(sub_text)
    })

    ## Return 1: edges
    edges <- do.call(rbind,lapply(hierarchy_levels,function(hl){
        hl <- gsub('\\)','',hl)
        parent <- strsplit(hl,'\\(')[[1]][1]
        children <- strsplit(strsplit(hl,'\\(')[[1]][2],'\\,')[[1]]
        data.frame(parent = parent,children = children)
    }))

    ## Return 2: layout
    edges <- edges %>% mutate(parent = gsub(paste0('\\',special_character),'',parent),
                              children = gsub(paste0('\\',special_character),'',children))
    graph <- graph_from_data_frame(edges %>%
                                   dplyr::rename(from = parent,to = children))
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

    ## plot
    if(plot){
        mygraph <- graph_from_data_frame(res$edges %>% dplyr::rename(from = parent, to = children))
        layout <- create_layout(mygraph,layout = 'dendrogram')

        (ggraph(layout, layout = 'dendrogram', circular = FALSE) +
            geom_node_point() +
            geom_edge_link(color = 'darkgrey') +
            theme_void(base_size = 15) +
            labs(y='') +
            coord_flip() +
            theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
                  panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
            scale_y_reverse(expand = expansion(mult = c(0.2, 0.5))) +
            geom_text(data = res$leaves_info %>% filter(leaf) %>% select(label,x,y) %>% unique,aes(x=x,y=y-0.1,label = label),hjust = 0,vjust = 0.5) +
            geom_text(data = res$leaves_info %>% filter(!leaf) %>% select(label,x,y) %>% unique,aes(x=x,y=y+0.1,label = label),hjust = 1,vjust = 0.5)) %>% print
    }

    return(res)

}
