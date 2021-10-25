#' Harmony Integration
#'
#' This function is developed using 'Seurat v3.2.2'. It takes in a raw count gene expression matrix and a sample meta data frame and performs harmony integration.
#' @param count             A raw count gene expression matrix with genes on rows and cells on columns. Note that cell barcode shall use ':' to separate sample name and barcode (i.e. "sample:barcode")
#' @param sample_meta       Sample metadata. Must contain a column named as 'sample'.
#' @param output_dir        Output directory
#' @param num_PCs           Number of PCs used in integration (by default: 20)
#' @param num_harmony       Number of harmony embedding used in integration (by default: 20)
#' @param num_features      Number of features used in integration (by default: 2000)
#' @param min_cells         Include features detected in at least this many cells (by default: 0). Same as \code{'min.cells'} parameter in \code{'CreateSeuratObject()'} function from 'Seurat' package.
#' @param min_features      Include cells where at least this many features are detected (by default: 0). Same as \code{'min.features'} parameter in \code{'CreateSeuratObject()'} function from 'Seurat' package.
#' @param pct_mito_cutoff   Include cells with less than this many percent of mitochondrial percent are detected (by default: 20). Ranges from 0 to 100. Will be used as a QC metric to subset the count matrix. Genes starting with 'MT-' are defined as a set of mitochondrial genes.
#' @param exclude_genes     Additional genes to be excluded from integration. Will subset the count matrix.
#' @param vars_to_regress   Variables to be regressed out during Harmony integration (by default: 'sample'). Same as \code{'group.by.vars'} in \code{'RunHarmony()'} function from 'harmony' package.
#' @param resolution        A clustering resolution (by default: 0.5). A higher (lower) value indicates larger (smaller) number of cell subclusters.
#' @param verbose           Show progress
#' @return A Seurat object
#' @export
#' @import Seurat harmony dplyr
#' @importFrom ape plot.phylo nodelabels
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#'
#' @examples
#' # default setting
#' treecor_harmony(count, sample_meta, output_dir)
#' # additionally regress out 'study' ID
#' treecor_harmony(count, sample_meta, output_dir,vars_to_regress = c('sample','study'))
#' # increase clustering resolution (with more refined cell clusters)
#' treecor_harmony(count, sample_meta, output_dir,resolution = 0.8)

treecor_harmony <- function(count, sample_meta, output_dir, num_PCs = 20, num_harmony = 20,num_features = 2000, min_cells = 0, min_features = 0, pct_mito_cutoff = 20, exclude_genes = NULL, vars_to_regress = c('sample'), resolution = 0.5,verbose = T){

    if(!dir.exists(output_dir)) dir.create(output_dir)
    set.seed(12345)
    ## 1. input data
    if(verbose){
        message('=== Read input dataset ===')
        message('Dimension of raw data:\n', nrow(count),' ', ncol(count))
    }
    ## 2. filter genes with zero expr across all cells
    if(isClass('dgCMatrix',count)){
        count <- count[Matrix::rowSums(count) >0,]
    }else{
        count <- count[rowSums(count) > 0,]
    }

    ## 3. harmony
    if(verbose){
        message('=== Harmony Integration ===')
    }
    u <- CreateSeuratObject(counts = count,min_cells = min_cells, min_features = min_features)
    u[["percent.mt"]] <- PercentageFeatureSet(u, pattern = "^MT-")
    u <- subset(u, subset = percent.mt < pct_mito_cutoff)
    u <- subset(u,features = setdiff(rownames(u),c(exclude_genes,rownames(u)[grepl('^MT-',rownames(u))])))
    u <- NormalizeData(u)
    if(verbose){
        message('Dimension of processed data after filtering:\n', nrow(u),' ', ncol(u))
    }
    u@meta.data$sample <- sub(paste0(':.*'),'',colnames(u)) # default format

    # add Meta data
    cell_meta <- suppressMessages(u[[]] %>% select(sample) %>% left_join(sample_meta) %>% select(-sample))
    rownames(cell_meta) <- rownames(u[[]])
    u <- AddMetaData(u, cell_meta)

    u.ls <- SplitObject(u,split.by = 'sample')
    if(verbose){
        message(paste0('Number of unique samples: ',length(u.ls)))
    }
    for(i in 1:length(u.ls)){

        u.ls[[i]] <- FindVariableFeatures(u.ls[[i]],selection.method = "vst", nfeatures = num_features)
        if(verbose) message(i)
    }

    integrated.features <- SelectIntegrationFeatures(u.ls, nfeatures = num_features)
    u <- ScaleData(u,features = integrated.features,verbose = FALSE)
    u <- RunPCA(u,npcs = num_PCs,features = integrated.features, verbose = FALSE)
    if(verbose){
        message('Variables to be regressed out: ', paste(vars_to_regress,collapse = ','))
        message(paste0('Clustering resolution: ',resolution))
    }

    u <- RunHarmony(u, group.by.vars = vars_to_regress)
    u <- RunUMAP(u, reduction = "harmony", dims = 1:num_harmony)
    u <- FindNeighbors(u,reduction = "harmony", dims = 1:num_harmony)
    u <- FindClusters(u,resolution = resolution)

    # Make cluster id start from 1 ###
    u$seurat_clusters <- as.factor(as.numeric(as.character(u$seurat_clusters)) + 1)
    Idents(u) <- u$seurat_clusters

    ## 4. build phylogenetic tree
    if(verbose){
        message('=== Build Tree ===')
    }
    u <- BuildClusterTree(u,dims = 1:num_PCs) # 3 choices: feature (hvg), PCs, graph
    saveRDS(u,paste0(output_dir,'integrate.rds'))

    pdf(paste0(output_dir,'phylo_tree.pdf'))
    phylo_tree <- Seurat::Tool(object = u, slot = 'BuildClusterTree')
    if(sum(phylo_tree$tip.label=='0')>0) phylo_tree$tip.label <- as.character(as.integer(phylo_tree$tip.label) + 1) # start from 1 (leaf node)
    ape::plot.phylo(phylo_tree)
    ape::nodelabels()
    dev.off()

    ## Differential FindAllMarkers
    if(verbose){
        message('=== Find gene markers for each cell cluster ===')
    }
    markers <- FindAllMarkers(u)
    write.csv(markers,paste0(output_dir,'markers.csv'),row.names = F)
    top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    write.csv(top10,paste0(output_dir,'markers_top10.csv'),row.names = F)

    return(u)
}
