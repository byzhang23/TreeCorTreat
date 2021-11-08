TreeCorTreat: Tree-based correlation screen for phenotype-associated transcriptomic features and cell types
====

## Overview
Single-cell RNA-seq experiments with multiple samples are increasingly used to discover cell types and their molecular features that may influence samples’ phenotype (e.g. disease). However, analyzing and visualizing the complex cell type-phenotype association remains nontrivial. TreeCorTreat is an open source R package that tackles this problem by using a tree-based correlation screen to analyze and visualize the association between phenotype and transcriptomic features and cell types at multiple cell type resolution levels. With TreeCorTreat, one can conveniently explore and compare different feature types, phenotypic traits, analysis protocols and datasets, and evaluate the impacts of potential confounders. 

TreeCorTreat takes a gene expression matrix (raw count), cell-level metadata and sample-level metadata as input. It provides a whole pipeline to integrate data across samples, identify cell clusters and their hierarchical structure, evaluate the association between sample phenotype and cell type at different resolution levels in terms of both cell type proportion and gene expression, and summarize and visualize the results in a tree structured TreeCorTreat plot. This pipeline consists of six functional modules: 

* Module 1: Data integration
* Module 2: Define cell types at multiple resolutions 
* Module 3: Identify association between cell type proportion and sample phenotype
* Module 4: Identify association between global gene expression and sample phenotype
* Module 5: Identify differentially expressed genes
* Module 6: Visualization via TreeCorTreat plot

The modular structure provides users with the flexibility to skip certain analysis steps and replace them by users’ own data or analysis functions. 

## TreeCorTreat Installation

TreeCorTreat software can be installed via Github. Users should have R installed on their computer before installing TreeCorTreat. R version needs to be at least 3.6.1 or higher. R can be downloaded here: http://www.r-project.org/.

Users shall also install the following R packages before installing TreeCorTreat:

* [Seurat v3](https://satijalab.org/seurat/index.html)
* [harmony](https://github.com/immunogenomics/harmony)


**TreeCorTreat package is developed based on Seurat v3 (e.g. v3.2.2)**, because there are some reported issues remained to be resolved in the latest Seurat v4 release. In order to implement TreeCorTreat successfully, please make sure to install Seurat v3 release with the following command:

```{r}
remotes::install_version("Seurat", version = "3.X.X") # version = '3.2.2'
```

To install harmony from CRAN:
```{r}
install.packages("harmony")
```

To install the latest version of TreeCorTreat package via Github, run following commands in R:
```{r}
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("byzhang23/TreeCorTreat")
```
## Tutorials

For the user manual of TreeCorTreat, please refer to:
https://byzhang23.github.io/TreeCorTreat/TreeCorTreat.html 

If you are not familiar with R programming and want to explore how to load in datasets from csv/txt files, please check here: https://byzhang23.github.io/TreeCorTreat/DataPreparation.html

If you are interested in applying the TreeCorTreat plot in other (non-genomics) fields/data types, please refer to: 
https://byzhang23.github.io/TreeCorTreat/GeneralizedTreeCorTreatPlot.html 

## Citation
Please cite the following paper:

Boyang Zhang, Zhicheng Ji and Hongkai Ji. Tree-based Correlation Screen and Visualization for Exploring Phenotype-Cell Type Association in Multiple Sample Single-Cell RNA-Sequencing Experiments. bioRxiv 2021.10.27.466024; doi: https://doi.org/10.1101/2021.10.27.466024

## Contact

If you encounter any bugs or have any suggestions, please feel free to contact Boyang Zhang bzhang34@jhu.edu, or open an issue on the Github page: https://github.com/byzhang23/TreeCorTreat/issues.
