#' Sample PCA plot
#'
#' @param sample_meta   Sample-level metadata. Must contain \code{sample} and variables to be plotted (e.g. \code{row_variable} and \code{col_variable}).
#' @param pca_matrix    PCA matrix with sample names on the rows and PC coordinates on the columns (e.g. PC1 and PC2).
#' @param response_variable A univariate phenotype. If multivariate phenotype is provided, either convert them to a univariate combined phenotype (e.g. using weight or PCA) or call this function separately for each outcome.
#' @param analysis_type Specify summary statistic, can be one of 'cancor' (default) or 'regression'.
#' @param num_permutations Number of permutations. Default is 1000.
#' @param alternative   Specify alternative hypothesis. Must be one of 'two.sided' (default) or 'less' or 'greater'.
#' @param row_variable  Specify variable for row. Must be categorical variable (e.g. `factor` or `character` class).
#' @param col_variable  Specify variable for column (optional). Must be categorical variable (e.g. `factor` or `character` class).
#' @param font_size     Font size. Default is 15.
#' @param point_size    Set point size. Default is 0.1.
#' @param point_shape   Set point shape. Default is 19.
#' @param point_alpha   Set point transparency. Default is 1.
#' @param line_type     Set line type as 'dashed' (default) or 'solid' or others.
#' @param line_color    Set optimal axis line color.
#' @return Sample pca plot.
#' @export
#' @import ggplot2 dplyr
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#' @examples
#' treecor_samplepcaplot(sample_meta, pca_matrix, response_variable = 'severity', row_variable = 'severity')

treecor_samplepcaplot <- function(sample_meta,pca_matrix,response_variable, analysis_type = 'cancor', num_permutations = 1000, alternative = 'two.sided', row_variable = NULL,col_variable = NULL,font_size = 15,point_size = 1,point_shape = 19,point_alpha = 1,line_type = 'dashed',line_color = 'black'){

    if(sum(colnames(pca_matrix) %in% c('PC1','PC2'))!=2){
        stop('PCA matrix does not contain PC1 and PC2 columns.')
    }
    if(analysis_type == 'regression'){
        analysis_type <- 'regression_concat'
    }

    if(length(response_variable)>1){
        stop('There are >1 phenotypes. Current PCA plot only supports univariate phenotype (can either use weight to convert multivariate phenotype to a univariate phenotype or call this function multiple times.)')

    }else if(analysis_type=='cancor'){
        tmp_meta <- sample_meta
        if(class(sample_meta[,response_variable]) %in% c('character','factor')){
            warning('Convert `response_variable` to factor.')
            tmp_meta[,response_variable] <- as.integer(as.factor(tmp_meta[,response_variable]))
        }

        mat <- pca_matrix[,c('PC1','PC2')]
        dat <- suppressMessages(inner_join(data.frame(sample = rownames(mat),PC1=mat[,1],PC2=mat[,2]),
                                           tmp_meta))
        cancor <- cancor(dat[,c('PC1','PC2'),drop=F],dat[,response_variable,drop=F])
        coef <- cancor$xcoef
        slope <- coef[2,1]/coef[1,1]
        Cancor <- cancor_permute(x = dat[,c('PC1','PC2'),drop=F],y = dat[,response_variable,drop=F],
                                method = 'cancor',num_permutations = num_permutations,alternative = alternative)

        dat2 <- suppressMessages(inner_join(data.frame(sample = rownames(mat),PC1=mat[,1],PC2=mat[,2]),
                                            sample_meta))

        p <- ggplot(dat2,aes_string(x='PC1',y='PC2',color = response_variable)) +
            geom_point(size = point_size, shape = point_shape, alpha = point_alpha) +
            geom_abline(slope = slope,linetype = line_type,color = line_color) +
            labs(title = paste0('Canonical correlation: ',round(Cancor[1],3),' (P-value: ',ifelse(Cancor[2]==0,paste0('<',1/num_permutations),format(Cancor[2],scientific =T,digits = 2)),')')) +
            theme_classic(base_size = font_size) +
            theme(legend.position = 'bottom',
                  strip.text.y.right = element_text(angle = 0),
                  strip.text = element_text(size = font_size))

    }else if(analysis_type == 'regression_concat'){
        mat <- pca_matrix[,c('PC1','PC2')]
        dat <- suppressMessages(inner_join(data.frame(sample = rownames(mat),PC1=mat[,1],PC2=mat[,2]),
                                           sample_meta))
        Fstat <- cancor_permute(x = as.matrix(dat[,c('PC1','PC2'),drop=F]),y = dat[,response_variable,drop=F],
                                method = 'regression_concat',num_permutations = num_permutations,alternative = alternative)

        p <- ggplot(dat,aes_string(x='PC1',y='PC2',color = response_variable)) +
            geom_point(size = point_size, shape = point_shape, alpha = point_alpha) +
            #  geom_abline(slope = slope,linetype = line_type,color = line_color) +
            labs(title = paste0('F-statistics: ',round(Fstat[1],3),' (P-value: ',ifelse(Fstat[2]==0,paste0('<',format(1/num_permutations,scentific = T, digits = 2)),format(Fstat[2],scientific =T,digits = 2)),')')) +
            theme_classic(base_size = font_size) +
            theme(legend.position = 'bottom',
                  strip.text.y.right = element_text(angle = 0),
                  strip.text = element_text(size = font_size))
    }


    if(!is.null(row_variable)){
        if(is.null(col_variable)){
            p <- p + facet_grid(as.formula(paste0('~',row_variable)))
        }else{
            p <- p + facet_grid(as.formula(paste0(row_variable,'~',col_variable)))
        }
    }
    return(p)
}
