#' Calculate summary statistic and perform permutation test
#'
#' Compute summary statistic (e.g. canonical correlation, F-statistic, likelihood ratio test statistic) and obtain p-value based on permutation test.
#' @param x                 A vector or data frame.
#' @param y                 A vector or data frame.
#' @param method            Specify analysis_type: cancor, pearson, spearman, regression
#' @param num_permutations  Number of permutations (by default: 100).
#' @param alternative       A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param componen_id       An index to extract canonical correlation component (by default: 1).
#' @return A data frame of summary statistic and p-value
#' @export
#' @import dplyr
#' @importFrom combinat permn
#' @importFrom lmtest lrtest
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji

cancor_permute <- function(x, y, method, num_permutations,alternative, component_id = 1){

    ## canonical_value
    if(is.vector(y)){
        y <- data.frame(y)
    }
    convert_column <- which(!sapply(y,is.numeric))
    response_unique <- unique(y)
    rank_zero <- ifelse(is.null(ncol(y)), length(response_unique)==1,nrow(response_unique)==1)

    ## Check if all values are numeric
    if(length(convert_column)>0 & method %in% c('pearson','spearman','cancor')){
        # warning(paste0(paste(colnames(y)[convert_column],collapse = ','),' is not numeric. Proceed with factorized numeric values.'))
        for(j in convert_column){
            y[,j] <- as.numeric(as.factor(y[,j]))
        }
    }

    ## Check if rank is 0
    if(rank_zero){
        warning('y has rank 0. Set canonical correlatioin = NA.')
        obs.cor <- NA
        p.value <- NA

    }else{
        ### permutation test
        # (1) List all possible permutations
        if(nrow(y)<7){
            new.idx <- do.call(rbind,combinat::permn(1:nrow(y)))[-1,] # 1st row is the original data
            if(nrow(new.idx)>num_permutations) new.idx <- new.idx[1:num_permutations,]

        }else{
        # (2) Sampling
            new.idx <- t(sapply(1:num_permutations,function(b){
                set.seed(b)
                vec = sample(1:nrow(y),replace = F)
                if(identical(vec,1:nrow(y))) vec = NULL
                vec
            })) %>% unique
        }

        ### Methods for canonical or simple correlation
        if(method == 'cancor'){
            # choose componen_id
            obs.cor <- cancor(x,y)$cor[component_id]

            new.cor <- sapply(1:nrow(new.idx),function(i){
                new.y <- y[new.idx[i,],]
                cancor(x,new.y)$cor[component_id]
            })
            p.value <- ifelse(alternative %in% c('two.sided','greater'),sum(new.cor>=obs.cor),sum(new.cor<=obs.cor))/nrow(new.idx)

        }else if(method == 'pearson'){
            ## pearson
            new.cor <- sapply(1:nrow(new.idx),function(i){
                new.y <- y[new.idx[i,],]
                cor(x,new.y,method = 'p')
            })

            obs.cor <- cor(x,y,method = 'p')[1]
            p.value <- ifelse(alternative =='greater',sum(new.cor>=obs.cor),
                             ifelse(alternative == 'less',sum(new.cor<=obs.cor),
                                    sum(abs(new.cor)>=abs(obs.cor))))/nrow(new.idx)

        }else if(method == 'spearman'){
            ## spearman
            new.cor <- sapply(1:nrow(new.idx),function(i){
                new.y <- y[new.idx[i,],]
                cor(x,new.y,method = 's')
            })

            obs.cor <- cor(x,y,method = 's')[1]
            p.value <- ifelse(alternative =='greater',sum(new.cor>=obs.cor),
                              ifelse(alternative == 'less',sum(new.cor<=obs.cor),
                                     sum(abs(new.cor)>=abs(obs.cor))))/nrow(new.idx)

        }else if(method == 'regression'){
            ## regression (handle non-linear and non-ordinal categorical variable)
            nested <- lm(formula = as.formula('x~1')) # x: cell prop; y: phenotype
            full <- lm(data = y,formula = as.formula(paste0('x~',paste(colnames(y),collapse = '+'))))

            new.cor <- sapply(1:nrow(new.idx),function(i){
                new.y <- y[new.idx[i,],,drop = F]
                new.full <- lm(data = new.y,formula = as.formula(paste0('x~',paste(colnames(new.y),collapse = '+'))))
                lmtest::lrtest(nested,new.full)$Chisq[2]
            })

            obs <- lmtest::lrtest(nested,full)
            obs.cor <- obs$Chisq[2]
            p.value <- ifelse(alternative =='greater',sum(new.cor>=obs.cor),
                              ifelse(alternative == 'less',sum(new.cor<=obs.cor),
                                     sum(abs(new.cor)>=abs(obs.cor))))/nrow(new.idx)

            # one column
            if(ncol(y)==1){
                obs.cor <- ifelse(coef(full)[2]>0,obs.cor,-obs.cor)
            }
        }else if(method == 'regression_concat'){

            ## MANOVA
            null <- lm(x~1)
            fit <- lm(x~.,data = y)
            type <- ifelse(sum(names(anova(null,fit)) %in% 'F')>0,'anova','manova')
            obs.cor <- ifelse(type == 'anova',anova(null,fit)$`F`[2], anova(null,fit)$`approx F`[2])  #`Pillai`[2]

            new.cor <- sapply(1:nrow(new.idx),function(i){
                new.y <- y[new.idx[i,],,drop = F]
                new.fit <- lm(x~.,data = new.y)
                ifelse(type == 'anova',anova(null,new.fit)$`F`[2], anova(null,new.fit)$`approx F`[2]) #`Pillai`[2]
            })
            p.value <- ifelse(alternative =='greater',sum(new.cor>=obs.cor),
                              ifelse(alternative == 'less',sum(new.cor<=obs.cor),
                                     sum(abs(new.cor)>=abs(obs.cor))))/nrow(new.idx)
        }
    }
    res <- data.frame(cancor = obs.cor, p = p.value)
}
