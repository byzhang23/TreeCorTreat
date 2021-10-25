#' Plot kernel density of cells in a UMAP embedding per group
#'
#' @param cell_meta     Cell-level metadata. Must contain `sample`,`UMAP_1` and `UMAP_2`columns.
#' @param sample_meta   Sample-level metadata. Must contain `sample` and variables to be plotted (e.g. row_variable and col_variable).
#' @param row_variable  Specify variable for stratifying on the row. Must be categorical variable (e.g. `factor` or `character` class).
#' @param col_variable  Specify variable for  stratifying on the column (optional). Must be categorical variable (e.g. `factor` or `character` class).
#' @param row_combined  An indicator for combining all subgroups on the rows and visualize the combined plot in the first row. Only works when both `row_variable` and `col_variable` exist.
#' @param font_size     Font size. Default is 15.
#' @param point_size    Set point size. Default is 0.1.
#' @param point_shape   Set point shape. Default is '.'.
#' @param point_alpha   Set point transparency. Default is 0.5.
#' @param breaks        A character vectors of breaks for color. Default is c(0,0.25,0.5,0.75,1).
#' @return A kernel density plot in UMAP embedding.
#' @export
#' @import ggplot2 dplyr
#' @importFrom MASS kde2d
#' @author Boyang Zhang <bzhang34@jhu.edu>, Hongkai Ji
#'
#' @examples
#' treecor_celldensityplot(cell_meta, sample_meta, row_variable = 'severity')


treecor_celldensityplot <- function(cell_meta,sample_meta,row_variable,col_variable = NULL,row_combined = T,font_size = 15,point_size = 0.1,point_shape = '.',point_alpha = 0.5,breaks = seq(0,1,0.25)){

    ## calculate density
    get_density <- function(x, y, ...) {
        dens <- MASS::kde2d(x, y, ...)
        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ii <- cbind(ix, iy)
        return(dens$z[ii])
    }

    # check UMAP_1 and UMAP_2 columns
    if(sum(colnames(cell_meta) %in% c('UMAP_1','UMAP_2'))==2){
        dat <- suppressMessages(inner_join(cell_meta,sample_meta))

        if(is.null(dat[,col_variable]) & class(dat[,row_variable]) %in% c('factor','character')){
            # one row
            for (i in unique(dat[,row_variable])) {
                id <- which(dat[,row_variable]==i)
                gd <- get_density(dat$UMAP_1[id], dat$UMAP_2[id], n = 100)
                dat[id,'Density'] <- gd/max(gd)
            }

            res <- ggplot(dat,aes(x = UMAP_1,y = UMAP_2, col = Density)) +
                geom_point(size = point_size, shape = point_shape, alpha = point_alpha) +
                theme_classic(base_size = font_size) +
                facet_grid(as.formula(paste0('~',row_variable))) +
                scale_colour_viridis_c(option = 'C',direction = 1,breaks = breaks) +
                theme(legend.position = 'bottom',
                      strip.text.y.right = element_text(angle = 0),
                      strip.text = element_text(size = font_size)) +
                guides(color = guide_colorbar(title.vjust = 1,
                                              barwidth = 12,
                                              barheight = 1.5))

        }else if(class(dat[,col_variable]) %in% c('factor','character') & class(dat[,row_variable]) %in% c('factor','character')){

            dat$row_col <- paste0(dat[,row_variable],':',dat[,col_variable])
            for (i in unique(dat$row_col)) {
                id <- which(dat$row_col==i)
                gd <- get_density(dat$UMAP_1[id], dat$UMAP_2[id], n = 100)
                dat[id,'Density'] <- gd/max(gd)
            }
            if(row_combined){
                dat2 <- dat
                for (i in unique(dat2[,col_variable])) {
                    id <- which(dat2[,col_variable]==i)
                    gd <- get_density(dat2$UMAP_1[id], dat2$UMAP_2[id], n = 100)
                    dat2[id,'Density'] <- gd/max(gd)
                }
                dat2[,row_variable] <- 'All'
                dat <- rbind(dat,dat2)

                res <- ggplot(dat,aes(x = UMAP_1,y = UMAP_2, col = Density)) +
                    geom_point(size = point_size, shape = point_shape, alpha = point_alpha) +
                    theme_classic(base_size = font_size) +
                    facet_wrap(as.formula(paste0('~',col_variable)),nrow = 1) +
                    facet_grid(as.formula(paste0(row_variable,'~',col_variable))) +
                    scale_colour_viridis_c(option = 'C',direction = 1,breaks = breaks) +
                    theme(legend.position = 'bottom',
                          strip.text.y.right = element_text(angle = 0),
                          strip.text = element_text(size = font_size)) +
                    guides(color = guide_colorbar(title.vjust = 1,
                                                  barwidth = 12,
                                                  barheight = 1.5))
            }else{
                res <- ggplot(dat,aes(x = UMAP_1,y = UMAP_2, col = Density)) +
                    geom_point(size = point_size, shape = point_shape, alpha = point_alpha) +
                    theme_classic(base_size = font_size) +
                    facet_grid(as.formula(paste0(row_variable,'~',col_variable))) +
                    scale_colour_viridis_c(option = 'C',direction = 1,breaks = breaks) +
                    theme(legend.position = 'bottom',
                          strip.text.y.right = element_text(angle = 0),
                          strip.text = element_text(size = font_size)) +
                    guides(color = guide_colorbar(title.vjust = 1,
                                                  barwidth = 12,
                                                  barheight = 1.5))
            }



        }else{
            stop('Current function only supports categorical variables (factor or character).')
        }

    }else{
        stop('Missing UMAP_1 and UMAP_2 columns.')
    }

    return(res)
}




