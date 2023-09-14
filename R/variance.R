# data: https://cellxgene.cziscience.com/collections/8e880741-bf9a-4c8e-9227-934204631d2a
# paper: https://doi.org/10.1016/j.isci.2022.104097

# title: High-resolution Slide-seqV2 spatial transcriptomics enables discovery 
# of disease-specific cell neighborhoods and pathways

#########################
library(Seurat)
library(ggplot2)
library(gghighlight)
library(patchwork)
library(dplyr)

library(gstat)
library(sf)

library(magrittr)
library(matrixStats)
library(scater)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(Voyager)

library(ggforce)


library(purrr)
library(furrr)
library(stringr)
library(tidyr)

library(rlang)

theme_set(theme_bw())
wd <- fs::path('/Users/kaylajackson/Library/CloudStorage/OneDrive-CaliforniaInstituteofTechnology/Projects/spatial_normalization/data/slide_seq/Marshall_2022/raw')
obj <- readRDS(fs::path(wd, 'puck_19112_05.rds'))

# create sfe object
m <- GetAssayData(obj, slot = "data", assay = 'RNA')
meta <- obj@meta.data


centroids <- FetchData(object = obj, vars = c("SPATIAL_1", "SPATIAL_2")) %>%
  tibble::rownames_to_column(var = 'barcode')


sfe <- SpatialFeatureExperiment(
  assays = list(counts = m),
  colData = meta,
  spatialCoords = as.matrix(centroids[,c("SPATIAL_1", "SPATIAL_2")])
)



ct = types[5]
for (ct in types[-c(1:3)]){
  # gene_inds <- rowSums(counts(sfe)[,colData(sfe)$cell_type == ct]) > 100
  sfe_sub <- sfe[, colData(sfe)$cell_type == ct]
  sfe_sub <- sfe_sub[rowSums(counts(sfe_sub)) > 100, ]
  
  
  dims <- dim(sfe_sub)
  
  if (any(map_lgl(dims, is_zero))) next 
  
  ct_sf <- colGeometry(sfe_sub, 'centroids')
  ct_sf  %<>% mutate(data.frame(t(counts(sfe_sub))))
  ct_sf <- st_set_crs(ct_sf, 3857)
  
  if (!exists(sym('vgm_'))) {
    print('computing variogram for boundary side-effect')
    print(ct)
    
    gene <- sym(colnames(ct_sf)[2])
    # vgm_ <- inject(gstat::variogram(!!gene ~ 1, ct_sf, cloud=FALSE))
    bnds <- attr(vgm_, 'boundaries')
    
    nb <- length(bnds)
    
    # distance matrix
    dis <- st_distance(ct_sf)
    dis <- as(dis, 'matrix')
    dis[lower.tri(dis)] <- 0
    
    # dis <- as(dis, "dgCMatrix")
    
    intervals <- purrr::map2(bnds[-nb], bnds[-1], c)
    where <- purrr::map(intervals, ~ which(dis>.x[1] & dis<=.x[2], arr.ind=TRUE))
    mean_dis <- purrr::map_dbl(intervals, ~ mean(dis[dis>.x[1] & dis<=.x[2]]))
  }
  
  # Generate plots for each gene
  pdf(paste0(ct,'_mixed_log.pdf'), width=15, height=7.5)
  for (gene in rownames(sfe_sub)) {
    
    gene <- sym(gene)
    
    spatial <- plotSpatialFeature(
      sfe_sub, features = as_string(gene), exprs_values = 'counts',
      ) +
      guides(
        color = guide_legend(
          override.aes = list(size = 2)
      )) + 
      labs(title = glue::glue("{as_string(gene)} log counts")) +
      expand_limits(
        y = c(1000,5000),
        x = c(1000, 5000)
      ) + 
      theme(
        #aspect.ratio = 1, 
        # legend.position = 'top',
        legend.title = element_blank(),
  
        legend.margin=margin(),
        legend.key = element_rect(colour = 'black', fill = 'white', size = 0.5),
        
        axis.title = element_text(size = rel(1), face = 'bold'),
        axis.text = element_text(size = rel(0.85), face = 'bold'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        
        panel.border = element_rect(size=1),
        panel.background = element_rect(fill = 'azure4'),
  
        plot.margin = margin(b=20),
        plot.title = element_text(vjust = -5, face = 'bold', size = rel(0.85))
    
      )
    
    spatial_ng <- plotSpatialFeature(
      sfe_sub, features = c('nGenes'), exprs_values = 'counts',
      ) +
      guides(
        color = guide_legend(
          override.aes = list(size = 2)
        )) + 
      expand_limits(
        y = c(1000,5000),
        x = c(1000, 5000)
      ) + 
      theme(
        plot.title = element_text(vjust=-10),
        #legend.title = element_blank(),
        
        legend.margin=margin(),
        legend.key = element_rect(colour = 'black', fill = 'white', size = 0.5),
        
        axis.title = element_text(size = rel(1), face = 'bold'),
        axis.text = element_text(size = rel(0.85), face = 'bold'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        
        panel.border = element_rect(size=1),
        panel.background = element_rect(fill = 'azure4'),
        
        plot.margin = margin(),
        
      )
    
    spatial_nc <- plotSpatialFeature(
      sfe_sub, features = c('nCounts_100'), exprs_values = 'counts',
    ) +
      guides(
        color = guide_legend(
          override.aes = list(size = 2)
        )) + 
      expand_limits(
        y = c(1000,5000),
        x = c(1000, 5000)
        ) + 
      theme(
        #aspect.ratio = 1, 
        # legend.position = 'top',
        plot.title = element_text(vjust=-10),
        #legend.title = element_blank(),
        
        legend.margin=margin(),
        legend.key = element_rect(colour = 'black', fill = 'white', size = 0.5),
        
        axis.title = element_text(size = rel(1), face = 'bold'),
        axis.text = element_text(size = rel(0.85), face = 'bold'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        
        panel.border = element_rect(size=1),
        panel.background = element_rect(fill = 'azure4'),
        
        plot.margin = margin(),
        
      )
    
    res <- purrr::map_dfr(
      where, function(x) {
        inds <- pull2(x, c(row,col))
        v1 <- pull(ct_sf, gene)[ inds[[1]] ]
        v2 <- pull(ct_sf, gene)[ inds[[2]] ]
                            
        mn <- (v1+v2)/2
        df2 <- (v1-v2)^2
                        
        rt <- mn/df2
        
        c(mean  = mean(mn), diff_sq = mean(df2), 
          gamma = mean(df2)/2, cv = mean(rt[!is.infinite(rt)], na.rm=TRUE))
                            
                            
      }
    )
    
    res <- res %>% 
      mutate(dist = mean_dis, ind = row_number()) %>%
      pivot_longer(cols = -c('ind','dist'), names_to = 'var')
    
    facet <- ggplot(res, aes(dist, value)) + 
      geom_point(
        color='black', 
        fill='#bdbdbd',
        size = rel(2),
        shape=21,
      ) +
      facet_wrap(~ var, scales = 'free_y') + 
      expand_limits(y = 0) + 
      labs(
        x = "Distance",
        y = element_blank(),
      ) +
      theme(
        #aspect.ratio = 1, 
        plot.title = element_text(size = rel(1.5), face = 'bold'),
        
        axis.title = element_text(size = rel(1), face = 'bold'),
        axis.text = element_text(size = rel(0.85), face = 'bold'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        
        plot.margin = margin(r=10, b=10),
        #plot.background = element_rect(fill = 'blue'),
        
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1),
        
        strip.text = element_text(size = rel(1), face = 'bold'),
        strip.background = element_rect(size=1)
        
      )
    
    
    print(
      wrap_plots(
        facet, spatial, spatial_nc, spatial_ng, nrow=2, guides = 'keep'
        ) +
        plot_annotation(title = ct)
    )
    
    
  }
 dev.off()
 
 rm(vgm_, intervals, bnds, where, dis, mean_dis, nb)
}

  
# row variance
rowData(sfe_sub)$gene_vc <- rowVars(counts(sfe_sub))

# KNN graph 
colGraph(sfe_sub, "knn10") <- findSpatialNeighbors(
  sfe_sub, method = "knearneigh", dist_type = "idw", 
  k = 10, style = "W")

rowData(sfe_sub)$moran <- calculateMoransI(sfe_sub, exprs_values = 'counts')
