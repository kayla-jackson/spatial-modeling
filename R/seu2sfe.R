library(fs)
library(SpatialExperiment)
library(SpatialFeatureExperiment)

library(AnnotationHub)
library(dplyr)

data_dir <- "../data/slide_seq/Marshall_2022/raw"
obj <- readRDS(fs::path(data_dir, 'puck_191112_05.rds'))

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

# Update gene names
ah = AnnotationHub()

query(ah, c("EnsDb.Mmusculus", "108"))
edb_108 <- ah[["AH109367"]]

annots <- select(edb_108, keys = rownames(sfe),
                 columns = "SYMBOL", keytype = 'GENEID')

matches <- match(rownames(sfe), annots$GENEID)
genes <- annots$SYMBOL[matches]

# some gene ids are retired in most up-to-date ENSEMBL database,
# so reset manually
na_genes <- which(is.na(genes))
genes[na_genes] <- rownames(sfe)[na_genes]

rowData(sfe)$ensembl_id <- rownames(sfe)
rownames(sfe) <- genes

nGenes <- dim(sfe)[1]

mt <- grepl('^mt-', rownames(sfe), ignore.case=T)
colData(sfe)$prop_mito <- colSums(counts(sfe)[mt,])/colSums(counts(sfe))
colData(sfe)$nCounts <- colSums(counts(sfe))

# number of genes expressed/cell
colData(sfe)$nGenes <- colSums(counts(sfe) > 0)

# number of cells that gene is expressed in
rowData(sfe)$nCellsExpr <- rowSums(counts(sfe) > 0)


saveRDS(sfe, fs::path(data_dir,'puck_191112_05-sfe.Rds'))
