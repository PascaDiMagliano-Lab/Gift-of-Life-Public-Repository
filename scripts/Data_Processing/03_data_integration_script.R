##################################################################
##                  loading libraries and data                  ##
##################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

merged_samples <- readRDS('../../results/merged_samples_fil_correctedcounts_w_extra_10.rds')
## Samples 4450-EC and 4451-EC will be discarded from the analysis due to quality issues
merged_samples <- merged_samples[,!(merged_samples$sample_name %in% c("4450-EC","4451-EC"))] 
DefaultAssay(merged_samples) <- 'CorrectedCounts'

##################################################################
##                       rPCA Integration                       ##
##################################################################

samples.list <- SplitObject(merged_samples, split.by = "sample_name")
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = samples.list)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = samples.list, anchor.features = features, reduction = "rpca")
rpca_integrated_by_samples <- IntegrateData(anchorset = anchors)
DefaultAssay(rpca_integrated_by_samples) <- "integrated"
rpca_integrated_by_samples <- ScaleData(rpca_integrated_by_samples, verbose = FALSE)
rpca_integrated_by_samples <- RunPCA(rpca_integrated_by_samples, verbose = FALSE)
rpca_integrated_by_samples <- RunUMAP(rpca_integrated_by_samples, reduction = "pca", dims = 1:30)
rpca_integrated_by_samples <- FindNeighbors(rpca_integrated_by_samples, reduction = "pca", dims = 1:30)
rpca_integrated_by_samples <- FindClusters(rpca_integrated_by_samples, resolution = 0.5)

saveRDS(rpca_integrated_by_samples, "../../results/rpca_integrated_samples.rds")