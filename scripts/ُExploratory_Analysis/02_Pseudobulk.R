
##################################################################
##                  Loading Data and Libraries                  ##
##################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(EnhancedVolcano)
  library(enrichplot)
  library(cowplot)
  source("../../src/pseudo_bulk.R")
})

all_samples <- readRDS("../../results/rpca_integrated_annotated.rds")
t_cells <- subset(all_samples, subset = cell_types == "T-cells")
fibroblasts <- subset(all_samples, subset = cell_types == "Fibroblast")
myeloid <- subset(all_samples, subset = cell_types %in% c("Myeloid","Granulocytes"))

##################################################################
##                Pseudobulk for each population                ##
##################################################################

all_cells_DGE <- pseudobulk_dge(seurat_object = all_samples, assay = "CorrectedCounts",
                              aggregate_by = "sample_name", comparison = "DiseaseState")
t_cells_DGE <- pseudobulk_dge(seurat_object = t_cells, assay = "CorrectedCounts",
                              aggregate_by = "sample_name", comparison = "DiseaseState")
fibroblasts_DGE <- pseudobulk_dge(seurat_object = fibroblasts, assay = "CorrectedCounts",
                              aggregate_by = "sample_name", comparison = "DiseaseState")
myeloid_DGE <- pseudobulk_dge(seurat_object = myeloid, assay = "CorrectedCounts",
                              aggregate_by = "sample_name", comparison = "DiseaseState")


#################################################################
##                          PCA Plots                          ##
#################################################################

all_pca <- pca_plot(all_cells_DGE, groups = "DiseaseState") + ggtitle("All cells")
tcell_pca <- pca_plot(t_cells_DGE, groups = "DiseaseState") + ggtitle("T-Cells Population")
fibro_pca <- pca_plot(fibroblasts_DGE, groups = "DiseaseState") + ggtitle("Fibroblasts Population")
myeloid_pca <- pca_plot(myeloid_DGE, groups = "DiseaseState") + ggtitle("Myeloid Population")
plot_grid(plotlist = list(all_pca, tcell_pca, fibro_pca, myeloid_pca), nrow = 2, ncol = 2)

##################################################################
##                     Correlation heatmaps                     ##
##################################################################

all_corr <- corr_heatmap(all_cells_DGE, groups = "DiseaseState", title = "Pseudobulk Correlation (All Cells)")$gtable
tcell_corr <- corr_heatmap(t_cells_DGE, groups = "DiseaseState", title = "Pseudobulk Correlation (T Cells)")$gtable
fibro_corr <- corr_heatmap(fibroblasts_DGE, groups = "DiseaseState", title = "Pseudobulk Correlation (Fibroblasts)")$gtable
myeloid_corr <- corr_heatmap(myeloid_DGE, groups = "DiseaseState", title = "Pseudobulk Correlation (Myeloid Cells)")$gtable
plot_grid(plotlist = list(all_corr, tcell_corr, fibro_corr, myeloid_corr))

##################################################################
##                Tumor vs. Healthy DGE Analysis                ##
##################################################################

TvH_Tcells_DEGs <- pairwise_dge(t_cells_DGE, group1 = "Tumor", group2 = "Healthy")
TvH_Tcells_GO <- GO_enrich(DGE_results = TvH_Tcells_DEGs, fcCutOff = 2, pvalCutOff = 0.05, GO_category = "BP", top = 20)

EnhancedVolcano(TvH_Tcells_DEGs, rownames(TvH_Tcells_DEGs), "log2FoldChange", "padj", title = "T-Cells", subtitle = NULL )
GOA_plot(TvH_Tcells_GO)
heatplot(TvH_Tcells_GO$upreg_GO,
         foldChange= TvH_Tcells_DEGs %>% pull(log2FoldChange) %>% `names<-`(rownames(TvH_Tcells_DEGs)),
         showCategory= 20)
termsim <- pairwise_termsim(TvH_Tcells_GO$upreg_GO)
treeplot(termsim, 50)

