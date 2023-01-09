
##################################################################
##                  Loading libraries and data                  ##
##################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(AUCell)
  library(msigdbr)
})

integrated_data <- readRDS('../../results/rpca_integrated_annotated.rds')
epi_subset <- subset(integrated_data, subset = cell_types %in% c("Acinar","Epithelial"))
h_gene_sets <- msigdbr(species = "human", category = "H")

## Hallmark KRAS UP and DOWN
KRAS_UP <- h_gene_sets %>% filter(gs_name == "HALLMARK_KRAS_SIGNALING_UP") %>% pull(human_gene_symbol)
KRAS_DN <- h_gene_sets %>% filter(gs_name == "HALLMARK_KRAS_SIGNALING_DN") %>% pull(human_gene_symbol)

# WashU paper gene sets (https://www.nature.com/articles/s41588-022-01157-1#Abs1)
washu_markers <- readxl::read_xlsx("../../data/washu_gene_list.xlsx") %>%
  as.list() %>%
  lapply(na.omit) %>%
  lapply(as.character)

# Panglao Db gene sets (https://panglaodb.se/markers.html?cell_type=%27all_cells%27)
markers <- read.table("../../data/PanglaoDB_markers_27_Mar_2020.tsv", sep = "\t", header = T)
PanglaoDB_acinar_markers <- markers[(markers$species == "Hs" | markers$species == "Mm Hs") &
                                      markers$cell.type == "Acinar cells", "official.gene.symbol"]
PanglaoDB_ductal_markers <- markers[(markers$species == "Hs" | markers$species == "Mm Hs") &
                                      markers$cell.type == "Ductal cells", "official.gene.symbol"]

GeneSets <- list("KRAS_UP" = KRAS_UP,
                 "KRAS_DN" = KRAS_DN,
                 "WashU_Acinar" = washu_markers$washu_acinar,
                 "WashU_Duct1" = washu_markers$washu_duct1,
                 "WashU_Duct2" = washu_markers$washu_duct2,
                 "WashU_PanIN" = washu_markers$washu_PanIN,
                 "WashU_Tumor" = washu_markers$washu_PDAC,
                 "PanglaoDB_Acinar" = PanglaoDB_acinar_markers,
                 "PanglaoDB_Ductal" = PanglaoDB_ductal_markers)

##################################################################
##                        Sub-clustering                        ##
##################################################################

DefaultAssay(epi_subset) <- "integrated"
epi_subset <- FindVariableFeatures(epi_subset, selection.method = "vst", nfeatures = 2000, assay = "CorrectedCounts")
epi_subset <- ScaleData(epi_subset, features = rownames(epi_subset),verbose = FALSE)
epi_subset <- RunPCA(epi_subset, verbose = TRUE)
ElbowPlot(epi_subset, 50)
epi_subset <- RunUMAP(epi_subset, reduction = "pca", dims = 1:20)
epi_subset <- FindNeighbors(epi_subset, reduction = "pca", dims = 1:20)
epi_subset <- FindClusters(epi_subset, resolution = 0.50)

##################################################################
##                       Plotting Results                       ##
##################################################################

DefaultAssay(epi_subset) <- "CorrectedCounts"

DimPlot(epi_subset, label = T, raster = T)
DimPlot(epi_subset, group.by = "DiseaseState",  raster = T)
DimPlot(epi_subset, group.by = "cell_types",  raster = T)
DimPlot(epi_subset, group.by = "lesion",  raster = T)
DimPlot(epi_subset, group.by = "sample_name", raster = T)

#################################################################
##                 Markers and Pathways Scores                 ##
#################################################################

expr <- GetAssayData(epi_subset,
                     assay = "CorrectedCounts",
                     slot = "counts")
geneRanking <- AUCell_buildRankings(expr, nCores = parallel::detectCores(), splitByBlocks = TRUE)
auc_score <- AUCell_calcAUC(rankings = geneRanking, geneSets = GeneSets, aucMaxRank = nrow(geneRanking)*0.1)
auc_score <- as.data.frame(t(assay(auc_score)))
auc_score_scaled <- as.data.frame(scale(auc_score))
epi_subset <- AddMetaData(epi_subset, metadata = auc_score)
FeaturePlot(epi_subset,features = names(GeneSets), cols = rev(brewer.pal(n = 11, name = "RdBu")))

saveRDS(epi_subset, "../../results/epi_subset.rds")