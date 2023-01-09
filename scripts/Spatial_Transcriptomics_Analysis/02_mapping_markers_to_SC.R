library(Seurat)
library(SeuratWrappers)
library(AUCell)
library(monocle3)
library(SummarizedExperiment)
library(GSVA)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(GeoMxWorkflows)

epi_subset <- readRDS("../../results/epi_subset.rds")
geneSets <- readRDS("../../results/lmm_main_labels_markers_list.rds")

#################################################################
##           Projecting Markers Gene Sets on SC Data           ##
#################################################################

##------------------------------------
##  Building cell ranking for AUCell
##------------------------------------

expr <- GetAssayData(epi_subset,
                     assay = "CorrectedCounts",
                     slot = "counts")
geneRanking <- AUCell_buildRankings(expr, nCores = 36, splitByBlocks = TRUE)

auc_score <- AUCell_calcAUC(rankings = geneRanking, geneSets = geneSets, aucMaxRank = nrow(geneRanking)*0.1)
auc_score <- as.data.frame(t(assay(auc_score)))
auc_score_scaled <- as.data.frame(scale(auc_score))
acinar_ductal <- AddMetaData(acinar_ductal, metadata = auc_score)
DimPlot(acinar_ductal, group.by = "DiseaseState", raster = T, cols = c("black","red","blue"), pt.size = 1) + ggtitle("Epithelial Population")
FeaturePlot(acinar_ductal, ncol = 3, pt.size = 0.5,
            features = c("panin_markers","adm_markers", "acinar_markers", "ductal_markers", "Glandular_Tumor_markers", "PoorlyDiffTumor_markers"),
            repel = TRUE, order = F, keep.scale = 'feature', col = rev(brewer.pal(n = 11, name = "RdBu")))




