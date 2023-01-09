##################################################################
##                  loading libraries and data                  ##
##################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
})

`%nin%` <- Negate(`%in%`)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

integrated_data <- readRDS("../../results/rpca_integrated_samples.rds")
DefaultAssay(integrated_data) <- "CorrectedCounts"

##-----------------------------------------
##  Finding Markers for different cluster
##-----------------------------------------

markers <- FindAllMarkers(integrated_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "markers.csv")

##-----------------------------------------
##  FeaturePlots for known markers
##-----------------------------------------

DimPlot(integrated_data, label = T, cols = getPalette(length(levels(integrated_data))),
        label.box = T, repel = T, raster = T) + ggtitle("Seurat Clusters")
DimPlot(integrated_data, label = F, group.by = "DiseaseState", raster = T, order = T) + ggtitle("DiseaseState")
DimPlot(integrated_data, label = F, group.by = "sample_name", repel = T, raster = T, order = T) + ggtitle("Sample Name")

DefaultAssay(integrated_data) <- "CorrectedCounts"
features <- c("SPINK1", "PRSS3", "CTRB1", "PRSS1", "AMY2A",
              "EPCAM","KRT19", "SOX9", "MUC1", "KRT18",
              "CDH5", "VWF",
              "PDGFRB", "PDGFRA", "PDPN", "ACTA2","TAGLN", "DPT", "RGS5",
              "PTPRC", "CD14", "ITGAM", "MARCO", "APOE", "C1QA","CD68",
              "CD3E", "CD4","CD8A",
              "CD19", "CD79A",'MS4A1',
              "HBA2", "HBA1",
              "MKI67",
              "INS",
              'NCAM1', 'NEGR1', 'NRN1')

DotPlot(integrated_data, features = features, assay = 'CorrectedCounts') + coord_flip()

#################################################################
##                      Renaming Clusters                      ##
#################################################################

cluster_annotation <- readxl::read_xlsx("../../data/cluster_annotation.xlsx")
cell_type <- cluster_annotation$cell_type
names(cell_type) <- cluster_annotation$cluster
integrated_data <- RenameIdents(integrated_data, cell_type)
integrated_data$cell_types <- Idents(integrated_data)

##-----------------------------
##  Plotting renamed clusters
##-----------------------------

DimPlot(integrated_data, group.by = "seurat_clusters",order = T, raster = T, label= T)
DimPlot(integrated_data, label = T, group.by = "cell_types", cols = getPalette(colourCount), order = T, raster = T)
DimPlot(integrated_data, group.by = "DiseaseState", order = T, raster = T)
DimPlot(integrated_data, group.by = "sample_name", order = T, raster = T)

integrated_data$cell_types <- factor(integrated_data$cell_types,
                                     levels = c("Acinar","Acinar_2","Ductal",
                                                "Endothelial","Fibroblast",
                                                "T-cells","B-cells","Myeloid","Granulocytes",
                                                "Mast","Cycling","Endocrine"))

DotPlot(integrated_data,
        assay = "CorrectedCounts",
        features = c("SPINK1", "CTRB1", "PRSS1", "AMY2A",
                     "KRT19", "SOX9", "MUC1", "KRT18","MUC2",
                     "CDH5", "VWF",
                     "PDGFRB", "PDGFRA", "PDPN", "ACTA2","TAGLN",
                     "PTPRC", "CD14", "APOE", "C1QA", "CD68", "HLA-DRA",
                     "CD3E", "CD4", "FOXP3", "CD8A", "GZMB",
                     "CD19", "CD79A", "IGJ", 'MS4A1',
                     "TPSAB1", "CPA3",
                     "S100A8", "CSF3R","G0S2",
                     "MKI67","INS"),
        group.by = "cell_types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "bold"), axis.text.y = element_text(size = 8, face = "bold")) +
  coord_flip() + ggtitle("Corrected Counts")

metadata <- integrated_data@meta.data[,c("DiseaseState","lesion")]
metadata$lesion <- gsub("PanIN","Tumor",metadata$lesion)
metadata[metadata$DiseaseState == "AdjNormal","lesion"] <- "AdjNorm"
metadata[is.na(metadata$lesion),"lesion"] <- "Lesion"
integrated_data$lesion <- metadata$lesion

# Saving results
saveRDS(integrated_data, "../../results/rpca_integrated_annotated.rds")

#################################################################
##         Plotting cell type abundance in each sample         ##
#################################################################

#splitting samples
samples.list <- unique(integrated_data$sample_name)
clusters <- lapply(samples.list, function(x){
  subset <- subset(integrated_data, subset = sample_name == x)
  dist <- data.frame(table(subset$cell_types))

  return(dist)
})

names(clusters) <- samples.list

#calculate relative freq (fractions) of each cell type
clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

#making things ggplot-friendly!
clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("cell_type","variable","value","sample")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("cell_type","variable","value","sample")

#calculating Shannon Entropy as index for diversity
# entropy <- reshape::melt(lapply(clusters, FUN = function(x){
#   Entropy(x$Freq)
# }))
# entropy$value <- scale(entropy$value)

#plotting
ggplot(clusters_dist, aes(fill=cell_type, y = value, x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("sample") + ggtitle("Cell Types Abundance") +
  scale_fill_manual(values = getPalette(colourCount))
ggplot(clusters_percent_dist, aes(fill=cell_type, y = value, x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("sample") + ggtitle("Relative Cell Types Abundance") +
  scale_fill_manual(values = getPalette(colourCount))

# ggplot(entropy, aes(x= L1, y = value)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_text(angle = 45, size = 7, face = "bold", hjust = 1)) +
#   ylab("count") + xlab("sample") + ggtitle("Shannon Diversity Index")
