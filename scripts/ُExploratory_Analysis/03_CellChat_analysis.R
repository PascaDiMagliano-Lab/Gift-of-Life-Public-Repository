
# loading libraries -------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(magrittr)
  library(patchwork)
  options(stringsAsFactors = FALSE)
  future::plan("multisession", workers = parallel::detectCores())
})


# loading data ------------------------------------------------------------

integrated_data <- readRDS("../../results/rpca_integrated_annotated.rds")
CellChatDB <- CellChatDB.human

# Creating tumor and healthy subsets --------------------------------------

integratedData.list <- SplitObject(integrated_data, split.by = "DiseaseState")
remove(integrated_data)

# Creating CellChat Object ------------------------------------------------

assayUsed <- 'CorrectedCounts'
CellChat.list <- lapply(integratedData.list, function(x){
  CellChat_obj <- createCellChat(object = GetAssayData(x, assay = assayUsed, slot = "data"),
                                 meta = x@meta.data,
                                 group.by = "cell_types")
})
remove(integratedData.list)

# Running CellChat --------------------------------------------------------

CellChat.list <- lapply(CellChat.list, function(x){
  x@DB <- CellChatDB
  cellchat <- subsetData(x)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat@meta$cell_types <- factor(cellchat@meta$cell_types)
  cellchat@meta$cell_types <- droplevels(cellchat@meta$cell_types,
                                         exclude = setdiff(levels(cellchat@meta$cell_types),unique(as.character(cellchat@meta$cell_types))))
  cellchat@idents <- cellchat@meta$cell_types
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
})

## Splitting objects
healthy <- CellChat.list$Healthy
tumor <- CellChat.list$Tumor
adjnorm <- CellChat.list$AdjNormal

# Interaction Weight / Strength -------------------------------------------

## cell-cell communication strength tables
healthy@net$weight %>% write.csv("../../results/healthy_cell_interaction_weights.csv")
tumor@net$weight %>% write.csv("../../results/tumor_cell_interaction_weights.csv")
adjnorm@net$weight %>% write.csv("../../results/adjnorm_cell_interaction_weights.csv")

## Plotting communication strength

par(mfrow = c(1,3), xpd=TRUE)
netVisual_circle(healthy@net$weight, vertex.weight = as.numeric(table(healthy@idents)), arrow.size = 0.1,
                 weight.scale = T, label.edge= F, title.name = "Healthy interaction strength")
netVisual_circle(tumor@net$weight, vertex.weight = as.numeric(table(tumor@idents)), arrow.size = 0.1,
                 weight.scale = T, label.edge= F, title.name = "Tumor interaction strength")
netVisual_circle(adjnorm@net$weight, vertex.weight = as.numeric(table(adjnorm@idents)), arrow.size = 0.1,
                 weight.scale = T, label.edge= F, title.name = "AdjNorm interaction strength")

netVisual_heatmap(healthy@net$weight, vertex.weight = as.numeric(table(healthy@idents)), arrow.size = 0.1,
                  weight.scale = T, label.edge= F, title.name = "Healthy interaction strength")
netVisual_circle(tumor@net$weight, vertex.weight = as.numeric(table(tumor@idents)), arrow.size = 0.1,
                 weight.scale = T, label.edge= F, title.name = "Tumor interaction strength")
netVisual_circle(adjnorm@net$weight, vertex.weight = as.numeric(table(adjnorm@idents)), arrow.size = 0.1,
                 weight.scale = T, label.edge= F, title.name = "AdjNorm interaction strength")


# Compare the total number of interactions and interaction strength -------

## Data Processing
CellChat.list <- lapply(CellChat.list, netAnalysis_computeCentrality)
group.new <- levels(CellChat.list[[2]]@idents)
CellChat.list[[1]] <- liftCellChat(CellChat.list[[1]], group.new)
CellChat <- mergeCellChat(CellChat.list[c(1,2)], add.names = names(CellChat.list)[c(1,2)], cell.prefix = TRUE)
CellChat@options$datatype <- "RNA"

## Plots
### Circos plots
gg1 <- compareInteractions(CellChat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(CellChat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

### Heatmaps
netVisual_heatmap(CellChat)
netVisual_heatmap(CellChat, measure = "weight", sources.use = seq(1:17)[-c(16)], targets.use = seq(1:17)[-c(16)], remove.isolate = T)

### Differential Circos plots
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(CellChat, sources.use = c(7,8,10), targets.use = seq(1:17)[-c(16)], remove.isolate = T, weight.scale = T, measure = "weight", title.name = '')
netVisual_diffInteraction(CellChat, sources.use = seq(1:17)[-c(7,8,10,16)], targets.use = seq(1:17)[-c(16)], weight.scale = T, measure = "weight", remove.isolate = T, title.name = '')

