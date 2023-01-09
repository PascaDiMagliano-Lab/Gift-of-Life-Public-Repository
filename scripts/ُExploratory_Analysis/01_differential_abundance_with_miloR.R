
##################################################################
##                  Loading libraries and data                  ##
##################################################################

suppressPackageStartupMessages({
  library(miloR)
  library(Seurat)
  library(tidyverse)
  library(SingleCellExperiment)
  library(RColorBrewer)
  source("../../src/miloR.R")
})


all_samples <- readRDS("../../results/rpca_integrated_annotated.rds")
colourCount <- length(unique(all_samples$cell_types))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

##################################################################
##                        MiloR Analysis                        ##
##################################################################


##-------------------------
##  Creating miloR object
##-------------------------

all_samples_sce <- as.SingleCellExperiment(all_samples)
reducedDim(all_samples_sce, "PCA") <- Embeddings(all_samples, "pca")
reducedDim(all_samples_sce, "UMAP") <- Embeddings(all_samples, "umap")
all_samples_milo <- Milo(all_samples_sce)


##--------------------
##  miloR processing
##--------------------

###  As per the authors supplementary notes here (https://www.nature.com/articles/s41587-021-01033-z),
###  choosing k is dependent on the dataset. They mentioned that it should be K >= S * 3-5 ,
###  where S is the number of samples. So we'll choose k = 30*4 = 120
###
###  They also state the follwoing "We recommend users to inspect the histogram of neighbourhood sizes after sampling of neighbourhoods
###  and to consider the number of cells that would be considered a “neighbourhood”in the dataset at hand.
###  As a heuristic for selecting a lower bound on K to increase the resolution of neighbourhoods for capturing rare sub-populations or states,
###  the user can select K such that the mean neighbourhood size is no more than 10% of the expected size of the rare population.
###  We provide the utility function plotNhoodSizeHist to visualize the neighbourhood size distribution as part of our R package."
###
###  As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples.
###
###  for the value of p, they state "we recommend initiating neighbourhood search withp= 0.05 for datasets with more than 100k cells
###  and p= 0.1 otherwise". so we'll pick up 0.1

k = 60 ## This is gave better mean nhood size value, around 80, which is around 10% of the endothelial cells (the rarest population here)
p =  0.1
traj_milo <- buildGraph(all_samples_milo, k = k, d = 30)
traj_milo <- makeNhoods(traj_milo, prop = p, k = k, d = 30, refined = TRUE)

traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample_name")
traj_design <- data.frame(colData(traj_milo))[,c("sample_name", "DiseaseState")] %>%
  `colnames<-`(c("Sample","Condition"))
traj_design <- distinct(traj_design) %>% `row.names<-`(.$Sample)
traj_milo <- calcNhoodDistance(traj_milo, d=30) ## this step takes AGES!
saveRDS(traj_milo, "../../results/traj_milo_k60.rds")

## checking mean nhood size
k60_hist <- plotNhoodSizeHist(traj_milo_k60)
k60_mean <- mean(k60_hist$data$nh_size)

##--------------------
##  Plotting Results
##--------------------

da_results60 <- DA_test(milo_obj = traj_milo_k60, samples_col = "sample_name", condition = "DiseaseState",
                         contrast = c("Tumor","Healthy"), annotation_col = "cell_types", plot = TRUE)
traj_milo_k60 <- buildNhoodGraph(traj_milo_k60)
DimPlot(all_samples, group.by = "DiseaseState", pt.size = 1, raster = T) + DimPlot(all_samples, group.by = "patient_id", raster = T)
plotNhoodGraphDA(traj_milo_k60, da_results60, alpha=0.1) + ggtitle("Differential Abudnance - K60")
plotDAbeeswarm(da_results60, group.by = "cell_types", alpha = 0.1) +
  ggtitle("Differential Abundance - K60")


