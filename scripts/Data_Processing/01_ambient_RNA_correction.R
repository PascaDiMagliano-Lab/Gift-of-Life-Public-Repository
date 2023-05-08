#################################################################
##                      loading libraries                      ##
#################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(parallel)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(gridExtra)
  library(tidyverse)
  source("../../src/correct_reads_function.r")
})


#################################################################
##                        loading files                        ##
#################################################################

## Organizing files into folders for each sample

dataDir <- "../../data/GSE229413_RAW/"
samples_names <- list.files(dataDir) %>%
  gsub(pattern = "_filtered_feature_bc_matrix.h5|_raw_feature_bc_matrix.h5|^GSM[0-9]*_", replacement = "") %>% unique()

sapply(samples_names, function(x){
  files <- list.files(dataDir, pattern = paste0(x,"_"), include.dirs = F)
  dir.create(file.path(dataDir,x), showWarnings = F, recursive = T)
  sapply(files, function(y){file.rename(from = file.path(dataDir,y), to = file.path(dataDir,x,y))})
})

samples_names <- list.files("../../data/GSE229413_RAW")
samples_dir <- list.files("../../data/GSE229413_RAW", full.names = T)

#################################################################
##           Ambient RNA Correction (w +10% and Cap)           ##
#################################################################

## Ambient RNA correction was done using soupX were the contamination fraction
## is estimated from empty droplets and negative markers. An extra 10% was added
## as per the authors recommendation to remove 95-98% of the soup (ambient RNA).
## https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
## while capping the contamination fraction to 20%

clust <- makeCluster(detectCores())
clusterExport(clust,
              c("correct_reads"),
              envir=environment())

corrected_samples <-
  parLapply(clust,
            samples_dir,
            fun = function(x){correct_reads(data_dir = x,  extra = 0.1, cap = TRUE, h5 = TRUE)}
  )

stopCluster(clust)

dir.create("objects", showWarnings = F)
names(corrected_samples) <- samples_names
saveRDS(object = corrected_samples, file = "../../results/corrected_samples_w_extra_10_cap.rds")
remove(corrected_samples)

