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
  source("../../src/correct_reads_function.r")
})


#################################################################
##                        loading files                        ##
#################################################################

samples_names <- list.files("../../data/samples/")
samples_dir <- list.files("../../data/samples/", full.names = T)


#################################################################
##           Ambient RNA Correction (w +10% and Cap)           ##
#################################################################

clust <- makeCluster(detectCores())
clusterExport(clust,
              c("correct_reads"),
              envir=environment())

corrected_samples <-
  parLapply(clust,
            samples_dir,
            fun = function(x){correct_reads(data_dir = x,  extra = 0.1, cap = TRUE)}
  )

stopCluster(clust)

dir.create("objects", showWarnings = F)
names(corrected_samples) <- samples_names
saveRDS(object = corrected_samples, file = "../../results/corrected_samples_w_extra_10_cap.rds")
remove(corrected_samples)

