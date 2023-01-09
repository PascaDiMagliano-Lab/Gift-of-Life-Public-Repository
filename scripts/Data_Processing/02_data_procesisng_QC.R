#################################################################
##                      Loading Libraries                      ##
#################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

#################################################################
##                 Loading and Processing Data                 ##
#################################################################

## Imporitng data and samples info
corrected_samples <- readRDS("../../results/corrected_samples_w_extra_10_cap.rds")
samples_info <- read.csv("../../data/samples_info.csv")

## Cleaning up some data
samples_info$patient_id <- gsub("AdjNormal$|Tumor$|TumorKeller$|^DS2019|^6135","",samples_info$patient_id)
rownames(samples_info) <- samples_info$folder_name

## Adding Metadata
for (sample in names(corrected_samples)){
  corrected_samples[[sample]]$sample_name <- samples_info[sample, "sample_name"]
  corrected_samples[[sample]]$patient_id <- samples_info[sample, "patient_id"]
  corrected_samples[[sample]]$DiseaseState <- samples_info[sample, "DiseaseState"]
  corrected_samples[[sample]]$runId <- samples_info[sample, "runId"]
  corrected_samples[[sample]]$location <- samples_info[sample, "location"]
  corrected_samples[[sample]]$lesion <- samples_info[sample, "lesion"]
}

## Merging all samples into one object
merged_samples <- corrected_samples[[1]]
for (i in 2:length(corrected_samples)){
  merged_samples <- merge(merged_samples, corrected_samples[[i]])
}

## Saving merged object
saveRDS(merged_samples, file = "../../results/merged_samples_corrected_w_extra_10.rds")

#################################################################
##                           Data QC                           ##
#################################################################

merged_samples[['percent.mt_RNA']] <- PercentageFeatureSet(merged_samples, pattern = "^MT-", assay = "RNA")
merged_samples[['percent.mt_CorrectedCounts']] <- PercentageFeatureSet(merged_samples, pattern = "^MT-", assay = "CorrectedCounts")
Idents(merged_samples) <- "all_samples"

VlnPlot(merged_samples, features = c("nFeature_RNA","nFeature_CorrectedCounts"), pt.size = 0)
VlnPlot(merged_samples, features = c("nCount_RNA","nCount_CorrectedCounts"), pt.size = 0)
VlnPlot(merged_samples, features = c("percent.mt_RNA","percent.mt_CorrectedCounts"), pt.size = 0)
VlnPlot(integrated_data, group.by = "sample_name", features = "percent.mt_CorrectedCounts", pt.size = 0)

merged_samples <- subset(merged_samples, subset = nFeature_CorrectedCounts > 200 & percent.mt_CorrectedCounts < 15)

## Saving merged object
saveRDS(merged_samples, file = "../../results/merged_samples_fil_correctedcounts_w_extra_10.rds")
