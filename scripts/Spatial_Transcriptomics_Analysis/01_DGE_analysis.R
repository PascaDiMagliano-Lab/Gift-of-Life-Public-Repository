##################################################################
##                  loading libraries and data                  ##
##################################################################

suppressPackageStartupMessages({
  library(NanoStringNCTools)
  library(GeomxTools)
  library(GeoMxWorkflows)
  library(ggthemes)
  library(ggiraph)
  library(pheatmap)
  library(tidyverse)
  library(knitr)
  library(dplyr)
  library(ggforce)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(reshape2)
  library(limma)
  library(umap)
  library(Rtsne)
  library(ggrepel)
  library(readxl)
  library(RColorBrewer)
  library(writexl)
  library(PCAtools)
  library(fgsea)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(VennDiagram)
  source("../../src/geomx_utils.R")
})

## Importing data
## This section will be re-written after uploading GeoMx to GEO
DCCFiles <- list.files("../../data/spatial_transcriptomics/dcc/", full.names = T, pattern = ".dcc$")
PKCFiles <- "../../data/spatial_transcriptomics/Hs_R_NGS_WTA_v1.0.pkc"
SampleAnnotationFile <- ("../../data/spatial_transcriptomics/pooled_annotation_eileen_10-5-22.xlsx")

ns_object <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = "sheet1",
                         phenoDataDccColName = "Sample_ID",
                         experimentDataColNames = c("panel"))

pkcs <- annotation(ns_object)
modules <- gsub(".pkc", "", pkcs)
pData(ns_object)$cell_type <- pData(ns_object)$`Cell Type_Updated`
pData(ns_object)$sample <- pData(ns_object)$patient_id
pData(ns_object)$tissue <- paste(pData(ns_object)$patient_id, pData(ns_object)$Location, sep = "_")
pData(ns_object)$main_label <-
  gsub(pattern = "MettoLN|Glandular_Tumor|PoorlyDiffTumor|TumorNerve",replacement = "Tumor", pData(ns_object)$cell_type)
count_mat <- dplyr::count(pData(ns_object), `scan name`, patient_id, main_label)
test_gr <- gather_set_data(count_mat, 1:3)
test_gr$x <- factor(test_gr$x,
                    levels = c("scan name","patient_id", "main_label"))
ns_object <- ns_object[,pData(ns_object)$cell_type %in% c("Acinar", "Duct", "PanIN", "ADM", "Tumor")]

#################################################################
##                  Data QC and Normalization                  ##
#################################################################

## QC
ns_object <- runQC(ns_object, stringent = TRUE)
ns_object <- ns_object[,!(pData(ns_object)$cell_type %in% c("MettoLN","TumorNerve"))]
# pData(ns_object)$cell_type <- factor(pData(ns_object)$cell_type, levels = c("Acinar","ADM","Duct","PanIN","Glandular_Tumor","PoorlyDiffTumor"))

## Data Summary
## Plot
count_mat <- dplyr::count(pData(ns_object), patient_id, cell_type)
test_gr <- gather_set_data(count_mat, 1:2)
test_gr$x <- factor(test_gr$x,
                    levels = c("patient_id","cell_type"))

ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = cell_type), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.3) +
  geom_parallel_sets_labels(color = "white", size = 4, angle = 0) +
  theme_classic(base_size = 17) +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  # axis.text.x = element_blank(),
  # axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = expansion(0)) +
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "")

## Summary
summary_table <- t(as.data.frame.matrix(table(pData(ns_object)[,c("patient_id","cell_type")]))) %>%
  as.data.frame() %>%
  mutate(total = rowSums(.))
summary_table$bio_rep <- apply(summary_table, 1, function(x){sum(x[1:6]!=0)})
write.csv(summary_table, "pooled_analysis_plots/summary_table.csv")

## Normalization
ns_object <- normalize(ns_object,
                     norm_method = "quant",
                     desiredQuantile = .75,
                     toElt = "q_norm")

assayDataElement(object = ns_object, elt = "log_q") <-
  assayDataApply(ns_object, 2, FUN = log, base = 2, elt = "q_norm")

saveRDS(ns_object, "../../results/processed_ns_object.RDS")

##################################################################
##                           PCA Plot                           ##
##################################################################

calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(ns_object,
                         elt = "log_q", MARGIN = 1, calc_CV)

CV_dat <- sort(CV_dat, decreasing = TRUE)

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]

# PCA
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

pre_pca <- pca(mat = assayDataElement(object = ns_object[GOI,], elt = "log_q"),
                   metadata = pData(ns_object))

pre <- biplot(pre_pca, showLoadings = F, boxedLoadingsNames = T,
       colby = 'cell_type', lab = pData(ns_object)$Sample, shape = 'patient_id',
       # colkey = c(Acinar = 'red', Duct = 'magenta', ADM = 'darkgreen', PanIN = 'navy'),
       title = "pre batch correction", subtitle = NULL,
       labSize = 0, pointSize = 2, drawConnectors = F,
       legendPosition = 'right',
       legendLabSize = 7, legendIconSize = 4, legendTitleSize = 10,
       encircle = T, encircleAlpha = 0.15)

assayDataElement(object = ns_object, elt = "limma_corrected_log_q") <-
  limma::removeBatchEffect(assayDataElement(object = ns_object, elt = "log_q"),
                           batch = pData(ns_object)$tissue,
                           design = model.matrix( ~ cell_type, data= pData(ns_object)))

limma_cor_pca <- pca(mat = assayDataElement(object = ns_object[GOI,], elt = "limma_corrected_log_q"),
                   metadata = pData(ns_object))

limma_post <- biplot(limma_cor_pca, showLoadings = F, boxedLoadingsNames = T,
               colby = 'cell_type', lab = pData(ns_object)$Sample, shape = 'patient_id',
               # colkey = c(Acinar = 'red', Duct = 'magenta', ADM = 'darkgreen', PanIN = 'navy'),
               title = "post batch correction", subtitle = NULL,
               labSize = 0, pointSize = 2, drawConnectors = F,
               legendPosition = 'right',
               legendLabSize = 7, legendIconSize = 4, legendTitleSize = 10,
               encircle = T, encircleAlpha = 0.15)

## FOR LATER MAKE IT INTERACTIVE ##
limma_post <- biplot(limma_cor_pca, showLoadings = F, boxedLoadingsNames = T,
                     colby = 'cell_type', lab = pData(ns_object)$roi, shape = "sample",
                     title = "post batch correction", subtitle = NULL,
                     labSize = 5, pointSize = 3, drawConnectors = T,
                     legendPosition = 'right', colkey = getPalette(colourCount),
                     legendLabSize = 7, legendIconSize = 4, legendTitleSize = 10,
                     encircle = T, encircleAlpha = 0.15)

cowplot::plot_grid(pre, limma_post)

#################################################################
##                             DGE                             ##
#################################################################

pan_markers <- extract_pan_markers(ns_object = ns_object, group = "cell_type",
                                   grouping_var = "tissue",
                                   pCutOff = 1, fcCutOff = -Inf)
## Top 20 markers
panin_pan_markers <- pan_markers %>% filter(Contrast == "PanIN" & p_adj < 0.05) %>% top_n(20, logFC) %>% pull(Gene)
adm_pan_markers <- pan_markers %>% filter(Contrast == "ADM" & p_adj < 0.05) %>% top_n(20, logFC) %>% pull(Gene)
acinar_pan_markers <- pan_markers %>% filter(Contrast == "Acinar" & p_adj < 0.05) %>% top_n(20, logFC) %>% pull(Gene)
ductal_pan_markers <- pan_markers %>% filter(Contrast == "Duct" & p_adj < 0.05) %>% top_n(20, logFC) %>% pull(Gene)
Glandular_Tumor_pan_markers <-  pan_markers %>% filter(Contrast == "Glandular_Tumor" & p_adj < 0.05) %>% top_n(20, logFC) %>% pull(Gene)
PoorlyDiffTumor_pan_markers <-  pan_markers %>% filter(Contrast == "PoorlyDiffTumor" & p_adj < 0.05) %>% top_n(20, logFC) %>% pull(Gene)

markers_list <- list('acinar_markers' = acinar_pan_markers,
                     'adm_markers' = adm_pan_markers,
                     'ductal_markers' = ductal_pan_markers,
                     'panin_markers' = panin_pan_markers,
                      'Glandular_Tumor_markers' = Glandular_Tumor_pan_markers,
                     'PoorlyDiffTumor_markers' = PoorlyDiffTumor_pan_markers
                      )

saveRDS(markers_list, "../../results/lmm_main_labels_markers_list.rds")

#################################################################
##                       Markers Heatmap                       ##
#################################################################

markers_df <- assayDataElement(ns_object, elt = "log_q")[unlist(markers_list),]
annot_df <- data.frame(cell_type = pData(ns_object)$cell_type, row.names = colnames(ns_object))
annot_df$cell_type <- factor(annot_df$cell_type, levels = c("Acinar","ADM", "Duct","PanIN","Glandular_Tumor","PoorlyDiffTumor"))
levels(annot_df$cell_type)

pheatmap(markers_df[,order(annot_df$cell_type)],
         color =  colorRampPalette(c("blue","white", "red"))(100) ,
         # color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         show_colnames = F,fontsize_row = 5,
         show_rownames = T, cluster_rows = F, cluster_cols = F,annotation_names_col = F,
         scale = "row", annotation_col = annot_df %>% arrange(cell_type))


#################################################################
##              PDAC Subtype Signature Enrichment              ##
#################################################################

##-------------------
##  Subtype Classes
##-------------------

## Collison geneset (https://pubmed.ncbi.nlm.nih.gov/21460848/)
collison_gs <- list(
  Collison_Exocrine_Like = c('REG1B','REG3A','REG1A','PNLIPRP2','CEL','PNLIP','PLA2G1B','CELA3A','CPB1','CELA3B','CTRB2','CLPS','CELA2B','PRSS2','PRSS1','GP2','SLC3A1','CFTR','SLC4A4','SPINK1'),
  Collison_Basal = c('AIM2','FAM26F','GPM6B','S100A2','KRT14','CAV1','LOX','SLC2A3','TWIST1','PAPPA','NT5E','CKS2','HMMR','SLC5A3','PMAIP1','PHLDA1','SLC16A1','FERMT1','HK2','AHNAK2'),
  Collison_Classical = c('TMEM45B','SDR16C5','GPRC5A','AGR2','S100P','FXYD3','ST6GALNAC1','CEACAM5','CEACAM6','TFF1','TFF3','CAPN8','FOXQ1','ELF3','ERBB3','TSPAN8','TOX3','LGALS4','PLS1','GPX2','ATP10B','MUC13'))

## Moffitt geneset (https://pubmed.ncbi.nlm.nih.gov/26343385/)
moffit_gs <- list(
  Moffit_Basal = c("VGLL1","UCA1","S100A2","LY6D","SPRR3","SPRR1B","LEMD1","KRT15","CTSL2","DHRS9","AREG","CST6","SERPINB3",
                   "SERPINB4","KRT6C","KRT6A","FAM83A","SCEL","FGFBP1","KRT7","KRT17","GPR87","TNS4","SLC2A1","ANXA8L2"),
  Moffit_Classical = c("BTNL8","FAM3D","ATAD4","AGR3","CTSE","LOC400573","LYZ","TFF2","TFF1","ANXA10","LGALS4","PLA2G10",
                       "CEACAM6","VSIG2","TSPAN8","ST6GALNAC1","AGR2","TFF3","CYP3A7","MYO1A","CLRN3","KRT20","CDH17","SPINK4","REG4")
)

## Baily geneset (https://www.nature.com/articles/nature16965)
baily_gs <- list()
baily_gs[['Baily_Basal']] <- c('PGAM1','ANGPTL4','ARPC2','ITGA3','PTPN1','GSDMC','DRAP1','ANXA8','LDHA','ARHGAP23','HCAR2','ULBP2','KIAA1609','PTHLH','CEBPB','SLC16A3','EGFR','PLXNA1','ZBED2','IL2A','TPD5212','DUSP14','FSCN1','RND3','GPR87','PANX1','S100A2','ANXA1','KRT6A','ADAM17','PTRF','EMLIN1','FBXL7','CERCAM','DDR2','GXYLT2','CHSY3','VSTM4','COL8A1','PRRX1','COL5A2','SPARC','FBN1','FAP','COL6A3','COL1A2','BNC2','EFEMP2','FSTL1','COL3A1','GSG2','KIF18A','CDCA5','CCNA2','MCM10','DEPDC1','UBE2C','KIF4A','CCNB1','KIF2C','SPAG5','TPX2','HIST1H2B0','BUB1','CKAP2L','PTTG1','CDC20','ERCC6L','NCAPG','NCAPH','PSMA7','CCT6A','CCT7','DKC1','HSPE1','EIF2S2','MRPL11','MRPL17','LYAR','PSMA1','EIF4A3','PPAT','HAT1','PAICS','BRIX1','RUVBL1','CKS1B','TOMM40','CCT5','ALYREF')
baily_gs[['Baily_Classical']] <- c('LRRC66','GJB1','LGALS4','PLA2G10','PLS1','BCAS1','CAPN5','HNF4G','ARHGEF38','CAPN8','AP001187.9','ABP1','ERBB3','C9orf152','SLC44A4','IHH','HNF4A','FAM83E','PRR15L','EPS8L3')
baily_gs[['Baily_ADEX']] <- c('REG3G','SYCN','REG1P','SERPINI2','CPB1','RP11-331F4.4','AQP8','RBPJL','CUZD1','PLA2G1B','CLPS','CEL','CELA3B','CELA3A','CTRB1','CTRC','CTRB2','PNLIPRP1','CPA1','CPA2','UNC79','GCK','GJD2','SYT4','KCNK16','NOL4','SCGN','INS','CABP7','CHGB','BEX1','SVOP','MR7-3HG','ABCC8','HMGCLL1','SLC30A8','SST','CELF3','PCSK2','SCG3')
# baily_gs[['Baily_Immunogenic']] <-  c('IGHV3-15','IGHV3-7','IGHV1-46','IGKV2-28','IGHV3-23','IGHV3-53','IGHV5-51','IGHA1','IGKV4-1','IGLV2-23','IGLV1-40','IGKV3-11','IGLL5','IGJ','IGLC3','IGKVLD-39','IGLV2-14','AC096579.7','IGKV3-20','IGKC','CSF1R','TYROBP','FCER1G','TLR8','C3AR1','CD86','WAS','SPI1','HAVCR2','SASH3','CYBB','MS4A4A','BTK','LAPTM5','PTPRC','MARCH1','CD53','AIF1','DOCK2','NCKAP1L','SIT1','GIMAP4','GPR174','SEPT1','CXCR6','ITK','TBC1D10C','GIMAP5','CD96','ZNF831','GZMK','PTPRCAP','SLAMF1','P2RY10','CD48','THEMIS','CD3E','CD3D','SH2D1A','CD2')

## Via GVSA
ns_object <- readRDS("pooled_analysis_objects/processed_ns_object.RDS")
subtypes_geneSets <- list("Collison_Exocrine_Like" = collison_gs$Collison_Exocrine_Like, "Baily_ADEX" = baily_gs$Baily_ADEX,
                          "Collison_Classical" = collison_gs$Collison_Classical, "Baily_Classical" = baily_gs$Baily_Classical,"Moffit_Classical" = moffit_gs$Moffit_Classical,
                          "Collison_Basal" = collison_gs$Collison_Basal, "Baily_Basal" = baily_gs$Baily_Basal, "Moffit_Basal" = moffit_gs$Moffit_Basal)
lesions_count_mt <- ns_object[,pData(ns_object)$main_label %in% c("Acinar", "Duct", "PanIN", "ADM", "Tumor")]
pdac_subtype_gvsa <- gsva(expr = assayDataElement(lesions_count_mt, elt = "q_norm"),
                          gset.idx.list = subtypes_geneSets)
pData(lesions_count_mt)$cell_type <- factor(pData(lesions_count_mt)$cell_type, levels = c("Acinar","ADM","PanIN","Duct","Glandular_Tumor","PoorlyDiffTumor"))
pdac_subtype_gvsa <- pdac_subtype_gvsa[,order(pData(lesions_count_mt)$cell_type)]
sample_type <- data.frame(pData(lesions_count_mt)$cell_type[order(pData(lesions_count_mt)$cell_type)], row.names = colnames(pdac_subtype_gvsa)) %>%
  `colnames<-`('cell_type')
sign_type <- data.frame(c("NA","NA","Classical","Classical","Classical","Basal","Basal","Basal"), row.names = names(subtypes_geneSets)) %>%
  `colnames<-`("Class")
ann_colors <- list(Class = c("NA" = 'black', "Classical" = "blue", "Basal" = "red"),
                   cell_type = c("Acinar" = 'yellow', "ADM" = "orange", "PanIN" = "darkred", "Duct" = "darkgreen", "Glandular_Tumor" = "darkblue", "PoorlyDiffTumor" = "purple"))

pheatmap(pdac_subtype_gvsa,
         annotation_col = sample_type, annotation_row = sign_type, annotation_colors = ann_colors,
         gaps_col = c(44,76,164,210,245), gaps_row = c(2,5),
         show_colnames = F, scale = "column", cluster_cols = F, cluster_rows = F,
         annotation_names_col = F)


