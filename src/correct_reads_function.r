correct_reads <- function(data_dir, rho = NULL , extra = NULL, cap = FALSE, h5 = TRUE){

  #' @title Correct Ambient RNA Contamination
  #' @description  Correct ambient RNA contamination of 10x Genomics scRNA seq of Pancreatic tissue with acinar transcipts contamination
  #'
  #' For this function, you can provide the data directory of cellranger output
  #' it should at least contain raw_feature_bc_matrix and filtered_feature_bc_matrix folders
  #' if you didn't provide (rho), the function will try to estimate the contamination fraction
  #' using the common genes that contaminate pancreatic tissue sample (acinar, RBC and endocrine markers)
  #' if you provide (rho) it will use it as given contamination fraction and won't estimate it
  #' you can also provide argument (extra) which is added on top of the manually estimated
  #' contamination fraction. This function doesn't include the automatic estimation of the
  #' as in the original paper contamination fraction as it was found inefficient in cleaning up the data

  #' @param data_dir directory with 10x Genomics output (i.e., raw_feature_bc_matrix and filtered_feature_bc_matrix folders)
  #' @param rho integer to set contamination fraction instead of estimating it
  #' @param extra integer to add extra percent to the contamination fraction
  #' @param cap bool to either cap the contamination fraction to 20\% or not
  #' @param h5 bool if input is '.h5' file
  #' @returns Returns Seurat object with extra assay of Corrected Counts

  # importing libraries
  library(Seurat)
  library(SoupX)
  library(dplyr)
  library(tidyr)

  # importing dataset
  if(h5){
    filt.matrix <- Read10X_h5(list.files(data_dir, full.names = T, pattern = '_filtered_feature_bc_matrix.h5'),use.names = T)
    raw.matrix  <- Read10X_h5(list.files(data_dir, full.names = T, pattern = '_raw_feature_bc_matrix.h5'),use.names = T)
    sc  <- SoupChannel(raw.matrix, filt.matrix)
  } else {
  sc = load10X(data_dir, verbose = F)
  }

  # Creating seurat clusters
  if(h5){
    seurat_obj <- Read10X_h5(list.files(data_dir, full.names = T, pattern = '_filtered_feature_bc_matrix.h5'),use.names = T)
  } else {
  seurat_obj <- Read10X(paste0(data_dir,"/filtered_feature_bc_matrix/"))
  }

  seurat_obj <- CreateSeuratObject(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

  # addin clusters and DR to soupx object
  clusters <- data.frame("clusters" = seurat_obj$seurat_clusters, row.names = colnames(seurat_obj))
  umap <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  sc = setClusters(sc, setNames(clusters$clusters, rownames(clusters)))
  sc = setDR(sc, umap)

if (missing(rho)){

  acinar_markers <- c('CTRB1','KLK1','RBPJL','PTF1A','CELA3A','PRSS1','SPINK1','ZG16','CEL','CELA2A','CPB1','CELA1','RNASE1','AMY2B','CPA2','CPA1','CELA3B','PNLIP','CTRB2','PLA2G1B','PRSS2','CLPS','REG1A','SYCN','PNLIPRP1','CTRC','REG3A','SERPINA3','PRSS3','REG1B','CFB','GDF15','MUC1','C15ORF48','DUOXA2','AKR1C3','OLFM4','GSTA1','LGALS2','PDZK1IP1','RARRES2','CXCL17','UBD','GSTA2','ANPEP','LYZ','ANGPTL4','ALDOB')
  acinar_markers <- acinar_markers[acinar_markers %in% row.names(sc$toc)]
  rbc_markers <- c('HEXA','AQP1','CCNA2','PKLR','LOX','ARG1','HBB','HBA1','HBA2','HBG1','HBG2','CEACAM1','CD36','GYPA','THBS1','ITGA4','GYPB','AHSP','HBM','FECH','EPOR','KIT','BPGM','DCAF12','HBD','HEMGN','HMBS','NCOA4','RBM38','UCP2','UBB','TMCC2','SLC4A1','SLC25A39','SLC25A37','FOSB','HTATSF1','EPB41','RSAD2','COX6B2','APOA1','ASNS','ALAD','S100A9','AHSG','S100A8','ERMAP','KCNN4','PAN3','PFAS','VAMP5','PNPO')
  rbc_markers <- rbc_markers[rbc_markers %in% row.names(sc$toc)]
  beta_markers <- c('GCGR','JPH3','CD40','HAMP','EZH1','NTRK1','PDX1','SLC2A2','NKX6-2','FXYD2','NPY','INS','RIMS1','MAFA','EFNA5','LMX1A','NKX2-2','NKX6-1','PAX4','IAPP','PCSK2','G6PC2','SLC30A8','PCSK1','GJD2','SCGN','IGF2','SYT13','FFAR2','NPTX2','PFKFB2','EDARADD','HOPX','SH3GL2','ADCYAP1','SCGB2A1','CASR','MAFB','PAX6','NEUROD1','ISL1','TGFBR3','SMAD9','SIX3','SIX2','BMP5','PIR','STXBP5','DLK1','MEG3','RGS16')
  beta_markers <- beta_markers[beta_markers %in% row.names(sc$toc)]

  nonExpressedGeneList = list(acinar_markers  = acinar_markers,
                              rbc_markers = rbc_markers,
                              beta_markers = beta_markers)

  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = nonExpressedGeneList)

  sc = calculateContaminationFraction(sc, nonExpressedGeneList = nonExpressedGeneList, useToEst = useToEst)
  print(paste0("The estimated contamination fraction is: ", round(sc$metaData$rho[1],3)))

  if (missing(extra)){

    note = "no extra value was added to the contamination fraction as it wasn't provided"
    if (cap){
    new_rho = min(sc$metaData$rho[1],0.2)
    sc = setContaminationFraction(sc, new_rho)
    out = adjustCounts(sc, roundToInt = T, verbose = F)
    } else {out = adjustCounts(sc, roundToInt = T, verbose = F)}

  } else {

    new_rho = sc$metaData$rho[1]+extra
    note = "added extra 10% to the contamination fraction as provided"
    if (cap){
      new_rho = min(new_rho,0.2)
      sc = setContaminationFraction(sc, new_rho)
      out = adjustCounts(sc, roundToInt = T, verbose = F)
    } else {
      sc = setContaminationFraction(sc, new_rho, forceAccept = T)
      out = adjustCounts(sc, roundToInt = T, verbose = F)
    }
  }

    # if (sc$metaData$rho[1] > 0.2) {
    #   note = "no extra value was added to the contamination fraction, Contamination fraction is already high (> 20%)"
    #   new_rho = min(sc$metaData$rho[1],0.2)
    #   sc = setContaminationFraction(sc, new_rho)
    #   out = adjustCounts(sc, roundToInt = T, verbose = F)

    # } else {
    # print(paste0("Adding ",extra," to the estimated fraction"))
    # new_rho = sc$metaData$rho[1]+extra
    # sc = setContaminationFraction(sc, new_rho)
    # note = "added extra 10% to the contamination fraction as provided"
    # out = adjustCounts(sc, roundToInt = T, verbose = F)
    # }


} else {

  sc = setContaminationFraction(sc, rho)
  out = adjustCounts(sc, roundToInt = T, verbose = F)
  note = "Contamination Fraction was provided by the user"

}

TopContGenes <- data.frame("Genes" = row.names(sc$soupProfile),
                           "counts" = sc$soupProfile$counts,
                           "fraction" = sc$soupProfile$est)

Misc(seurat_obj, "TopContGenes") <- TopContGenes[order(TopContGenes$counts, decreasing = T),]
Misc(seurat_obj, "ContaminationFraction") <- round(sc$metaData$rho[1],3)
Misc(seurat_obj, "notes") <- note

seurat_obj[['CorrectedCounts']] <- CreateAssayObject(counts = out)
DefaultAssay(seurat_obj) <- 'CorrectedCounts'

print(paste0("Returning seurat object with new default assay (CorrectedCounts) and top contaminating genes in Misc slot"))

return(seurat_obj)
}
