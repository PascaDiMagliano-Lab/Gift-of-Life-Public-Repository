pseudobulk_dge <- function(seurat_object, assay = "RNA", aggregate_by = "sample_name", comparison = "DiseaseState", filter = FALSE){
  
  #' @title PseudoBulk DGE
  #' @description  this function performs pseudobulk DGE on single cell seurat object
  #'
  #' @param seurat_object seurat object
  #' @param assay assay to be used
  #' @param aggregate_by column name with sample names to aggregate counts
  #' @param groups column name with groups for comparisons
  
  ## Loading Libraries
  suppressPackageStartupMessages({
    require(Seurat)
    require(SingleCellExperiment)
    require(Matrix.utils)
    require(tidyverse)
    require(magrittr)
    require(DESeq2)
    require(PCAtools)
    require(pheatmap)
  })

  message("Creating and aggregating count matrix")
  ## Extracting count matrix
  counts <- seurat_object@assays[[assay]]@counts
  metadata <- seurat_object@meta.data
  ## Converting to SCE
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  ## Creating Aggregate Matrix 
  groups <- colData(sce)[, aggregate_by]
  pb <- t(aggregate.Matrix(t(counts(sce)), 
                           groupings = groups, fun = "sum"))
  pb <- pb[rowSums(pb) > 0,] ## Removing genes with zero count in all samples
  
  pb_melt <- reshape2::melt(as.data.frame(pb), measure.vars = colnames(pb))
  colnames(pb_melt) <- c("sample","gene")
  
  ggplot(pb_melt, aes(x = sample, y = log2(gene)))+
    geom_violin() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    scale_y_continuous(breaks=c(0,seq(from = 0, to = 20, by = 0.5))) +
    ylab("log2(genes count)")
  
hist(log2(rowMeans(pb)))
  
  message("Processing metadata")
  ## Processing Metadata
  metadata <- metadata[,c(aggregate_by, comparison)] %>%
    distinct()
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[colnames(pb),]
  
  message("Creating DESeq object")
  ## DESeq object
  Design <- as.formula(paste0("~ ", comparison))
  dds <- DESeqDataSetFromMatrix(pb, 
                                colData = metadata, 
                                design = Design)
  dds <- DESeq(dds, quiet = T)
  return(dds)
}

pairwise_dge <- function(dds, group1 = NULL, group2 = NULL){
  
  #' @title pairwise_DGE
  #' @description  this function performs pairwise DGE on DDS object
  #'
  #' @param dds DESeq2 object
  #' @param group1 group1 of the comparison
  #' @param group2 group2 of the comparison
  
  suppressPackageStartupMessages({
    require(DESeq2)   
  })
  
  design <- dds@design[[2]]
  comp <- results(dds, 
                  contrast = as.character(c(design, group1, group2)),
                  alpha = 0.05)
  
  comp <- lfcShrink(dds, 
                    contrast = as.character(c(design, group1, group2)),
                    res=comp, 
                    type = "ashr") %>%
    as.data.frame()
  
  comp <- comp[!is.na(comp$padj),]
  
  return(comp)
  }

pca_plot <- function(dds, groups){
  #' @title Plot PCA from Pseudobulk DDS object
  #'
  #' @param dds DESeq2 object created from Seurat Object
  #' @param groups column name with groups for comparisons
  
  suppressPackageStartupMessages({
    require(DESeq2)   
  })
  
  vst <- assay(vst(dds))
  rv <- rowVars(vst)
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(dds)))]
  pca <- pca(vst[select,], metadata = colData(dds))
  pca_plot <- biplot(pca, showLoadings = F, boxedLoadingsNames = T,
                     colby = groups, colLegendTitle = groups,
                     labSize = 2, pointSize = 5, drawConnectors = T,
                     legendPosition = 'right',
                     legendLabSize = 7, legendIconSize = 4, legendTitleSize = 10,
                     encircle = T, encircleAlpha = 0.15)
  
  return(pca_plot)
}
  
corr_heatmap <- function(dds, groups, title){
  #' @title Plot correlation heatmap from Pseudobulk DDS object
  #'
  #' @param dds DESeq2 object created from Seurat Object
  #' @param groups column name with groups for comparisons
  #' @param title title of the plot
  
  suppressPackageStartupMessages({
    require(DESeq2)   
  })
  
  vst <- assay(vst(dds))
  vst_cor <- cor(vst)
  ann <- as.data.frame(colData(dds)[,groups])
  rownames(ann) <- colnames(dds)
  heatmap <- pheatmap(vst_cor, main = title,
                      annotation_col = ann,
                      annotation_row = ann,
                      fontsize = 7, border_color = NA, 
                      annotation_names_row = F, annotation_names_col = F)
  
  return(heatmap)
}

GO_enrich <- function(DGE_results, fcCutOff, pvalCutOff, GO_category = "BP", top = 10){
  
  suppressPackageStartupMessages({
    require(org.Hs.eg.db)
    require(clusterProfiler)
  })
  
  
  message("Extracting DEGs")  
  upreg_genes <- DGE_results %>%
      filter(padj < pvalCutOff & log2FoldChange > fcCutOff) %>%
      rownames() %>%
      mapIds(x = org.Hs.eg.db, column = "ENTREZID", keytype = "SYMBOL") %>%
      na.omit()
    
    dnreg_genes <- DGE_results %>%
      filter(padj < pvalCutOff & log2FoldChange < -fcCutOff) %>%
      rownames() %>%
      mapIds(x = org.Hs.eg.db, column = "ENTREZID", keytype = "SYMBOL") %>%
      na.omit()
    
    bg_genes <- mapIds(org.Hs.eg.db, rownames(DGE_results), "ENTREZID", "SYMBOL", multiVals = 'first') %>%
      na.omit()
    
    message("Performing GO Enrichment")
    up_GO <- enrichGO(gene          = upreg_genes,
                      universe      = bg_genes,
                      OrgDb         = org.Hs.eg.db,
                      ont           = GO_category,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable      = TRUE)
    
    dn_GO <- enrichGO(gene          = dnreg_genes,
                      universe      = bg_genes,
                      OrgDb         = org.Hs.eg.db,
                      ont           = GO_category,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable      = TRUE)
    
    GOA <- list('upreg_GO' = up_GO, 'dnreg_GO' = dn_GO)
}

GOA_plot <- function(x){
  GOA <- lapply(x, function(y){y <- y %>% as.data.frame() %>% head(10)})
  GOA[[1]]$p.adjust<- -log10(GOA[[1]]$p.adjust)
  GOA[[2]]$p.adjust<- log10(GOA[[2]]$p.adjust)
  final <- rbind(GOA[[1]],GOA[[2]]) %>% arrange(p.adjust)
  final$Description <- factor(final$Description, levels = final$Description)
  ggplot(final, aes(x =  Description, y = p.adjust, fill = p.adjust)) +
    geom_bar(stat = "identity") + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
    xlab(NULL) + ylab(NULL) +
    coord_flip()
}

GSEA_enrich <- function(x){
  
  
}
