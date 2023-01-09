library(cowplot)

plot_correction <- function(object, genes){
  
  DefaultAssay(object) <- "RNA"
  original <- FeaturePlot(object, features = genes, ncol = length(genes), slot = "data", order = T)
  DefaultAssay(object) <- "CorrectedCounts"
  corrected <- FeaturePlot(object, features = genes, ncol = length(genes), slot = "data", order = T)
  
  plots_list <- plot_grid(original, corrected, labels = c("original","corrected"), nrow = 2, ncol = 1)
  
  return(plots_list)
}
