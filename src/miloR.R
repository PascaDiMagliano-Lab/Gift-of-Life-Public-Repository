prepare_milo <- function(seurat_object){
  
  require(miloR)
  require(Seurat)
  require(SingleCellExperiment)
  
  message("Creating SCE object")
  sce <- as.SingleCellExperiment(seurat_object)
  reducedDim(sce, "PCA") <- Embeddings(seurat_object, "pca")
  reducedDim(sce, "UMAP") <- Embeddings(seurat_object, "umap")
  message("Creating Milo object")
  milo <- Milo(sce)
  
  return(milo)
}


pick_k <- function(milo, k_range = c(5,10,15,20,25,30,35,40,45,50), clusters_col){
  
  require(tidyverse)
  require(miloR)
  
  results <- sapply(k_range, function(k){
    
    message("K = ", k)
    
    traj_milo <- buildGraph(milo, k = k, d = 30)
    traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = k, d=30, refined = TRUE)
    
    hist <- plotNhoodSizeHist(traj_milo)
    mean <- mean(hist$data$nh_size)
    rare_pop_size <- min(table(as.character(colData(traj_milo)[,clusters_col])))
    ratio <- mean/rare_pop_size
    
    return(c(k, mean, rare_pop_size, ratio*100)) 
  }) %>%
    t() %>% as.data.frame() %>%
    `colnames<-`(c("K", "Mean Nhood Size", "Rarest Population Size", "Ratio"))
}

run_milo <- function(milo_obj, k, p, d = 30, samples_col, graph = "PCA"){
  
  require(miloR)
  
  message("::> Building Graph")
  traj_milo <- buildGraph(milo_obj, k = k, d = d, reduced.dim = graph)
  message("::> Making Nhoods")
  traj_milo <- makeNhoods(traj_milo, prop = p, k = k, d = 30, refined = TRUE)
  message("::> Counting Cells")
  traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample = samples_col)
  message("::> Calculating Nhood Distance (This might take a while ..)")
  traj_milo <- calcNhoodDistance(traj_milo, d = 30) 
  traj_milo <- buildNhoodGraph(traj_milo)
  
  return(traj_milo)
}

DA_test <- function(milo_obj, samples_col, condition, contrast, annotation_col, plot = FALSE){
  
  require(miloR)
  
  message("::> Creating Design Matrix")
  traj_design <- data.frame(colData(milo_obj))[,c(samples_col, condition)] %>%
    `colnames<-`(c("Sample","Condition"))
  traj_design <- distinct(traj_design) %>% `row.names<-`(.$Sample)
  # the syntax is <VariableName><ConditionLevel> - <VariableName><ControlLevel>
  contrast_mtx <- paste0("Condition",contrast[1]," - Condition", contrast[2])
  
  message("::> Testing Nhoods")
  da_results <- testNhoods(milo_obj, design = ~ 0 + Condition, design.df = traj_design,
                           model.contrasts = contrast_mtx, fdr.weighting = "graph-overlap")
  
  message("::> Annotating Nhoods")
  da_results <- annotateNhoods(milo_obj, da_results, coldata_col = annotation_col)
  
  
  if(plot){
    milo_obj <- buildNhoodGraph(milo_obj)
    plotNhoodGraphDA(milo_obj, da_results, alpha=0.05) + 
      ggtitle(paste0("Differential Abudnance - ", contrast[1], " vs. ", contrast[2])) &
    plotDAbeeswarm(da_results, group.by = annotation_col, alpha = 0.05) + 
      ggtitle(paste0("Differential Abudnance - ", contrast[1], " vs. ", contrast[2]))}

  return(da_results)
}

