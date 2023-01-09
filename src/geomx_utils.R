library(ggplot2)
library(EnhancedVolcano)
library(cowplot)
library(stringr)

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         split_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + 
    # guides(fill = "none") +
    facet_wrap(as.formula(paste("~", split_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

# running QC on object
runQC <- function(x, stringent = TRUE){
  
  ns_data <- shiftCountsOne(x, useDALogic = TRUE)
  
  QC_params <-
    list(minSegmentReads = 1000, 
         percentTrimmed = 80,    
         percentStitched = 80,   
         percentAligned = 75,    
         percentSaturation = 0, 
         minNegativeCount = 1,   
         maxNTCCount = 9000,     
         minNuclei = 20,         
         minArea = 1000)   
  
  cat("Segment QC \n")
  ns_data <-
    setSegmentQCFlags(ns_data, 
                      qcCutoffs = QC_params)      
  
  QCResults <- protocolData(ns_data)[["QCFlags"]]
  flag_columns <- colnames(QCResults)
  QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                           Warning = colSums(QCResults[, flag_columns]))
  QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
  })
  QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))
  
  # calculate the negative geometric means for each module
  negativeGeoMeans <- 
    esBy(negativeControlSubset(ns_data), 
         GROUP = "Module", 
         FUN = function(x) { 
           assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
  protocolData(ns_data)[["NegGeoMean"]] <- negativeGeoMeans
  
  pkcs <- annotation(ns_data)
  modules <- gsub(".pkc", "", pkcs)
  
  # explicitly copy the Negative geoMeans from sData to pData
  negCols <- paste0("NegGeoMean_", modules)
  pData(ns_data)[, negCols] <- sData(ns_data)[["NegGeoMean"]]
  pData(ns_data)$sample <- as.factor(pData(ns_data)$sample)
  pData(ns_data)$segment <- as.factor(pData(ns_data)$segment)
  pData(ns_data) <- pData(ns_data)[, !colnames(pData(ns_data)) %in% negCols]
  
  ns_data <- ns_data[, QCResults$QCStatus == "PASS"]
  
  ns_data <- setBioProbeQCFlags(ns_data, 
                                  qcCutoffs = list(minProbeRatio = 0.1,
                                                   percentFailGrubbs = 20), 
                                  removeLocalOutliers = TRUE)
  cat("Probe QC \n")
  ProbeQCResults <- fData(ns_data)[["QCFlags"]]
  
  ProbeQCPassed <- 
    subset(ns_data, 
           fData(ns_data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
             fData(ns_data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
  
  ns_data <- ProbeQCPassed 
  
  # collapse to targets
  target_ns_data <- aggregateCounts(ns_data)
  
  cutoff <- 2
  minLOQ <- 2
  cat("Calculating LOQ \n")
  # Calculate LOQ per module tested
  LOQ <- data.frame(row.names = colnames(target_ns_data))
  for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(target_ns_data)))) {
      LOQ[, module] <-
        pmax(minLOQ,
             pData(target_ns_data)[, vars[1]] * 
               pData(target_ns_data)[, vars[2]] ^ cutoff)
    }
  }
  pData(target_ns_data)$LOQ <- LOQ
  
  LOQ_Mat <- c()
  for(module in modules) {
    ind <- fData(target_ns_data)$Module == module
    Mat_i <- t(esApply(target_ns_data[ind, ], MARGIN = 1,
                       FUN = function(x) {
                         x > LOQ[, module]
                       }))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
  }
  
  # ensure ordering since this is stored outside of the geomxSet
  LOQ_Mat <- LOQ_Mat[fData(target_ns_data)$TargetName, ]
  
  pData(target_ns_data)$GenesDetected <- 
    colSums(LOQ_Mat, na.rm = TRUE)
  pData(target_ns_data)$GeneDetectionRate <-
    pData(target_ns_data)$GenesDetected / nrow(target_ns_data)
  
  # Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
  pData(target_ns_data)$DetectionThreshold <- 
    cut(pData(target_ns_data)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
  
  if (stringent){seg_threshold = 0.1} else {seg_threshold=0.01}
  target_ns_data <-
    target_ns_data[, pData(target_ns_data)$GeneDetectionRate >= seg_threshold]
  
  LOQ_Mat <- LOQ_Mat[, colnames(target_ns_data)]
  fData(target_ns_data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
  fData(target_ns_data)$DetectionRate <-
    fData(target_ns_data)$DetectedSegments / nrow(pData(target_ns_data))
  
  negativeProbefData <- subset(fData(target_ns_data), CodeClass == "Negative")
  neg_probes <- unique(negativeProbefData$TargetName)
  if (stringent){gene_threshold = 0.1} else {gene_threshold=0.01}
  target_ns_data <- 
    target_ns_data[fData(target_ns_data)$DetectionRate >= gene_threshold |
                       fData(target_ns_data)$TargetName %in% neg_probes, ]
  
  return(target_ns_data)
}

calc_CV <- function(x) {sd(x) / mean(x)}

volcano_plot <- function(results, title = NULL){
  title <- ggdraw() + 
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 7)
    )
  
  plot <- plot_grid(
    
    EnhancedVolcano(toptable = results,
                    x = "logFC",
                    y = "p_adj",
                    lab = rownames(results),
                    FCcutoff = 0.2,
                    pCutoff = 0.05,
                    labSize = 3, labFace = "bold",
                    drawConnectors = T, max.overlaps = 20,
                    # ylim = c(0,max(results$p_adj)+3),
                    xlim = c(min(results$logFC)-0.2, max(results$logFC)+0.2),
                    title = NULL,
                    subtitle = "p adjusted",
    ) + theme(legend.position = "None"),
    
    EnhancedVolcano(toptable = results,
                    x = "logFC",
                    y = "p_value",
                    lab = rownames(results),
                    FCcutoff = 0.2,
                    pCutoff = 0.05,
                    labSize = 3, labFace = "bold",
                    drawConnectors = T, max.overlaps = 20,
                    # ylim = c(0,max(results$p_value)+3),
                    xlim = c(min(results$logFC)-0.2, max(results$logFC)+0.2),
                    title = NULL,
                    subtitle = "p nominal",
                    legendLabels = F
    ) + theme(legend.position = "None"),
    
    nrow = 1, ncol = 2
  )
  
  plot_grid(title, plot, ncol = 1, rel_heights = c(0.1, 1))
  
}

violin_plot <- function(data, gene){
  ggplot(pData(data),
         aes(x = Progression_Status, fill = Progression_Status,
             y = assayDataElement(data[gene, ],
                                  elt = "q_norm"))) +
    geom_violin() +
    # geom_jitter(width = 0.2, aes(color = Patient_ID)) +
    geom_jitter(width = 0.1) +
    labs(y = paste0(gene," Expression")) +
    scale_y_continuous(trans = "log2") +
    theme_bw()
}

create_cluster_col <- function(ns_object, clusters){
  clusters_names <- as.character(unique(pData(ns_object)[[clusters]]))
  for (cluster in clusters_names){
    pData(ns_object)[[cluster]] <- paste0("Non",cluster)
    pData(ns_object)[pData(ns_object)[[clusters]] == cluster,][[cluster]] <- paste(cluster)
    pData(ns_object)[[cluster]] <- factor(pData(ns_object)[[cluster]], levels = c(paste0(cluster), paste0("Non",cluster)))
  }
  return(ns_object)
}

process_lmm_results <- function(lmm_results){
  r_test <- do.call(rbind, lmm_results["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  r_test$Gene <- 
    unlist(lapply(colnames(lmm_results),
                  rep, nrow(lmm_results["lsmeans", ][[1]])))
  
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  colnames(r_test) <- c("Gene", "Contrast", "logFC","p_value", "p_adj")
  
  return(r_test)
}

extract_pairwise_markers <- function(lmm_results, cell_type, FC_cutoff, only_genes = FALSE){
  genes1 <- lmm_results %>% filter(p_adj < 0.05 & logFC > FC_cutoff & grepl(paste0(cell_type," -"), .$Contrast))
  genes2 <- lmm_results %>% filter(p_adj < 0.05 & logFC < - FC_cutoff & grepl(paste0("- ",cell_type), .$Contrast))
  genes <- rbind(genes1, genes2)
  if(! only_genes){
  return(genes)
  } else {
  return(unique(genes$Gene))
  }
}

extract_pan_markers <- function(ns_object, group, grouping_var, pCutOff, fcCutOff){
  
  require(tidyr)
  
  group_names <- as.character(unique(pData(ns_object)[[group]]))
  cores <- parallel::detectCores()
  results <- matrix(ncol = 5)
  
  for (sub_group in group_names){
    message(stringr::str_glue("Extracting {sub_group} markers"))
    pData(ns_object)[[sub_group]] <- paste0("Non",sub_group)
    pData(ns_object)[pData(ns_object)[[group]] == sub_group,][[sub_group]] <- paste(sub_group)
    pData(ns_object)[[sub_group]] <- factor(pData(ns_object)[[sub_group]], levels = c(paste0(sub_group), paste0("Non",sub_group)))
    
    formula <- as.formula(stringr::str_glue("~ {sub_group} + (1 + {sub_group} | {grouping_var})"))
    
    markers <-
      mixedModelDE(ns_object,
                   elt = "log_q",
                   modelFormula = formula,
                   groupVar = paste0(sub_group), nCores = cores)
    
    markers <- process_lmm_results(markers)
    markers$Contrast <- sub_group
    
    results <- rbind(results, as.matrix(markers))
  }
  
  results <- results %>% 
    as.data.frame() %>% 
    drop_na()  %>%
    mutate_at(c("logFC","p_value","p_adj"), as.numeric) %>%
    `rownames<-`(NULL) %>%
    filter(p_adj < pCutOff & logFC > fcCutOff)
}
