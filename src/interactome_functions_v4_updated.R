# Interactome functions v4 updated 
# 2/22/21
# DO NOT CHANGE

check_genes <- function(genes, database, object_genes) {
  'Check to make sure that wanted genes are in reference you provide and in object
  
   Args:
   genes (chr vector): list of potential wanted genes
   database (data.frame): table with reference ligand/receptor pairs
   object_genes (chr vector): list of genes present in Seurat object
   
   Returns:
   not_found_list (chr list): list of genes not found in database and/or object
  '
  
  database_genes <- as.vector(unlist(database))
  not_found_database <- c()
  not_found_object <- c()
  
  for (i in 1:length(genes)) {
    if (!(genes[i] %in% database_genes)) {
      not_found_database <- c(not_found_database, genes[i])
    }
    
    if (!(genes[i] %in% object_genes)) {
      not_found_object <- c(not_found_object, genes[i])
    }
  }
  
  not_found_list <- list(database=not_found_database, object=not_found_object)
  return(not_found_list)
}

make_LR_pairs <- function(ligands, receptors, database) {
  'Make all LR pairs based on wanted ligands and/or receptors and database provided
  
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   database (data.frame): table with reference ligand/receptor pairs
   
   Returns:
   wanted_LR (data.frame): data.frame with ligand/receptor pairs from wanted ligands and receptors
  '
  
  wanted_LR <- data.frame(Ligand = character(0), Receptor = character(0))
  
  for (ligand in ligands){
    # list of corresponding receptors
    corresponding_receptors <- unique(as.character(database[,2][database[,1] == ligand]))
    
    for (receptor in corresponding_receptors) {
      LR_row <- data.frame(Ligand = ligand, Receptor = receptor)
      wanted_LR <- rbind(wanted_LR, LR_row)
    }
  }
  
  # filter out unwanted receptors
  wanted_LR <- wanted_LR[which(wanted_LR$Receptor %in% wanted_receptors),]
  
  return(wanted_LR)
}

create_LR_table <- function(ligands, receptors, cell_types, LRs, avg0, avg1) {
  'Create table w/ ligand, receptor, source and target cells, average expressions, and IDs
   
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   cell_types (chr vector): list of common cell types between two objects
   LRs (data.frame): table with with wanted ligand/receptor pairs (from make_LR_pairs function)
   avg0 (list of num data.frame): average expression table of object 0
   avg1 (list of num data.frame): average expression table of object 1
   
   Returns:
   LR_table (data.frame): table of potential ligand/receptor pairs
    Contains ligand, receptor, source, target, average expression of ligand and receptors from 
    source and target cells. Also contains IDs (important for Cytoscape). 
  '
  
  LR_table <- data.frame()
  count <- 0
  
  for (ligand in ligands) {
    known_receptors <- LRs$Receptor[which(LRs$Ligand == ligand)]
    
    for (receptor in known_receptors) {
      for (i in c(1:length(cell_types))) {
        for (j in c(1:length(cell_types))) {
          LR_table <- rbind(LR_table,data.frame(ligands = ligand, receptors = receptor, 
                                                source = cell_types[i], target = cell_types[j], 
                                                avg_lig_0 = avg0[ligand, cell_types[i]],
                                                avg_rec_0 = avg0[receptor, cell_types[j]],
                                                avg_lig_1 = avg1[ligand, cell_types[i]],
                                                avg_rec_1 = avg1[receptor, cell_types[j]]))
        }
      }
    }
    
    cat(count, ' ')
    count <- count+1
  }
  
  return(LR_table)
}

avg_LR_filt <- function(table, threshold) {
  'Calculates states (ON/OFF --> 1/0) and filter if there is no expression in both groups
  
   Args:
   table (data.frame): table with potential ligand/receptors (from create_LR_table)
   threshold (num): average expression threshold 
   
   Returns:
   table (data.frame): filtered table based on average expression
  '
  
  # Find states of pairs in each group
  LR_states <- data.frame(lig_0 = character(0), rec_0 = character(0),
                          lig_1 = character(0), rec_1 = character(0),
                          is_on = logical(0))
  
  for (i in c(1:length(table$avg_lig_0))) {
    row_states <- data.frame(lig_0 = table$avg_lig_0[i] > threshold, 
                             rec_0 = table$avg_rec_0[i] > threshold,
                             lig_1 = table$avg_lig_1[i] > threshold, 
                             rec_1 = table$avg_rec_1[i] > threshold)
    row_states$is_on <- (row_states$lig_0 & row_states$rec_0) | (row_states$lig_1 & row_states$rec_1)
    
    LR_states <- rbind(LR_states, row_states)
  }
  
  table <- cbind(table, ifelse((LR_states$lig_0 & LR_states$rec_0), 1, 0))
  table <- cbind(table, ifelse((LR_states$lig_1 & LR_states$rec_1), 1, 0))
  
  colnames(table)[9] <- 'is_on_0'
  colnames(table)[10] <- 'is_on_1'
  
  # Filter out pairs if pairs in both group are less than threshold 
  table <- table[which(LR_states$is_on == TRUE),]
  
  return(table)
}

LR_diff <- function(table, data_0, data_1, genes, label, alpha = 0.05) {
  'Calculate Wilcoxon-rank test on ligands and receptors (separately) between groups
    Order of test: data_0 --> data_1
  
    Args:
    table (data.frame): table with potential ligand/receptors (from create_LR_table)
    data_0 (data.frame): table of gene expression from object 0
    data_1 (data.frame): table of gene expression from object 1
    genes (chr vector): list of wanted genes
    label (chr): name of meta.data slot in Seurat object
    alpha (num): alpha level for Wilcoxon-rank test
    
    Returns:
    table (data.frame): filtered table based on Wilcoxon-rank test
  '
  
  table$lig_diff <- rep(0, length(table[,1]))
  table$rec_diff <- rep(0, length(table[,1]))
  
  table$lig_diff_p <- rep(0, length(table[,1]))
  table$rec_diff_p <- rep(0, length(table[,1]))
  
  for (i in 1:length(table[,1])) {
    ligand <- table$ligands[i]
    receptor <- table$receptors[i]
    source <- table$source[i]
    target <- table$target[i]
    
    lig_0_data <- data_0[which(data_0[,label] == source), ligand]
    rec_0_data <- data_0[which(data_0[,label] == target), receptor]
    
    lig_1_data <- data_1[which(data_1[,label] == source), ligand]
    rec_1_data <- data_1[which(data_1[,label] == target), receptor]
    
    lig_wilcox <- wilcox.test(lig_0_data, lig_1_data, exact = F, paired = F)
    rec_wilcox <- wilcox.test(rec_0_data, rec_1_data, exact = F, paired = F)
    
    table$lig_diff_p[i] <- lig_wilcox$p.value
    table$rec_diff_p[i] <- rec_wilcox$p.value
  }
  
  table$lig_diff_p_adj <- p.adjust(table$lig_diff_p, method = 'bonferroni')
  table$rec_diff_p_adj <- p.adjust(table$rec_diff_p, method = 'bonferroni')
  
  # If not significant, then 0
  # If significant, then 1
  for (i in 1:length(table[,1])) {
    table$lig_diff[i] <- ifelse(table$lig_diff_p_adj[i] < alpha, 1, 0)
    table$rec_diff[i] <- ifelse(table$rec_diff_p_adj[i] < alpha, 1, 0)
  }
  
  # If there is difference, then find if increase/decrease
  # for ligands
  for (i in 1:length(table[,1])) {
    if (table$lig_diff[i] == 1) {
      table$lig_diff[i] <- ifelse(table$avg_lig_0[i] > table$avg_lig_1[i], yes = -1, no = 1)
    }
    
    if (table$rec_diff[i] == 1) {
      table$rec_diff[i] <- ifelse(table$avg_rec_0[i] > table$avg_rec_1[i], yes = -1, no = 1)
    }
  }
  
  return(table)
}

