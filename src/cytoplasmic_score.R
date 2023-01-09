#' QC score based on marker genes
#' 
#' sum of marker genes is used as the qc.score
#'
#' @param data : count matrix
#' @param qc.gene.FN : filename holding QC gene list 
#' @return qc.score : qc score
#' @examples
#' qc.score <- cal_qc_score(data, qc.gene.FN="AIBS_qc_genes_10X.csv") 

cal_qc_score <- function (data, qc.gene.FN="AIBS_qc_genes_10X.csv") {
  
  library("scrattch.hicat")
  
  qc.genes = read.csv(qc.gene.FN)[, "Gene"]
  common.qc.genes = intersect(qc.genes, row.names(data))
  norm.dat = logCPM(data[common.qc.genes,])
  qc.score = apply(norm.dat, 2, sum)
  return(qc.score)
}