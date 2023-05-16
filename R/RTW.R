#'@title A method for differential expression analysis of single-cell RNA-seq 
#'data based on the discrete generalized beta distribution
#'@description  The function based on the discrete generalized beta distribution (DGBD), 
#'introduced the T-test and Wilcoxon test to further explain the variability of genes between different
#'cell types, and combined the obtained results using Fisher's method for comparing gene expression 
#'between single-cell groups with different phenotypes
#'@param object  Take SingleCellExperiment object as input 
#'@return Gene name, P-values, and adjusted P-values of BH
#'@export ROSeq
RTW<- function(object){
  result_Wilcoxon <- Execute_Wilcoxon(object)
  result_Ttest <- Execute_Ttest(object)
  result_DGBD <- Execute_DGBD(object)
  object <- data.frame(result_Wilcoxon[[2]],result_Ttest[[2]],result_DGBD[[2]])
  #object <- data.frame(result_Wilcoxon[[3]],result_Ttest[[3]],result_DGBD[[2]])
  cpm <- object
  idx <- 1:nrow(cpm)
  names(idx) <- rownames(cpm)
  p <- sapply(idx, function(i){
    sumlog(cpm[i,])$p
  })
  FDR <- p.adjust(p, method = "BH")
  result_p<- list(gene_names = names(p),
                  pvalue = p,
                  FDR = FDR)
  return(result_p)
}

# perform Wilcoxon test
Execute_Wilcoxon <- function(object){
  cpm <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))
  groups <- factor(object$label)
  idx <- 1:nrow(cpm)
  names(idx) <- rownames(cpm)
  wilcox_p <- sapply(idx, function(i){
    wilcox.test(cpm[i,]~ groups)$p.value
  })
  FDR <- p.adjust(wilcox_p, method = "BH")
  result_Wilcoxon <- list(gene_names = names(wilcox_p),
                          pvalue = wilcox_p,
                          FDR = FDR)
  return(result_Wilcoxon)
}


# perform T-test
Execute_Ttest <- function(object){
  cpm <- as.matrix(SummarizedExperiment::assay(object, "normcounts"))
  logcpm <- log2(cpm + 1)
  groups <- factor(object$label)
  idx <- seq_len(nrow(logcpm))
  names(idx) <- rownames(logcpm)
  ttest_p <- sapply(idx, function(i){
    t.test(logcpm[i,]~ groups)$p.value
  })
  FDR <- p.adjust(ttest_p, method = "BH")
  result_Ttest <- list(gene_names = names(ttest_p),
                       pvalue = ttest_p,
                       FDR = FDR)
  return(result_Ttest)
}

# perform DGBD
Execute_DGBD <- function(object){
  cpm <- as.matrix(SummarizedExperiment::assay(object, "counts"))
  count<-limma::voom(ROSeq::TMMnormalization(cpm))
  group<-colData(object,"label")
  groups<-group$label
  res<-ROSeq(countData=count, condition = groups)
  result_DGBD <- list(gene_names = row.names(res),
                      pvalue = res[,1],
                      FDR = res[,2])
  return(result_DGBD)
}






