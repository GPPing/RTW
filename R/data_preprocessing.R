preprocessing <- function(Data, group, norm.form = "CPM",  is.normalized = FALSE){
  options(warn = -1)
  
  if(is.normalized){
    
    normcounts <- Data
    gene_df <- data.frame(Gene = rownames(Data))
    cell_df <- data.frame(cell = colnames(Data))
    # pd <- new("AnnotatedDataFrame", data = cell_df)
    # fd <- new("AnnotatedDataFrame", data = gene_df)
    # transfer <- new("CellDataSet", exprs = as.matrix(Data))
    cds <- monocle::newCellDataSet(cellData = Data, expressionFamily = tobit())
    counts_relative <- monocle::relative2abs(cds)
    counts_relative <- floor(counts_relative)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_relative, normcounts = normcounts))
    gene_df <- DataFrame(Gene = rownames(sce))
    cell_df <- DataFrame(label = group, cell = colnames(sce))
    rownames(gene_df) <- gene_df$Gene
    rownames(cell_df) <- cell_df$cell
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_relative, normcounts = normcounts),
                                                      colData = cell_df,
                                                      rowData = gene_df)
    # sce <- scater::calculateQCMetrics(sce)
  } else
  {
    normcounts <- normalized(Data, method = norm.form)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Data, normcounts = normcounts))
    gene_df <- DataFrame(Gene = rownames(sce))
    cell_df <- DataFrame(label = group, cell = colnames(sce))
    rownames(gene_df) <- gene_df$Gene
    rownames(cell_df) <- cell_df$cell
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Data, normcounts = normcounts),
                                                      colData = cell_df,
                                                      rowData = gene_df)
    # sce <- scater::calculateQCMetrics(sce)
  }
  
  return(sce)
}

#' Normalized process
#' The function provide several normalized methods
#' @importFrom edgeR calcNormFactors
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats median
#' @param counts.matrix count expression matrix
#' @param method character. "TMM", "RLE", "CPM". The default value is "CPM".
normalized <- function(counts.matrix, method = "CPM"){
  if(method == "TMM"){
    norm_factor <- edgeR::calcNormFactors(counts.matrix, method = method)
    norm.item <- t(t(counts.matrix)/norm_factor)
    return(norm.item)
  }
  if(method == "RLE"){
    geomeans <- exp(rowMeans(log(counts.matrix)))
    SF <-function(cnts){
      stats::median((cnts/geomeans)[is.finite(geomeans) & geomeans >0])
    }
    norm_factor <- apply(counts.matrix, 2, SF)
    norm.item <- t(t(counts.matrix)/norm_factor)
    return(norm.item)
  }
  if(method == "CPM"){
    # sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix))
    # gene_df <- DataFrame(Gene = rownames(sce))
    # cell_df <- DataFrame(cell = colnames(sce))
    # rownames(gene_df) <- gene_df$Gene
    # rownames(cell_df) <- cell_df$cell
    # sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix),
    #                                                   colData = cell_df,
    #                                                   rowData = gene_df)
    # norm.item <- scater::calculateCPM(sce)
    norm_factor <- colSums(counts.matrix)
    norm.item <- t(t(counts.matrix)/norm_factor) * 1e6
    return(norm.item)
  }
  # if(method == "TPM"){
  #   sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix))
  #   gene_df <- DataFrame(Gene = rownames(sce))
  #   cell_df <- DataFrame(cell = colnames(sce))
  #   rownames(gene_df) <- gene_df$Gene
  #   rownames(cell_df) <- cell_df$cell
  #   sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix),
  #                                                     colData = cell_df,
  #                                                     rowData = gene_df)
  #   norm.item <- scater::calculateTPM(sce, exprs_values = "counts")
  #   return(norm.item)
  # }
}