#' @objective to compute the correlation range values for all genes in the gene-gene correlation matrix
#' @param correlation_matrix gene-gene correlation matrix
#' @param log.data log-transformed gene-expression matrix
#' @return list of genes with their z-transformed correlation range values 
getGGC <- function(log.data, is_noisy = FALSE) {
    
    # Bin genes by mean expression
    num.bins = 20
    gene.mean <- rowMeans(log.data)
    
    gene.bins <- cut(x = gene.mean, breaks = num.bins)
    names(gene.bins) <- names(gene.mean)
    
    rangeObj <- getCorrelationRange(log.data = log.data)
    
    logRange <- log(1+rangeObj$range)
    
    # Z-transform correlation range
    zRangeList <- tapply(X = logRange, INDEX = gene.bins, FUN = function(x){
        (x - mean(x))/sd(x)
    })
    
    zRangeList <- lapply(zRangeList, function(x){
        x[is.na(x)] <- 0
        return(x)
    })
    
    zRange <- unlist(zRangeList)
    names(zRange) <- unlist(lapply(strsplit(x = names(zRange), split = "\\]\\."), "[[", 2))
    
    geneset <- names(which(zRange > 0.7))
    
    topGenes <- geneset
    
    if(is_noisy) {
        seurat_obj <- CreateSeuratObject(counts = filt.data)
        seurat_obj <- NormalizeData(object = seurat_obj)
        seurat_obj <- FindVariableFeatures(object = seurat_obj, selection.method = "dispersion")
        hvg.info = seurat_obj@assays$RNA@meta.features
        hvg.info <- hvg.info[order(hvg.info$mvp.dispersion.scaled, decreasing = T), ]
        
        topHVGenes <- rownames(subset(hvg.info, mvp.dispersion.scaled > 1))
        topGenes <- union(topGenes, topHVGenes)
    }
    
    # Reduce the GGC to use only these genes
    ggc <- fastCor(xt = t(as.matrix(log.data[topGenes, ])), nSplit = 5, upperTri = F, optBLAS = T)
    
    # Return as matrix
    return(as.matrix(ggc))
}

#' @objective to compute the correlation range values for all genes in the gene-gene correlation matrix
#' @param correlation_matrix gene-gene correlation matrix
#' @return list of p-values, adjusted p-values and correlation ranges for each gene 
getCorrelationRange <- function(log.data) {
    
    # Run fastCor for quicker correlations
    system.time({
        correlation_matrix <- fastCor(xt = t(as.matrix(log.data)), nSplit = 10, upperTri = T, optBLAS = T)
    })
    correlation_matrix[which(is.na(correlation_matrix))] <- 0
    correlation_matrix <- correlation_matrix + t(correlation_matrix)
    diag(correlation_matrix) <- 1
    
    
    # Select the second largest correlation for each gene
    max_corr = apply(correlation_matrix, 2, Rfast::nth, 3, descending = T)
    
    # Select the smallest correlation for each gene
    min_corr = apply(correlation_matrix,2,Rfast::nth, 1)
    
    # Compute the correlation range
    diff_corr = max_corr-(0.75*min_corr)
    # names(diff_corr) <- names(rank_2nd_corr)
    
    # Sort genes in descending order of correlation range
    diff_corr <- sort(diff_corr, decreasing = T)
    
    # Return as list
    return(list("range"=diff_corr))
}