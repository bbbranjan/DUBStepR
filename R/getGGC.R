#' @title Compute the correlation range values for all genes in the gene-gene correlation matrix
#' @param data log-transformed gene-expression matrix
#' @return list of genes with their z-transformed correlation range values
#'
#' @export
#'
getGGC <- function(data) {

    # Bin genes by mean expression
    num.bins = min(20, (nrow(data)-1))
    gene.mean <- Matrix::rowMeans(data)

    gene.bins <- cut(x = gene.mean, breaks = num.bins)
    names(gene.bins) <- names(gene.mean)

    # Run fastCor for quicker correlation matrix computation
    if(require(HiClimR)) {
        correlation_matrix <- HiClimR::fastCor(xt = t(as.matrix(data)), nSplit = 10, upperTri = T, optBLAS = T, verbose = F)
        correlation_matrix[which(is.na(correlation_matrix))] <- 0
        correlation_matrix <- correlation_matrix + t(correlation_matrix)
        diag(correlation_matrix) <- 1
    } else {
        correlation_matrix <- cor(x = as.matrix(t(data)), method = "pearson")
    }



    # Obtain correlation range
    rangeObj <- getCorrelationRange(correlation_matrix = correlation_matrix)

    logRange <- log(1+rangeObj$range)

    # Z-transform correlation range
    zRangeList <- tapply(X = logRange, INDEX = gene.bins, FUN = function(x){
        (x - mean(x))/stats::sd(x)
    })

    if(all(is.na(zRangeList))) {
        stop("Feature correlation range could not be obtained.")
    }

    # Set NAs to 0
    zRangeList <- lapply(zRangeList, function(x){
        x[is.na(x)] <- 0
        return(x)
    })
    zRange <- unlist(zRangeList)
    names(zRange) <- unlist(lapply(strsplit(x = names(zRange), split = "\\]\\."), "[[", 2))

    # Obtain geneset with zRange > 0.7
    topGenes <- names(which(zRange > 0.7))

    # Reduce the GGC to use only these genes
    ggc <- correlation_matrix[topGenes, topGenes]

    # Return as matrix
    return(list("corr.range" = zRange, "ggc" = as.matrix(ggc)))
}

#' @title Compute the correlation range values for all genes in the gene-gene correlation matrix
#' @param correlation_matrix gene-gene correlation matrix
#' @return list of p-values, adjusted p-values and correlation ranges for each gene
getCorrelationRange <- function(correlation_matrix) {

    # Sort each column of the correlation matrix
    sorted_corr = apply(correlation_matrix, 2, sort, decreasing = T)

    # Select the second largest correlation for each gene
    max_corr = sorted_corr[3, ]

    # Select the smallest correlation for each gene
    min_corr = sorted_corr[nrow(sorted_corr), ]

    # Compute the correlation range
    diff_corr = max_corr-(0.75*min_corr)
    # names(diff_corr) <- names(rank_2nd_corr)

    # Sort genes in descending order of correlation range
    diff_corr <- sort(diff_corr, decreasing = T)

    # Return as list
    return(list("range"=diff_corr))
}
