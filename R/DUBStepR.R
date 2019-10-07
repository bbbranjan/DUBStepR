#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param log.data normalised log-transformed genes x cells single-cell RNA-seq data
#' @param min.cells minimum number of cells to filter genes out and smooth data over
#' @return Returns a list of all feature genes and the optimal feature set
#'
#' @export
#'

DUBStepR <- function(log.data, min.cells = 10) {

    # Filter genes
    log.filt.data <- getfilteredData(data = log.data, min.cells = min.cells)

    # Smooth data
    smooth.log.filt.data <- kNNSmoothing(log.filt.data = log.filt.data, k = min.cells)

    # Compute gene-gene correlation matrix
    ggc <- getGGC(log.data = smooth.log.filt.data)

    # Order genes using stepwise regression
    ordered.genes <- runStepwiseReg(ggc = ggc)

    # Use Compactness Index (CI) to obtain optimal feature set
    feature.genes <- getOptimalFeatureSet(log.data = log.filt.data, ordered.genes = ordered.genes)


    # Return filtered data
    return(list("optimal.feature.genes" = feature.genes, "all.feature.genes" = ordered.genes))
}
