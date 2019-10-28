#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param raw.data raw gene expression matrix (genes x cells)
#' @param log.normalize boolean indicating whether input matrix needs to be log-normalized. Default is TRUE.
#' @param min.cells minimum number of cells to filter genes out and smooth data over
#' @param k number of nearest neighbours
#' @return Returns optimal feature set
#'
#' @export
#'
DUBStepR <- function(raw.data, log.normalize = T, min.cells = 0.05 * ncol(raw.data), k = 10) {
    # Log-normalize data
    log.data <- logNormalize(raw.data = raw.data)

    # Filter genes
    log.filt.data <- getfilteredData(data = log.data, min.cells = min.cells)

    # Smooth data using k-NN
    smooth.log.filt.data <- kNNSmoothing(log.filt.data = log.filt.data, k = k)

    # Compute gene-gene correlation matrix
    ggc <- getGGC(log.data = smooth.log.filt.data)

    # Obtain optimal feature set using stepwise regression
    featureSet <- runStepwiseReg(ggc = ggc, log.filt.data = log.filt.data)

    # Return feature genes
    return(featureSet)
}
