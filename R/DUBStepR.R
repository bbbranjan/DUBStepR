#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param input.data input gene expression matrix (genes x cells)
#' @param min.cells minimum number of cells to filter genes out and smooth data over
#' @param k number of nearest neighbours. Default is 10
#' @param num.pcs number of principal components to represent sc data. Default is 15.
#' @return Returns optimal feature set
#'
#' @export
#'
DUBStepR <- function(input.data, min.cells = 0.05 * ncol(raw.data), k = 10, num.pcs = 15) {

    # Filter genes
    log.filt.data <- getfilteredData(data = input.data, min.cells = min.cells)

    # Smooth data using k-NN
    smooth.log.filt.data <- kNNSmoothing(log.filt.data = log.filt.data, k = k, num.pcs = num.pcs)

    # Compute gene-gene correlation matrix
    ggc <- getGGC(log.data = smooth.log.filt.data)

    # Obtain optimal feature set using stepwise regression
    featureSet <- runStepwiseReg(ggc = ggc, log.filt.data = log.filt.data, num.pcs = num.pcs)

    # Return feature genes
    return(featureSet)
}
