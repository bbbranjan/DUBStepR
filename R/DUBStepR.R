#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param log.data normalised log-transformed genes x cells single-cell RNA-seq data
#' @param min.cells minimum number of cells to filter genes out and smooth data over
#' @return Returns a list of all feature genes and the optimal feature set
#'
#' @export
#'

DUBStepR <- function(log.data, min.cells = 0.01*ncol(log.data)) {

    # Filter genes
    log.filt.data <- getfilteredData(data = log.data, min.cells = min.cells)

    # Smooth data using k-NN
    smooth.log.filt.data <- kNNSmoothing(log.filt.data = log.filt.data)

    # Compute gene-gene correlation matrix
    ggc <- getGGC(log.data = smooth.log.filt.data)

    # Order genes and obtain optimal feature set using stepwise regression
    geneList <- runStepwiseReg(ggc = ggc, smooth.log.filt.data = smooth.log.filt.data)

    # Return list of feature genes
    return(geneList)
}
