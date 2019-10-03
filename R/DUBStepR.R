#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param log.data normalised log-transformed genes x cells single-cell RNA-seq data
#' @param is_noisy TRUE if data is noisy, FALSE otherwise.
#' @return Returns a list of all feature genes and the optimal feature set
#'
#' @export
#'

DUBStepR <- function(log.data, is_noisy = FALSE) {

    # Filter genes
    log.filt.data <- getfilteredData(data = log.data)

    # Smooth data
    smooth.log.filt.data <- kNNSmoothing(log.filt.data = log.filt.data)

    # Compute gene-gene correlation matrix
    ggc <- getGGC(log.data = smooth.log.filt.data, is_noisy = is_noisy)

    # Order genes using stepwise regression
    ordered.genes <- runStepwiseReg(ggc = ggc)

    # Use Compactness Index (CI) to obtain optimal feature set
    feature.genes <- getOptimalFeatureSet(log.data = log.filt.data, ordered.genes = ordered.genes)


    # Return filtered data
    return(log.data[feature.genes, ])
}
