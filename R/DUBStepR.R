#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param input.data input gene expression matrix (genes x cells)
#' @param min.cells minimum number of cells to filter genes out and smooth data over
#' @param optimise.features Determine optimal feature set using density index. (Time-consuming step).
#' @param k number of nearest neighbours. Default is 10.
#' @param num.pcs number of principal components to represent sc data. Default is 15.
#' @param error Acceptable error margin for kNN computation. Default is 0.
#' @return Returns optimal feature set
#'
#' @export
#'
DUBStepR <- function(input.data, min.cells = 100, optimise.features = F, k = 10, num.pcs = 15, error = 0) {

    # Filter genes
    filt.data <- getfilteredData(data = input.data, min.cells = min.cells)

    # Compute gene-gene correlation matrix
    ggc.out <- getGGC(data = filt.data)

    # Run stepwise regression
    swreg.out <- runStepwiseReg(ggc = ggc.out$ggc, filt.data = filt.data)
    corr.info <- data.frame(feature.genes = swreg.out$feature.genes, corr.range = ggc.out$corr.range[swreg.out$feature.genes])

    if(optimise.features) {
        # Obtain optimal feature set
        density.out <- getOptimalFeatureSet(filt.data = filt.data, ordered.genes = swreg.out$feature.genes, elbow.pt = swreg.out$elbow.pt)

        dubStepR.out <- list("corr.info" = corr.info, "optimal.feature.genes" = density.out$optimal.feature.genes, "density.index" = density.out$density.index)
    } else {
        # Return all correlated genes
        dubStepR.out <- list("corr.info" = corr.info)
    }


    # Return DUBStepR output list
    return(dubStepR.out)
}
