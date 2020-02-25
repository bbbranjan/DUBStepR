#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param input.data input gene expression matrix (genes x cells)
#' @param min.cells minimum number of cells to filter genes out and smooth data over
#' @param k number of nearest neighbours. Default is 10 - no smoothing.
#' @param smooth Boolean determining whether smoothing is necessary or not.
#' @param num.pcs number of principal components to represent sc data. Default is 15.
#' @param error Acceptable error margin for kNN computation. Default is 0.
#' @return Returns optimal feature set
#'
#' @export
#'
DUBStepR <- function(input.data, min.cells = 100, k = 10, smooth = F, num.pcs = 15, error = 0) {

    # Filter genes
    filt.data <- getfilteredData(data = input.data, min.cells = min.cells)

    # Smooth data using k-NN
    if(smooth)
        smooth.filt.data <- kNNSmoothing(filt.data = filt.data, k = k, num.pcs = num.pcs)
    else {
        print("No kNN Smoothing required.")
        smooth.filt.data <- filt.data
    }


    # Compute gene-gene correlation matrix
    ggc.out <- getGGC(data = smooth.filt.data)

    # Obtain optimal feature set using stepwise regression
    swreg.out <- runStepwiseReg(ggc = ggc.out$ggc, filt.data = filt.data, k = k, num.pcs = num.pcs, error = error)

    corr.info <- data.frame(feature.genes = swreg.out$feature.genes, corr.range = ggc.out$corr.range[swreg.out$feature.genes])
    dubStepR.out <- list("corr.info" = corr.info, "optimal.feature.genes" = swreg.out$optimal.feature.genes, "density.index" = swreg.out$density.index)

    # Return feature genes
    return(dubStepR.out)
}
