#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param log.data normalised log-transformed genes x cells single-cell RNA-seq data
#' @param is_noisy TRUE if data is noisy, FALSE otherwise.
#' @return
#'
#' @export
#'

DUBStepR <- function(log.data, is_noisy = FALSE) {

    filt.data <- getfilteredData(data = log.data)
    ggc <- getGGC(log.data = filt.data, is_noisy = is_noisy)
    ordered.genes <- runStepwiseReg(ggc = ggc)
    feature.genes <- getOptimalFeatureSet(log.data = log.data, ordered.genes = ordered.genes)


    # Return filtered data
    return(log.data[feature.genes, ])
}


### Load packages
if (!require(WGCNA)) {
    source("http://bioconductor.org/biocLite.R")

    biocLite(c("impute", "GO.db", "preprocessCore"))

    install.packages("WGCNA")

}
require(WGCNA)
if (!require(tibble))
    install.packages("tibble", repos = "http://cran.us.r-project.org")
require(tibble)
if (!require(Matrix))
    install.packages("Matrix", repos = "http://cran.us.r-project.org")
require(Matrix)
if (!require(MASS))
    install.packages("MASS", repos = "http://cran.us.r-project.org")
require(matrixcalc)
if (!require(matrixcalc))
    install.packages("matrixcalc", repos = "http://cran.us.r-project.org")
require(matrixcalc)
if (!require(Rfast))
    install.packages("Rfast", repos = "http://cran.us.r-project.org")
require(Rfast)
if (!require(flashClust))
    install.packages("flashClust", repos = "http://cran.us.r-project.org")
require(flashClust)
library(Hmisc)
if (!require(Seurat))
    install.packages("Seurat", repos = "http://cran.us.r-project.org")
require(Seurat)
if (!require(rtracklayer)) {
    BiocManager::install(pkgs = "rtracklayer")
}
require(rtracklayer)
if (!require(RANN))
    install.packages("RANN", repos = "http://cran.us.r-project.org")
require(RANN)
if (!require(HiClimR))
    install.packages("HiClimR", repos = "http://cran.us.r-project.org")
require(HiClimR)
if (!require(qlcMatrix))
    install.packages("qlcMatrix", repos = "http://cran.us.r-project.org")
require(qlcMatrix)

require(parallel)
