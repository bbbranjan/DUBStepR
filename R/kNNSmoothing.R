#' @title k-NN Smoothing of gene expression data
#' @param log.filt.data filtered and normalised log-transformed genes x cells single-cell RNA-seq data matrix
#' @param k number of nearest neighbours the data is smoothed over
#' @return k-NN Smoothed data matrix
#'
#' @export
#'

kNNSmoothing <- function(log.filt.data, k = 11) {

    pca.data <-
        irlba::prcomp_irlba(
            x = t(as.matrix(log.filt.data)),
            n = 15,
            center = TRUE,
            scale. = FALSE
        )$x
    rownames(pca.data) <- colnames(log.filt.data)

    my.knn <-
        nn2(
            data = pca.data,
            k = k,
            searchtype = "priority",
            eps = 0
        )
    nn.idx <- my.knn$nn.idx
    nn.idx <- nn.idx[, -1]
    rownames(nn.idx) <- rownames(pca.data)

    smooth.log.filt.data <- apply(nn.idx, 1, function(neighbours) {
        Matrix::rowMeans(log.filt.data[, neighbours])
    })

    return(smooth.log.filt.data)
}
