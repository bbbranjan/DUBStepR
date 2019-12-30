#' @title k-NN Smoothing of gene expression data
#' @param filt.data filtered and normalised log-transformed genes x cells single-cell RNA-seq data matrix
#' @param k number of nearest neighbours the data is smoothed over
#' @param num.pcs number of principal components to represent sc data. Default is 15.
#' @return k-NN smoothed data matrix
#'
#' @export
#'

kNNSmoothing <- function(filt.data, k = 1, num.pcs = 15) {

    # Sanity checks
    if(k > ncol(filt.data)) {
        stop("k is too large for the input data.")
    }

    if(k == 1) {
        print("No kNN Smoothing required.")
        return(filt.data)
    }

    if(num.pcs > nrow(filt.data)) {
        stop("num.pcs is too large for the input data.")
    }

    # Reduce expression data to PC space for k-NN computation
    pca.data <-
        irlba::prcomp_irlba(
            x = t(as.matrix(filt.data)),
            n = num.pcs,
            center = TRUE,
            scale. = FALSE
        )$x
    rownames(pca.data) <- colnames(filt.data)

    # Determine k-nearest neighbours
    my.knn <-
        RANN::nn2(
            data = pca.data,
            k = (k+1),
            searchtype = "standard",
            eps = 0
        )

    # Extract indices of k-nearest neighbours
    nn.idx <- my.knn$nn.idx

    # Discard self as nearest neighbour and set cell names
    nn.idx <- nn.idx[, -1]
    rownames(nn.idx) <- rownames(pca.data)

    # Replace nearest neighbour indices with cell names
    nn.mat <- t(apply(nn.idx, 1, function(neighbours){
        rownames(nn.idx)[neighbours]
    }))
    colnames(nn.mat) <- paste0(seq(1,ncol(nn.mat)))

    # Melt nearest neighbour matrix
    nn.melt <- reshape2::melt(data = nn.mat)
    nn.melt$Var2 <- NULL

    # Construct adjacency matrix
    adj.mat <- igraph::get.adjacency(graph = igraph::graph.edgelist(as.matrix(nn.melt), directed = TRUE))

    # Compute smoothed data based on adjacency matrix
    smooth.filt.data <- Matrix::t(Matrix::tcrossprod(x = adj.mat, y = filt.data)/k)

    print("kNN Smoothing - Done")

    # Return smoothed data
    return(smooth.filt.data)
}
