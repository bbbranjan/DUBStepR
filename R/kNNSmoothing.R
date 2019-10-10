#' @title k-NN Smoothing of gene expression data
#' @param log.filt.data filtered and normalised log-transformed genes x cells single-cell RNA-seq data matrix
#' @param k number of nearest neighbours the data is smoothed over
#' @return k-NN smoothed data matrix
#'
#' @export
#'

kNNSmoothing <- function(log.filt.data, k = 11) {

    # Reduce expression data to PC space for k-NN computation
    pca.data <-
        irlba::prcomp_irlba(
            x = t(as.matrix(log.filt.data)),
            n = 15,
            center = TRUE,
            scale. = FALSE
        )$x
    rownames(pca.data) <- colnames(log.filt.data)

    # Determine k-nearest neighbours
    my.knn <-
        RANN::nn2(
            data = pca.data,
            k = k,
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
    smooth.log.filt.data <- t(tcrossprod(x = adj.mat, y = log.filt.data)/(k-1))

    print("kNN Smoothing - Done")

    # Return smoothed data
    return(smooth.log.filt.data)
}
