getOptimalFeatureSet <- function(log.data, ordered.genes) {
    
    # Initialise variables
    mean_knn_vec <- c()
    minNumGenes = ""
    numStepsUnchangedMin = 0
    
    # For each neighbour
    for(num_genes in seq(from = 100, to = length(ordered.genes), by = 25)) {
        
        # Run PCA on the feature data
        log.feature.data <-
            log.data[ordered.genes[1:num_genes],]
        pca.data <- irlba::prcomp_irlba(x = t(log.feature.data), n = min(nrow(log.feature.data) - 1, 15), center = TRUE, scale. = FALSE)$x
        rownames(pca.data) <- colnames(log.feature.data)
        
        # Compute k-NN distance
        my.knn <-
            nn2(
                data = pca.data,
                k = k,
                searchtype = "priority",
                eps = 0
            )
        nn.dists <- my.knn$nn.dists
        rownames(nn.dists) <- rownames(pca.data)
        
        # Remove first column as it consists of zeroes
        nn.dists <- nn.dists[, -1]
        
        # Calculate length scale to normalise distances
        sdVec <- apply(X = pca.data,
                       MARGIN = 2,
                       FUN = sd)
        length_scale <- sqrt(sum(sdVec ^ 2))
        
        # Scale k-NN distances by length scale
        mean_nn_dist <- mean(x = nn.dists)
        scaled_mean_nn_dist <- mean_nn_dist / length_scale
        names(scaled_mean_nn_dist) <- num_genes
        
        mean_knn_vec <- append(mean_knn_vec, scaled_mean_nn_dist)
        
        # Check if the minima has been updated
        if(which.min(mean_knn_vec) != minNumGenes) {
            minNumGenes = which.min(mean_knn_vec)
            numStepsUnchangedMin = 0
        } else {
            numStepsUnchangedMin = numStepsUnchangedMin + 1
        }
        
        print(minNumGenes)
        
        # If no updates in the last 5 rounds, break
        if(numStepsUnchangedMin > 5)
            break
    }
    
    # Determine optimal feature set
    optimal_feature_genes <- ordered.genes[1:as.numeric(names(minNumGenes))]
    
    return(optimal_feature_genes)
}