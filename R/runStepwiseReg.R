#' @title Run step-wise regression to order the features
#' @param ggc gene-gene correlation matrix
#' @param log.filt.data filtered and normalised log-transformed genes x cells single-cell RNA-seq data matrix
#' @param k number of nearest neighbours for CI computation
#' @return optimal feature set
#'
#' @export
#'
runStepwiseReg <- function(ggc, log.filt.data, k = 10) {

    # Initialize variables
    ggc_centered <- scale(x = ggc, center = TRUE, scale = FALSE)
    step = 1
    num_steps = 100
    step_seq = seq(from = step, to = num_steps, by = step)
    scree_values <- c(matrixcalc::frobenius.norm(ggc_centered))
    names(scree_values) <- c("0")
    feature_genes <- c()

    # For each step in the stepwise regression process
    for(i in step_seq) {

        # Print statement
        print(paste("Num genes:", i))

        # Compute GGC'*GGC
        ggc_ggc <- Rfast::mat.mult(t(ggc_centered), ggc_centered)
        dimnames(ggc_ggc) <- dimnames(ggc_centered)

        # Compute variance explained
        gcNormVec <- apply(ggc_ggc, 1, function(x){
            matrixcalc::frobenius.norm(x)
        })

        # Compute norm of gene vectors
        gSeq = 1:ncol(ggc_ggc)

        gNormVec <- sapply(gSeq, function(idx){
            sqrt(ggc_ggc[idx,idx])
        })

        # Select gene to regress out
        varExpVec <- gcNormVec/gNormVec
        regressed.genes <- names(sort(varExpVec, decreasing = T))[1:step]

        # Add regressed gene to feature set
        feature_genes <- union(feature_genes, regressed.genes)

        # Obtain the variance explained by the regressed genes
        regressed.g <- ggc_centered[, regressed.genes]
        gtg <- as.numeric(t(regressed.g) %*% regressed.g)
        explained <- (regressed.g %*% t(ggc_ggc[regressed.genes,]))/gtg
        eps <- ggc_centered - explained

        # Append variance explained value to vector for scree plot
        scree_values[paste0(i)] <- sum(varExpVec[regressed.genes])

        # Update the GGC to be the residual after regressing out these genes
        ggc_centered = eps

    }

    # Plot scree plot to see change in Frobenius norm over genes
    scree_values <- scree_values[which(scree_values != 0)]

    # Find elbow point
    elbow_id <- findElbow(y = log(scree_values)[1:100], ylab = "Log Variance Explained")

    # Initialise variables to add neighbours
    elbow_feature_genes <- feature_genes[1:elbow_id]
    neighbour_feature_genes <- elbow_feature_genes

    neighbour_fg_list <- as.list(elbow_feature_genes)
    names(neighbour_fg_list) <- elbow_feature_genes

    # Get list of GGC
    ggc_list <- apply(ggc, 2, function(x) {
        list(x)
    })
    ggc_list <- lapply(X = ggc_list, unlist)
    ggc_list <- lapply(X = ggc_list,
                       FUN = sort,
                       decreasing = TRUE)
    ggc_list <- lapply(ggc_list, function(x) {
        x[!(names(x) %in% neighbour_feature_genes)]
    })

    # Select potential candidates for next neighbour
    candidateGGCList <-
        lapply(ggc_list[neighbour_feature_genes], function(x) {
            x[which.max(x)]
        })

    # Use Compactness Index (CI) as stopping criterion
    mean_knn_vec <- c()
    numStepsUnchangedMin = 0
    minNumGenes = ""

    # Adding neighbours of each gene
    for (i in 1:(nrow(ggc) - length(neighbour_feature_genes))) {

        # Select potential candidates for next neighbour
        candidateGGCNames <-
            unlist(lapply(candidateGGCList, names))

        candidateGGCVec <-
            unlist(candidateGGCList, use.names = F)
        names(candidateGGCVec) <- candidateGGCNames

        # Select candidate with largest correlation to feature set
        nearest.index <- which.max(candidateGGCVec)
        whoseCandidate <- names(candidateGGCList[nearest.index])
        nearestNeighbour <- candidateGGCVec[nearest.index]

        # Add new neighbour to feature set, and remove from next neighbour candidate list
        neighbour_feature_genes <-
            append(neighbour_feature_genes, names(nearestNeighbour))

        # Update GGC list
        ggc_list <- lapply(ggc_list, function(x) {
            x[!(names(x) == names(nearestNeighbour))]
        })

        # Update list of nearest neighbour candidates
        candidateGGCList <-
            lapply(ggc_list[neighbour_feature_genes], function(x) {
                x[which.max(x)]
            })

        # For every 25 genes added
        if((i%%25 == 0 | length(neighbour_feature_genes) == nrow(ggc))) {

            # Initialise number of genes
            num_genes = length(neighbour_feature_genes)

            # Run PCA on the feature data
            log.feature.data <-
                log.filt.data[neighbour_feature_genes,]
            pca.data <- irlba::prcomp_irlba(x = Matrix::t(log.feature.data), n = 15, center = TRUE, scale. = FALSE)$x
            rownames(pca.data) <- colnames(log.feature.data)

            # Compute k-NN distance
            system.time(
                my.knn <- RANN::nn2(
                    data = pca.data,
                    k = (k+1),
                    treetype = "kd",
                    searchtype = "standard",
                    eps = 0
                )
            )

            # system.time(
            #     nn.dists <- FNN::knn.dist(
            #         data = pca.data,
            #         k = 11,
            #         algorithm = "CR"
            #     )
            # )
            nn.dists <- my.knn$nn.dists
            rownames(nn.dists) <- rownames(pca.data)

            # Remove first column as it consists of zeroes
            nn.dists <- nn.dists[, -1]

            # Calculate length scale to normalise distances
            sdVec <- apply(X = pca.data,
                           MARGIN = 2,
                           FUN = stats::sd)
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

        }
    }

    # Determine optimal feature set
    optimal_feature_genes <- neighbour_feature_genes[1:as.numeric(names(minNumGenes))]

    # Return
    return(optimal_feature_genes)
}
