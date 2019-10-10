#' @title Run step-wise regression to order the features
#' @param ggc gene-gene correlation matrix
#' @return list of ordered genes
#'
#' @export
#'
runStepwiseReg <- function(ggc) {

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

        varExpVec <- gcNormVec/gNormVec

        regressed.genes <- names(sort(varExpVec, decreasing = T))[1:step]
        # print(regressed.genes)
        feature_genes <- union(feature_genes, regressed.genes)

        # Obtain the variance explained by the regressed genes
        regressed.g <- ggc_centered[, regressed.genes]
        gtg <- as.numeric(t(regressed.g) %*% regressed.g)

        explained <- (regressed.g %*% t(ggc_ggc[regressed.genes,]))/gtg
        eps <- ggc_centered - explained

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

    # Average 6-NN distance as a stopping criterion
    k = 7
    mean_knn_vec <- c()


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
    }

    # Return ordered genes
    return(neighbour_feature_genes)
}
