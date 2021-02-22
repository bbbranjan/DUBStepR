#' @title Filter the dataset by removing lowly expressed genes and mitochondrial, spike-in and ribosomal genes
#' @param data gene expression matrix
#' @param min.cells gene expression matrix
#' @return filtered gene-expression matrix
#'
#' @export
#'
getFilteredData <- function(data, min.cells = 0.05*ncol(data)) {
    message(" ")
    message("Running DUBStepR...")
    
    # Filter data by gene expression
    filt.data <- data[Matrix::rowSums(data > 0) > min.cells, ]
    
    message(" ")
    message("Expression Filtering Done ✓")
    
    # Remove mitochondrial, spike-in and ribosomal genes as they do not assist in cell type separation
    filt.data <- filt.data[!grepl(pattern = "^MT-|^ERCC-|^RPS|^RPL", x = rownames(filt.data)), ]
    
    # Remove pseudogenes
    filt.data <- filt.data[!(rownames(filt.data) %in% pseudo_genes),]
    
    message("Mitochondrial, Ribosomal and Pseudo Genes Filtering Done ✓")
    
    return(filt.data)
}
