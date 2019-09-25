#' @title Filter the dataset by removing lowly expressed genes and mitochondrial, spike-in and ribosomal genes
#' @param data gene expression matrix
#' @param min.cells gene expression matrix
#' @return filtered gene-expression matrix
#'
#' @export
#'
getfilteredData <- function(data, min.cells = 0.05*ncol(data)) {

    # Filter data by gene expression
    filt.data <- data[apply(data, 1, function(x){
        sum(x>0) > 0.05*ncol(data)
    }), ]


    print("Expression Filtering - Done")

    # Remove mitochondrial, spike-in and ribosomal genes as they do not assist in cell type separation
    filt.data <- filt.data[!grepl(pattern = "^MT-|^ERCC-|^RPS|^RPL", x = rownames(filt.data)), ]

    # Remove pseudogenes
    data(pseudo_genes, envir=environment())
    filt.data <- filt.data[!(rownames(filt.data) %in% pseudo_genes),]

    print("Mitochondrial, Ribosomal and Pseudo Genes Filtering - Done")

    return(filt.data)
}
