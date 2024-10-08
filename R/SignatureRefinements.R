#' Rename Cell Types of Gene Expression Signature
#'
#' Rename columns of a gene expression signature matrix.
#'
#' @param signature_mat gene expression signature
#' @param cell_type cell type to rename
#' @param new_celltype_name new cell type name
#'
#' @returns gene expression signature with updated cell type names
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#'
#' # rename "B" to "B.cells"
#' signature <- renameCellType(signature, "B", "B.cells")
renameCellType <- function(signature_mat, cell_type, new_celltype_name) {
  if (is.null(signature_mat)) {
    stop("Please provide a signature_mat")
  }

  if (is.null(cell_type) | is.null(new_celltype_name)) {
    stop("cell_type or new_celltype_name is NULL, cannot rename")
  }

  if (cell_type == "" | new_celltype_name == "") {
    stop("cell_type or new_celltype_name empty! Cannot rename")
  }

  if (!(cell_type %in% colnames(signature_mat))) {
    stop("Cannot rename cell type: cell type does not exist in signature_mat")
  }

  newSignature <- as.data.frame(signature_mat)

  names(newSignature)[names(newSignature) == cell_type] <- new_celltype_name

  return(as.matrix(newSignature))
}

#' Remove Rows from base expression matrix with a specified amount of zeros in a row
#'
#' Remove Rows (Genes) from an expression matrix with more zeros than specified by
#' the percentage threshold parameter. This function aims to reduce sequencing artifacts by
#' removing genes that are only detected in a few cells/samples.
#'
#' @param signature_mat GenesXcelltype Matrix with expression values
#' @param max_percentage_zeroes maximum percentage of row values allowed to be 0
#' @returns A signature which matches the criteria above
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#' dim(signature)
#' signature <- removePercentZeros(signature, max_percentage_zeroes = 0.5)
#' dim(signature)
removePercentZeros <- function(signature_mat, max_percentage_zeroes = 0.5) {
  if (is.null(signature_mat)) {
    stop("Please provide a signature")
  }

  if (max_percentage_zeroes > 1 | max_percentage_zeroes < 0.00001) {
    stop("Please provide a valid percentage between 0 and 1")
  }

  threshold <- ncol(signature_mat) * max_percentage_zeroes # max number of zeroes allowed
  signature <- signature_mat[rowSums(signature_mat == 0) <= threshold, ]

  return(signature)
}


#' Remove unspecific Genes of a Gene Expression Signature
#'
#' Remove genes expressed in an unspecific manner. The expression range is divided into
#' a user selected number of bins. Only genes expressed high in <max_count> cell types are returned.
#' Genes expressed high in more than <max_count> cell types are discarded.
#'
#' @param signature_mat gene Expression Signature
#' @param number_of_bins number of bins to categorize the data into
#' @param max_count number of Cell Types allowed to be in the highest bin,
#' all other cells are required to be in lower expressed bins
#' @param labels vector of bin names, required if number_of_bins != 3
#'
#' @returns a gene expression signature containing only genes matching the passed requirements
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#' dim(signature)
#'
#' signature <- removeUnspecificGenes(signature, number_of_bins = 3, max_count = 1)
#' dim(signature)
removeUnspecificGenes <- function(signature_mat,
                                  number_of_bins = 3,
                                  max_count = 2,
                                  labels = c("low", "medium", "high")) {
  if (is.null(signature_mat)) {
    stop("Please provide a signature_mat")
  }

  if (number_of_bins < 2) {
    stop("number_of_bins has to be >= 2!")
  }

  if (max_count <= 0) {
    stop("max_count has to be a positive integer")
  }

  if (length(labels) != number_of_bins) {
    stop("number_of_bins does not match label length")
  }

  signature_mat <- as.matrix(signature_mat)

  to_keep <- sapply(1:nrow(signature_mat), function(i){
    row <- signature_mat[i, ] # has colnames! drop FALSE is mandatory !!!!!
    
    # calculate bins to prevent error
    breaks <- seq(floor(min(row)), ceiling(max(row)), length.out = number_of_bins + 1)
    
    # cut into bins, seperate for each gene
    bins <- cut(row, breaks = breaks, labels = labels, include.lowest = TRUE)
    
    nHighBins <- sum(bins == "high") # not working when labels is something else
    
    # this value needs to be greater than one, depending of the step in the pipeline there arent
    # any rows producing zeros left but that is not the case for all signatures
    if (nHighBins <= max_count & nHighBins > 0) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  

  refinedSignature <- signature_mat[to_keep, ]

  # # turn back to a matrix
  # refinedSignature <- as.matrix(refinedSignature)

  return(refinedSignature)
}


#' Select a specified amount of genes for each cell type based on a score, discard all other
#'
#' Reduce the amount of signature genes by selecting the best-scored genes for each cell type.
#' As scoring methods "Entropy" and "Gini" can be applied.
#'
#' @param signature_mat gene expression matrix
#' @param scoring_method method to score the genes ("entropy", "gini")
#' @param select_celltype method to select the cell type the gene is contributing to, used to balance the number of genes between cell types
#' @param genes_per_cell_type maximum of genes selected for each cell type
#'
#' @return A data frame with the compacted signatures
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#' dim(signature)
#'
#' signature <- selectGenesByScore(signature, "gini", genes_per_cell_type = 50)
#' dim(signature)
selectGenesByScore <- function(signature_mat,
                               scoring_method = "entropy",
                               select_celltype = "max",
                               genes_per_cell_type = 20) {
  # TODO Checks #####

  # SCORE THE MATRIX
  scoresByCellType <- NULL
  for (celltype in colnames(signature_mat)) {
    scoresByCellType[[celltype]] <- list()
  }


  # not using apply because i need the return values, did not work otherwise
  # scoring is happening here
  for (i in 1:nrow(signature_mat)) {
    row <- signature_mat[i, , drop = FALSE] # has colnames!
    gene <- rownames(row)

    maxCelltype <- colnames(signature_mat)[max.col(row)] # this might be problematic

    score <- list()
    if (scoring_method == "entropy") {
      # score[gene] <- scoreEntropy(row) # calculate score and save named result
      score[gene] <- BioQC::entropySpecificity(rbind(row, row))[1]
    } else if (scoring_method == "gini") {
      score[gene] <- 1 - BioQC::gini(row) # need to flip the value since lower scores schould be better (entropy!)
    }

    # append score to most expressed celltype
    scoresByCellType[[maxCelltype]] <- append(scoresByCellType[[maxCelltype]], score)
  }

  # initialize new refined signature with colnames, rows stay empty
  refinedSignature <- matrix(nrow = 0, ncol = length(colnames(signature_mat)), dimnames = list(NULL, colnames(signature_mat)))

  # for each celltype get the lowest scores
  for (celltype in names(scoresByCellType)) {
    # get best genes for each score, sorted ascending
    # reminder: the lower the score the better in entropy case, might need a change to fit all
    scores <- sort(unlist(scoresByCellType[[celltype]]))

    # check number of genes!!!
    if (genes_per_cell_type > length(names(scores))) {
      genes_per_cell_type <- length(names(scores))
    }

    # iterate genes
    for (gene in names(scores)[1:genes_per_cell_type]) {
      row <- signature_mat[gene, , drop = FALSE]

      refinedSignature <- rbind(refinedSignature, row)
    }
  }

  return(refinedSignature)
}
