#' Rename Cell Types of Gene Expression Signature
#'
#' Rename columns of a gene expression signature matrix.
#'
#' @param signature gene expression signature
#' @param cell_type cell type to rename
#' @param newName new cell type name
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
renameCellType <- function(signature, cell_type, newName) {
  if (is.null(signature)) {
    stop("Please provide a signature")
  }

  if (is.null(cell_type) | is.null(newName)) {
    stop("cell_type or newName is NULL, cannot rename")
  }

  if (cell_type == "" | newName == "") {
    stop("cell_type or newName empty! Cannot rename")
  }

  if (!(cell_type %in% colnames(signature))) {
    stop("Cannot rename cell type: cell type does not exist in signature")
  }

  newSignature <- as.data.frame(signature)

  names(newSignature)[names(newSignature) == cell_type] <- newName

  return(as.matrix(newSignature))
}

#' Remove Rows from base expression matrix with a specified amount of zeros in a row
#'
#' Remove Rows (Genes) from an expression matrix with more zeros than specified by
#' the percentage threshold parameter. This function aims to reduce sequencing artifacts by
#' removing genes that are only detected in a few cells/samples.
#'
#' @param signature_mat GenesXcelltype Matrix with expression values
#' @param percentage maximum percentage of row values allowed to be 0
#' @returns A signature which matches the criteria above
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#' dim(signature)
#' signature <- removePercentZeros(signature, percentage = 0.5)
#' dim(signature)
removePercentZeros <- function(signature_mat, percentage = 0.5) {
  if (is.null(signature_mat)) {
    stop("Please provide a signature")
  }

  if (percentage > 1 | percentage < 0.00001) {
    stop("Please provide a valid percentage between 0 and 1")
  }

  threshold <- ncol(signature_mat) * percentage # max number of zeroes allowed
  signature <- signature_mat[rowSums(signature_mat == 0) <= threshold, ]

  return(signature)
}


#' Remove unspecific Genes of a Gene Expression Signature
#'
#' Remove genes expressed in an unspecific manner. The expression range is devided into
#' a user selected number of bins. Only genes expressed high in <max_count> celltypes are returned.
#' Genes expressed high in more than <max_count> cell types are discarded.
#'
#' @param signature gene Expression Signature
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
removeUnspecificGenes <- function(signature,
                                  number_of_bins = 3,
                                  max_count = 2,
                                  labels = c("low", "medium", "high")) {
  if (is.null(signature)) {
    stop("Please provide a signature")
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

  signature <- as.matrix(signature)

  to_keep <- vector(length = nrow(signature))

  for (i in 1:nrow(signature)) {
    row <- signature[i, ] # has colnames! drop FALSE is mandatory !!!!!

    # calculate bins to prevent error
    breaks <- seq(floor(min(row)), ceiling(max(row)), length.out = number_of_bins + 1)

    # cut into bins, seperate for each gene
    bins <- cut(row, breaks = breaks, labels = labels, include.lowest = TRUE)

    nHighBins <- sum(bins == "high") # not working when labels is something else

    # this value needs to be greater than one, depending of the step in the pipeline there arent
    # any rows producing zeros left but that is not the case for all  signatures
    if (nHighBins <= max_count & nHighBins > 0) {
      to_keep[i] <- TRUE
    }
  }

  refinedSignature <- signature[to_keep, ]

  # # turn back to a matrix
  # refinedSignature <- as.matrix(refinedSignature)

  return(refinedSignature)
}


#' Select a specified amount of genes for each cell type based on a score, discard all other
#'
#' Reduce the amount of signature genes by selecting the best-scored genes for each cell type.
#' As scoring methods "Entropy" and "Gini" can be applied.
#'
#' @param signature gene expression matrix
#' @param scoring_method method to score the genes ("entropy", "gini")
#' @param selectCellType method to select the cell type the gene is contributing to, used to balance the number of genes between cell types
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
selectGenesByScore <- function(signature,
                               scoring_method = "entropy",
                               selectCellType = "max",
                               genes_per_cell_type = 20) {
  # TODO Checks #####

  # SCORE THE MATRIX
  scoresByCellType <- NULL
  for (celltype in colnames(signature)) {
    scoresByCellType[[celltype]] <- list()
  }


  # not using apply because i need the return values, did not work otherwise
  # scoring is happening here
  for (i in 1:nrow(signature)) {
    row <- signature[i, , drop = FALSE] # has colnames!
    gene <- rownames(row)

    maxCelltype <- colnames(signature)[max.col(row)] # this might be problematic

    score <- list()
    if (scoring_method == "entropy") {
      score[gene] <- scoreEntropy(row) # calculate score and save named result
    } else if (scoring_method == "gini") {
      score[gene] <- 1 - BioQC::gini(row) # need to flip the value since lower scores schould be better (entropy!)
    }

    # append score to most expressed celltype
    scoresByCellType[[maxCelltype]] <- append(scoresByCellType[[maxCelltype]], score)
  }

  # initialize new refined signature with colnames, rows stay empty
  refinedSignature <- matrix(nrow = 0, ncol = length(colnames(signature)), dimnames = list(NULL, colnames(signature)))

  # for each celltype get the lowest scores
  for (celltype in names(scoresByCellType)) {
    # get best genes for each score, sorted ascending
    # reminder: the lower the score the better in entropy case, might need a change to fit all
    scores <- scoresByCellType[[celltype]] %>%
      unlist() %>%
      sort()

    # check number of genes!!!
    if (genes_per_cell_type > length(names(scores))) {
      genes_per_cell_type <- length(names(scores))
    }

    # iterate genes
    for (gene in names(scores)[1:genes_per_cell_type]) {
      row <- signature[gene, , drop = FALSE]

      refinedSignature <- rbind(refinedSignature, row)
    }
  }

  return(refinedSignature)
}

#' Score Gene Expression of a single Gene based on information entropy
#'
#' Score Genes Expression of a single gene across celltypes. The function returns
#' the calculated entropy of the expression value distribution.
#'
#' @param expression_feature row from Gene Expression Matrix = Expression Data for a single Gene
#' @returns Score for the given gene based on information entropy
#' Here: The lower the better
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#'
#' entropy <- scoreEntropy(signature[1, ]) # scoring the first gene
scoreEntropy <- function(expression_feature) {
  # TODO add parameter checks ####
  probs <- list()

  # turn expression data to a list of probabilities
  for (val in expression_feature) {
    if (val == 0) {
      next
    }
    probs <- append(probs, val / sum(expression_feature)) # turn in to probabilities
  }

  entropy <- -sum(unlist(lapply(probs, function(x) log(x) * x)))

  return(entropy)
}
