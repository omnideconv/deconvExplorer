

#' Remove Rows from base expression matrix with a specified amount of zeros in a row
#' 
#' @param baseSignatue GenesXcelltype Matrix with expression values
#' @param percentage maximum percentage of row values allowed to be 0
#' @returns A signature which matches the criteria above
removePercentZeros <- function (baseSignature, percentage = 0.5){
  if(is.null(baseSignature)){
    stop("Please provide a signature")
  }
  
  if (percentage > 1 | percentage < 0.00001){
    stop("Please provide a valid percentage between 0 and 1")
  }
  
  shiny::showNotification(paste0("Removing genes with more than ", percentage*100, "% zeroes in a row"))
  
  threshold = ncol(baseSignature)*percentage # max number of zeroes allowed
  signature = baseSignature[rowSums(baseSignature==0) <= threshold, ] 
  
  shiny::showNotification(paste0("Removed a total of ", nrow(baseSignature) - nrow(signature), " genes")) 
  
  return (signature)
} 


#' Remove unspecific Genes of a Gene Expression Signature
#' 
#' @param signature gene Expression Signature
#' @param numberOfBins number of bins to categorize the data into 
#' @param maxCount number of Cell Types allowed to be in the highest bin, all other cells are required to be in lower expressed bins
#' @param labels vector of bin names, required if numberOfBins != 3
#' 
#' @returns a gene expression signature containing only genes matching the passed requirements
removeUnspecificGenes = function (signature, numberOfBins = 3, maxCount = 2, labels = c("low", "medium", "high")){
  if (is.null(signature)){
    stop("Please provide a signature")
  }
  
  if (numberOfBins < 2){
    stop("numberOfBins has to be >= 2!")
  }
  
  if(maxCount<=0){
    stop("maxCount has to be a positive integer")
  }
  
  if (length(labels) != numberOfBins){
    stop("numberOfBins does not match label length")
  }
  
  shiny::showNotification("removing unspecific genes from signature")
  
  # initialize new refined signature with colnames, rows stay empty
  refinedSignature = matrix(nrow = 0, ncol = length(colnames(signature)), dimnames = list(NULL, colnames(signature)))
  
  # calculate bins for each gene and keep only the genes where <= maxCount genes are in the highest bin 
  for (i in 1:nrow(signature)){
    row <- signature[i, , drop=FALSE] # has colnames! drop FALSE is mandatory !!!!!
    
    # calculate bins to prevent error
    breaks <- seq(floor(min(row)), ceiling(max(row)), length.out=numberOfBins+1)
    
    # cut into bins, seperate for each gene
    bins <- cut(row, breaks = breaks, labels = labels, include.lowest = TRUE)
    
    nHighBins = sum(bins=="high") # not working when labels is something else
    
    # this value needs to be greater than one, depending of the step in the pipeline there arent 
    # any rows producing zeros left but that is not the case for all  signatures
    if (nHighBins<= maxCount & nHighBins>0){
      refinedSignature = rbind(refinedSignature, row)
    }
  }
  
  shiny::showNotification(paste("removed", nrow(signature) - nrow(refinedSignature), "unspecific genes from the signature."))
  
  # turn back to a matrix
  refinedSignature <- as.matrix(refinedSignature)
  
  return (refinedSignature)
}


#' Select a specified amount of genes for each cell type based on a score, discard all other 
#' 
#' @param signature gene expression matrix
#' @param method method to score the genes
#' @param selectCellType method to select the cell type the gene is contributing to, used to balance the number of genes between cell types
#' @param genesPerCellType maximum of genes selected for each cell type
selectGenesByScore <- function (signature, method = "entropy", selectCellType = "max", genesPerCellType = 20){
  # TODO Checks #####
  
  shiny::showNotification(paste0("Refining Signature by score: ", method))
  
  # SCORE THE MATRIX
  scoresByCellType = NULL
  for (celltype in colnames(signature)){
    scoresByCellType[[celltype]] <-  list()
  }
  
  
  # not using apply because i need the return values, did not work otherwise
  # scoring is happening here
  for (i in 1:nrow(signature)){
    row <- signature[i, , drop=FALSE] # has colnames!
    gene <- rownames(row)
    
    maxCelltype = colnames(signature)[max.col(row)] # this might be problematic
    
    score = list()
    if (method == "entropy"){
      score[gene] = scoreEntropy(row) # calculate score and save named result
    } else if (method == "gini"){
      score[gene] <- 1- BioQC::gini(row) # need to flip the value since lower scores schould be better (entropy!)
    }
    
    
    
    # append score to most expressed celltype
    scoresByCellType[[maxCelltype]] <- append(scoresByCellType[[maxCelltype]], score)
  }
  
  
  
  # initialize new refined signature with colnames, rows stay empty
  refinedSignature = matrix(nrow = 0, ncol = length(colnames(signature)), dimnames = list(NULL, colnames(signature)))
  
  # for each celltype get the lowest scores
  for (celltype in names(scoresByCellType)){
    # get best genes for each score, sorted ascending
    # reminder: the lower the score the better in entropy case, might need a change to fit all
    scores <- scoresByCellType[[celltype]] %>% unlist() %>% sort()
    
    # check number of genes!!!
    if (genesPerCellType>length(names(scores))){
      genesPerCellType <- length(names(scores))
    }
    
    # iterate genes 
    for (gene in names(scores)[1:genesPerCellType]){
      row <- signature[gene, , drop=FALSE] 
      
      refinedSignature <- rbind(refinedSignature, row)
    }
  }
  
  shiny::showNotification(paste0("removed a total of ", nrow(signature) - nrow(refinedSignature), " genes"))
  
  return (refinedSignature)
  
}

#' Score Gene Expression of a single Gene based on information entropy
#' 
#' @param geneExpression row from Gene Expression Matrix = Expression Data for a single Gene
#' @returns Score for the given gene based on information entropy
#' Here: The lower the better
scoreEntropy <- function (geneExpression){
  # TODO add parameter checks ####
  probs <- list()
  
  # turn expression data to a list of probabilities 
  for (val in geneExpression){
    if (val == 0){
      next
    }
    probs <- append(probs, val/sum(geneExpression)) # turn in to propabilities
  }
  
  entropy <- - sum (unlist(lapply(probs, function (x) log(x)*x)))
  
  return (entropy)
} 