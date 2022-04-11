#' Calculate Barplot of Signature Genes per Method
#'
#' @param signatures Named List of sigantures, names are the calculation methods
#' @param palette RColorBrewer palette name, standard = Set1
#'
#' @returns A Barplot
#'

plot_signatureGenesPerMethod <- function(signatures, palette="Set1") {
  df <- data.frame(method = character(), number_of_genes = numeric())

  # calculate number of genes per method
  for (name in names(signatures)) {
    number_of_genes <- dim(signatures[[name]])[1]
    df[nrow(df) + 1, ] <- list(name, number_of_genes)
  }

  # plot
  plot <- ggplot(data = df, aes(
    x = method, y = number_of_genes,
    fill = method, text = paste0(
      "Method: ",
      method, "\nGene Count: ",
      number_of_genes
    )
  )) +
    geom_col() +
    # ggtitle("Number of Signature Genes per Method") +
    labs(x = "Method", y = "Number of Genes", fill = "Method") +
    geom_text(aes(label = number_of_genes),
      fontface = "bold", size = 7,
      nudge_y = -1000, family = "Helvetica",
      color = "white"
    ) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    bbc_style() +
    theme(legend.position = "none") + 
    ggplot2::scale_fill_manual(values=RColorBrewer::brewer.pal(8, palette)[1:length(names(signatures))])

  plot
}

#' Calculate Barplot of Signature Genes per Method
#'
#' @param signatures Named List of sigantures, names are the calculation methods
#' @param palette RColorBrewer Palette name, standard = Set1
#'
#' @returns A Barplot

plot_conditionNumberPerMethod <- function(signatures, palette="Set1") {
  df <- data.frame(method = character(), kappa = numeric())

  # calculate condition number for each method
  for (name in names(signatures)) {
    kappa <- kappa(signatures[[name]][, -1], exact = TRUE)
    df[nrow(df) + 1, ] <- list(name, kappa)
  }

  # plot
  plot <- ggplot(data = df, aes(
    x = method, y = kappa, fill = method,
    text = paste0("Method: ", method, "\nKappa: ", kappa)
  )) +
    geom_col() +
    # ggtitle("5. Condition Number per Method") +
    geom_text(aes(label = round(kappa, 2)),
      fontface = "bold", nudge_y = -7,
      color = "white", size = 7, family = "Helvetica"
    ) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    bbc_style() +
    labs(x = "Method", y = "Kappa") +
    theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values=RColorBrewer::brewer.pal(8, palette)[1:length(names(signatures))])

  plot
}

#' Plot Mean Entropy for a set of singatures
#' 
#' @param signatures named List of Signatures
#' @param palette RColorBrewerPalette
#' 
#' @returns a plot
plot_meanEntropyPerMethod <- function(signatures, palette = "Set1"){
  
  entropies <- data.frame(method=character(), meanEntropy=numeric())
  
  # calculate Mean Entropy for each signature
  for (name in names(signatures)) {
    meanEntropy <- mean(apply(signatures[[name]], 1, scoreEntropy))
    entropies[nrow(entropies) + 1, ] <- list(name, meanEntropy)
  }
  
  plot <- ggplot(data = entropies, aes(
    x = method, y = meanEntropy, fill = method,
    text = paste0("Method: ", method, "\nEntropy: ", meanEntropy),
  )) +
    geom_col() +
    # ggtitle("5. Condition Number per Method") +
    geom_text(aes(label = round(meanEntropy, 2)),
              fontface = "bold", nudge_y = 0.1,
              color = "black", size = 7, family = "Helvetica"
    ) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    bbc_style() +
    labs(x = "Method", y = "Entropy") +
    theme(legend.position = "none") +
    # ggplot2::ylim(0, 5)+ # could be changed
    ggplot2::scale_fill_manual(values=RColorBrewer::brewer.pal(8, palette)[1:length(names(signatures))])
  
  plot
}


#' Calculate Clustered Heatmap of Signature Genes
#'
#' @param signature One Signature to plot
#' @param palette RColorBrewer Palette name, standard = Spectral
#' @param score The score used to annotate the genes
#' @param annotation_type How the score is rendered
#'
#' @returns A Heatmap
plot_signatureClustered <- function(signature, score="entropy", annotation_type="line", palette="Spectral") {
  if (is.null(signature)){
    stop("Please provide a signature")
  }
  
  if (!(score %in% c("entropy"))){
    stop("Score Method not supported")
  }
  
  if (!(annotation_type %in% c("line", "bar"))){
    stop("annotation_type not supported")
  }
  
  df <- data.frame(signature)

  df <- cbind("X" = rownames(df), df) # add gene names as column

  # log first
  df[, -1] <- log10(df[, -1] + 1)

  # calc mean and sd
  df$mean <- rowMeans(df[, -1])
  df$sd <- apply(df[, 2:(ncol(df) - 1)], 1, sd)

  # pivot for z-score calc
  df <- tidyr::pivot_longer(df, !c("X", "mean", "sd"), names_to = "cell_type", values_to = "value")
  df$z <- (df$value - df$mean) / df$sd # z-score

  # delete unused columns
  df$sd <- NULL
  df$value <- NULL
  df$mean <- NULL

  # pivot longer and convert to matrix
  df <- tidyr::pivot_wider(df, names_from = "cell_type", values_from = "z")
  mat <- as.matrix(df[, -1]) # without gene names
  rownames(mat) <- df$X # set gene names
  
  #mat <- stats::na.omit(mat) #####
  
  # calculate color palette
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c(RColorBrewer::brewer.pal(8, palette)[8:8], # first color of palette
                                                "white", # middle color
                                                RColorBrewer::brewer.pal(8, palette)[1:1] # last color of palette
                                                )
                                 )
  
  
  # render the signature annotation, this might also render multiple annotations
  # -> iterate over a list
  
  
  annotation <- NULL
  
  if (score == "entropy"){
    if (annotation_type == "line"){
      annotation <- ComplexHeatmap::columnAnnotation(entropy = ComplexHeatmap::anno_lines(apply(signature, 1, scoreEntropy), which = "row"))
      
    } else if (annotation_type == "bar"){
      annotation <- ComplexHeatmap::columnAnnotation(entropy = ComplexHeatmap::anno_barplot(apply(signature, 1, scoreEntropy), which = "row"))
    }
    
  }  else if (score == "gini"){
    annotation = NULL
  }

  
  # Plot with complex heatmap
  heatmap <- ComplexHeatmap::Heatmap(t(mat),
    name = "z-score", show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    row_title = NULL, row_names_side = "left",
    border = TRUE, col=col_fun, 
    #cluster_columns = agnes(mat), cluster_rows = diana(t(mat))
    cluster_columns = TRUE, cluster_rows = TRUE,  # clustering_method_columns = "euclidean",
    top_annotation = annotation
    
  )


  heatmap <- ComplexHeatmap::draw(heatmap)

  return(heatmap)

}


#' Calculate UpSet Plot Signature Genes
#'
#' @param signatures named list of deconvolution signatures
#' @param mode upSet Mode (distinct, intersect, union)
#' @param minDegree minimal set degree to display in the plot 
#' @param maxDegree maxmiaml set degree to display in the plot, NULL to display all sets
#' @param order order Sets by Size or Degree (size, degree)
#' @param invert invert the order of the Sets, standard = FALSE
#' @param colorDegrees color sets according to their degree, standard = TRUE
#' @param palette Name of a RColorBrewer palette, standard = Set1
#'
#' @returns UpSet Plot

plot_signatureUpset <- function(signatures, mode = "distinct", minDegree = 1,
                                maxDegree = NULL, order = "size", invert = FALSE,
                                colorDegrees = TRUE, palette = "Set1") {
  # takes list of signatures
  sets <- list()

  for (name in names(signatures)) {
    sets[[name]] <- rownames(signatures[[name]])
  }

  # modes available: distinct, intersect and union
  mat <- ComplexHeatmap::make_comb_mat(sets, mode = mode)

  # subset plot according to minDegree and maxDegree
  # if maxDegree is NULL, get maxDegree from data
  if (is.null(maxDegree)) {
    maxDegree <- max(ComplexHeatmap::comb_degree(mat))
  }

  mat <- mat[ComplexHeatmap::comb_degree(mat) >= minDegree] # lower
  mat <- mat[ComplexHeatmap::comb_degree(mat) <= maxDegree] # upper

  # calculate order: size, degree
  if (order == "size") {
    combOrder <- order(ComplexHeatmap::comb_size(mat), decreasing = !invert) # invert = FALSE -> will sort decreasing
  } else { # order=="degree"
    combOrder <- order(ComplexHeatmap::comb_degree(mat), decreasing = !invert)
  }

  # calculate colors
  if (colorDegrees == TRUE) {
    upSetColors <- RColorBrewer::brewer.pal(8, palette)[ComplexHeatmap::comb_degree(mat)] # max five different right now
  } else { # =FALSE
    upSetColors <- c("black")
  }

  plot <- ComplexHeatmap::UpSet(mat,
    comb_order = combOrder,
    top_annotation = ComplexHeatmap::upset_top_annotation(mat,
      add_numbers = TRUE,
      numbers_gp = grid::gpar(
        fontsize = "14",
        fontface = "bold"
      )
    ),
    pt_size = grid::unit(8, "mm"), lwd = 6,
    comb_col = upSetColors
  )

  return(list(plot, mat))
  # here is something missing, should evaluate the data here....
}

# helper, migth go elsewere in the code but belongs to upset plot function
download_signatureUpset <- function(signatures, combination, mode = "distinct") {
  # in case no set is selected return NULL
  if (is.null(combination)) {
    return(NULL)
  } else {
    # case: minimum of 1 Set selected
    sets <- list()
    token <- ""

    for (name in names(signatures)) {
      # add genes to set list
      sets[[name]] <- rownames(signatures[[name]])

      # check if name in combination and construct token
      if (name %in% combination) {
        token <- paste0(token, "1")
      } else {
        token <- paste0(token, "0")
      }
    }

    # modes available: distinct, intersect and union
    mat <- ComplexHeatmap::make_comb_mat(sets, mode = mode)

    # construct subset string from given signatures, if "000" return NULL
    return(ComplexHeatmap::extract_comb(mat, token))
  }
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
