#' Calculate Barplot of Signature Genes per Method
#'
#' This Barplot allows the comparison of the size of differnt signatures by plotting the
#' number of genes for each signature as a barplot
#'
#' @param signatures Named List of signatures, names are the calculation methods
#' @param palette RColorBrewer palette name, standard = Set1
#'
#' @returns A Barplot
#' 
#' @export
#'
#' @examples
#' library(DeconvExplorer)
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package="DeconvExplorer"))
#' 
#' # list containing deconvolution results
#' signatureList = list("bisque" =signature, "momf" = signature)
#' 
#' plot_signatureGenesPerMethod(signatureList)
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
    ggplot2::geom_text(aes(label = number_of_genes),
      fontface = "bold", size = 7,
      vjust = -1, family = "Helvetica",
      color = "black"
    ) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    bbc_style() +
    theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values=RColorBrewer::brewer.pal(8, palette)[1:length(names(signatures))])+
    ggplot2::ylim (0, max(df$number_of_genes)*1.1) # scale y axis to contain bar label


  plot
}

#' Calculate Condition Number per Method
#'
#' This plot focuses on the condition number for each signature and simplifies
#' the comparison by providion a barplot
#'
#' @param signatures Named List of signatures, names are the calculation methods
#' @param palette RColorBrewer Palette name, standard = Set1
#'
#' @returns A Barplot
#' 
#' @export
#' 
#' @examples
#' library(DeconvExplorer)
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package="DeconvExplorer"))
#' 
#' # list containting deconvolution results
#' signatureList = list("bisque" =signature, "momf" = signature)
#' 
#' plot_conditionNumberPerMethod(signatureList)
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
      fontface = "bold", vjust = -1,
      color = "black", size = 7, family = "Helvetica"
    ) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    bbc_style() +
    labs(x = "Method", y = "Kappa") +
    theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values=RColorBrewer::brewer.pal(8, palette)[1:length(names(signatures))]) +
    ggplot2::ylim(0, max(df$kappa)*1.1) # scale y axis to contain bar label

  plot
}

#' Plot Mean Entropy for multiple signatures
#'
#' This plot focuses on the mean entropy of a signature and simplifies the comparison
#' to other signature with a barplot
#'
#' @param signatures named List of Signatures
#' @param palette RColorBrewerPalette
#'
#' @returns a barplot
#' 
#' @export
#'
#' @examples
#' library(DeconvExplorer)
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package="DeconvExplorer"))
#' 
#' # list containting deconvolution results
#' signatureList = list("bisque" =signature, "momf" = signature)
#' 
#' plot_meanEntropyPerMethod(signatureList)
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
              fontface = "bold", vjust = -1,
              color = "black", size = 7, family = "Helvetica"
    ) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    bbc_style() +
    labs(x = "Method", y = "Entropy") +
    theme(legend.position = "none") +
    # ggplot2::ylim(0, 5)+ # could be changed
    ggplot2::scale_fill_manual(values=RColorBrewer::brewer.pal(8, palette)[1:length(names(signatures))])+
    ggplot2::ylim(0, max(entropies$meanEntropy)*1.1) # scale y axis to contain bar label

  plot
}


#' Calculate Clustered Heatmap of Signature Genes
#'
#' This Heatmap displayes a z-scored signature in heatmap form. The plot is annotated
#' by a gene scores ranking the distinctness of a gene in the signature.
#'
#' @param signature One Signature to plot
#' @param palette RColorBrewer Palette name, standard = Spectral
#' @param score The score used to annotate the genes (entropy, gini)
#' @param annotation_type How the score is rendered (line, bar)
#'
#' @returns A Heatmap
#' @export
#' @examples
#' library(DeconvExplorer)
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package="DeconvExplorer"))
#' plot_signatureClustered(signature, score="gini", annotation_type="bar")
plot_signatureClustered <- function(signature, score="entropy", annotation_type="line", palette="Spectral") {
  if (is.null(signature)){
    stop("Please provide a signature")
  }

  if (!(score %in% c("entropy", "gini"))){
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
    if (annotation_type == "line"){
      annotation <- ComplexHeatmap::columnAnnotation(gini_index = ComplexHeatmap::anno_lines(apply(signature, 1, BioQC::gini), which = "row"))
    } else if(annotation_type == "bar"){
      annotation <- ComplexHeatmap::columnAnnotation(gini_index = ComplexHeatmap::anno_barplot(apply(signature, 1, BioQC::gini), which = "row"))
    }
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


#' Calculate UpSet Plot for Signature Gene Sets
#'
#' This plot allows the comparison of multiple signature gene sets by utilizing
#' UpSet Plots.
#'
#' @param signatures named list of deconvolution signatures
#' @param mode upSet Mode (distinct, intersect, union)
#' @param minDegree minimal set degree to display in the plot
#' @param maxDegree maximal set degree to display in the plot, NULL to display all sets
#' @param order order Sets by Size or Degree (size, degree)
#' @param invert invert the order of the Sets, standard = FALSE
#' @param colorDegrees color sets according to their degree, standard = TRUE
#' @param palette Name of a RColorBrewer palette, standard = Set1
#'
#' @returns UpSet Plot
#' @export
#'
#' @examples
#' library(DeconvExplorer)
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package="DeconvExplorer"))
#' signatures <- list("dwls" = signature, "momf" = signature, "bisque" = signature)
#' plot_signatureUpset(signatures, mode="union")
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
    top_annotation = upset_top_annotation(mat,
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

#' Download a gene subset of multiple signature
#'
#' Returns gene sets of signatures according to the selected combination. As
#' intersection mode "distinct", "intersect" and "union" are available. 
#'
#' @param signatures list of named signatures
#' @param combination vector of signaature names that should be intersected
#' @param mode intersection type c("distinct", "intersect", "union")
#'
#' @return List of genes
#' @export
#'
#' @examples
#' library(DeconvExplorer)
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package="DeconvExplorer"))
#' 
#' signatures = list("dwls" = signature, "momf" = signature, "bisque" = signature)
#' download_signatureUpset(signatures, c("dwls", "bisque"), "intersect")
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
