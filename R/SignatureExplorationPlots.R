#' Calculate Barplot of Signature Genes per Method
#'
#' This Barplot allows the comparison of the size of different signatures by plotting the
#' number of genes for each signature as a barplot
#'
#' @param signature_list Named List of signatures, names are the calculation methods
#' @param color_palette RColorBrewer palette name, standard = Set1
#'
#' @returns A Barplot
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#'
#' # list containing deconvolution results
#' signatureList <- list("bisque" = signature, "momf" = signature)
#'
#' plot_signatureGenesPerMethod(signatureList)
plot_signatureGenesPerMethod <- function(signature_list,
                                         color_palette = "Set1") {
  df <- data.frame(method = character(), number_of_genes = numeric())

  # calculate number of genes per method
  for (name in names(signature_list)) {
    number_of_genes <- dim(signature_list[[name]])[1]
    df[nrow(df) + 1, ] <- list(name, number_of_genes)
  }

  # plot
  p <- ggplot(data = df, aes(
    x = .data$method, y = .data$number_of_genes,
    fill = .data$method, text = paste0(
      "Method: ",
      .data$method, "\nGene Count: ",
      .data$number_of_genes
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
    geom_hline(yintercept = 0, linewidth = 1, colour = "#333333") +
    theme_minimal() +
    theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(8, color_palette)[1:length(names(signature_list))]) +
    ggplot2::ylim(0, max(df$number_of_genes) * 1.1) # scale y axis to contain bar label

  return(p)
}

#' Calculate Condition Number per Method
#'
#' This plot focuses on the condition number for each signature and simplifies
#' the comparison by providion a barplot
#'
#' @param signature_list Named List of signatures, names are the calculation methods
#' @param color_palette RColorBrewer Palette name, standard = Set1
#'
#' @returns A Barplot
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#'
#' # list containting deconvolution results
#' signatureList <- list("bisque" = signature, "momf" = signature)
#'
#' plot_conditionNumberPerMethod(signatureList)
plot_conditionNumberPerMethod <- function(signature_list,
                                          color_palette = "Set1") {
  df <- data.frame(method = character(), kappa = numeric())

  # calculate condition number for each method
  for (name in names(signature_list)) {
    kappa <- kappa(signature_list[[name]][, -1], exact = TRUE)
    df[nrow(df) + 1, ] <- list(name, kappa)
  }

  # plot
  p <- ggplot(data = df, aes(
    x = .data$method, y = .data$kappa, fill = .data$method,
    text = paste0("Method: ", .data$method, "\nKappa: ", .data$kappa)
  )) +
    geom_col() +
    # ggtitle("5. Condition Number per Method") +
    geom_text(aes(label = round(.data$kappa, 2)),
      fontface = "bold", vjust = -1,
      color = "black", size = 7, family = "Helvetica"
    ) +
    geom_hline(yintercept = 0, linewidth = 1, colour = "#333333") +
    theme_minimal() +
    labs(x = "Method", y = "Kappa") +
    theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(8, color_palette)[1:length(names(signature_list))]) +
    ggplot2::ylim(0, max(df$kappa) * 1.1) # scale y axis to contain bar label

  return(p)
}

#' Plot Mean Entropy for multiple signatures
#'
#' This plot focuses on the mean entropy of a signature and simplifies the comparison
#' to other signature with a barplot
#'
#' @param signature_list named List of Signatures
#' @param color_palette RColorBrewerPalette
#'
#' @returns a barplot
#'
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#'
#' # list containting deconvolution results
#' signatureList <- list("bisque" = signature, "momf" = signature)
#'
#' plot_meanEntropyPerMethod(signatureList)
plot_meanEntropyPerMethod <- function(signature_list,
                                      color_palette = "Set1") {
  entropies <- data.frame(method = character(), meanEntropy = numeric())

  # calculate Mean Entropy for each signature
  for (name in names(signature_list)) {
    meanEntropy <- mean(apply(signature_list[[name]], 1, scoreEntropy))
    entropies[nrow(entropies) + 1, ] <- list(name, meanEntropy)
  }

  p <- ggplot(data = entropies, aes(
    x = .data$method, y = .data$meanEntropy, fill = .data$method,
    text = paste0("Method: ", .data$method, "\nEntropy: ", .data$meanEntropy),
  )) +
    geom_col() +
    # ggtitle("5. Condition Number per Method") +
    geom_text(aes(label = round(.data$meanEntropy, 2)),
      fontface = "bold", vjust = -1,
      color = "black", size = 7, family = "Helvetica"
    ) +
    geom_hline(yintercept = 0, linewidth = 1, colour = "#333333") +
    theme_minimal() +
    labs(x = "Method", y = "Entropy") +
    theme(legend.position = "none") +
    # ggplot2::ylim(0, 5)+ # could be changed
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(8, color_palette)[1:length(names(signature_list))]) +
    ggplot2::ylim(0, max(entropies$meanEntropy) * 1.1) # scale y axis to contain bar label

  return(p)
}


#' Calculate Clustered Heatmap of Signature Genes
#'
#' This Heatmap displays a z-scored signature in heatmap form. The plot is annotated
#' by a gene scores ranking the distinctness of a gene in the signature.
#'
#' @param signature_mat One Signature to plot
#' @param color_palette RColorBrewer Palette name, standard = Spectral
#' @param scoring_method The score used to annotate the genes (entropy, gini)
#' @param annotation_type How the score is rendered (line, bar)
#'
#' @returns A Heatmap
#' @export
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#' plot_signatureClustered(signature, scoring_method = "gini", annotation_type = "bar")
plot_signatureClustered <- function(signature_mat,
                                    scoring_method = "entropy",
                                    annotation_type = "line",
                                    color_palette = "Spectral",
                                    order_rows = 'cluster') {
  if (is.null(signature_mat)) {
    stop("Please provide a signature")
  }

  if (!(scoring_method %in% c("entropy", "gini"))) {
    stop("Score Method not supported")
  }

  if (!(annotation_type %in% c("line", "bar"))) {
    stop("annotation_type not supported")
  }

  df <- data.frame(signature_mat)

  df <- cbind("X" = rownames(df), df) # add gene names as column

  # log first
  df[, -1] <- log10(df[, -1] + 1)

  # calc mean and sd
  df$mean <- rowMeans(df[, -1])
  df$sd <- apply(df[, 2:(ncol(df) - 1)], 1, sd) + 0.0001 # add pseudocount to not divide by 0 in case of SD=0

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

  # mat <- stats::na.omit(mat) #####

  # calculate color palette
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c(
    RColorBrewer::brewer.pal(8, color_palette)[8:8], # first color of palette
    "white", # middle color
    RColorBrewer::brewer.pal(8, color_palette)[1:1] # last color of palette
  ))


  # render the signature annotation, this might also render multiple annotations
  # -> iterate over a list

  annotation <- NULL

  if (scoring_method == "entropy") {
    if (annotation_type == "line") {
      annotation <- ComplexHeatmap::columnAnnotation(entropy = ComplexHeatmap::anno_lines(apply(signature_mat, 1, scoreEntropy), which = "row"))
    } else if (annotation_type == "bar") {
      annotation <- ComplexHeatmap::columnAnnotation(entropy = ComplexHeatmap::anno_barplot(apply(signature_mat, 1, scoreEntropy), which = "row"))
    }
  } else if (scoring_method == "gini") {
    if (annotation_type == "line") {
      annotation <- ComplexHeatmap::columnAnnotation(gini_index = ComplexHeatmap::anno_lines(apply(signature_mat, 1, BioQC::gini), which = "row"))
    } else if (annotation_type == "bar") {
      annotation <- ComplexHeatmap::columnAnnotation(gini_index = ComplexHeatmap::anno_barplot(apply(signature_mat, 1, BioQC::gini), which = "row"))
    }
  }

  if(order_rows == 'cluster'){
    cluster_rows <- TRUE
  }else if(order_rows == 'no_cluster'){
    mat <- mat[,order(colnames(mat))]
    cluster_rows = FALSE
  }

  # Plot with complex heatmap
  heatmap <- ComplexHeatmap::Heatmap(t(mat),
    name = "z-score", show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    row_title = NULL, row_names_side = "left",
    border = TRUE, col = col_fun,
    # cluster_columns = agnes(mat), cluster_rows = diana(t(mat))
    cluster_columns = TRUE, cluster_rows = cluster_rows, # clustering_method_columns = "euclidean",
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
#' @param signature_list named list of deconvolution signatures
#' @param upset_mode upSet Mode (distinct, intersect, union)
#' @param min_degree minimal set degree to display in the plot
#' @param max_degree maximal set degree to display in the plot, NULL to display all sets
#' @param order_sets order sets by Size or Degree (size, degree)
#' @param invert_sets Logical value. Inverts the order of the sets, defaults to FALSE
#' @param color_by_degrees Logical value. Whether to color sets according to their
#' degree, defaulting to TRUE
#' @param color_palette Name of a RColorBrewer palette, standard Set1
#'
#' @returns UpSet Plot
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#' signatures <- list("dwls" = signature, "momf" = signature, "bisque" = signature)
#' plot_signatureUpset(signatures, upset_mode = "union")
plot_signatureUpset <- function(signature_list,
                                upset_mode = "distinct",
                                min_degree = 1,
                                max_degree = NULL,
                                order_sets = "size",
                                invert_sets = FALSE,
                                color_by_degrees = TRUE,
                                color_palette = "Set1") {
  # takes list of signatures
  sets <- list()

  for (name in names(signature_list)) {
    sets[[name]] <- rownames(signature_list[[name]])
  }

  # modes available: distinct, intersect and union
  mat <- ComplexHeatmap::make_comb_mat(sets, mode = upset_mode)

  # subset plot according to min_degree and max_degree
  # if max_degree is NULL, get max_degree from data
  if (is.null(max_degree)) {
    max_degree <- max(ComplexHeatmap::comb_degree(mat))
  }

  mat <- mat[ComplexHeatmap::comb_degree(mat) >= min_degree] # lower
  mat <- mat[ComplexHeatmap::comb_degree(mat) <= max_degree] # upper

  # calculate order: size, degree
  if (order_sets == "size") {
    combOrder <- order(ComplexHeatmap::comb_size(mat), decreasing = !invert_sets) # invert_sets = FALSE -> will sort decreasing
  } else {
    # order_sets=="degree"
    combOrder <- order(ComplexHeatmap::comb_degree(mat), decreasing = !invert_sets)
  }

  # calculate colors
  if (color_by_degrees == TRUE) {
    upSetColors <- RColorBrewer::brewer.pal(8, color_palette)[ComplexHeatmap::comb_degree(mat)] # max five different right now
  } else {
    # =FALSE
    upSetColors <- c("black")
  }

  top_annotation <- ComplexHeatmap::upset_top_annotation(
    mat,
    add_numbers = TRUE,
    numbers_gp = grid::gpar(
     fontsize = "14",
     fontface = "bold"
    )
  )
  
  p <- ComplexHeatmap::UpSet(mat,
    comb_order = combOrder,
    top_annotation = top_annotation,
    pt_size = grid::unit(8, "mm"), lwd = 6,
    comb_col = upSetColors
  )

  return(list(p, mat))
  # here is something missing, should evaluate the data here....
}

#' Download a gene subset of multiple signature
#'
#' Returns gene sets of signatures according to the selected combination. As
#' intersection mode "distinct", "intersect" and "union" are available.
#'
#' @param signature_list list of named signatures
#' @param combination_to_include vector of signature names that should be intersected
#' @param upset_mode intersection type c("distinct", "intersect", "union")
#'
#' @return List of genes
#' @export
#'
#' @examples
#' signature <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
#'
#' signatures <- list("dwls" = signature, "momf" = signature, "bisque" = signature)
#' download_signatureUpset(signatures, c("dwls", "bisque"), "intersect")
download_signatureUpset <- function(signature_list,
                                    combination_to_include,
                                    upset_mode = "distinct") {
  # in case no set is selected return NULL
  if (is.null(combination_to_include)) {
    return(NULL)
  } else {
    # case: minimum of 1 Set selected
    sets <- list()
    token <- ""

    for (name in names(signature_list)) {
      # add genes to set list
      sets[[name]] <- rownames(signature_list[[name]])

      # check if name in combination and construct token
      if (name %in% combination_to_include) {
        token <- paste0(token, "1")
      } else {
        token <- paste0(token, "0")
      }
    }

    # modes available: distinct, intersect and union
    mat <- ComplexHeatmap::make_comb_mat(sets, mode = upset_mode)

    # construct subset string from given signatures, if "000" return NULL
    return(ComplexHeatmap::extract_comb(mat, token))
  }
}
