#' Calculate Barplot of Signature Genes per Method
#'
#' @param signatures Named List of sigantures, names are the calculation methods
#'
#' @returns A Barplot
#' 

plot_signatureGenesPerMethod <- function(signatures) {
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
  #ggtitle("Number of Signature Genes per Method") +
  labs(x = "Method", y = "Number of Genes", fill = "Method") +
  geom_text(aes(label = number_of_genes),
    fontface = "bold", size = 7,
    nudge_y = -1000, family = "Helvetica",
    color = "white"
  ) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  bbc_style() +
  theme(legend.position = "none")

  plot
}

#' Calculate Barplot of Signature Genes per Method 
#' 
#' @param signatures Named List of sigantures, names are the calculation methods
#'
#' @returns A Barplot

plot_conditionNumberPerMethod <- function(signatures) {
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
    #ggtitle("5. Condition Number per Method") +
    geom_text(aes(label = round(kappa, 2)),
      fontface = "bold", nudge_y = -7,
      color = "white", size = 7, family = "Helvetica"
    ) +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    bbc_style() +
    labs(x = "Method", y = "Kappa") +
    theme(legend.position = "none")

  plot
}

#' Calculate Clustered Heatmap of Signature Genes 
#' 
#' @param signature One Signature to plot 
#' 
#' @returns A Heatmap 
plot_signatureClustered <- function(signature){
  df = data.frame(signature)
  
  df <- cbind ("X" = rownames(df), df) # add gene names as column 
  
  # log first 
  df[,-1] <-  log10(df[,-1] + 1)
  
  # calc mean and sd
  df$mean <- rowMeans(df[,-1])
  df$sd <- apply(df[,2:(ncol(df)-1)], 1, sd)
  
  # pivot for z-score calc 
  df <-  tidyr::pivot_longer(df, !c("X", "mean", "sd"), names_to="cell_type", values_to="value")
  df$z <- (df$value - df$mean) / df$sd # z-score
  
  # delete unused columns 
  df$sd <- NULL
  df$value <- NULL
  df$mean <- NULL
  
  # pivot longer and convert to matrix
  df <- tidyr::pivot_wider(df, names_from = "cell_type", values_from="z")
  mat <- as.matrix(df[,-1]) # without gene names
  rownames(mat) <- df$X # set gene names 
  
  # Plot with complex heatmap 
  ComplexHeatmap::Heatmap(t(mat), name="z-score", show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE, 
                          # column_title = "10. Clustered Heatmap of Signature \nLog10 and z-scored, k-means partitioned", 
                          # column_title_gp = gpar(fontsize=20, fontface="bold"), 
                          row_title = NULL, row_split = ncol(mat), row_names_side = "left", 
                          cluster_columns = TRUE, column_km = ncol(mat), 
                          border=TRUE)
  
  #TODO Make Column Order deterministic!
  
}

plot_signatureUpset <- function(signatures, mode="distinct"){
  # takes list of signatures
  sets <- list()
  
  for (name in names(signatures)){
    sets[[name]] <- rownames(signatures[[name]])
  }
  
  # modes available: distinct, intersect and union 
  mat <- ComplexHeatmap::make_comb_mat(sets, mode = mode)
  
  # optional subset if intersect then remove the single ones (full set)
  #mat <- mat[comb_degree(mat)>=2]
  
  ComplexHeatmap::UpSet(mat, comb_order = order(
    ComplexHeatmap::comb_size(mat), decreasing = TRUE),
    top_annotation = ComplexHeatmap::upset_top_annotation(mat, add_numbers = TRUE))
  
  # here is something missing, should evaluate the data here.... 
}

# helper, migth go elsewere in the code but belongs to upset plot function 
download_signatureUpset <- function(signatures, combination){
  # takes list of signatures
  sets <- list()
  
  for (name in names(signatures)){
    sets[[name]] <- rownames(signatures[[name]])
  }
  
  # modes available: distinct, intersect and union 
  mat <- ComplexHeatmap::make_comb_mat(sets, mode = "distinct")
  
  # hier: return intersection genes
}


