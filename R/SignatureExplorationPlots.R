#' Calculate Barplot of Signature Genes per Method
#'
#' @param signatures Named List of sigantures, names are the calculation methods
#'
#' @returns A Barplot

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
  ggtitle("Number of Signature Genes per Method") +
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
  for (name in names(results)) {
    kappa <- kappa(results[[name]][, -1], exact = TRUE)
    df[nrow(df) + 1, ] <- list(name, kappa)
  }

  # plot
  plot <- ggplot(data = df, aes(
    x = method, y = kappa, fill = method,
    text = paste0("Method: ", method, "\nKappa: ", kappa)
  )) +
    geom_col() +
    ggtitle("5. Condition Number per Method") +
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