
#' plot benchmark
#' @param gtruth dataframe of gtruth/simulation fractions
#' @param estimate deconvolution result list (named)
#' @param palette RColorBrewer Palette
plot_benchmark_scatter <- function(gtruth, estimate, palette = "Spectral") {
  
  ref <- gtruth %>% as.data.frame()
  ref$sample <- rownames(ref)
  ref <- tidyr::pivot_longer(ref, !sample, names_to = "cell_type", values_to = "truth")
  
  df <- NULL 
  # prepare data and add method and rownames as column
  for (method in names(estimate)) {
    deconvolution <- as.data.frame(estimate[[method]])
    deconvolution$sample <- rownames(deconvolution)
    deconvolution$method <- rep(method, nrow(deconvolution))
    
    df <- rbind(df, tidyr::pivot_longer(deconvolution, !c("sample", "method"), names_to = "cell_type", values_to = "estimate"))
  }
  
  # merge
  merged.df <- merge(df, ref, all.x = TRUE) # keep all deconvolution results
  
  # build plot
  plot <- ggplot(merged.df, aes(x = truth, y = estimate, color = cell_type)) +
    geom_point(show.legend = FALSE) +
    ggpubr::stat_cor(label.sep = "\n", size = 3, color = "black", label.x.npc = 0.01) +
    ggforce::facet_grid_paginate(method ~ cell_type,
                        margins = c("cell_type"), scales = "free"
    ) +
    ggplot2::theme_bw() +
    ggpubr::rotate_x_text(angle = 60) +
    labs(x = "true cellular fractions", y = "cell type estimates", title = "") +
    theme(legend.position = "none", text = element_text(size = 15))
  
  plot <- plot + geom_abline() +
    ggforce::facet_grid_paginate(method ~ cell_type,
                        margins = c("cell_type"), scales = "free"
    ) +
    ggpubr::rotate_x_text(angle = 60)
  
  
  # get palette
  max_colors <- RColorBrewer::brewer.pal.info[palette, ]$maxcolors # for brewer.pal()
  n_cell_types <- length(unique(merged.df$cell_type)) # number of needed colors
  getPalette <- colorRampPalette(brewer.pal(max_colors, palette)) # function to return custom interpolated palettes
  
  # color
  plot <- plot + ggplot2::scale_color_manual(values = getPalette(n_cell_types + 1))
  
  plot
}


#' plot benchmark correlation
#' @param gtruth dataframe of gtruth/simulation fractions
#' @param estimate deconvolution result list (named)
#' @param plot_method method to plot, one of c("circle", "square", "ellipse", "number", "shade", "color", "pie")
#' @param pValueColor color of p value annotation, "white" or "black"
#' @param pValueType one of the following c("p-value", "label_sig", "n"), see corrplot package for further info
plot_benchmark_correlation <- function(gtruth, estimate, pValueType = "label_sig", pValueColor="black", plot_method = "number") {
  if (!plot_method %in% c("circle", "square", "ellipse", "number", "shade", "color", "pie")) {
    stop("correlation plot method not supported")
  }
  
  if (!pValueType %in% c("p-value", "label_sig", "n")){
    stop("P Value Method not supported")
  }
  
  if (!pValueColor %in% c("black", "white")){
    stop("P Value Annotation color must be white or black")
  }
  
  ref <- gtruth %>% as.data.frame()
  ref$sample <- rownames(ref)
  
  ref <- tidyr::pivot_longer(ref, !sample, names_to = "cell_type", values_to = "truth")
  
  df <- NULL
  
  # prepare data and add method and rownames as column
  for (method in names(estimate)) {
    deconvolution <- as.data.frame(estimate[[method]])
    deconvolution$sample <- rownames(deconvolution)
    deconvolution$method <- rep(method, nrow(deconvolution))
    
    df <- rbind(df, tidyr::pivot_longer(deconvolution, !c("sample", "method"), names_to = "cell_type", values_to = "estimate"))
  }
  
  # merge
  merged.df <- merge(df, ref, all.x = TRUE) # keep all deconvolution results
  
  cor.df <- data.frame("method" = character(), "cell_type" = character(), "correlation" = numeric())
  p.df <- data.frame("method" = character(), "cell_type" = character(), "p" = numeric())
  
  # for each cell type and each method calculate correlation
  for (method in names(estimate)) {
    subset <- merged.df[merged.df$method == method, ]
    for (cellType in unique(subset$cell_type)) {
      subsubset <- subset[subset$cell_type == cellType, ]
      cor <- cor(subsubset$truth, subsubset$estimate)
      
      p <- tryCatch(
        {
          cor.test(subsubset$truth, subsubset$estimate)$p.value
        },
        error = function(e) {
          NA
        }
      )
      
      row <- c(method, cellType, cor)
      cor.df[nrow(cor.df) + 1, ] <- row
      
      row <- c(method, cellType, p)
      p.df[nrow(p.df) + 1, ] <- row
    }
  }
  
  # round
  cor.df$correlation <- round(as.numeric(cor.df$correlation), 2)
  p.df$p <- round(as.numeric(p.df$p), 2)
  
  
  # pivot longer and turn to matrix for corrplot
  
  cor.df <- tidyr::pivot_wider(cor.df, names_from = "cell_type", values_from = "correlation") %>% as.data.frame()
  rownames(cor.df) <- cor.df$method
  cor.df$method <- NULL
  
  
  p.df <- tidyr::pivot_wider(p.df, names_from = "cell_type", values_from = "p") %>% as.data.frame()
  rownames(p.df) <- p.df$method
  p.df$method <- NULL
  
  cor.df <- as.matrix(cor.df)
  p.df <- as.matrix(p.df)
  
  # return plot
  return(corrplot::corrplot(cor.df,
                            p.mat = p.df, insig = pValueType, sig.level = c(0.05, 0.1, 0.2), pch.cex = 4, 
                            pch.col = pValueColor, method = plot_method,
                            na.label = "NA", tl.col = "black", tl.srt = 60,
                            cl.pos = "r", cl.align.text = "r", tl.cex = 2, cl.cex = 1.5,
                            number.cex = 1.5, na.label.col = "#7F7F7F"
  ))
}
