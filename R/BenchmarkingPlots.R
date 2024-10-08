#' plot benchmarking scatterplot
#'
#' Plot a deconvolution results against it's corresponding ground truth.
#'
#' @param gtruth_df dataframe of gtruth/simulation fractions
#' @param deconv_list deconvolution result list (named)
#' @param color_palette RColorBrewer Palette
#'
#' @return A `ggplot` object
#'
#' @export
#'
#' @examples
#' data("RefData", package = "omnideconv")
#' RefData <- as.data.frame(RefData)
#' deconv <- readRDS(system.file("extdata", "deconvolution_example.rds",
#'   package = "DeconvExplorer"
#' ))
#' deconvList <- list("momf" = deconv, "bisque" = deconv)
#' plot_benchmark_scatter(RefData, deconvList)
plot_benchmark_scatter <- function(gtruth_df,
                                   deconv_list,
                                   color_palette = "Spectral") {
  stopifnot(is.data.frame(gtruth_df))

  # is a list AND not a data frame
  stopifnot(is.list(deconv_list))
  if (is.data.frame(deconv_list)) {
    stop("You should provide the set of estimates as a list, not as a single data frame")
  }
  # and contains data frames as expected
  stopifnot(all(unlist(lapply(deconv_list, is.data.frame))))

  ref <- as.data.frame(gtruth_df)
  ref$sample <- rownames(ref)
  ref <- tidyr::pivot_longer(ref, !sample, names_to = "cell_type", values_to = "truth")

  df <- NULL
  # prepare data and add method and rownames as column
  for (method in names(deconv_list)) {
    deconvolution <- as.data.frame(deconv_list[[method]])
    deconvolution$sample <- rownames(deconvolution)
    deconvolution$method <- rep(method, nrow(deconvolution))

    df <- rbind(df, tidyr::pivot_longer(deconvolution, !c("sample", "method"), names_to = "cell_type", values_to = "estimate"))
  }

  # merge
  merged.df <- merge(df, ref, all.x = TRUE) # keep all deconvolution results
  merged.df <- merged.df[complete.cases(merged.df), ]

  # build plot
  plot <- ggplot(merged.df, aes(x = .data$truth, y = .data$estimate, color = .data$cell_type)) +
    geom_point(show.legend = FALSE, size = 3, alpha = .8) +
    ggpubr::stat_cor(label.sep = "\n", size = 3, color = "black", label.x.npc = 0.05, label.y = max(merged.df$estimate) + 0.05, vjust = 1) +
    ggforce::facet_grid_paginate(method ~ .data$cell_type,
      margins = c("cell_type")
    ) +
    ggplot2::theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), strip.background = ggplot2::element_rect(fill = "white")) +
    labs(x = "true cellular fractions", y = "cell type estimates", title = "") +
    theme(legend.position = "none", text = element_text(size = 15)) +
    ggplot::geom_abline(linetype = "dashed")

  # get palette
  max_colors <- RColorBrewer::brewer.pal.info[color_palette, ]$maxcolors # for brewer.pal()
  n_cell_types <- length(unique(merged.df$cell_type)) # number of needed colors
  getPalette <- colorRampPalette(brewer.pal(max_colors, color_palette)) # function to return custom interpolated palettes

  # color
  plot <- plot + ggplot2::scale_color_manual(values = getPalette(n_cell_types + 1))

  return(plot)
}


#' Plot benchmark correlation
#'
#' Plot the correlation of a deconvolution results and it's corresponding ground truth
#'
#' @param gtruth_df dataframe of gtruth/simulation fractions
#' @param deconv_list deconvolution result list (named)
#' @param plot_method method to plot, one of c("circle", "square", "ellipse", "number", "shade", "color", "pie")
#' @param pvalue_color color of p value annotation, "white" or "black"
#' @param pvalue_type one of the following c("p-value", "label_sig", "n"), see corrplot package for further info
#'
#' @return A list, as returned by the `corrplot` function
#'
#' @export
#'
#' @examples
#' data("RefData", package = "omnideconv")
#' RefData <- as.data.frame(RefData)
#' deconv <- readRDS(system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer"))
#' deconvList <- list("momf" = deconv, "bisque" = deconv)
#' plot_benchmark_correlation(RefData, deconvList)
plot_benchmark_correlation <- function(gtruth_df,
                                       deconv_list,
                                       pvalue_type = "label_sig",
                                       pvalue_color = "black",
                                       plot_method = "number") {
  stopifnot(is.data.frame(gtruth_df))

  # is a list AND not a data frame
  stopifnot(is.list(deconv_list))
  if (is.data.frame(deconv_list)) {
    stop("You should provide the set of estimates as a list, not as a single data frame")
  }
  # and contains data frames as expected
  stopifnot(all(unlist(lapply(deconv_list, is.data.frame))))

  if (!plot_method %in% c("circle", "square", "ellipse", "number", "shade", "color", "pie")) {
    stop("correlation plot method not supported")
  }

  if (!pvalue_type %in% c("p-value", "label_sig", "n")) {
    stop("P Value Method not supported")
  }

  if (!pvalue_color %in% c("black", "white")) {
    stop("P Value Annotation color must be white or black")
  }

  ref <- as.data.frame(gtruth_df)
  ref$sample <- rownames(ref)

  ref <- tidyr::pivot_longer(ref, !sample, names_to = "cell_type", values_to = "truth")

  df <- NULL

  # prepare data and add method and rownames as column
  for (method in names(deconv_list)) {
    deconvolution <- as.data.frame(deconv_list[[method]])
    deconvolution$sample <- rownames(deconvolution)
    deconvolution$method <- rep(method, nrow(deconvolution))

    df <- rbind(df, tidyr::pivot_longer(deconvolution, !c("sample", "method"), names_to = "cell_type", values_to = "estimate"))
  }

  # merge
  merged.df <- merge(df, ref, all.x = TRUE) # keep all deconvolution results

  cor.df <- data.frame("method" = character(), "cell_type" = character(), "correlation" = numeric())
  p.df <- data.frame("method" = character(), "cell_type" = character(), "p" = numeric())

  # for each cell type and each method calculate correlation
  for (method in names(deconv_list)) {
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

  cor.df <- tidyr::pivot_wider(cor.df, names_from = "cell_type", values_from = "correlation") |>
    as.data.frame()
  rownames(cor.df) <- cor.df$method
  cor.df$method <- NULL


  p.df <- tidyr::pivot_wider(p.df, names_from = "cell_type", values_from = "p") |>
    as.data.frame()
  rownames(p.df) <- p.df$method
  p.df$method <- NULL

  cor.df <- as.matrix(cor.df)
  p.df <- as.matrix(p.df)

  # return plot
  return(corrplot::corrplot(cor.df,
    p.mat = p.df, insig = pvalue_type, sig.level = c(0.05, 0.1, 0.2), pch.cex = 4,
    pch.col = pvalue_color, method = plot_method,
    na.label = "NA", tl.col = "black", tl.srt = 60,
    cl.pos = "r", cl.align.text = "r", tl.cex = 2, cl.cex = 1.5,
    number.cex = 1.5, na.label.col = "#7F7F7F"
  ))
}

#' plot benchmark RMSE
#'
#' Plot RMSE (root mean squared error) for a list of deconvolution results and the
#' corresponding ground truth
#'
#' @param gtruth_df dataframe of gtruth/simulation fractions
#' @param deconv_list deconvolution result list (named)
#' @param hm_method method to plot, one of c("circle", "square", "ellipse", "number", "shade", "color", "pie")
#' @param plot_type "heatmap" or "boxplot"
#' @param color_palette RColorBrewer Palette
#'
#' @return A `ggplot` object, or a list as returned by `corrplot()`
#'
#' @export
#'
#' @examples
#' data("RefData", package = "omnideconv")
#' RefData <- as.data.frame(RefData)
#' deconv <- readRDS(system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer"))
#' deconvList <- list("momf" = deconv, "bisque" = deconv)
#' plot_benchmark_rmse(RefData, deconvList)
plot_benchmark_rmse <- function(gtruth_df,
                                deconv_list,
                                plot_type = "heatmap",
                                hm_method = "color",
                                color_palette = "Spectral") {
  stopifnot(is.data.frame(gtruth_df))

  # is a list AND not a data frame
  stopifnot(is.list(deconv_list))
  if (is.data.frame(deconv_list)) {
    stop("You should provide the set of estimates as a list, not as a single data frame")
  }
  # and contains data frames as expected
  stopifnot(all(unlist(lapply(deconv_list, is.data.frame))))

  if (!(plot_type %in% c("heatmap", "boxplot"))) {
    stop("plot_type not supported")
  }

  if (!(hm_method %in% c("circle", "square", "ellipse", "number", "shade", "color", "pie"))) {
    stop("hm_method not supported")
  }

  if (is.null(gtruth_df) | is.null(deconv_list)) {
    stop("gtruth_df or deconv_list not provided")
  }

  ref <- as.data.frame(gtruth_df)
  ref$sample <- rownames(ref)

  ref <- tidyr::pivot_longer(ref, !sample, names_to = "cell_type", values_to = "truth")

  df <- NULL

  # prepare data and add method and rownames as column
  for (method in names(deconv_list)) {
    deconvolution <- as.data.frame(deconv_list[[method]])
    deconvolution$sample <- rownames(deconvolution)
    deconvolution$method <- rep(method, nrow(deconvolution))

    df <- rbind(df, tidyr::pivot_longer(deconvolution, !c("sample", "method"), names_to = "cell_type", values_to = "estimate"))
  }

  # merge
  merged.df <- merge(df, ref, all.x = TRUE) # keep all deconvolution results

  rmse.df <- data.frame("method" = character(), "cell_type" = character(), "rmse" = numeric())

  # for each cell type and each method calculate correlation
  for (method in names(deconv_list)) {
    subset <- merged.df[merged.df$method == method, ]
    for (cellType in unique(subset$cell_type)) {
      subsubset <- subset[subset$cell_type == cellType, ]
      rmse <- sqrt(mean((subsubset$truth - subsubset$estimate)^2))

      row <- c(method, cellType, rmse)
      rmse.df[nrow(rmse.df) + 1, ] <- row
    }
  }

  # round
  rmse.df$rmse <- round(as.numeric(rmse.df$rmse), 2)

  # return plot
  if (plot_type == "heatmap") {
    # pivot wider for corrplot
    rmse.df <- tidyr::pivot_wider(rmse.df, names_from = "cell_type", values_from = "rmse") |>
      as.data.frame()
    rownames(rmse.df) <- rmse.df$method
    rmse.df$method <- NULL
    rmse.df <- as.matrix(rmse.df)

    # plot
    return(corrplot::corrplot(rmse.df,
      method = hm_method,
      na.label = "NA", tl.col = "black", tl.srt = 60,
      cl.pos = "r", cl.align.text = "r", tl.cex = 2, cl.cex = 1.5,
      number.cex = 1.5, na.label.col = "#7F7F7F"
    ))
  } else if (plot_type == "boxplot") {
    # get color first
    max_colors <- RColorBrewer::brewer.pal.info[color_palette, ]$maxcolors # for brewer.pal()
    n_methods <- length(unique(merged.df$method)) # number of needed colors
    getPalette <- colorRampPalette(brewer.pal(max_colors, color_palette)) # function to return custom interpolated palettes


    plot <- ggplot(rmse.df, aes(x = method, y = rmse, fill = method)) +
      geom_boxplot() +
      ggplot2::theme_bw() +
      labs(x = "Method", y = "RMSE", fill = "Method") +
      ggplot2::scale_color_manual(values = getPalette(n_methods))

    # some theming
    plot <- plot + ggplot2::theme(
      # panel.grid.major = ggplot2::element_blank()
      panel.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#bdbebd"),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16)
    )

    return(plot)
  }
}
