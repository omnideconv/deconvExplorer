# library(ggplot2)

plot_deconvolution <- function(to_plot_list, plotMethod, all_deconvolutions) {

  # load all deconvolutions from all_deconvolutions into a list
  deconvolution_list <- list()
  for (deconvolution in to_plot_list) {
    deconvolution_list[length(deconvolution_list) + 1] <- all_deconvolutions[[deconvolution]][1]
  }

  # preformat data
  deconvolution_list <- lapply(deconvolution_list, function(deconvolution) {
    cbind(deconvolution, samples = rownames(deconvolution)) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(!samples, names_to = "cell_type", values_to = "fraction")
  })

  # add all different plots together and add facet information
  data <- do.call("rbind", deconvolution_list)
  data$facets <- rep(to_plot_list, each = nrow(data) / length(to_plot_list))

  plot <- ggplot(data, aes(x = reorder(cell_type, desc(cell_type)), y = as.numeric(fraction), fill = cell_type, text = ""))

  if (plotMethod == "bar") {
    plot <- plot + geom_col(aes(
      y = samples, x = as.numeric(fraction), fill = cell_type,
      text = paste0(
        "Cell Type: ", cell_type, "\nFraction: ",
        sprintf("%1.2f%%", 100 * as.numeric(fraction))
      )
    )) +
      # facet_wrap(~data$facets) +
      labs(x = "estimated fraction", y = "sample", fill = "cell type")
  } else if (plotMethod == "jitter") {
    plot <- plot + geom_jitter(aes(color = cell_type, text = paste0(
      "Sample: ", samples, "\nFraction: ",
      sprintf("%1.2f%%", 100 * as.numeric(fraction))
    ))) +
      labs(x = "cell type", y = "estimated fraction", color = "cell type", fill = "")
  } else if (plotMethod == "scatter") {
    plot <- plot + geom_point(aes(color = cell_type, text = paste0(
      "Sample: ", samples, "\nFraction: ",
      sprintf("%1.2f%%", 100 * as.numeric(fraction))
    ))) +
      labs(x = "cell type", y = "estimated fraction", color = "cell type", fill = "")
  } else if (plotMethod == "box") {
    plot <- plot + geom_boxplot(aes(text = "")) +
      labs(x = "cell type", y = "estimated fraction", fill = "cell type")
  } else if (plotMethod == "sina") {
    plot <- plot + geom_violin(colour = "grey", fill = "grey") +
      ggforce::geom_sina() +
      labs(x = "cell type", y = "estimated fraction", alpha = "", fill = "cell type")
  } else if (plotMethod == "heatmap") {
    plot <- plot + geom_tile(aes(y = samples, fill = as.numeric(fraction), text = sprintf("%1.2f%%", 100 * as.numeric(fraction)))) +
      # geom_text(aes(y = samples, x = cell_type, label =  sprintf("%1.2f%%", 100*as.numeric(fraction))))+
      labs(x = "cell type", y = "sample", fill = "estimated fraction") +
      scale_fill_gradient(low = "white", high = "blue") +
      guides(fill = guide_colorbar(barwith = 0.5, barheight = 20))
  }

  plot <- plot + facet_wrap(~ data$facets)

  # render
  plotly::ggplotly(plot, tooltip = c("text")) %>%
    plotly::config(
      displaylogo = FALSE, showTips = FALSE, toImageButtonOptions = list(filename = paste0(plotMethod, "_plot")),
      modeBarButtonsToRemove = list(
        "hoverCLosestCartesian",
        "hoverCompareCartesian",
        "zoomIn2d", "zoomOut2d",
        "lasso2d", "zoom2d",
        "pan2d", "autoScale2d", "select2d"
      )
    ) %>%
    plotly::layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE))
}
