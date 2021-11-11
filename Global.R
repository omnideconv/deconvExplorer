
# axis definition ---------------------------------------------------------

axis <- list(
  "method" = list( # x, y, fill/color
    "bar" = list("fraction", "sample", "cell_type"),
    "jitter" = list("fraction", "cell_type", "cell_type"),
    "scatter" = list("fraction", "cell_type", "cell_type"),
    "box" = list("cell_type", "fraction", "cell_type"),
    "sina" = list("cell_type", "fraction", "cell_type"),
    "heatmap" = list("cell_type", "sample", "fraction")
  ),
  "cell_type" = list(
    "bar" = list("fraction", "sample", "method"),
    "jitter" = list("fraction", "method", "sample"),
    "scatter" = list("fraction", "method", "sample"),
    "box" = list("method", "fraction", "method"),
    "sina" = list("method", "fraction", "cell_type"),
    "heatmap" = list("sample", "method", "fraction")
  ),
  "sample" = list(
    "bar" = list("cell_type", "fraction", "method"),
    "jitter" = list("cell_type", "fraction", "method"),
    "scatter" = list("cell_type", "fraction", "method"),
    "box" = list("method", "fraction", "method"),
    "sina" = list("method", "fraction", "method"),
    "heatmap" = list("cell_type", "method", "fraction")
  )
)

# aes definition  ---------------------------------------------------------
getAes <- function(facets, plotMethod) {
  x <- axis[[facets]][[plotMethod]][[1]]
  y <- axis[[facets]][[plotMethod]][[2]]
  col <- axis[[facets]][[plotMethod]][[3]]

  if (plotMethod %in% c("jitter", "scatter")) {
    return(ggplot2::aes_string(x = x, y = y, col = col))
  } else if (plotMethod == "sina") { # does not work with aes_string, using aes_
    return(ggplot2::aes_(x = as.name(x), y = as.name(y), col = as.name(col)))
  } else {
    return(ggplot2::aes_string(x = x, y = y, fill = col))
  }
}

getLabs <- function(facets, plotMethod) {
  x <- axis[[facets]][[plotMethod]][[1]]
  y <- axis[[facets]][[plotMethod]][[2]]
  col <- axis[[facets]][[plotMethod]][[3]]
  if (plotMethod %in% c("jitter", "scatter")) {
    return(ggplot2::labs(x = x, y = y, col = col))
  } else {
    return(ggplot2::labs(x = x, y = y, fill = col))
  }
}

# functions ---------------------------------------------------------------

plot_deconvolution <- function(to_plot_list, plotMethod, facets, all_deconvolutions) {
  # load deconvolutions from all_deconvolutions into a list
  deconvolution_list <- list()
  for (deconvolution in to_plot_list) {
    deconvolution_list[length(deconvolution_list) + 1] <- all_deconvolutions[[deconvolution]][1]
  }

  # preformat data into a dataframe
  deconvolution_list <- lapply(deconvolution_list, function(deconvolution) {
    cbind(deconvolution, sample = rownames(deconvolution)) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(!sample, names_to = "cell_type", values_to = "fraction")
  })

  # add all different plots together and add facet information
  data <- do.call("rbind", deconvolution_list)
  data$fraction <- as.numeric(data$fraction) # change datatype of column
  data$method <- rep(to_plot_list, each = nrow(data) / length(to_plot_list)) # add computation method as column

  # calculate tooltip based on chosen "group by"
  tooltip <- switch(facets,
    "method" = aes(
      text = paste0("Cell Type: ", cell_type, "\nFraction: ", sprintf("%1.2f%%", 100 * fraction), "\nSample: ", sample)
    ),
    "cell_type" = aes(
      text = paste0("Sample: ", sample, "\nFraction: ", sprintf("%1.2f%%", 100 * fraction), "\nMethod: ", method)
    ),
    "sample" = aes(
      text = paste0("Cell Type: ", cell_type, "\nFraction: ", sprintf("%1.2f%%", 100 * fraction), "\nMethod: ", method)
    )
  )

  # put plot together
  plot <- ggplot(data, getAes(facets, plotMethod))
  plot <- plot + facet_wrap(~ data[[facets]])

  if (plotMethod == "bar") {
    plot <- plot + geom_col(tooltip)+
      getLabs(facets, plotMethod)
  } else if (plotMethod == "jitter") {
    plot <- plot + geom_jitter(tooltip) +
      getLabs(facets, plotMethod)
  } else if (plotMethod == "scatter") {
    plot <- plot + geom_point(tooltip) +
      getLabs(facets, plotMethod)
  } else if (plotMethod == "box") {
    plot <- plot + geom_boxplot(tooltip) +
      coord_flip() +
      getLabs(facets, plotMethod)
  } else if (plotMethod == "sina") {
    plot <- plot + ggforce::geom_sina() +
      coord_flip() +
      getLabs(facets, plotMethod)
  } else if (plotMethod == "heatmap") {
    plot <- plot + geom_tile(tooltip) +
      # geom_text(aes(label =  sprintf("%1.2f%%", 100*as.numeric(fraction))))+
      getLabs(facets, plotMethod) +
      theme(axis.text.x = element_text(angle = 90)) +
      scale_fill_gradient(low = "white", high = "blue") +
      guides(fill = guide_colorbar(barwith = 0.5, barheight = 20))
  }

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
