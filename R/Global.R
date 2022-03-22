

# functions ---------------------------------------------------------------

# named deconvolution list
returnSelectedDeconvolutions <- function(to_plot_list, all_deconvolutions){
  
  deconvolutions <- list()
  
  for (deconvolution in to_plot_list) {
    deconvolutions[deconvolution] <- all_deconvolutions[deconvolution]
  }
  
  return (deconvolutions)
}





#' Plot Deconvolution results
#'
#' @param deconvolutions A named list of deconvolution results
#' @param plotMethod Type of plot to be rendered  ("bar", "jitter", "scatter", "box", "heatmap")
#' @param facets Variable for grouping the plots ("method", "cell_type", "sample")
#' @param palette RColorBrewer palette name (optional), standard = "Set1"
#'
#' @returns ggplot rendered by plotly for interactivity
#' 
#' @examples 
#' 
#' @export
plot_deconvolution <- function(deconvolutions, plotMethod, facets, palette = "Set1") {
  # data needs to be a named deconvolution list
  if (is.null(names(deconvolutions))){
    stop("Please supply a NAMED list, names(deconvolution) returns NULL")
  }
  
  # check plot Method
  if (!(plotMethod %in% c("bar", "jitter", "scatter", "box", "heatmap"))){
    stop("plot_method not supported")
  }
  
  # check facet parameter
  if (!(facets %in% c("method", "cell_type", "sample"))){
    stop("facet not supported. Please provide one of the following: 'method', 'cell_type', 'sample'")
  }
  
  

  # preformat data into a dataframe
  deconvolutions <- lapply(deconvolutions, function(deconvolution) {
    cbind(deconvolution, sample = rownames(deconvolution)) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(!sample, names_to = "cell_type", values_to = "fraction")
  })
  
  # combine list to one dataframe and add calculation method
  data <- do.call("rbind", deconvolutions)
  data$fraction <- as.numeric(data$fraction) # change datatype of column
  data$method <- rep(names(deconvolutions), each = nrow(data) / length(names(deconvolutions))) # add computation method as column

  # calculate tooltip based on chosen facet
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
  
  # axis information
  axis <- list(
    "method" = list( # x, y, fill/color
      "bar" = list("fraction", "sample", "cell_type"),
      "jitter" = list("fraction", "cell_type", "cell_type"),
      "scatter" = list("fraction", "cell_type", "cell_type"),
      "box" = list("cell_type", "fraction", "cell_type"),
      "heatmap" = list("cell_type", "sample", "fraction")
    ),
    "cell_type" = list(
      "bar" = list("fraction", "sample", "method"),
      "jitter" = list("fraction", "method", "sample"),
      "scatter" = list("fraction", "method", "sample"),
      "box" = list("method", "fraction", "method"),
      "heatmap" = list("sample", "method", "fraction")
    ),
    "sample" = list(
      "bar" = list("fraction", "cell_type", "method"),
      "jitter" = list("fraction", "cell_type", "method"),
      "scatter" = list("fraction", "cell_type", "method"),
      "box" = list("method", "fraction", "method"),
      "heatmap" = list("cell_type", "method", "fraction")
    )
  )
  
  # extract ggplot aes
  x <- axis[[facets]][[plotMethod]][[1]] # x axis
  y <- axis[[facets]][[plotMethod]][[2]] # y axis
  col <- axis[[facets]][[plotMethod]][[3]] # color / fill, depending on plot type
  
  aes <- NULL
  
  if (plotMethod %in% c("jitter", "scatter")) {
    aes <- aes_string(x = x, y = y, col = col)
  } else {
    aes <- aes_string(x = x, y = y, fill = col)
  }

  # put plot together
  plot <- ggplot(data, aes)
  plot <- plot + facet_wrap(~ data[[facets]])

  # general theme
  plot <- plot + bbplot::bbc_style() +
    theme(
      legend.title = element_text(size = 16), # legend title font size
      legend.text = element_text(size = 14), # legend element font size
      axis.text.x = element_text(size = 14), # x axis font size
      axis.text.y = element_text(size = 14), # y axis font size
      axis.ticks.x = ggplot2::element_line(colour = "#333333"), # vertical ticks for fractions
      axis.ticks.length = grid::unit(0.26, "cm"), # tick length in cm
      strip.text = element_text(size = 16)
    )
  
  # add plot content
  if (plotMethod == "bar") {
    if (facets != "method") {
      plot <- plot + geom_col(tooltip, position = "dodge") # not stacked
    } else { 
      plot <- plot + geom_col(tooltip) +
        theme(panel.grid.major.y = ggplot2::element_blank()) # remove vertical lines
    }
  } else if (plotMethod == "jitter") {
    plot <- plot + geom_jitter(tooltip) 
  } else if (plotMethod == "scatter") {
    plot <- plot + geom_point(tooltip) 
  } else if (plotMethod == "box") {
    plot <- plot + geom_boxplot(tooltip) +
      coord_flip() # this is mandatory here
  } else if (plotMethod == "heatmap") {
    plot <- plot + geom_tile(tooltip) +
      theme(axis.text.x = element_text(angle = 90)) +
      guides(fill = guide_colorbar(barwith = 0.5, barheight = 20)) 
  }
  
  # add color theme based on plot method
  if (plotMethod %in% c("jitter", "scatter")) {
    plot <- plot + ggplot2::scale_colour_brewer(palette=palette)
  } else if (plotMethod =="heatmap"){
    scale_fill_gradient(low = "white", high = RColorBrewer::brewer.pal(3, palette)[1:1])
  } else {
    plot <- plot + ggplot2::scale_fill_brewer(palette=palette)
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

plot_benchmark <- function(deconvolutions) {
  # import and preformat data
  
  #View(deconvolutions)
  
  #print ("BENCHMARK")

  # deconvolution_list <- list()
  # for (deconvolution in to_plot_list) {
  #   # deconvolution_list[length(deconvolution_list) + 1] <- all_deconvolutions[[deconvolution]][1]
  #   deconvolution_list[deconvolution] <- all_deconvolutions[[deconvolution]][1]
  # }

  # add samples and deconvolution method
  deconvolutions <- lapply(deconvolutions, function(x) cbind(x, sample = rownames(x)))
  deconvolutions <- lapply(names(deconvolutions), function(x) {
    cbind(deconvolutions[[x]], method = rep(x, nrow(deconvolutions[[x]])))
  })

  deconvolutions <- lapply(deconvolutions, function(x) {
    tidyr::pivot_longer(data.frame(x), !c("sample", "method"),
      names_to = "cell_type", values_to = "predicted_fraction"
    )
  })

  # combine to one dataframe
  data <- do.call("rbind", deconvolutions)

  # preformat reference data
  ref <- omnideconv::RefData
  ref$sample <- rownames(ref)
  ref <- tidyr::pivot_longer(ref, !sample, names_to = "cell_type", values_to = "true_fraction")

  # merge the reference data with the deconvolution results
  data <- merge(ref, data, by = c("sample", "cell_type"))

  # change datatype to numeric
  data$predicted_fraction <- as.numeric(data$predicted_fraction)
  data$true_fraction <- as.numeric(data$true_fraction)

  # calculate max width/heigth -> plot symmetric and line @ 45 Degrees
  max_value <- max(max(data$true_fraction), max(data$predicted_fraction)) + 0.1

  # create plot
  plot <- ggplot(data, aes(
    x = true_fraction, y = predicted_fraction, color = cell_type,
    text = paste0("Sample: ", sample, "\nTrue: ", true_fraction, "\nPredicted: ", predicted_fraction)
  )) +
    geom_point(size = 4) +
    facet_wrap(~method) +
    geom_abline(color = "black") +
    labs(x = "True Fraction", y = "predicted Fraction", color = "cell type") +
    coord_cartesian(xlim = c(0, max_value), ylim = c(0, max_value))

  # render
  plotly::ggplotly(plot, tooltip = c("text")) %>%
    plotly::config(
      displaylogo = FALSE, showTips = FALSE, toImageButtonOptions = list(filename = paste0("plotMethod", "_plot")),
      modeBarButtonsToRemove = list(
        "hoverClosestCartesian",
        "hoverCompareCartesian",
        "zoomIn2d", "zoomOut2d",
        "lasso2d", "zoom2d",
        "pan2d", "autoScale2d", "select2d"
      )
    ) %>%
    plotly::layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE))
}
