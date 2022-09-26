
#' Get subset of deconvolution results
#'
#' This function returns the requested deconvolution results as a list
#'
#' @param to_plot_list list of deconvolution result names which should be plotted
#' @param all_deconvolutions List containing all available deconvolution results
#'
#' @return list of deconvolution results
#' @export
#'
#' @examples
#' deconv <- readRDS(system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer"))
#' 
#' # list containting deconvolution results
#' deconvList = list("bisque" =deconv, "momf" = deconv)
#' 
#' returnSelectedDeconvolutions(c("momf"), deconvList)
returnSelectedDeconvolutions <- function(to_plot_list, all_deconvolutions){

  deconvolutions <- list()

  for (deconvolution in to_plot_list) {
    deconvolutions[deconvolution] <- all_deconvolutions[deconvolution]
  }

  return (deconvolutions)
}


###########################################
## Deprecated , but commented for testing #
###########################################

#' Title
#'
#' CosTODO
#'
#' @param deconvolutions CosTODO
#'
#' @return CosTODO
#' @export
#'
#' @examples
#' # CosTODO
# plot_benchmark <- function(deconvolutions) {
#   # import and preformat data
#
#   # deconvolution_list <- list()
#   # for (deconvolution in to_plot_list) {
#   #   # deconvolution_list[length(deconvolution_list) + 1] <- all_deconvolutions[[deconvolution]][1]
#   #   deconvolution_list[deconvolution] <- all_deconvolutions[[deconvolution]][1]
#   # }
#
#   # add samples and deconvolution method
#   deconvolutions <- lapply(deconvolutions, function(x) cbind(x, sample = rownames(x)))
#   deconvolutions <- lapply(names(deconvolutions), function(x) {
#     cbind(deconvolutions[[x]], method = rep(x, nrow(deconvolutions[[x]])))
#   })
#
#   deconvolutions <- lapply(deconvolutions, function(x) {
#     tidyr::pivot_longer(data.frame(x), !c("sample", "method"),
#       names_to = "cell_type", values_to = "predicted_fraction"
#     )
#   })
#
#   # combine to one dataframe
#   data <- do.call("rbind", deconvolutions)
#
#   # preformat reference data
#   ref <- omnideconv::RefData
#   ref$sample <- rownames(ref)
#   ref <- tidyr::pivot_longer(ref, !sample, names_to = "cell_type", values_to = "true_fraction")
#
#   # merge the reference data with the deconvolution results
#   data <- merge(ref, data, by = c("sample", "cell_type"))
#
#   # change datatype to numeric
#   data$predicted_fraction <- as.numeric(data$predicted_fraction)
#   data$true_fraction <- as.numeric(data$true_fraction)
#
#   # calculate max width/heigth -> plot symmetric and line @ 45 Degrees
#   max_value <- max(max(data$true_fraction), max(data$predicted_fraction)) + 0.1
#
#   # create plot
#   plot <- ggplot(data, aes(
#     x = .data$true_fraction, y = predicted_fraction, color = cell_type,
#     text = paste0("Sample: ", sample, "\nTrue: ", true_fraction, "\nPredicted: ", predicted_fraction)
#   )) +
#     geom_point(size = 4) +
#     facet_wrap(~method) +
#     geom_abline(color = "black") +
#     labs(x = "True Fraction", y = "predicted Fraction", color = "cell type") +
#     coord_cartesian(xlim = c(0, max_value), ylim = c(0, max_value))
#
#   # render
#   plotly::ggplotly(plot, tooltip = c("text")) %>%
#     plotly::config(
#       displaylogo = FALSE, showTips = FALSE, toImageButtonOptions = list(filename = paste0("plotMethod", "_plot")),
#       modeBarButtonsToRemove = list(
#         "hoverClosestCartesian",
#         "hoverCompareCartesian",
#         "zoomIn2d", "zoomOut2d",
#         "lasso2d", "zoom2d",
#         "pan2d", "autoScale2d", "select2d"
#       )
#     ) %>%
#     plotly::layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE))
# }
