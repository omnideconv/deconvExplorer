
library(dplyr)
library(ggplot2)


#### Helper #####
plot_deconv_barplot = function(deconvolution){
  data = cbind(deconvolution, samples = rownames(deconvolution))%>%
    as.data.frame() %>%
    tidyr::pivot_longer(!samples, names_to = "cell_type", values_to = "fraction")
  
  plot = ggplot(data, aes(y = samples, x = as.numeric(fraction), fill=cell_type, text=paste0("Cell Type: ", cell_type, "\nFraction: ", sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
    geom_bar( stat="identity") + 
    labs(x = "predicted fraction", y = "sample", fill="cell type") 
  # theme_fivethirtyeight()
  
  plotly::ggplotly(plot, tooltip = c("text"))%>%
  plotly::config(displaylogo = FALSE, showTips = FALSE, toImageButtonOptions = list(filename = "plot.png"), 
    modeBarButtonsToRemove = list("hoverCLosestCartesian", 
                                  "hoverCompareCartesian", 
                                  "zoomIn2d", "zoomOut2d", 
                                  "lasso2d", "zoom2d", 
                                  "pan2d", "autoScale2d", "select2d" )) %>%
  plotly::layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE))
  
  #return (final_plot)
  
  
  
}


function(input, output) {
  
  single_cell_data = omnideconv::single_cell_data_1
  cell_type_annotations <- omnideconv::cell_type_annotations_1
  batch_ids <- omnideconv::batch_ids_1
  bulk <- omnideconv::bulk
  
  # reactive Signature Matrix
  signatureMatrix = reactive(omnideconv::build_model(
    single_cell_object = single_cell_data,
    cell_type_annotations = cell_type_annotations,
    method=input$deconvMethod,
    batch_ids = batch_ids, 
    bulk_gene_expression = bulk))
  
  # Onclick Deconvolution 
  
  # trying to load an example before
  t = readRDS("deconvolution_example.rds")
  output$distPlot = plotly::renderPlotly(plot_deconv_barplot(t))
  
  # deconvolute when button is clicked
  deconvolution = eventReactive(input$deconvolute, omnideconv::deconvolute(bulk, signatureMatrix(), input$deconvMethod, single_cell_data, cell_type_annotations, batch_ids))
  
  # render the deconvolution
  output$distPlot = plotly::renderPlotly(plot_deconv_barplot(deconvolution()))
}