library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(ggplot2)
library(omnideconv)


deconvExplorer = function(usr_bulk = NULL, usr_singleCell = NULL, usr_cellAnnotation = NULL, usr_batch = NULL){
  
  # possible Methods for signature interchangability
  methods_reduced = c("Bisque" = "bisque", "BSeq-sc" = "bseqsc", "CDSeq" = "cdseq", "CIBERSORTx" = "cibersortx", "CPM" = "cpm", "DWLS" = "dwls", "MOMF" = "momf")
  
  
  # box definitions
  data_upload_box = box(title = "Upload your Data", status="primary", solidHeader = TRUE, 
                        helpText("If no file is provided the analysis will be run with a sample dataset"), 
                        fileInput("userBulk", "Upload Bulk RNAseq Data"),
                        div(style="margin-top: -20px"), 
                        fileInput("userSingleCell", "Upload Single Cell RNASeq Data"), 
                        div(style="margin-top: -20px"), 
                        fileInput("userCellTypes", "Upload Cell Type Annotations"), 
                        div(style="margin-top: -20px"), 
                        fileInput("userBatchId", "Upload Batch IDs"))
  
  settings_box = box(title="Deconvolution Settings", status="primary", solidHeader = TRUE, 
                     img(src="logo.jpg", width = "100%"), br(),
                     selectInput("deconvMethod", "Deconvolution Method", choices = omnideconv::deconvolution_methods),
                     conditionalPanel(
                       condition = "input.deconvMethod == 'bisque'|| input.deconvMethod == 'cibersortx' || input.deconvMethod == 'dwls' || input.deconvMethod == 'momf'", 
                       selectInput("sigMethod", "Signature Calculation Method", choices = methods_reduced)
                     ),
                     actionButton("deconvolute", "Deconvolute"))
  
  deconv_box = box(title="Deconvolution Result", status = "warning", solidHeader = TRUE, width = 12, 
                   plotly::plotlyOutput("distPlot") %>% withSpinner())
  
  
  ##### BEGIN UI DEFINITION ####´
  ui =  dashboardPage(
    dashboardHeader(title = "Omnideconv"),
    dashboardSidebar(sidebarMenu(
      menuItem("Deconvolution", tabName = "deconv"), 
      menuItem("Further Information", tabName= "fInfo")
    )),
    dashboardBody(tabItems(
      tabItem(tabName="deconv", fluidRow(data_upload_box, settings_box), fluidRow( deconv_box)),
      tabItem(tabName="fInfo", fluidRow(h2("Test")))
    ))
  )
  
  #### END UI DEFINITION ####

  #####BEGIN SERVER ##########
  server = shinyServer(function(input, output) {
    values = reactiveValues() # storing everything
    
    ### ersetzten mit dem laden der user Uploads!
    values$single_cell = omnideconv::single_cell_data_1
    values$cell_annotations = omnideconv::cell_type_annotations_1
    values$batch_ids = omnideconv::batch_ids_1
    values$bulk = omnideconv::bulk
    
    values$deconvolution_result = readRDS("deconvolution_example.rds")
    
    
    
    
    
    signature = reactive(omnideconv::build_model(
      single_cell_object = values$single_cell, 
      cell_type_annotations = values$cell_annotations, 
      method = input$deconvMethod, 
      batch_ids = values$batch_ids, 
      bulk_gene_expression = values$bulk
    ))
    
    
    
    # deconvolute when button is clicked
    observeEvent(input$deconvolute, {
      values$deconvolution_result = 
        omnideconv::deconvolute(
          bulk_gene_expression = values$bulk,
          signature = signature(),
          method = input$deconvMethod,
          single_cell_object = values$single_cell,
          cell_type_annotations = values$cell_annotations,
          batch_ids = values$batch_ids
        )
      }
    )
    
    # render the deconvolution
    output$distPlot = plotly::renderPlotly({
      data = cbind(values$deconvolution_result, samples = rownames(values$deconvolution_result))%>%
        as.data.frame() %>%
        tidyr::pivot_longer(!samples, names_to = "cell_type", values_to = "fraction")
      
      plot = ggplot(data, aes(y = samples, x = as.numeric(fraction), fill=cell_type, 
                              text=paste0("Cell Type: ", cell_type, "\nFraction: ", sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
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
    })
  }) ####END OF Shiny Server
  
  #######END SERVER ########
  
  shiny::shinyApp(ui = ui, server = server)
}


deconvExplorer()