library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(knitr)
library(ggplot2)
library(omnideconv)
library(RColorBrewer)


deconvExplorer = function(usr_bulk = NULL, usr_singleCell = NULL, usr_cellAnnotation = NULL, usr_batch = NULL){
  
  
  # possible Methods for signature interchangability
  methods_reduced = c("Bisque" = "bisque", "BSeq-sc" = "bseqsc", "CDSeq" = "cdseq", "CIBERSORTx" = "cibersortx", "CPM" = "cpm", "DWLS" = "dwls", "MOMF" = "momf")
  methods_interchangeable = c("Bisque" = "bisque", "CIBERSORTx" = "cibersortx", "DWLS" = "dwls", "MOMF" = "momf")
  
  # box definitions ---------------------------------------------------------
  data_upload_box = box(title = "Upload your Data", status="primary", solidHeader = TRUE, height = "31.5em", 
                        helpText("If no file is provided the analysis will be run with a sample dataset"), 
                        fileInput("userBulk", "Upload Bulk RNAseq Data"),
                        div(style="margin-top: -20px"), 
                        fileInput("userSingleCell", "Upload Single Cell RNASeq Data"), 
                        div(style="margin-top: -20px"), 
                        fileInput("userCellTypes", "Upload Cell Type Annotations"), 
                        div(style="margin-top: -20px"), 
                        fileInput("userBatchId", "Upload Batch IDs"))



  
  settings_box = box(title="Deconvolution Settings", status="primary", solidHeader = TRUE, height = "31.5em", 
                     img(src="logo.jpg", width = "100%"), br(),
                     selectInput("deconvMethod", "Deconvolution Method", choices = omnideconv::deconvolution_methods),
                     conditionalPanel(
                       condition = "input.deconvMethod == 'bisque'|| input.deconvMethod == 'cibersortx' || input.deconvMethod == 'dwls' || input.deconvMethod == 'momf'",
                       selectInput("sigMethod", "Signature Calculation Method", choices = methods_reduced)
                     ),
                     actionButton("deconvolute", "Deconvolute"))
  
  deconv_plot_box = box(title="Deconvolution Plot", status = "warning", solidHeader = TRUE, width = 12, 
                        selectInput("plotMethod", "Plot as: ", choices = c("Bar Plot" = "bar", "Scatter" = "scatter", "Jitter Plot" = "jitter", "Box Plot" = "box", "Sina Plot" = "sina", "Heatmap" = "heatmap")),
                        plotly::plotlyOutput("plotBox") %>% withSpinner())
  
  deconv_table_box = box(title="Deconvolution Table", status = "warning", solidHeader = TRUE, width = 12, 
                         DT::dataTableOutput("tableBox") %>%  withSpinner())
  

  # ui definition  ----------------------------------------------------------


  ui =  dashboardPage(
    dashboardHeader(title = "Omnideconv"),
    dashboardSidebar(sidebarMenu(
      menuItem("Deconvolution", tabName = "deconv"), 
      menuItem("Further Information", tabName= "fInfo")
    )),
    dashboardBody(tabItems(
      tabItem(tabName="deconv", fluidRow(data_upload_box, settings_box), fluidRow(deconv_plot_box), fluidRow(deconv_table_box)),
      tabItem(tabName="fInfo")
    ))
  )
  

  # server definition  ------------------------------------------------------
  
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
      method = input$sigMethod, 
      batch_ids = values$batch_ids, 
      bulk_gene_expression = values$bulk
    ))
    
    
    
    # deconvolute when button is clicked
    observeEvent(input$deconvolute, {
      signature_Method = input$sigMethod
      if (!(input$deconvMethod %in% methods_interchangeable)){
        signature_Method = input$deconvMethod
      }
      message(paste0("Starting Deconvolution. Deconvolution: ", input$deconvMethod, " Signature: ", signature_Method))
      
      #### At this Point the correct Variables for the Deconvolution are: 
      #### Deconvolution Method: input$deconvMethod
      #### Signature Calculation Method: signature_Method
      
      
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
    

    # Plots -------------------------------------------------------------------
    
    output$plotBox  = plotly::renderPlotly({
      data = cbind(values$deconvolution_result, samples = rownames(values$deconvolution_result)) %>%
        as.data.frame() %>%
        tidyr::pivot_longer(!samples, names_to = "cell_type", values_to = "fraction")
    
      #### potential calc of tooltip text
      
      plot = ggplot(data, aes(x = reorder(cell_type, desc(cell_type)), y = as.numeric(fraction), fill=cell_type)) 
    
      if (input$plotMethod == "bar"){
        plot = plot + geom_col(aes(y = samples, x = as.numeric(fraction), fill = cell_type, 
                                   text=paste0("Cell Type: ", cell_type, "\nFraction: ", 
                                   sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
                      labs(x = "predicted fraction", y = "sample", fill = "cell type")
      } else if (input$plotMethod == "jitter"){
        plot = plot + geom_jitter(aes(color = cell_type, text = paste0("Sample: ", samples, "\nFraction: ", 
                                      sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
                      labs(x = "cell type", y = "predicted fraction", color="cell type", fill= "")
      } else if (input$plotMethod == "scatter"){
        plot = plot + geom_point(aes(color = cell_type, text = paste0("Sample: ", samples, "\nFraction: ", 
                                                                       sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
          labs(x = "cell type", y = "predicted fraction", color="cell type", fill= "")
      } else if (input$plotMethod == "box"){
        plot = plot + geom_boxplot(aes(text="")) + # optional aes(color= cell_type), black bars get invisible
                      labs(x = "cell type", y = "predicted fraction", fill="cell type")
      } else if (input$plotMethod == "sina"){
        plot = plot + geom_violin(aes(alpha = 0.1, hoverinfo="none"), colour = "grey", fill="grey") + 
                      ggforce::geom_sina(aes(text="")) +
                      labs(x = "cell type", y = "predicted fraction", alpha = "", fill = "cell type")
      } else if (input$plotMethod == "heatmap"){
        plot = plot + geom_tile(aes(y = samples, fill = as.numeric(fraction), text = sprintf("%1.2f%%", 100*as.numeric(fraction)))) + 
                      labs(x = "cell type", y = "sample", fill = "predicted fraction")
      }
      
      # render
      plotly::ggplotly(plot, tooltip = c("text"))%>%
        plotly::config(displaylogo = FALSE, showTips = FALSE, toImageButtonOptions = list(filename = paste0(input$plotMethod, "_plot")),
                           modeBarButtonsToRemove = list("hoverCLosestCartesian",
                                                         "hoverCompareCartesian",
                                                         "zoomIn2d", "zoomOut2d",
                                                         "lasso2d", "zoom2d",
                                                         "pan2d", "autoScale2d", "select2d" )) %>%
        plotly::layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE))
    })
    
    output$tableBox = DT::renderDataTable({
                      DT::datatable(values$deconvolution_result, 
                                    extensions = "Buttons", 
                                    options = list(dom = "Bfrtip", 
                                                   buttons = c("copy", "csv", "excel", "pdf"))) %>%
                        DT::formatPercentage(c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK"), 2)
                      
      })
   
  })
  
  shiny::shinyApp(ui = ui, server = server)
}


deconvExplorer()