# library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(ggplot2)
library(omnideconv)
library(RColorBrewer)
library(waiter)

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
                       selectInput("signatureMethod", "Signature Calculation Method", choices = methods_reduced)
                     ),
                     conditionalPanel(
                       condition = "input.deconvMethod == 'bseqsc'",
                       fileInput("userMarker", "Marker Genes")
                     ),
                     actionButton("deconvolute", "Deconvolute"), 
                     waiter::useWaitress(),)
  
  deconv_plot_box = box(title="Deconvolution Plot", status = "warning", solidHeader = TRUE, width = 12, 
                        selectInput("plotMethod", "Plot as: ", choices = c("Bar Plot" = "bar", "Scatter" = "scatter", "Jitter Plot" = "jitter", "Box Plot" = "box", "Sina Plot" = "sina", "Heatmap" = "heatmap")),
                        plotly::plotlyOutput("plotBox") %>% withSpinner())
  
  deconv_table_box = box(title="Deconvolution Table", status = "warning", solidHeader = TRUE, width = 12, 
                         DT::dataTableOutput("tableBox") %>%  withSpinner())
  
  deconv_signature_box = box(title="Deconvolution Signature", status = "info", solidHeader = TRUE, width = 12, 
                         downloadButton("signatureDownload", "Download Signature"),
                         DT::dataTableOutput("signatureBox") %>% withSpinner())
  deconv_all_results = box(title = "All Deconvolutions", status = "info", solidHeader = TRUE, width = 12, 
                           selectInput ("computedDeconvMethod", "Deconvolution Method", choices = NULL), 
                           selectInput ("computedSignatureMethod", "Signature Method", choices = NULL), 
                           actionButton("loadDeconvolution", "Load Deconvolution Result"))
  

  # ui definition  ----------------------------------------------------------


  ui =  dashboardPage(
    dashboardHeader(title = "Omnideconv"),
    dashboardSidebar(sidebarMenu(
      menuItem("Deconvolution", tabName = "deconv"), 
      menuItem("Benchmarking", tabName = "benchmark"),
      menuItem("Further Information", tabName= "fInfo")
    )),
    dashboardBody(tabItems(
      tabItem(tabName="deconv", fluidPage(fluidRow(deconv_all_results), fluidRow(data_upload_box, settings_box), fluidRow(deconv_plot_box, deconv_table_box, deconv_signature_box))),
      tabItem(tabName="benchmark", fluidPage(textOutput("test"))), 
      tabItem(tabName="fInfo", fluidPage(includeMarkdown("omnideconv_vignette.md")))
    ))
  )
  

  
  # server definition  ------------------------------------------------------
  
  server = shinyServer(function(input, output, session) {
    ### background datastructure to store several 
    
    # idea: just use a reactive value, access: all_deconvolutions[[deconvMethod_signatureMethod]] = c(DeconvResult, Signature)
    # only storing the results! Not the input files
    # vectors must be the same data structure! should be a problem as both are double matrices 
    # indexing starts at 1 
    
    
    
    all_deconvolutions = reactiveValues()
    
    
    
    
    # methods= unname(omnideconv::deconvolution_methods)
    # 
    # df = data.frame(matrix(nrow = length(methods), ncol = length (methods)))
    # rownames(df) = methods
    # colnames(df) = methods
    
    # for later: run all possible combinations! i and j cover all valid
    # deconvolution / signature combinations 
    
    # for (i in methods){
    #   if (i %in% methods_interchangeable){
    #     for (j in methods_reduced){
    #       df[j, i] = "X" # insert as df[signature, deconvMethod]
    #     }
    #   }
    #   df[i, i] = "X" # todo: skip the ones already inserted in for loop above
    # }
    # 
    
    
    
    
    
    waitress = Waitress$new("#deconvolute", infinite=TRUE)
    
    values = reactiveValues() # storing the current deconvolution 
    
    ### ersetzten mit dem laden der user Uploads!
    values$single_cell = omnideconv::single_cell_data_1
    values$cell_annotations = omnideconv::cell_type_annotations_1
    values$batch_ids = omnideconv::batch_ids_1
    values$bulk = omnideconv::bulk
    
    values$deconvolution_result = readRDS("deconvolution_example.rds")
    
    
    
    signature = reactive(omnideconv::build_model(
      single_cell_object = values$single_cell, 
      cell_type_annotations = values$cell_annotations, 
      method = input$signatureMethod, 
      batch_ids = values$batch_ids, 
      bulk_gene_expression = values$bulk
    ))
    
    
    
    # deconvolute when button is clicked
    observeEvent(input$deconvolute, {
      waitress$start()
      signature_Method = input$signatureMethod
      if (!(input$deconvMethod %in% methods_interchangeable)){
        signature_Method = input$deconvMethod
        
      }
      message(paste0("Starting Deconvolution. Deconvolution: ", input$deconvMethod, ", Signature: ", signature_Method))
     
      # save reactive siganture to values
      values$signature = signature()
      
      #### At this Point the correct Variables for the Deconvolution are: 
      #### Deconvolution Method: input$deconvMethod
      #### Signature Calculation Method: signature_Method
      #### Signature: values$signature
      
      values$deconvolution_result =
        omnideconv::deconvolute(
          bulk_gene_expression = values$bulk,
          signature = values$signature,
          method = input$deconvMethod,
          single_cell_object = values$single_cell,
          cell_type_annotations = values$cell_annotations,
          batch_ids = values$batch_ids
        )
      
      # insert result into the all_deconvolutions reactive Value
      all_deconvolutions[[paste0(input$deconvMethod, "_", signature_Method)]] = list(values$deconvolution_result, values$signature)
      
      # update select input for selection, just the deconv Method
      deconv_choices = names(all_deconvolutions) %>%
                       strsplit("_") %>%
                       unlist()
      deconv_choices = deconv_choices[seq(1, length(deconv_choices), 2)] # jede 2, startend von 1
      updateSelectInput(session, inputId = "computedDeconvMethod" , choices = deconv_choices)

      waitress$close()
      }
    )
    
    # load Deconvolution result
    observeEvent(input$loadDeconvolution, {
      result = all_deconvolutions[[paste0(input$computedDeconvMethod, "_", input$computedSignatureMethod)]]
      message(typeof(result))
      # values$deconvolution_result = result[1]
      # values$signature = result[2]
    })
    
    # update signature Method choices when selecting deconvolution to load 
    observe({
      signature_choices = names(all_deconvolutions) %>%
        stringr::str_subset(pattern = paste0(input$computedDeconvMethod, "_"))  %>%
        strsplit("_") %>%
        unlist()
      
      if (length(signature_choices)>1){
        signature_choices = signature_choices[seq(2, length(signature_choices), 2)] # jede zweite ab dem zweiten 
      } else {
        signature_choices = signature_choices[2] # nur das zweite
      }

      updateSelectInput(session, inputId = "computedSignatureMethod", choices = signature_choices)
    })
    

    # Plots -------------------------------------------------------------------
    
    output$test = renderText({
      names(all_deconvolutions)
    })
    
    output$plotBox  = plotly::renderPlotly({
      data = cbind(values$deconvolution_result, samples = rownames(values$deconvolution_result)) %>%
        as.data.frame() %>%
        tidyr::pivot_longer(!samples, names_to = "cell_type", values_to = "fraction")
    
      #### potential calc of tooltip text
      
      plot = ggplot(data, aes(x = reorder(cell_type, desc(cell_type)), y = as.numeric(fraction), fill=cell_type, text=""))
    
      if (input$plotMethod == "bar"){
        plot = plot + geom_col(aes(y = samples, x = as.numeric(fraction), fill = cell_type, 
                              text=paste0("Cell Type: ", cell_type, "\nFraction: ", 
                              sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
                      labs(x = "estimated fraction", y = "sample", fill = "cell type")
        
      } else if (input$plotMethod == "jitter"){
        plot = plot + geom_jitter(aes(color = cell_type, text = paste0("Sample: ", samples, "\nFraction: ", 
                                  sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
                      labs(x = "cell type", y = "estimated fraction", color="cell type", fill= "")
        
      } else if (input$plotMethod == "scatter"){
        plot = plot + geom_point(aes(color = cell_type, text = paste0("Sample: ", samples, "\nFraction: ", 
                                 sprintf("%1.2f%%", 100*as.numeric(fraction))))) +
                      labs(x = "cell type", y = "estimated fraction", color="cell type", fill= "")
        
      } else if (input$plotMethod == "box"){
        plot = plot + geom_boxplot(aes(text="")) + 
                      labs(x = "cell type", y = "estimated fraction", fill="cell type")
        
      } else if (input$plotMethod == "sina"){
        plot = plot + geom_violin(colour = "grey", fill="grey") + 
                      ggforce::geom_sina() +
                      labs(x = "cell type", y = "estimated fraction", alpha = "", fill = "cell type")
        
      } else if (input$plotMethod == "heatmap"){
        plot = plot + geom_tile(aes(y = samples, fill = as.numeric(fraction), text =  sprintf("%1.2f%%", 100*as.numeric(fraction)))) + 
                      #geom_text(aes(y = samples, x = cell_type, label =  sprintf("%1.2f%%", 100*as.numeric(fraction))))+
                      labs(x = "cell type", y = "sample", fill = "estimated fraction") +
                      scale_fill_gradient(low = "white", high = "blue") + 
                      guides(fill = guide_colorbar(barwith = 0.5, barheight = 20))
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
    
    output$signatureBox = DT::renderDataTable({
                          DT::datatable(values$signature)# %>%
                         # DT::formatRound(c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK"), 2)
      })
    
    
    output$signatureDownload = downloadHandler(
                              filename = function(){
                                paste("signature", ".csv", sep="")
                              }, 
                              content = function(file){
                                write.csv(values$signature, file)
                              }
    )
    
  })
  
  shiny::shinyApp(ui = ui, server = server)
}


deconvExplorer()