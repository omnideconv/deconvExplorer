# library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)
library(ggplot2)
library(omnideconv)
library(RColorBrewer)
library(waiter)

source("Global.R")



deconvExplorer <- function(usr_bulk = NULL, usr_singleCell = NULL, usr_cellAnnotation = NULL, usr_batch = NULL) {


  # possible Methods for signature interchangability
  methods_reduced <- c("Bisque" = "bisque", "BSeq-sc" = "bseqsc", "CDSeq" = "cdseq", "CIBERSORTx" = "cibersortx", "CPM" = "cpm", "DWLS" = "dwls", "MOMF" = "momf")
  methods_interchangeable <- c("Bisque" = "bisque", "CIBERSORTx" = "cibersortx", "DWLS" = "dwls", "MOMF" = "momf")

  # box definitions ---------------------------------------------------------
  data_upload_box <- box(
    title = "Upload your Data", status = "primary", solidHeader = TRUE, height = "31.5em",
    helpText("If no file is provided the analysis will be run with a sample dataset"),
    fileInput("userBulk", "Upload Bulk RNAseq Data"),
    div(style = "margin-top: -20px"),
    fileInput("userSingleCell", "Upload Single Cell RNASeq Data"),
    div(style = "margin-top: -20px"),
    fileInput("userCellTypes", "Upload Cell Type Annotations"),
    div(style = "margin-top: -20px"),
    fileInput("userBatchId", "Upload Batch IDs")
  )




  settings_box <- box(
    title = "Deconvolution Settings", status = "primary", solidHeader = TRUE, height = "31.5em",
    img(src = "logo.jpg", width = "100%"), br(),
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
    waiter::useWaitress(),
  )

  deconv_plot_box <- box(
    title = "Deconvolution Plot", status = "warning", solidHeader = TRUE, width = 12,
    selectInput("plotMethod", "Plot as: ", choices = c("Bar Plot" = "bar", "Scatter" = "scatter", "Jitter Plot" = "jitter", "Box Plot" = "box", "Sina Plot" = "sina", "Heatmap" = "heatmap")),
    plotly::plotlyOutput("plotBox") %>% withSpinner()
  )

  deconv_table_box <- box(
    title = "Deconvolution Table", status = "warning", solidHeader = TRUE, width = 12,
    DT::dataTableOutput("tableBox") %>% withSpinner()
  )

  deconv_signature_box <- box(
    title = "Deconvolution Signature", status = "info", solidHeader = TRUE, width = 12,
    downloadButton("signatureDownload", "Download Signature"),
    DT::dataTableOutput("signatureBox") %>% withSpinner()
  )
  deconv_all_results <- box(
    title = "All Deconvolutions", status = "info", solidHeader = TRUE, width = 12,
    selectInput("computedDeconvMethod", "Deconvolution Method", choices = NULL),
    selectInput("computedSignatureMethod", "Signature Method", choices = NULL),
    actionButton("loadDeconvolution", "Load Deconvolution Result"),
    actionButton("addToPlot", "Compare: Add to Plot"),
    actionButton("removeFromPlot", "Compare: Remove from Plot") # ,
    # downloadButton("downloadSession", "Download Session"),
    # fileInput("uploadSession", "Upload Session File")
  )


  # ui definition  ----------------------------------------------------------


  ui <- dashboardPage(
    dashboardHeader(
      title = "Omnideconv",
      dropdownMenu(
        type = "task",
        icon = icon("download"),
        headerText = "Use session to proceed your work later",
        badgeStatus = NULL,
        notificationItem(text = downloadButton("downloadSession", "Download Session")),
        notificationItem(text = fileInput("uploadSession", "Upload Session File"))
      )
    ),
    dashboardSidebar(sidebarMenu(
      shinyjs::useShinyjs(),
      menuItem("Deconvolution", tabName = "deconv"),
      menuItem("Benchmark", tabName = "benchmark"),
      menuItem("Further Information", tabName = "fInfo")
    )),
    dashboardBody(tabItems(
      tabItem(tabName = "deconv", fluidPage(fluidRow(deconv_all_results), fluidRow(data_upload_box, settings_box), fluidRow(deconv_plot_box))), # remove 3 ), # deconv_table_box, deconv_signature_box))),
      tabItem(tabName = "benchmark", fluidPage()),
      tabItem(tabName = "fInfo", fluidPage(includeMarkdown("omnideconv_vignette.md")))
    ))
  )



  # server definition  ------------------------------------------------------

  server <- shinyServer(function(input, output, session) {
    ### background datastructure to store several deconvolution results


    all_deconvolutions <- reactiveValues()

    ### es gibt eine reactive values to list funktion für den export!


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


    # idea: change values$deconv/signature to a named list (name= deconv/signature), the plot function would need to be changed to take a named list
    # (in values$deconvolution, etc) as an input, so everytime the list updates the plots are recalced,
    # running a new deconv: list is overwritten and only the newest result is displayed
    # adding deconvs to the plot with conditional buttons (add / remove), condition checks list for selected deconvolution
    # idea: values only contains uploaded data and a list of deconvolution results, not doubled! The plot function then takes the list of
    # "to Plot" deconvolutions as an input and collects the datasets from all_deconvolutions

    # TODO:
    # - Plot Funktion für mehrere Ergebnisse Optimieren
    # - what to do with the signatures? right now, just save them, but is it possible to recycle them?


    waitress <- Waitress$new("#deconvolute", infinite = TRUE)

    values <- reactiveValues() # storing the current deconvolution

    ### ersetzten mit dem laden der user Uploads!
    values$single_cell <- omnideconv::single_cell_data_1
    values$cell_annotations <- omnideconv::cell_type_annotations_1
    values$batch_ids <- omnideconv::batch_ids_1
    values$bulk <- omnideconv::bulk

    # values$deconvolution_result <- readRDS("deconvolution_example.rds")
    values$deconvolution_result <- c("bisque_bisque")
    all_deconvolutions[["bisque_bisque"]] <- list(readRDS("deconvolution_example.rds"), readRDS("signature_example.rds"))


    # Observers and Eventhandling ---------------------------------------------



    # restore session with file upload
    observeEvent(input$uploadSession, {
      sessionFile <- readRDS(input$uploadSession$datapath)

      for (deconvolution in names(sessionFile)) {
        all_deconvolutions[[deconvolution]] <- sessionFile[[deconvolution]]
        message("Loaded Deconvolution: ", deconvolution)
      }
    })

    signature <- reactive(omnideconv::build_model(
      single_cell_object = values$single_cell,
      cell_type_annotations = values$cell_annotations,
      method = input$signatureMethod,
      batch_ids = values$batch_ids,
      bulk_gene_expression = values$bulk
    ))

    # deconvolute when button is clicked
    observeEvent(input$deconvolute, {
      ### todo: add deconvolution to the "to plot" list
      waitress$start()

      signature_Method <- input$signatureMethod
      if (!(input$deconvMethod %in% methods_interchangeable)) {
        signature_Method <- input$deconvMethod
      }
      message(
        paste0("Starting Deconvolution. Deconvolution: ", input$deconvMethod, ", Signature: ", signature_Method)
      )

      signature <- signature()

      deconvolution_result <-
        omnideconv::deconvolute(
          bulk_gene_expression = values$bulk,
          signature = signature,
          method = input$deconvMethod,
          single_cell_object = values$single_cell,
          cell_type_annotations = values$cell_annotations,
          batch_ids = values$batch_ids
        )

      # insert result into the all_deconvolutions reactive Value
      all_deconvolutions[[paste0(input$deconvMethod, "_", signature_Method)]] <- list(deconvolution_result, signature)

      waitress$close()
    })

    # update Deconvolution Method choices when all_deconvolution changes
    observe({
      deconv_choices <- names(all_deconvolutions) %>%
        strsplit("_") %>%
        unlist()
      deconv_choices <- deconv_choices[seq(1, length(deconv_choices), 2)] # jede 2, startend von 1
      updateSelectInput(session, inputId = "computedDeconvMethod", choices = deconv_choices)
    })

    # update signature Method choices when selecting deconvolution to load
    observeEvent(input$computedDeconvMethod, {
      signature_choices <- names(all_deconvolutions) %>%
        stringr::str_subset(pattern = paste0(input$computedDeconvMethod, "_")) %>%
        strsplit("_") %>%
        unlist()

      if (length(signature_choices) > 1) {
        signature_choices <- signature_choices[seq(2, length(signature_choices), 2)] # jede zweite ab dem zweiten
      } else {
        signature_choices <- signature_choices[2] # nur das zweite
      }

      updateSelectInput(session, inputId = "computedSignatureMethod", choices = signature_choices)
    })


    # add Deconvolution to To Plot list
    observeEvent(input$addToPlot, {
      tmp <- values$deconvolution_result
      values$deconvolution_result <- c(tmp, getSelectionToPlot())
    })

    # remove deconvolution from To Plot list
    observeEvent(input$removeFromPlot, {
      tmp <- values$deconvolution_result
      tmp <- tmp[!tmp %in% getSelectionToPlot()]
      values$deconvolution_result <- tmp
    })

    # load Deconvolution result
    observeEvent(input$loadDeconvolution, {
      values$deconvolution_result <- c(getSelectionToPlot())
      values$signature <- c(getSelectionToPlot())
    })

    # observe the selection to plot and show buttons if conditions match
    observe({
      if (getSelectionToPlot() %in% values$deconvolution_result) {
        shinyjs::hide("addToPlot")
        shinyjs::show("removeFromPlot")
      } else {
        shinyjs::hide("removeFromPlot")
        shinyjs::show("addToPlot")
      }
    })

    # Plots -------------------------------------------------------------------

    output$plotBox <- plotly::renderPlotly(plot_deconvolution(values$deconvolution_result, input$plotMethod, all_deconvolutions))

    output$tableBox <- DT::renderDataTable({
      # update: work with a list of  deconvoltutions from values$deconvolution

      DT::datatable(values$deconvolution_result,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel", "pdf")
        )
      ) %>%
        DT::formatPercentage(c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK"), 2)
    })

    output$signatureBox <- DT::renderDataTable({
      # update: work with a list of  sigantures from values$signature

      DT::datatable(values$signature) # %>%
      # DT::formatRound(c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK"), 2)
    })


    output$signatureDownload <- downloadHandler(
      filename = function() {
        paste("signature", ".csv", sep = "")
      },
      content = function(file) {
        write.csv(values$signature, file)
      }
    )

    output$downloadSession <- downloadHandler(
      filename = function() {
        paste0("session.rds")
      },
      content = function(file) {
        saveRDS(reactiveValuesToList(all_deconvolutions), file)
      }
    )

    # functions ---------------------------------------------------------------

    getSelectionToPlot <- function() {
      return(paste0(input$computedDeconvMethod, "_", input$computedSignatureMethod))
    }
  })

  shiny::shinyApp(ui = ui, server = server)
}


deconvExplorer()
