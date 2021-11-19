# library(shiny)
# library(shinydashboard)
# library(shinycssloaders)
# library(dplyr)
# library(ggplot2)
# #library(omnideconv)
# library(RColorBrewer)
# library(waiter)
# library(rintrojs)

# source("./R/Global.R")

#' Run DeconvExplorer
#'
#' @param usr_bulk Bulk Sequencing data which will be deconvoluted
#' @param usr_singleCell Single Cell Data which is used to calculate the signature matrix
#' @param usr_cellAnnotation Cell Type annotations for the single cell data
#' @param usr_batch Batch IDs, only for some deconvolution methods
#'
#' @export
deconvExplorer <- function(usr_bulk = NULL,
                           usr_singleCell = NULL,
                           usr_cellAnnotation = NULL,
                           usr_batch = NULL) {

  # possible Methods for signature interchangability
  methods_reduced <- c(
    "Bisque" = "bisque", "BSeq-sc" = "bseqsc",
    "CDSeq" = "cdseq", "CIBERSORTx" = "cibersortx",
    "CPM" = "cpm", "DWLS" = "dwls", "MOMF" = "momf"
  )
  methods_interchangeable <- c(
    "Bisque" = "bisque", "CIBERSORTx" = "cibersortx",
    "DWLS" = "dwls", "MOMF" = "momf"
  )

  # box definitions ---------------------------------------------------------
  data_upload_box <- shinydashboard::box(
    title = "Upload your Data", status = "primary",
    solidHeader = TRUE, height = "31.5em",
    introBox(
      helpText("If no file is provided the analysis will be run with a sample dataset"),
      fileInput("userBulk", "Upload Bulk RNAseq Data"),
      div(style = "margin-top: -20px"),
      fileInput("userSingleCell", "Upload Single Cell RNASeq Data"),
      div(style = "margin-top: -20px"),
      fileInput("userCellTypes", "Upload Cell Type Annotations"),
      div(style = "margin-top: -20px"),
      fileInput("userBatchId", "Upload Batch IDs"),
      data.step = 1, data.intro = "Upload your Data"
    )
  )

  settings_box <- shinydashboard::box(
    title = "Deconvolution Settings", status = "primary",
    solidHeader = TRUE, height = "31.5em",
    introBox(
      img(src = system.file("www", "logo.jpg", package = "DeconvExplorer"), width = "100%"), br(),
      selectInput("deconvMethod", "Deconvolution Method",
        choices = omnideconv::deconvolution_methods
      ),
      conditionalPanel(
        condition = "input.deconvMethod == 'bisque'||
                    input.deconvMethod == 'cibersortx' ||
                    input.deconvMethod == 'dwls' ||
                      input.deconvMethod == 'momf'",
        selectInput("signatureMethod", "Signature Calculation Method",
          choices = methods_reduced
        )
      ),
      conditionalPanel(
        condition = "input.deconvMethod == 'bseqsc'",
        fileInput("userMarker", "Marker Genes")
      ),
      actionButton("deconvolute", "Deconvolute"),
      waiter::useWaitress(),
      actionButton("deconvoluteAll", "Deconvolute All"),
      data.step = 2, data.intro = "Select preferred Deconvolution and Signature calculation Method"
    )
  )

  deconv_plot_box <- shinydashboard::box(
    title = span("Deconvolution Plot ", icon("tasks", lib = "glyphicon")),
    status = "warning", solidHeader = TRUE, width = 12,
    introBox(
      selectInput("plotMethod", "Plot as: ",
        choices = c(
          "Bar Plot" = "bar", "Scatter" = "scatter",
          "Jitter Plot" = "jitter", "Box Plot" = "box",
          "Sina Plot" = "sina", "Heatmap" = "heatmap"
        )
      ),
      selectInput("facets", "Group Plots By",
        choices = c(
          "Deconvolution Method" = "method",
          "Cell Type" = "cell_type", "Sample" = "sample"
        )
      ),
      shinycssloaders::withSpinner(
        plotly::plotlyOutput("plotBox")
      ),
      data.step = 4, data.intro = "View the deconvolution results and compare "
    )
  )

  deconv_table_box <- shinydashboard::box(
    title = span("Deconvolution Table ", icon("th", lib = "glyphicon")),
    status = "warning", solidHeader = TRUE, width = 12,
    selectInput("deconvolutionToTable", "Deconvolution Result", choices = NULL),
    shinycssloaders::withSpinner(
      DT::dataTableOutput("tableBox")
    )
  )

  deconv_signature_box <- shinydashboard::box(
    title = span("Deconvolution Signature ", icon("fingerprint")),
    status = "info", solidHeader = TRUE, width = 12,
    downloadButton("signatureDownload", "Download Signature"),
    selectInput("signatureToTable", "Signature", choices = NULL),
    shinycssloaders::withSpinner(
      DT::dataTableOutput("signatureBox")
    )
  )
  deconv_all_results <- shinydashboard::box(
    title = "All Deconvolutions", status = "info", solidHeader = TRUE, width = 12,
    introBox(
      selectInput("computedDeconvMethod", "Deconvolution Method", choices = NULL),
      selectInput("computedSignatureMethod", "Signature Method", choices = NULL),
      actionButton("loadDeconvolution", "Load Deconvolution Result"),
      actionButton("addToPlot", "Compare: Add to Plot"),
      actionButton("removeFromPlot", "Compare: Remove from Plot"),
      data.step = 3, data.intro = "Deconvolution results are stored and can be reloaded for visualization and comparison"
    )
  )

  benchmark_all_results <- shinydashboard::box(
    title = "All Deconvolutions", status = "info", solidHeader = TRUE, width = 12,
    selectInput("computedDeconvMethod", "Deconvolution Method", choices = NULL),
    selectInput("computedSignatureMethod", "Signature Method", choices = NULL),
    fileInput("groundTruth", "Reference Data (Ground Truth)"),
    actionButton("benchmark", "Benchmark")
  )

  benchmark_plot_box <- shinydashboard::box(
    title = "Benchmark", status = "info", solidHeader = TRUE, width = 12,
    shinycssloaders::withSpinner(
      plotly::plotlyOutput("benchmarkPlot")
    )
  )

  # ui definition  ----------------------------------------------------------


  de_ui <- dashboardPage(
    dashboardHeader(
      title = "Omnideconv",
      dropdownMenu(
        type = "task",
        icon = icon("question-circle"),
        headerText = "View a tour or look up the source code",
        badgeStatus = NULL,
        notificationItem(
          text = actionButton("startTour", "Start Tour",
            icon = icon("directions")
          ),
          icon = icon("", verify_fa = FALSE)
        ),
        notificationItem(text = actionButton("githubLink", "View the Code",
          onclick = "window.open('https://github.com/omnideconv', '_blank')",
          icon = icon("github")
        ), icon = icon("", verify_fa = FALSE))
      ),
      dropdownMenu(
        type = "task",
        icon = icon("bookmark", lib = "glyphicon"),
        headerText = "Use session to proceed your work later",
        badgeStatus = NULL,
        notificationItem(
          text = downloadButton("downloadSession", "Download Session"),
          icon = icon("", verify_fa = FALSE)
        ),
        notificationItem(
          text = fileInput("uploadSession", "Upload Session File"),
          icon = icon("", verfiy_fa = FALSE), status = "primary"
        )
      )
    ),
    dashboardSidebar(sidebarMenu(
      shinyjs::useShinyjs(),
      introjsUI(),
      menuItem("Deconvolution", tabName = "deconv"),
      menuItem("Benchmark", tabName = "benchmark"),
      menuItem("Further Information", tabName = "fInfo")
    )),
    dashboardBody(tabItems(
      tabItem(tabName = "deconv", fluidPage(
        fluidRow(deconv_all_results),
        fluidRow(data_upload_box, settings_box),
        fluidRow(deconv_plot_box, deconv_table_box, deconv_signature_box)
      )),
      tabItem(tabName = "benchmark", fluidPage(fluidRow(benchmark_plot_box))),
      tabItem(tabName = "fInfo", fluidPage(
        includeMarkdown(
          system.file("extdata", "omnideconv_vignette.md", package = "DeconvExplorer")
        )
      ))
    ))
  )

  # server definition  ------------------------------------------------------

  de_server <- shinyServer(function(input, output, session) {
    ### background datastructure to store several deconvolution results


    all_deconvolutions <- reactiveValues()


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
    waitress <- Waitress$new("#deconvolute", infinite = TRUE)

    values <- reactiveValues() # storing the current deconvolution

    ### ersetzten mit dem laden der user Uploads!
    values$single_cell <- omnideconv::single_cell_data_1
    values$cell_annotations <- omnideconv::cell_type_annotations_1
    values$batch_ids <- omnideconv::batch_ids_1
    values$bulk <- omnideconv::bulk

    # values$deconvolution_result <- readRDS("deconvolution_example.rds")
    values$deconvolution_result <- c("bisque_bisque")

    all_deconvolutions[["bisque_bisque"]] <- list(
      readRDS(
        system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer")
      ),
      readRDS(
        system.file("extdata", "signature_example.rds", package = "DeconvExplorer")
      )
    )
    # updateTableSelection()

    # Observers and Eventhandling ---------------------------------------------

    # start the tour
    observeEvent(input$startTour, {
      introjs(session, options = list(
        "nextLabel" = ">",
        "prevLabel" = "<",
        "skipLabel" = "X"
      ))
    })
    
    # handle user file upload 
    # ids: userBulk, userSingleCell, userCellTypes,userBatchId
    # from the App function: usr_bulk, usr_singleCell, usr_cellAnnotation, usr_batch
    values$bulk <- {
      req(input$userBulk || usr_bulk) # workd for both!
      # todo
    }

    # restore session with file upload
    observeEvent(input$uploadSession, {
      sessionFile <- readRDS(input$uploadSession$datapath)

      for (deconvolution in names(sessionFile)) {
        all_deconvolutions[[deconvolution]] <- sessionFile[[deconvolution]]
        showNotification(paste0("Loaded Deconvolution: ", deconvolution))
        # message("Loaded Deconvolution: ", deconvolution)
      }
      updateTableSelection()
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
      showNotification("Deconvolution started", type = "warning")

      signature_Method <- input$signatureMethod
      if (!(input$deconvMethod %in% methods_interchangeable)) {
        signature_Method <- input$deconvMethod
      }

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
      showNotification("Deconvolution finished", type = "message")
    })

    # update Deconvolution Method and signature method choices when new deconvolution result is calculated
    observe({
      # deconvolution to load
      deconv_choices <- unlist(strsplit(names(all_deconvolutions), "_"))
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

      # update selection choices for the tables
      updateTableSelection()
    })

    # remove deconvolution from To Plot list
    observeEvent(input$removeFromPlot, {
      tmp <- values$deconvolution_result
      tmp <- tmp[!tmp %in% getSelectionToPlot()]
      values$deconvolution_result <- tmp

      # update selection choices for the tables
      updateTableSelection()
    })

    # load Deconvolution result
    observeEvent(input$loadDeconvolution, {
      values$deconvolution_result <- c(getSelectionToPlot())
      # values$signature <- c(getSelectionToPlot())

      # update selection choices for the tables
      updateTableSelection()
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

    output$plotBox <- plotly::renderPlotly(
      plot_deconvolution(
        values$deconvolution_result,
        input$plotMethod,
        input$facets,
        all_deconvolutions
      )
    )

    output$benchmarkPlot <- plotly::renderPlotly(
      plot_benchmark(
        values$deconvolution_result,
        all_deconvolutions
      )
    )


    # Tables ------------------------------------------------------------------


    output$tableBox <- DT::renderDataTable({
      # update: work with a list of  deconvoltutions from values$deconvolution
      req(input$deconvolutionToTable)
      deconvolution <- all_deconvolutions[[input$deconvolutionToTable]][[1]]

      DT::datatable(deconvolution,
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
      req(input$signatureToTable)

      signature <- all_deconvolutions[[input$signatureToTable]][[2]]

      DT::datatable(signature) # %>%
      # DT::formatRound(c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK"), 2)
    })


    output$signatureDownload <- downloadHandler(
      filename = function() {
        paste("signature", ".csv", sep = "")
      },
      content = function(file) {
        data <- all_deconvolutions[[input$signatureToTable]][[2]]
        write.csv(data, file)
      }
    )

    output$downloadSession <- downloadHandler(
      filename = function() {
        # paste0("session.rds")
        paste0("omnideconv_", Sys.Date(), ".rds")
      },
      content = function(file) {
        saveRDS(reactiveValuesToList(all_deconvolutions), file)
      }
    )

    # functions ---------------------------------------------------------------

    getSelectionToPlot <- function() {
      return(paste0(input$computedDeconvMethod, "_", input$computedSignatureMethod))
    }

    updateTableSelection <- function() {
      updateSelectInput(session, inputId = "deconvolutionToTable", choices = values$deconvolution_result)
      updateSelectInput(session, inputId = "signatureToTable", choices = values$deconvolution_result)
    }
  })

  shiny::shinyApp(ui = de_ui, server = de_server)
}
