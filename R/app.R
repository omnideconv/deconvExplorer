
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
    solidHeader = TRUE, height = "30em", #collapsible = TRUE,
    introBox(
      helpText("If no file is provided the analysis will be run with a sample dataset"),
      fileInput("userBulk", "Upload Bulk RNAseq Data"),
      div(style = "margin-top: -20px"),
      fileInput("userSingleCell", "Upload Single Cell RNASeq Data"),
      div(style = "margin-top: -20px"),
      fileInput("userCellTypeAnnotations", "Upload Cell Type Annotations"),
      div(style = "margin-top: -20px"),
      fileInput("userBatchIDs", "Upload Batch IDs"),
      data.step = 1, data.intro = "Upload your Data. Allowed formats: txt, csv, tsv"
    )
  )

  settings_box <- shinydashboard::box(
    title = "Deconvolution Settings", status = "primary",
    solidHeader = TRUE, height = "30em", #collapsible = TRUE, 
    introBox(
      imageOutput("logo", height = "auto"), br(), 
      column(6,
      selectInput("deconvMethod", "Deconvolution Method",
        choices = omnideconv::deconvolution_methods
      )),
      column(6,
      conditionalPanel(
        condition = "input.deconvMethod == 'bisque'||
                    input.deconvMethod == 'cibersortx' ||
                    input.deconvMethod == 'dwls' ||
                      input.deconvMethod == 'momf'",
        selectInput("signatureMethod", "Signature Calculation Method",
          choices = methods_reduced
        )
      )),
      column(6,
      conditionalPanel(
        condition = "input.deconvMethod == 'bseqsc'",
        fileInput("userMarker", "Marker Genes")
      )),
      column(12,
      actionButton("deconvolute", "Deconvolute"),
      waiter::useWaitress(),
      actionButton("deconvoluteAll", "Deconvolute All")),
      data.step = 2, data.intro = "Select preferred Deconvolution and Signature calculation Method"
    )
  )

  deconv_plot_box <- shinydashboard::box(
    title = span("Deconvolution Plot ", icon("tasks", lib = "glyphicon")),
    status = "warning", solidHeader = TRUE, width = 12,
    introBox(
      column(3,
      selectInput("plotMethod", "Plot as: ",
        choices = c(
          "Bar Plot" = "bar", "Scatter" = "scatter",
          "Jitter Plot" = "jitter", "Box Plot" = "box",
          "Sina Plot" = "sina", "Heatmap" = "heatmap"
        )
      )
      ),
      column(3,
      selectInput("facets", "Group Plots By",
        choices = c(
          "Deconvolution Method" = "method",
          "Cell Type" = "cell_type", "Sample" = "sample"
        )
      )),
      column(12,
      shinycssloaders::withSpinner(
        plotly::plotlyOutput("plotBox")
      )),
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
    title = "Plotting Settings", status = "info", solidHeader = TRUE, width = 12,
    introBox(
      column(3, selectInput("computedDeconvMethod", "Deconvolution Method", choices = NULL)),
      column(3, selectInput("computedSignatureMethod", "Signature Method", choices = NULL)),
      column(4, 
      actionButton("loadDeconvolution", "Load Deconvolution Result", style="margin-top: 1.7em"),
      actionButton("addToPlot", "Compare: Add to Plot", style="margin-top: 1.7em"),
      actionButton("removeFromPlot", "Compare: Remove from Plot", style="margin-top: 1.7em")),
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
        icon = icon("cog", lib = "glyphicon"),
        headerText = "Set your CIBERSORTx Credentials",
        badgeStatus = NULL,
        notificationItem(
          text = textInput("csxEmail", "Email Adress"),
          icon = icon("", verify_fa = FALSE)
        ),
        notificationItem(
          text = textInput("csxToken", "Token"),
          icon = icon("", verify_fa = FALSE)
        ),
        notificationItem(
          text = actionButton("setCSX", "Set CIBERSORTx Credentials"),
          icon = icon("", verify_fa = FALSE)
        )
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
    dashboardBody(
      tags$head(tags$style(
        HTML(".wrapper {height: auto !important; 
             position:relative; overflow-x:hidden; overflow-y:hidden}"))),
      tabItems(
        tabItem(tabName = "deconv", fluidPage(
          fluidRow(data_upload_box, settings_box),
          fluidRow(deconv_all_results),
          fluidRow(deconv_plot_box, deconv_table_box, deconv_signature_box)
        )),
        tabItem(tabName = "benchmark", fluidPage(fluidRow(benchmark_plot_box))),
        tabItem(tabName = "fInfo", fluidPage(
          includeMarkdown(
            system.file("extdata", "omnideconv_vignette.md", package = "DeconvExplorer")
          )
        ))
      )
    )
  )

  # server definition  ------------------------------------------------------

  de_server <- shinyServer(function(input, output, session) {
    ### background datastructure to store several deconvolution results


    # storing all calculated deconvolution
    all_deconvolutions <- reactiveValues()

    userData <- reactiveValues()


    # options
    options(shiny.maxRequestSize = 10 * 1024^2 * 100) # 1GB

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

    ### ersetzten mit dem laden der user Uploads!
    userData$singleCell <- omnideconv::single_cell_data_1
    userData$cellTypeAnnotations <- omnideconv::cell_type_annotations_1
    userData$batchIDs <- omnideconv::batch_ids_1
    userData$bulk <- omnideconv::bulk
    # updateTableSelection()

    # userData$deconvolution_result <- readRDS("deconvolution_example.rds")
    userData$deconvolution_result <- c("bisque_bisque")

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


    # User Upload: Bulk Expression Data, preuploaded files will be overwritten
    observeEvent(input$userBulk, {
      userData$bulk <- loadFile(input$userBulk)
    })

    observeEvent(input$userSingleCell, {
      userData$singleCell <- loadFile(input$userSingleCell)
    })

    observeEvent(input$userCellTypeAnnotations, {
      userData$cellTypeAnnotations <- loadFile(input$userCellTypeAnnotations, type = "vector")
    })

    observeEvent(input$userBatchIDs, {
      userData$batchIDs <- loadFile(input$userBatchIDs, type = "vector")
    })

    observeEvent(input$userMarker, {
      userData$marker <- loadFile(input$userMarker)
    })

    # set CIBERSORTx Credentials from User Input
    observeEvent(input$setCSX, {
      req(input$csxEmail, input$csxToken) # email and token have to be set
      omnideconv::set_cibersortx_credentials(input$csxEmail, input$csxToken)
      showNotification("CIBERSORTx Credentials set")
    })

    # restore session with file upload
    observeEvent(input$uploadSession, {
      sessionFile <- readRDS(input$uploadSession$datapath)
      message(input$uploadSession$datapath)
      for (deconvolution in names(sessionFile)) {
        all_deconvolutions[[deconvolution]] <- sessionFile[[deconvolution]]
        showNotification(paste0("Loaded Deconvolution: ", deconvolution))
        # message("Loaded Deconvolution: ", deconvolution)
      }
      # updateTableSelection()
    })

    # signature <- reactive(omnideconv::build_model(
    #   single_cell_object = userData$singleCell,
    #   cell_type_annotations = userData$cellTypeAnnotations,
    #   method = input$signatureMethod,
    #   batch_ids = userData$batchIDs,
    #   bulk_gene_expression = userData$bulk,
    #   markers = userData$marker
    # ))

    # deconvolute when button is clicked
    observeEvent(input$deconvolute, {
      ### todo: add deconvolution to the "to plot" list

      waitress$start()

      # check signature method interchangeability
      # when not interchangeable set signatureMethod = DeconvMethod
      signature_Method <- input$signatureMethod
      if (!(input$deconvMethod %in% methods_interchangeable)) {
        signature_Method <- input$deconvMethod
      }

      showNotification(paste0("Building Signature: ", signature_Method), type = "warning")

      # get signature or calculate new one
      # signature <- signature()
      signature <- omnideconv::build_model(
        single_cell_object = userData$singleCell,
        bulk_gene_expression = userData$bulk,
        method = signature_Method,
        batch_ids = userData$batchIDs,
        cell_type_annotations = userData$cellTypeAnnotations,
        markers = userData$marker
      )

      # deconvolute
      showNotification(paste0("Deconvolution started: ", input$deconvMethod), type = "warning")
      deconvolution_result <-
        omnideconv::deconvolute(
          bulk_gene_expression = userData$bulk,
          signature = signature,
          method = input$deconvMethod,
          single_cell_object = userData$singleCell,
          cell_type_annotations = userData$cellTypeAnnotations,
          batch_ids = userData$batchIDs
        )

      # insert result into the all_deconvolutions reactive Value
      all_deconvolutions[[paste0(input$deconvMethod, "_", signature_Method)]] <- list(deconvolution_result, signature)

      waitress$close()
      showNotification("Deconvolution finished", type = "message")
    }) # end deconvolution´

    # update Deconvolution Method and signature method choices when new deconvolution result is calculated
    observe({
      deconv_choices <- unlist(strsplit(names(all_deconvolutions), "_"))
      deconv_choices <- deconv_choices[seq(1, length(deconv_choices), 2)] # jede 2, startend von 1
      updateSelectInput(session, inputId = "computedDeconvMethod", choices = deconv_choices)
    })

    # update signature Method choices when selecting deconvolution to load
    observe({
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
      tmp <- userData$deconvolution_result
      userData$deconvolution_result <- c(tmp, getSelectionToPlot())

      # updateTableSelection()
    })

    # remove deconvolution from To Plot list
    observeEvent(input$removeFromPlot, {
      tmp <- userData$deconvolution_result
      tmp <- tmp[!tmp %in% getSelectionToPlot()]
      userData$deconvolution_result <- tmp

      # updateTableSelection()
    })

    # load Deconvolution result
    observeEvent(input$loadDeconvolution, {
      userData$deconvolution_result <- c(getSelectionToPlot())
      # userData$signature <- c(getSelectionToPlot())

      # updateTableSelection()
    })

    # observe the selection to plot and show buttons if conditions match
    observe({
      if (getSelectionToPlot() %in% userData$deconvolution_result) {
        shinyjs::hide("addToPlot")
        shinyjs::show("removeFromPlot")
      } else {
        shinyjs::hide("removeFromPlot")
        shinyjs::show("addToPlot")
      }
    })

    # update selection inputs if deconvolution result list gets changed
    observeEvent(userData$deconvolution_result, {
      updateSelectInput(session, inputId = "deconvolutionToTable", choices = userData$deconvolution_result)
      updateSelectInput(session, inputId = "signatureToTable", choices = userData$deconvolution_result)
    })

    # Plots -------------------------------------------------------------------

    output$plotBox <- plotly::renderPlotly(
      plot_deconvolution(
        userData$deconvolution_result,
        input$plotMethod,
        input$facets,
        all_deconvolutions
      )
    )

    output$benchmarkPlot <- plotly::renderPlotly(
      plot_benchmark(
        userData$deconvolution_result,
        all_deconvolutions
      )
    )


    # Tables ------------------------------------------------------------------


    output$tableBox <- DT::renderDataTable({

      # needs a deconvolution to be selected
      req(input$deconvolutionToTable)

      # load deconvolution
      deconvolution <- all_deconvolutions[[input$deconvolutionToTable]][[1]]

      # render table
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
      # TODO: update: work with a list of  sigantures from userData$signature

      # run only if variable contains correct signature
      req(
        !is.null(input$signatureToTable),
        input$signatureToTable != "autogenes_autogenes", # will store a link to tmp file
        input$signatureToTable != "scaden_scaden" # will store a link to tmp file
      )

      # load signature
      signature <- all_deconvolutions[[input$signatureToTable]][[2]]

      # render table
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
    

    # Images ------------------------------------------------------------------
    output$logo <- renderImage({
      list(src=system.file("www", "logo.jpg", package = "DeconvExplorer"),
           contentType="image/jpeg", 
           width="100%")
    }, deleteFile = TRUE
    )

    # functions ---------------------------------------------------------------

    getSelectionToPlot <- function() {
      return(paste0(input$computedDeconvMethod, "_", input$computedSignatureMethod))
    }

    # updateTableSelection <- function() {
    #   updateSelectInput(session, inputId = "deconvolutionToTable", choices = userData$deconvolution_result)
    #   updateSelectInput(session, inputId = "signatureToTable", choices = userData$deconvolution_result)
    # }

    # load user file, file information from fileInput()
    loadFile <- function(file, type = "") {
      # get file extension and path
      path <- file$datapath
      ext <- tools::file_ext(path)
      content <- NULL

      # load file, depending on extension
      if (ext == "txt") {
        content <- utils::read.table(path)
      } else if (ext == "csv") {
        content <- vroom::vroom(path, delim = ",")
      } else if (ext == "tsv") {
        content <- vroom::vroom(path, delim = "\t")
      } else {
        showNotification(paste("File extension ", ext, " not supported.
                               Please view documentation for further information."), type = "error")
      }

      # file is loaded, perform checks
      if (!is.null(content)) {
        # check data ... (content, structure, gene names, etc.)

        # convert to dataframe and set colnames
        content <- as.data.frame(content)
        rownames(content) <- content[, 1] # first column
        content[, 1] <- NULL # remove first column

        # if requested a vector turn into character vector
        if (type == "vector") {
          if (dim(content)[2] != 1) {
            message("Wrong file dimensions for file ", file$name)
          }
          # turn into vector
          content <- as.vector(t(content))
        }

        showNotification(paste("Successfully Loaded File: ", file$name), type = "default")
      }
      content # case NULL = File not loaded, error already displayed to user
    }
  })

  shiny::shinyApp(ui = de_ui, server = de_server)
}
