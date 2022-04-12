
# library(shiny)
# library(shinydashboard)
# library(shinycssloaders)
# library(dplyr)
# library(ggplot2)
# library(omnideconv)
# library(RColorBrewer)
# library(waiter)
# library(rintrojs)
# 
# source("Global.R")

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
  #library(Biobase) # as long as music isn't fixed

  # methods that produce a signature
  produces_signature <- c(
    # "BSeq-sc" = "bseqsc", # markers!!!
    "CIBERSORTx" = "cibersortx",
   "DWLS" = "dwls", "MOMF" = "momf",
   "MuSiC" = "music" # basicly a one step method but allow calculating a signature
  )
  
  # methods that allow the input of a custom signature
  two_step_methods <- c(
    "CIBERSORTx" = "cibersortx",
    "DWLS" = "dwls", "MOMF" = "momf"
  )

  # box definitions ---------------------------------------------------------
  data_upload_box <- shinydashboard::box(
    title = "Upload your Data", status = "primary",
    solidHeader = TRUE, height = "34em", # collapsible = TRUE, # used to be 30em
    introBox(
      helpText("If no file is provided the analysis will be run with a sample dataset"),
      fileInput("userBulk", "Upload Bulk RNAseq Data"),
      div(style = "margin-top: -20px"),
      fileInput("userSingleCell", "Upload Single Cell RNASeq Data"),
      div(style = "margin-top: -20px"),
      fileInput("userCellTypeAnnotations", "Upload Cell Type Annotations"),
      div(style = "margin-top: -20px"),
      fileInput("userBatchIDs", "Upload Batch IDs"),
      div(style = "margin-top: -20px"), 
      fileInput("userSignature", "Upload your own Signature", multiple = TRUE),
      data.step = 1, data.intro = "Upload your Data. Allowed formats: txt, csv, tsv"
    )
  )

  settings_box <- shinydashboard::box(
    title = "Deconvolution Settings", status = "primary",
    solidHeader = TRUE, height = "34em", # collapsible = TRUE, # used to be 30em 
    introBox(
      imageOutput("logo", height = "auto"), br(),
      column(
        6,
        selectInput("deconvMethod", "Deconvolution Method",
          choices = omnideconv::deconvolution_methods
        )
      ),
      column(
        6,
        conditionalPanel(
          
          # all methods that take another signature as input
          condition = "input.deconvMethod == 'cibersortx' ||
                    input.deconvMethod == 'dwls' ||
                      input.deconvMethod == 'momf'",
          selectInput("signatureMethod", "Signature Calculation Method",
            choices = produces_signature
          )
        )
      ),
      column(
        6,
        conditionalPanel(
          condition = "input.deconvMethod == 'bseqsc'",
          fileInput("userMarker", "Marker Genes")
        )
      ),
      column(
        12,
        actionButton("deconvolute", "Deconvolute"),
        waiter::useWaitress(),
        actionButton("deconvoluteAll", "Deconvolute All")
      ),
      data.step = 2, data.intro = "Select preferred Deconvolution and Signature calculation Method"
    )
  )

  deconv_plot_box <- shinydashboard::box(
    title = span("Deconvolution Plot ", icon("tasks", lib = "glyphicon")),
    status = "warning", solidHeader = TRUE, width = 12,
    introBox(
      column(
        3,
        selectInput("plotMethod", "Plot as: ",
          choices = c(
            "Bar Plot" = "bar", "Scatter" = "scatter",
            "Jitter Plot" = "jitter", "Box Plot" = "box",
            "Heatmap" = "heatmap"
          )
        )
      ),
      column(
        3,
        selectInput("facets", "Group Plots By",
          choices = c(
            "Deconvolution Method" = "method",
            "Cell Type" = "cell_type", "Sample" = "sample"
          )
        )
      ),
      column(
        12,
        shinycssloaders::withSpinner(
          plotly::plotlyOutput("plotBox", height = "500px") # standard is 400
        )
      ),
      data.step = 4, data.intro = "View the deconvolution results and compare "
    )
  )

  deconv_table_box <- shinydashboard::box(
    title = span("Deconvolution Table ", icon("th", lib = "glyphicon")),
    status = "warning", solidHeader = TRUE, width = 12,
    column(3,
    selectInput("deconvolutionToTable", "Deconvolution Result", choices = NULL)),
    column(3, div(downloadButton("deconvolutionDownload", "Download Deconvolution"), style="margin-top:1.9em")), 
    column(12,
    shinycssloaders::withSpinner(
      DT::dataTableOutput("tableBox")
    ))
  )

  deconv_signature_box <- shinydashboard::box(
    title = span("Deconvolution Signature ", icon("fingerprint")),
    status = "info", solidHeader = TRUE, width = 12,
    column(3, 
    selectInput("signatureToTable", "Signature", choices = NULL)),
    column(3, 
    div(downloadButton("signatureDownload", "Download Signature"), style="margin-top:1.9em")),
    column(12, 
    shinycssloaders::withSpinner(
      DT::dataTableOutput("signatureBox")
    ))
  )
  deconv_all_results <- shinydashboard::box(
    title = "Plotting Settings", status = "info", solidHeader = TRUE, width = 12,
    introBox(
      column(3, selectInput("computedDeconvMethod", "Deconvolution Method", choices = NULL)),
      column(3, selectInput("computedSignatureMethod", "Signature Method", choices = NULL)),
      column(
        4,
        actionButton("loadDeconvolution", "Load Deconvolution Result", style = "margin-top: 1.7em"),
        actionButton("addToPlot", "Compare: Add to Plot", style = "margin-top: 1.7em"),
        actionButton("removeFromPlot", "Compare: Remove from Plot", style = "margin-top: 1.7em")
      ),
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


  # Signature Exploration Boxes ---------------------------------------------
  signature_genesPerMethod <- shinydashboard::box(
    title = "Genes per Method", status = "info", solidHeader = TRUE, width = 4,
    shinycssloaders::withSpinner(plotOutput("signatureGenesPerMethod"))
  )

  signature_kappaPerMethod <- shinydashboard::box(
    title = "Condition Number per Method", status = "info", solidHeader = TRUE, width = 4,
    shinycssloaders::withSpinner(plotOutput("kappaPerMethod"))
  )
  
  signature_entropyPerMethod <- shinydashboard::box(
    title = "Mean Entropy per Method", status = "info", solidHeader = TRUE, width = 4, 
    shinycssloaders::withSpinner(plotOutput("signatureEntropyPerMethod"))
  )

  signature_clusteredHeatmap <- shinydashboard::box(
    title = "Clustered Signature", status = "info", solidHeader = TRUE, width = 12,
    column(4,
    selectInput("signatureToHeatmap", "Select a Signature", choices = NULL)),
    column(2, selectInput("signatureAnnotationScore", "Select an annotation score", choices = c("Entropy" = "entropy", "Gini Index" = "gini"))),
    column(2, selectInput("signatureAnnotationPlotType", "Annotation Type", choices = c("Bars" = "bar", "Lines" = "line"))),
    column(4,
    div(downloadButton("signatureSelectedGenesDownloadButton", "Download selected Genes"),style="margin-top:1.9em")),
    column(12,
    InteractiveComplexHeatmap::originalHeatmapOutput("clusteredHeatmapOneSignature",
      width="1500px", height="450px", containment = TRUE
    ))
  )
  
  signature_clusteredHeatmapSubPlot <- shinydashboard::box(
    title = "Sub Heatmap", status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
    column(
      8,
      InteractiveComplexHeatmap::subHeatmapOutput("clusteredHeatmapOneSignature", width = "1000px")
    ),
    column(
      4, shinycssloaders::withSpinner(
        DT::dataTableOutput("signatureHeatmap_SelectedGenesTable")
      )
    ), 
    conditionalPanel(condition = "false",
                     InteractiveComplexHeatmap::HeatmapInfoOutput("clusteredHeatmapOneSignature")) # necessary, will not display if function not used
  )

  signature_upsetPlot <- shinydashboard::box(
    title = "UpSet Plot", status = "info", solidHeader = TRUE, width = 8, height = "33em",
    shinycssloaders::withSpinner(plotOutput("signatureUpset"))
  )
  signature_upsetPlotSettings <- shinydashboard::box(
    title = "UpSet Plot Settings", status = "info", solidHeader = TRUE, width = 4, height = "33em",
    column(
      11,
      selectInput("upsetMode", "Upset Plot Mode", choices = c(
        "Distinct" = "distinct",
        "Intersect" = "intersect",
        "Union" = "union"
      ))
    ),
    column(1,
      # link to help
      tags$a(href = "https://jokergoo.github.io/ComplexHeatmap-reference/book/08-upset_files/figure-html/unnamed-chunk-7-1.png", target = "_blank", icon("question-circle")),
      style = "margin-top:2em"
    ),

    # plot settings
    sliderInput("upSetDegree",
      label = "Intersection Sizes to display", min = 1, max = 5,
      value = c(1, 5), round = TRUE, step = 1, ticks = FALSE
    ),
    column(
      5,
      selectInput("upSetOrder", label = "Order Sets by", choices = c(
        "Set Size" = "size",
        "Set Degree" = "degree"
      ))
    ),
    column(
      3,
      div(checkboxInput("upSetInvert", label = "Invert Order", value = FALSE), style = "margin-top:2em")
    ),
    column(
      4,
      div(checkboxInput("upSetColorDegrees", label = "Color Degrees", value = TRUE), style = "margin-top:2em")
    ),

    # download of results
    checkboxGroupInput("upSetDownloadSelection", h3("Download Genes of a specific subset"),
      choices = NULL, inline = TRUE
    ),
    downloadButton("upSetDownloadButton", label = "Download Subset Genes")
  )
  

  # Signature Refinement Boxes ----------------------------------------------

  refinementHeatmapBox <- shinydashboard::box(
    title = "Signature", solidHeader = TRUE, width = 12, status = "info", 
    column(2, selectInput("refinementHeatmapScore", "Gene Score", choices = c("Entropy" = "entropy", "Gini Index" = "gini"))),
    column(2, selectInput("refinementHeatmapScorePlotType", "Score Plot Type", choices = c("Bars" = "bar", "Line" = "line"))),
    column(8, NULL),
    column(12, shinycssloaders::withSpinner(plotOutput("refinementHeatmapPlot")))
  )
  
  refinementSettingsBox <- shinydashboard::box(
    title = "Settings", solidHeader = TRUE, width = 4, status = "info",
    column(8, selectInput("signatureToRefine", "Choose a signature to refine", choices = NULL)),
    column(4, actionButton("loadRefinementSignature", "Load")), 
    # column with text/instructions
    column(8, textInput("refinementNewName", "New Signature Name")), 
    column(4, actionButton("saveRefinedSignature", "Save"))
  )
  
  refinementFunctionsBox <- shinydashboard::box(
    title = "Refine your signature", solidHeader = TRUE, width = 8, status = "info",
    column(8, NULL)
  )
  

  # ui definition  ----------------------------------------------------------


  de_ui <- dashboardPage(
    dashboardHeader(
      title = "DeconvExplorer",
      dropdownMenu(
        type = "task",
        icon = icon("question-circle"),
        headerText = "View a tour or look up the source code",
        badgeStatus = NULL,
        notificationItem(
          text = actionButton("startTour", "Start Tour",
            icon = icon("directions")
          ),
          icon = icon("info", verify_fa = FALSE)
        ),
        notificationItem(text = actionButton("githubLink", "View the Code",
          onclick = "window.open('https://github.com/omnideconv', '_blank')",
          icon = icon("github")
        ), icon = icon("info", verify_fa = FALSE))
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
          icon = icon("info", verify_fa = FALSE)
        ),
        notificationItem(
          text = fileInput("uploadSession", "Upload Session File"),
          icon = icon("info", verfiy_fa = FALSE), status = "primary"
        )
      )
    ),
    dashboardSidebar(sidebarMenu(
      shinyjs::useShinyjs(),
      introjsUI(),
      menuItem("Deconvolution", tabName = "deconv"),
      menuItem("Signature Exploration", tabName = "signatureExploration"),
      menuItem("Signature Refinement", tabName = "signatureRefinement"),
      menuItem("Benchmark", tabName = "benchmark"),
      menuItem("Further Information", tabName = "fInfo"), 
      selectInput("globalColor", "Select Plot Color Palette", 
                  choices = c("Set1", "Set2", "Set3", "Paired", "Dark2", "Spectral", "Accent"), 
                  selected="Spectral")
    )),
    dashboardBody(
      tags$head(tags$style(
        HTML(".wrapper {height: auto !important;
             position:relative; overflow-x:hidden; overflow-y:hidden}")
      )),
      tabItems(
        tabItem(tabName = "deconv", fluidPage(
          fluidRow(data_upload_box, settings_box),
          fluidRow(deconv_all_results),
          fluidRow(deconv_plot_box, deconv_table_box, deconv_signature_box)
        )),
        tabItem(tabName = "signatureExploration", fluidPage(
          fluidRow(signature_genesPerMethod, signature_kappaPerMethod, signature_entropyPerMethod),
          fluidRow(signature_clusteredHeatmap),
          fluidRow(signature_clusteredHeatmapSubPlot),
          fluidRow(signature_upsetPlot, signature_upsetPlotSettings)
        )),
        tabItem(tabName = "signatureRefinement", fluidPage(
          fluidRow(refinementHeatmapBox), 
          fluidRow(refinementFunctionsBox, refinementSettingsBox)
        )),
        tabItem(tabName = "benchmark", fluidPage(
          fluidRow(benchmark_plot_box)
        )),
        tabItem(tabName = "fInfo", fluidPage(
          includeMarkdown(
            system.file("www", "vignette.md", package = "DeconvExplorer")
          )
        ))
      )
    )
  )

  # server definition  ------------------------------------------------------

  de_server <- shinyServer(function(input, output, session) {
    ### background datastructure to store several deconvolution results

    # functions
    getSelectionToPlot <- function() {
      req(input$computedDeconvMethod != "", 
          input$computedSignatureMethod != "")
      return(paste0(input$computedDeconvMethod, "_", input$computedSignatureMethod))
    }


    # General Setup -----------------------------------------------------------

    
    # storing all calculated deconvolutions and signatures
    all_deconvolutions <- reactiveValues()
    all_signatures <- reactiveValues() 

    userData <- reactiveValues() # whatever this does

    # options
    options(shiny.maxRequestSize = 10 * 1024^2 * 100) # 1GB

    # for later: run all possible combinations! i and j cover all valid
    # deconvolution / signature combinations

    # for (i in methods){
    #   if (i %in% two_step_methods){
    #     for (j in produces_signature){
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

    
    
    # SAMPLE DATA
    
    userData$deconvolution_result <- c("momf_momf")
    
    all_deconvolutions[["momf_momf"]] <- readRDS(system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer"))
    all_signatures[["momf"]] <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))




    # Reactives ---------------------------------------------------------------

    # collect all available signature option (Calculate new one and already available)
    allSignatureOptions = reactive({
      # collect all available precalculated signatures
      precalcSignatures = NULL
      
      # add token to make clear that these represent precalculated signatures
      for (name in names(all_signatures)){
        token = stringr::str_to_title(name)
        if (grepl("precalculated_", name)){
          token  = stringr::str_split(token, "_")[[1]][2] # remove precalc and 
          precalcSignatures[token] = paste0("precalculated_", name)
        }
        precalcSignatures[token] <- paste0("precalculated_", name)
      }
      
      list("Calculate New Signature" = produces_signature,
           "Available Signatures" = precalcSignatures)
    })
    
    
    # reactiveVal of refinable signature, separated from the rest of signatures
    signatureRefined <- reactiveVal("")
    
    # init for later
    signatureSelectedGenesDownloadContent <- reactiveVal("") # set empty reactiveVal

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
    
    observeEvent(input$userSignature, {
      # checks and filename
      filename <- input$userSignature$name
      
      # add to all_signatures
      all_signatures[[filename]] <- loadFile(input$userSignature)
    })
    
    # update Signature Select Options
    observeEvent(allSignatureOptions(), {
      updateSelectInput(session, "signatureMethod", choices = allSignatureOptions())
    })
    
    # for selecting a siganture to refine, users should be able to select from all already available signatures
    observe({
      updateSelectInput(session, "signatureToRefine", choices = names(all_signatures))
    })
    
    # when "load Refinement" is clickes, load siganture in reactive Value
    observeEvent(input$loadRefinementSignature, {
      req(input$signatureToRefine)
      showNotification(paste0("Loading Signature for Refinement: ", input$signatureToRefine))
      signatureRefined(all_signatures[[input$signatureToRefine]])
    })
    

    # set CIBERSORTx Credentials from User Input
    observeEvent(input$setCSX, {
      req(input$csxEmail, input$csxToken)
      omnideconv::set_cibersortx_credentials(input$csxEmail, input$csxToken)
      showNotification("CIBERSORTx Credentials set")
    })

    # restore session with file upload
    observeEvent(input$uploadSession, {
      sessionFile <- readRDS(input$uploadSession$datapath)
      
      # separate signatures and deconvolutions from session file
      session_deconvolutions <- sessionFile[["deconvolutions"]] 
      session_signatures <- sessionFile[["signatures"]] 
      
      nDeconvolutions = length(session_deconvolutions)
      nSignatures = length(session_signatures)

      # works, i checked that
      for (name in names(session_deconvolutions)){
        all_deconvolutions[[name]] <- session_deconvolutions[[name]]
      }
      
      for (name in names(session_signatures)){
        all_signatures[[name]] <- session_signatures[[name]]
      }
      
      showNotification(paste0("Loaded ", nDeconvolutions, " deconvolutions and ", nSignatures, " signatures"))
    })

    # deconvolute when button is clicked
    observeEvent(input$deconvolute, {

      waitress$start()

      # check signature method interchangeability
      # when not interchangeable set signatureMethod = DeconvMethod
      signature_Method <- input$signatureMethod
      if (!(input$deconvMethod %in% two_step_methods)) {
        signature_Method <- input$deconvMethod
      }
      
      # check if signature needs to be calculated or loaded
      if (grepl("precalculated", signature_Method)){
        # load signature
        token = stringr::str_split(signature_Method, "_")[[1]][2] # get the signature name
        signature_Method = token
        signature <- all_signatures[[token]]
        showNotification(paste0("Using Available Signature ", signature_Method, " for deconvolution"))
        
      } else {
        # calculate signature from signature method
        showNotification(paste0("Building Signature: ", signature_Method), type = "warning")
        
        signature <- omnideconv::build_model(
          single_cell_object = userData$singleCell,
          bulk_gene_expression = userData$bulk,
          method = signature_Method,
          batch_ids = userData$batchIDs,
          cell_type_annotations = userData$cellTypeAnnotations,
          markers = userData$marker,
          verbose = TRUE
        )
      }
      
      # deconvolute
      showNotification(paste0("Deconvolution started: ", input$deconvMethod), type = "warning")
      deconvolution_result <-
        omnideconv::deconvolute(
          bulk_gene_expression = userData$bulk,
          signature = signature,
          method = input$deconvMethod,
          single_cell_object = userData$singleCell,
          cell_type_annotations = userData$cellTypeAnnotations,
          batch_ids = userData$batchIDs,
          verbose = TRUE
        )
    
      # insert result into the all_deconvolutions reactive Value
      all_deconvolutions[[paste0(input$deconvMethod, "_", signature_Method)]] <- deconvolution_result
      
      # only add signature if not null
      if (!is.null(signature)&& signature_Method != "autogenes" && signature_Method != "scaden"){
        all_signatures[[signature_Method]] <- signature
      } 
      
      waitress$close()
      showNotification("Deconvolution finished", type = "message")
      print ("Finished Deconvolution") # debug reasons
    }) 

    # update Deconvolution Method and signature method choices when new deconvolution result is calculated
    observe({
      deconv_choices <- unlist(strsplit(names(all_deconvolutions), "_"))
      deconv_choices <- deconv_choices[seq(1, length(deconv_choices), 2)] # jede 2, startend von 1
      updateSelectInput(session,
        inputId = "computedDeconvMethod",
        choices = deconv_choices
      )
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
      
      updateSelectInput(session,
        inputId = "computedSignatureMethod",
        choices = signature_choices
      )
    })

    # update Signature Tab Choices when new Deconvolution Added

    observe({
      updateSelectInput(session, inputId = "signatureToHeatmap", choices = names(all_signatures)) # used to be allSignatures()
    })

    # add Deconvolution to ToPlot list
    observeEvent(input$addToPlot, {
      tmp <- userData$deconvolution_result
      userData$deconvolution_result <- c(tmp, getSelectionToPlot())
    })

    # remove deconvolution from To Plot list
    observeEvent(input$removeFromPlot, {
      tmp <- userData$deconvolution_result
      tmp <- tmp[!tmp %in% getSelectionToPlot()]
      userData$deconvolution_result <- tmp
    })

    # load Deconvolution result
    observeEvent(input$loadDeconvolution, {
      userData$deconvolution_result <- c(getSelectionToPlot())
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

    # update selection inputs if deconvolution gets added
    observe({
      updateSelectInput(session, inputId = "deconvolutionToTable", choices = names(all_deconvolutions))
      updateSelectInput(session, inputId = "signatureToTable", choices = names(all_signatures))
    })

    # Plots -------------------------------------------------------------------

    output$plotBox <- plotly::renderPlotly({
      req(userData$deconvolution_result)
      omnideconv::plot_deconvolution(
        returnSelectedDeconvolutions(userData$deconvolution_result, shiny::reactiveValuesToList(all_deconvolutions)),
        input$plotMethod,
        input$facets,
        input$globalColor
      )
    })

    output$benchmarkPlot <- plotly::renderPlotly({
      plot_benchmark(returnSelectedDeconvolutions(userData$deconvolution_result, 
                                                  shiny::reactiveValuesToList(all_deconvolutions)))
    })

    # Number Of Genes Barplot
    output$signatureGenesPerMethod <- renderPlot({
      req(all_signatures)
      signatures <- shiny::reactiveValuesToList(all_signatures)
      plot_signatureGenesPerMethod(signatures, input$globalColor)
    })

    #Condition Number Plot
    output$kappaPerMethod <- renderPlot({
      req(all_signatures)
      signatures <- shiny::reactiveValuesToList(all_signatures)
      plot_conditionNumberPerMethod(signatures, input$globalColor)
    })
    
    output$signatureEntropyPerMethod <- renderPlot({
      req(all_signatures)
      signatures <- shiny::reactiveValuesToList(all_signatures)
      plot_meanEntropyPerMethod(signatures, input$globalColor)
    }
    )

    # plot interactive heatmap
    observe({
      req(input$signatureToHeatmap,
          input$signatureAnnotationScore,
          input$signatureAnnotationPlotType)
      signature <- all_signatures[[input$signatureToHeatmap]]
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input,
                                                               output,
                                                               session,
                                                               plot_signatureClustered(signature, 
                                                                                       score = input$signatureAnnotationScore, 
                                                                                       annotation_type = input$signatureAnnotationPlotType, 
                                                                                       palette = input$globalColor),
                                                               "clusteredHeatmapOneSignature",
                                                               brush_action = brush_action)
    })

    # UpSet Plot
    output$signatureUpset <- renderPlot({
      req(all_signatures, input$upSetDegree, input$upSetOrder)

      # update checkbox of setting box before rendering the plot
      # needs to be done with every plot rerendering, data could have been changed!
      updateCheckboxGroupInput(session, "upSetDownloadSelection", choices = names(all_signatures), inline = TRUE)

      # get upset Degree Choices from slider Input
      minDegree <- input$upSetDegree[[1]]
      maxDegree <- input$upSetDegree[[2]]

      # calculate the plot
      result <- plot_signatureUpset(shiny::reactiveValuesToList(all_signatures),
        mode = input$upsetMode,
        minDegree = minDegree,
        maxDegree = maxDegree,
        order = input$upSetOrder,
        invert = input$upSetInvert,
        colorDegrees = input$upSetColorDegrees,
        palette = input$globalColor
      )

      # update settings
      # probably going with a preselected range of values
      # might also be possible to update to min=1, max=numberofsamples, and if maxDegree now higher than selected: value=c(min, newmax)
      # updateSliderInput(session, inputId = "upSetDegree", max=max(ComplexHeatmap::comb_degree(result[[2]])))

      # show the plot
      result[[1]]
    })
    
    output$refinementHeatmapPlot <- renderPlot({
      req(input$refinementHeatmapScore, input$refinementHeatmapScorePlotType, signatureRefined()) # und die signature
      
      plot_signatureClustered(signatureRefined(), score=input$refinementHeatmapScore, 
                              annotation_type = input$refinementHeatmapScorePlotType, 
                              palette=input$globalColor)
      
      
      
    })

    # Tables ------------------------------------------------------------------


    output$tableBox <- DT::renderDataTable({
      req(input$deconvolutionToTable != "",
          all_deconvolutions[[input$deconvolutionToTable]])

      # load deconvolution
      deconvolution <- all_deconvolutions[[input$deconvolutionToTable]]
      
      # turn rownames to column to enable DT search
      deconvolution <- data.frame("Gene" = rownames(deconvolution), deconvolution, check.names = FALSE) # check.names prevents cell type names from beeing changed
      rownames(deconvolution) <- NULL
      
      columns <- colnames(deconvolution)[-1]

      # render table
      DT::datatable(deconvolution, filter="top",
        options = list(
          dom = "tip"
        )
      ) %>%
        DT::formatPercentage(columns, 2)
    })

    output$signatureBox <- DT::renderDataTable({
      # run only if variable contains correct signature
      req(
        input$signatureToTable != "",
        input$signatureToTable != "autogenes", 
        input$signatureToTable != "scaden",
        all_signatures[[input$signatureToTable]]
      )

      # load signature
      signature <- all_signatures[[input$signatureToTable]]
      
      # turn rownames to column to enable DT Search
      signature <- data.frame("Gene" = rownames(signature), signature, check.names = FALSE) # check.names prevents Cell Type names to be changed
      rownames(signature) <- NULL
      
      columns <- colnames(signature)[-1]
      
      # render table
      DT::datatable(signature, filter="top", options=list(dom="tip")) %>%
        DT::formatRound(columns, 2)
    })


    # Downloads ---------------------------------------------------------------

    output$signatureDownload <- downloadHandler(
      filename = function() {
        paste("signature_", input$signatureToTable ,".csv", sep = "")
      },
      content = function(file) {
        #data <- all_deconvolutions[[input$signatureToTable]][[2]]
        data <- all_signatures[[input$signatureToTable]]
        write.csv(data, file)
      }
    )
    
    output$deconvolutionDownload <- downloadHandler(
      filename=function(){
        paste("deconvolution_", input$deconvolutionToTable, ".csv", sep="")
      },
      content = function(file){
        data <- all_deconvolutions[[input$deconvolutionToTable]] # removed [[1]]
        write.csv(data, file)
      }
    )

    # save all deconvolutions and signatures to .RDS 
    output$downloadSession <- downloadHandler(
      filename = function() {
        paste0("omnideconv_", Sys.Date(), ".rds")
      },
      content = function(file) {
        data <- list()
        
        # save separate for later distinction
        data[["deconvolutions"]] <- shiny::reactiveValuesToList(all_deconvolutions)
        data[["signatures"]] <- shiny::reactiveValuesToList(all_signatures)
        
        # save data
        saveRDS(data, file)
      }
    )

    # TODO: UPDATE TO NEW DATATYPE ####
    output$upSetDownloadButton <- downloadHandler(
      filename = function() {
        paste0("subset_", paste0(input$upSetDownloadSelection, collapse = "_"), ".txt")
      },
      content = function(file) {
        # get subset selection from checkbox
        # Variable which contains the info: input$upSetDownloadSelection
        signatures <- shiny::reactiveValuesToList(all_signatures)
        
        data <- download_signatureUpset(signatures,
          combination = input$upSetDownloadSelection,
          mode = input$upsetMode
        )

        # get genes from function
        write.table(data, file)
      }
    )
    
    # download selected Genes from the Interactive Signature Heatmap
    output$signatureSelectedGenesDownloadButton <- downloadHandler(
      filename = function(){
        paste0("Selection_", input$signatureToHeatmap, ".txt")
      }, 
      content = function (file){
        data <- signatureSelectedGenesDownloadContent()

        # write file
        write.table(data, file)
      }
    )

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
        #content <- as.matrix(content)
        rownames(content) <- content[, 1] # first column
        content[, 1] <- NULL # remove first column

        # if requested a vector turn into character vector
        if (type == "vector") {
          if (dim(content)[2] != 1) {
            message("Wrong file dimensions for file ", file$name)
          }
          # turn into vector
          content <- as.vector(t(content))
        } else { # all other cases: matrix?
          content <- as.matrix(content)
        }

        showNotification(paste("Successfully Loaded File: ", file$name), type = "default")
      }
      content # case NULL = File not loaded, error already displayed to user
    }
    
    # Images ------------------------------------------------------------------
    output$logo <- renderImage(
      {
        list(
          src = system.file("www", "omnideconv_logo.svg", package = "DeconvExplorer"),
          contentType = "image/svg+xml",
          width = "100%"
        )
      },
      deleteFile = TRUE
    )
    
    # functions ---------------------------------------------------------------
    brush_action <- function(df, input, output, session) {
      req(all_signatures, input$signatureToHeatmap) # used to contain all_deconvolutions
      
      #ClusteredHeatmapSelectedGenes(Table)
      
      # get index of selected columns 
      column_index <- unique(unlist(df$column_index))
      
      # get full dataset
      #signature <- allSignatures()[[input$signatureToHeatmap]]
      signature <- all_signatures[[input$signatureToHeatmap]]
      
      # get selected subset
      selected <- signature[column_index,]   
      
      # Output Table of selected Genes
      output$signatureHeatmap_SelectedGenesTable <- DT::renderDataTable(DT::formatRound(DT::datatable(selected), columns=1:ncol(selected), digits=2))
      
      # Output List of Gene Names for Download
      signatureSelectedGenesDownloadContent(paste(rownames(selected), sep = "\n"))
    }
  })

  shiny::shinyApp(ui = de_ui, server = de_server)
}

#deconvExplorer()
