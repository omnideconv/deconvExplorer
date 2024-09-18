#' Run DeconvExplorer
#'
#' This function launches a Shiny app to facilitate cell type deconvolution using both bulk
#' and single-cell RNA sequencing data. It provides a comprehensive interface for data upload,
#' deconvolution execution, and result visualization. The app supports various deconvolution
#' methods and offers tools for signature matrix refinement.
#'
#' @param deconvexp_bulk Optional; a matrix or data frame containing bulk sequencing data to be deconvoluted.
#'   Rows should represent genes, and columns should represent samples. The data can also be
#'   uploaded directly in the app.
#'
#' @param deconvexp_singlecelldata Optional; a matrix, data frame, or SingleCellExperiment object
#'   containing single-cell data used to calculate the signature matrix. Rows should represent genes,
#'   and columns should represent single cells.
#'
#' @param deconvexp_cell_annotation Optional; a vector providing cell type annotations
#'   for the single-cell data. Each entry corresponds to the cell type of the respective column
#'   in `deconvexp_singlecelldata`.
#'
#' @param deconvexp_batch Optional; a vector indicating the batch ID for each sample or cell
#'   in `deconvexp_singlecelldata`. This is relevant for methods that can adjust for batch effects.
#'
#' @param maxsize_upload Numeric; specifies the maximum file size in MB acceptable for upload
#'   during runtime. This is particularly important when files are uploaded directly through the
#'   app interface. Defaults to 50 MB.
#'
#' @return Starts a shiny app
#'
#' @export
DeconvExplorer <- function(deconvexp_bulk = NULL,
                           deconvexp_singlecelldata = NULL,
                           deconvexp_cell_annotation = NULL,
                           deconvexp_batch = NULL,
                           maxsize_upload = 50) {
  # options management
  oopt <- options(
    spinner.type = 6,
    spinner.color = "#0092AC",
    shiny.maxRequestSize = maxsize_upload * 1024^2
  )
  # play nice with other previously chosen options
  on.exit(options(oopt))

  # methods that produce a signature
  produces_signature <- c(
    # "BSeq-sc" = "bseqsc", # markers!!!
    "CIBERSORTx" = "cibersortx",
    "DWLS" = "dwls",
    "MOMF" = "momf",
    "MuSiC" = "music" # basicly a one step method but allow calculating a signature
  )

  # methods that allow the input of a custom signature
  two_step_methods <- c(
    "CIBERSORTx" = "cibersortx",
    "DWLS" = "dwls",
    "MOMF" = "momf"
  )

  overlay_color <- "rgb(51, 62, 72, .5)"

  # Data Upload Boxes -------------------------------------------------------

  data_simbu_box <- shinydashboard::box(
    id = "tour_simbu",
    title = "SimBu", solidHeader = TRUE, status = "primary", width = 12,
    column(
      width = 8,
      fileInput("data_simbu_simulation", "SimBu Simulation", accept = c(".rds"))
    ),
    column(
      width = 4,
      helpText("Upload a SimBu simulation result as .rds to load a bulk sample with corresponding cell fractions")
    )
  )

  data_deconvolution <- shinydashboard::box(
    id = "tour_upload",
    title = span("Input files for Deconvolution", icon("question-circle"), id = "uploadDeconvolutionQ"),
    solidHeader = TRUE, status = "primary", width = 12,
    fileInput("userBulkUpload", "Upload Bulk Data"),
    div(style = "margin-top: -20px"),
    fileInput("userSingleCellUpload", "Upload Single Cell Data"),
    div(style = "margin-top: -20px"),
    fileInput("userAnnotationUpload", "Upload Cell Type Annotations"),
    div(style = "margin-top: -20px"),
    fileInput("userBatchUpload", "Upload Batch IDs"),
    div(style = "margin-top: -20px"),
    fileInput("userMarkerUpload", "Upload Marker Genes"),
    div(style = "margin-top: -20px"), collapsible = T, collapsed = T,
    shinyWidgets::actionBttn("selectDeconvolution", "Perform deconvolution", icon = icon("arrow-right"), color = "success", style = "simple")
  )

  deconvUploadPopover <-
    shinyBS::bsPopover(
      id = "uploadDeconvolutionQ",
      title = "Upload Data for deconvolution",
      content = "Provide bulk and single cell data as well as the respective cell type annotation. This data will be used as input for deconvolution. Some methods might require additional batch labels or marker genes.",
      trigger = "hover"
    )

  data_load_sample <- shinydashboard::box(
    id = "tour_sample",
    title = span("Load Example Data", icon("question-circle", id = "exampleDataQ")),
    solidHeader = TRUE, status = "primary", width = 12,
    column(
      width = 3,
      div(actionButton("loadSample", "Load Example Files"), style = "margin-top:0.5em")
    ),
    column(
      width = 8,
      helpText("Ground truth data will be loaded as 'Example Ground-truth'")
    )
  )

  exampleDataPopover <-
    shinyBS::bsPopover(
      id = "exampleDataQ",
      title = "Example Data",
      content = "Load a sample dataset that can be used to showcase deconvExplorers features."
    )

  data_load_signature <- shinydashboard::box(
    id = "tour_signatureUpload",
    title = span("Upload Signature", icon("question-circle", id = "uploadSignatureQ")),
    solidHeader = TRUE, status = "primary",
    width = 12,
    fileInput("userSignatureUpload", "Upload Signature"),
    div(style = "margin-top: -25px"),
    p("You can upload a previsouly generated signature matrix of a deconvolution method and analyse it with DeconvExplorer. Multiple uploads are possible."),
    fluidRow(
      column(4, shinyWidgets::actionBttn("selectSigExploration", "Explore the signature", icon = icon("arrow-right"), color = "success", style = "simple")),
      column(4, shinyWidgets::actionBttn("selectSigRefinement", "Refine the signature", icon = icon("arrow-right"), color = "success", style = "simple"))
    )
  )

  signatureUploadPopover <-
    shinyBS::bsPopover(
      id = "uploadSignatureQ",
      title = "Gene Expression Signature",
      content = "Upload a gene expression signature. The signature can further be analyzed in the Signature Exploration and Refinement modules or used as input in deconvolution."
    )

  data_load_fractions <- shinydashboard::box(
    title = span("Upload cell-type fractions", icon("question-circle", id = "uploadFractionsQ")),
    solidHeader = TRUE, status = "primary",
    width = 12,
    fileInput("userFractionsUpload", "Upload table with cell-type fractions"),
    div(style = "margin-top: -25px"),
    p("You can upload a table containing cell-type fractions, either coming from a deconvolution method or a ground-truth dataset with which you want to compare your deconvolution result. Multiple uploads are possible."),
    shinyWidgets::actionBttn("selectBenchmark", "Compare fractions", icon = icon("arrow-right"), color = "success", style = "simple")
  )

  fractionsUploadPopover <-
    shinyBS::bsPopover(
      id = "uploadFractionsQ",
      title = "Deconvolution results or Ground Truth",
      content = "Upload cell fractions or ground truth to be used in comparisons or benchmarking."
    )

  data_info <- shinydashboard::box(
    title = span(icon("info-circle"), "Input data formats and requirements"),
    solidHeader = FALSE, width = 12,
    collapsible = TRUE, collapsed = TRUE,
    includeMarkdown(system.file("extdata", "data_info.md", package = "DeconvExplorer"))
  )

  # Deconvolution Boxes -------------------------------------------------------
  data_upload_box <- shinydashboard::box(
    title = span("Select your Data", icon("question-circle", id = "deconvSelectDataQ")),
    status = "primary",
    solidHeader = TRUE, height = "30em", # collapsible = TRUE, # used to be 30em
    selectInput("bulkSelection", "Select a bulk dataset", choices = NULL),
    div(style = "margin-top: -10px"),
    selectInput("singleCellSelection", "Select a single cell dataset", choices = NULL),
    div(style = "margin-top: -10px"),
    selectInput("annotationSelection", "Select Cell Type annotations", choices = NULL),
    div(style = "margin-top: -10px"),
    selectInput("batchSelection", "Select Batch IDs", choices = NULL),
    div(style = "margin-top: -10px"),
    selectInput("markerSelection", "Select Marker Genes", choices = NULL)
  )

  deconvSelectDataPopover <-
    shinyBS::bsPopover(
      id = "deconvSelectDataQ",
      title = "",
      content = "Select the datasets you want to use as input for deconvolution. Upload your Data in the Data Upload module. Your first uploaded dataset will be selected automatically. "
    )

  settings_box <- shinydashboard::box(
    id = "tour_deconvSettings",
    title = span("Deconvolution Settings", icon("question-circle", id = "deconvSettingsQ")),
    status = "primary",
    solidHeader = TRUE, height = "30em", # collapsible = TRUE, # used to be 30em
    imageOutput("logo", height = "auto"),
    column(
      width = 12,
      h2("Robust deconvolution of cell types from any tissue", style = "text-align: center; font-weight: bold; color:#003F5C; margin-top: 0.5em; margin-bottom: 0.5em")
    ),
    column(
      width = 4,
      selectInput("deconvMethod", "Deconvolution Method",
        choices = c('MuSiC'='music', omnideconv::deconvolution_methods[-10])
      )
    ),
    column(
      width = 5,
      conditionalPanel(

        # all methods that take another signature as input
        condition = "input.deconvMethod == 'cibersortx' ||
                  input.deconvMethod == 'dwls' ||
                    input.deconvMethod == 'momf'",
        selectInput("signatureMethod", "Signature",
          choices = produces_signature
        )
      )
    ),
    column(
      width = 3,
      div(
        shinyBS::popify(shinyWidgets::actionBttn("deconvolute", "Deconvolute", style = 'simple', icon = icon('triangle-exclamation'), color = 'warning'),
                        "Attention", "Some methods are considerably slower than others; please keep this in mind when using DeconvExplorer for deconvolution."),
        style = "margin-top:1.7em"
      )
    ),
    waiter::useWaitress()
  )

  deconvSettingsPopover <-
    shinyBS::bsPopover(
      id = "deconvSettingsQ",
      title = "",
      content = "Select a deconvolution method to run. If required and supported by the deconvolution method you can additionally select a custom signature to be used in computation. Please note this is an advanced feature and should be used with caution. "
    )
  

  deconv_plot_box <- shinydashboard::box(
    id = "tour_deconvPlot",
    title = span("Deconvolution Plot ", icon("tasks", lib = "glyphicon"), icon("question-circle", id = "deconvPlotQ")),
    status = "warning", solidHeader = TRUE, width = 12,
    column(
      width = 3,
      selectInput("plotMethod", "Plot as: ",
        choices = c(
          "Bar Plot" = "bar", "Scatter" = "scatter",
          "Jitter Plot" = "jitter", "Box Plot" = "box",
          "Heatmap" = "heatmap"
        )
      )
    ),
    column(
      width = 3,
      selectInput("facets", "Group Plots By",
        choices = c(
          "Deconvolution Method" = "method",
          "Cell Type" = "cell_type", "Sample" = "sample"
        )
      )
    ),
    column(
      width = 12,
      withSpinner(
        plotly::plotlyOutput("plotBox", height = "500px")
      )
    )
  )

  deconvPlotPopover <-
    shinyBS::bsPopover(
      id = "deconvPlotQ",
      title = "",
      content = "Customize plot type and grouping of the results. You can change the selected deconvolution results above. "
    )

  deconv_table_box <- shinydashboard::box(
    title = span("Deconvolution Table ", icon("th", lib = "glyphicon")),
    status = "warning", solidHeader = TRUE, width = 12,
    column(
      width = 3,
      selectInput("deconvolutionToTable", "Deconvolution Result", choices = NULL)
    ),
    column(
      width = 3,
      div(downloadButton("deconvolutionDownload", "Download Deconvolution"),
        actionButton("deconvolutionToTableDelete", icon("trash")),
        style = "margin-top:1.9em"
      )
    ),
    column(
      width = 12,
      withSpinner(
        DT::dataTableOutput("tableBox")
      )
    )
  )

  deconv_signature_box <- shinydashboard::box(
    title = span("Deconvolution Signature ", icon("fingerprint")),
    status = "info", solidHeader = TRUE, width = 12,
    column(
      width = 3,
      selectInput("signatureToTable", "Signature", choices = NULL)
    ),
    column(
      width = 2,
      div(downloadButton("signatureDownload", "Download Signature"),
        actionButton("signatureToTableDelete", icon("trash")),
        style = "margin-top:1.9em"
      )
    ),
    column(
      width = 12,
      withSpinner(
        DT::dataTableOutput("signatureBox")
      )
    )
  )

  deconv_all_results <- shinydashboard::box(
    id = "tour_deconvPlotSettings",
    title = span("Results", icon("question-circle", id = "deconvResultQ")), status = NULL, solidHeader = FALSE, width = 12,
    column(
      width = 5,
      selectInput("deconvolutionToPlot", "Select Deconvolution results",
        choices = c(""), multiple = TRUE
      )
    ),
    column(
      width = 4,
      helpText("Select the deconvolution results to be plotted on the left side."),
      helpText("Deconvolution results get identified by the selected method and signature: ", shiny::tags$b("DeconvolutionMethod_Signature"))
    ),
    column(
      width = 2,
      selectInput("deconvolutionToDelete", "Delete a deconvolution result", choices = NULL)
    ),
    column(
      width = 1,
      div(
        actionButton("deconvolutionToDeleteButton", icon("trash")),
        style = "margin-top:1.7em"
      )
    )
  )

  deconvResultPopover <- shinyBS::bsPopover(
    id = "deconvResultsQ",
    title = "",
    content = "Select one or multiple deconvolution results to be plotted below. "
  )


  # Benchmarking Boxes ------------------------------------------------------
  benchmark_deconvolutionSelection <- shinydashboard::box(
    title = span("Benchmarking Selection", icon("question-circle", id = "benchSettingsQ")),
    status = "info", solidHeader = TRUE, width = 12,
    selectInput("benchmark_reference", "Reference", choices = NULL),
    selectInput("benchmark_ToPlot", "Select Deconvolution to benchmark", choices = NULL, multiple = TRUE)
  )

  benchSettingsPopover <-
    shinyBS::bsPopover(
      id = "benchSettingsQ",
      title = "",
      content = "Compare deconvolution results to a ground truth or between each other. The ground truth can be uploaded as cell fractions in the data module."
    )

  benchmark_plot_box <- shinydashboard::tabBox(
    title = "Benchmark",
    width = 12,
    tabPanel(
      "Scatter Plot",
      withSpinner(
        shiny::plotOutput("benchmark_scatter")
      )
    ),
    tabPanel(
      "Correlation",
      column(
        width = 2,
        selectInput("correlationPlotType", "Plot Type",
          choices = c(
            "Circle" = "circle",
            "Square" = "square",
            "Ellipse" = "ellipse",
            "Number" = "number",
            "Shade" = "shade",
            "Color" = "color",
            "Pie" = "pie"
          ),
          selected = "color"
        ),
      ),
      column(
        width = 2,
        selectInput("correlationAnnotationType", "P Value Annotation Type",
          choices = c(
            "None" = "n",
            "Value" = "p-value",
            "Significance" = "label_sig"
          ),
          selected = "label_sig"
        ),
      ),
      column(
        width = 2,
        selectInput("correlationAnntotationColor", "Annotation Color",
          choices = c(
            "Black" = "black",
            "White" = "white"
          ),
          selected = "white"
        ),
      ),
      withSpinner(
        shiny::plotOutput("benchmark_correlation")
      )
    ),
    tabPanel(
      "RMSE",
      column(
        width = 2,
        selectInput("rmsePlotType", "RMSE Plot Type",
          choices = c(
            "Heatmap" = "heatmap",
            "Boxplot" = "boxplot"
          ),
          selected = "heatmap"
        )
      ),
      column(
        width = 2,
        conditionalPanel(
          condition = "input.rmsePlotType == 'heatmap'",
          selectInput("rmseHeatmapMethod", "RMSE Heatmap Method",
            choices = c(
              "circle",
              "square",
              "ellipse",
              "number",
              "shade",
              "color",
              "pie"
            ),
            selected = "color"
          )
        ),
      ),
      withSpinner(
        shiny::plotOutput("benchmark_rmse")
      )
    )
  )

  # Signature Exploration Boxes ---------------------------------------------
  signature_genesPerMethod <- shinydashboard::box(
    id = "tour_genesPlot",
    title = span("Genes per Method", icon("question-circle", id = "sigGenesQ")),
    status = "info", solidHeader = TRUE,
    width = 4,
    withSpinner(
      plotOutput("signatureGenesPerMethod")
    ),
    downloadButton("downloadSignatureGenesPerMethod", label = "Download as PDF")
  )

  sigGenesPopover <-
    shinyBS::bsPopover(
      id = "sigGenesQ",
      title = "",
      content = "Compare the size of available signatures. The plot displays the number of genes contained in each signature. Smaller signatures result in faster computation time."
    )

  signature_kappaPerMethod <- shinydashboard::box(
    title = span("Condition Number per Method", icon("question-circle", id = "sigConditionNumberQ")),
    status = "info", solidHeader = TRUE,
    width = 4,
    withSpinner(
      plotOutput("kappaPerMethod")
    ),
    downloadButton("downloadKappaPerMethod", label = "Download as PDF")
  )

  sigConditionNumberPopover <-
    shinyBS::bsPopover(
      id = "sigConditionNumberQ",
      title = "",
      content = "The condition number quantifies how robust the signature is to noise in the input bulk data. In general, a lower number is considered to provide more stable deconvolution results."
    )

  signature_entropyPerMethod <- shinydashboard::box(
    title = span("Mean Entropy per Method", icon("question-circle", id = "sigEntropyQ")),
    status = "info", solidHeader = TRUE,
    width = 4,
    withSpinner(
      plotOutput("signatureEntropyPerMethod")
    ),
    downloadButton("downloadSignatureEntropyPerMethod", label = "Download as PDF")
  )

  sigEntropyPopover <-
    shinyBS::bsPopover(
      id = "sigEntropyQ",
      title = "",
      content = "The entropy quantifies how informative a gene's expression pattern is across cell-types. Lower values indicate a more cell-type specific expression pattern. Please note that very low values could also hint sequencing artifacts."
    )

  signature_clusteredHeatmap <- shinydashboard::box(
    title = span("Clustered Signature", icon("question-circle", id = "sigHeatmapQ")),
    status = "info", solidHeader = TRUE,
    width = 12,
    fluidRow(    
      column(
        width = 4,
        selectInput("signatureToHeatmap", "Select a Signature", choices = NULL)
      ),
      column(
        width = 2,
        selectInput("signatureAnnotationScore", "Select an annotation score",
                    choices = c("Entropy" = "entropy", "Gini Index" = "gini")
        )
      ),
      column(
        width = 2,
        selectInput("signatureAnnotationPlotType", "Annotation Type",
                    choices = c("Bars" = "bar", "Lines" = "line")
        )
      ),
      column(
        width = 2,
        selectInput("clusterCelltypes", "Order rows (cell types)",
                    choices = c(".. by cell-type similarity" = "cluster", ".. alphabetically" = "no_cluster")
        )
      ),
      column(
        width = 2,
        selectInput("clusterGenes", "Order columns (genes)",
                    choices = c(".. by maximal z-score per cell type" = "z-score cutoff",
                                ".. hierarchically based on euclidean distance" = "hierarchical clustering",
                                ".. alphabetically" = "alphabetical")
        )
      )
    ),
    fluidRow(
      column(
        width = 12,
        InteractiveComplexHeatmap::originalHeatmapOutput("clusteredHeatmapOneSignature",
                                                         width = "1250px", height = "450px", containment = TRUE
        )
      )
    ),
    fluidRow(
      column(
        width = 2,
        div(downloadButton("signatureSelectedGenesDownloadButton", "Download selected Genes"), style = "margin-top:1.9em")
      )
    )
  )

  sigHeatmapPopover <-
    shinyBS::bsPopover(
      id = "sigHeatmapQ",
      title = "",
      content = "This heatmap displays the z-scored expression profile of all genes in the signature. By selecting a subsection of the plot a more in-depth analysis is possible in the next plot below."
    )

  signature_clusteredHeatmapSubPlot <- shinydashboard::box(
    title = "Sub Selection Heatmap", status = "info", solidHeader = TRUE,
    width = 12, collapsible = TRUE, collapsed = TRUE,
    InteractiveComplexHeatmap::subHeatmapOutput("clusteredHeatmapOneSignature", width = "1250px"),
    conditionalPanel(
      condition = "false",
      InteractiveComplexHeatmap::HeatmapInfoOutput("clusteredHeatmapOneSignature")
    ) # necessary, will not display if function not used
  )

  signature_clusteredHeatmapSubTable <- shinydashboard::box(
    title = "Sub Selection Table", status = "info", solidHeader = TRUE,
    width = 12, collapsible = TRUE, collapsed = TRUE,
    withSpinner(
      DT::dataTableOutput("signatureHeatmap_SelectedGenesTable")
    )
  )

  signature_upsetPlot <- shinydashboard::box(
    title = span("UpSet Plot", icon("question-circle", id = "sigUpsetQ")),
    status = "info", solidHeader = TRUE, width = 8, height = "33em",
    withSpinner(
      plotOutput("signatureUpset")
    )
  )

  sigUpsetPopover <-
    shinyBS::bsPopover(
      id = "sigUpsetQ",
      title = "",
      content = "Compare and analyze the available genes and composition from all provided signatures."
    )

  signature_upsetPlotSettings <- shinydashboard::box(
    title = "UpSet Plot Settings", status = "info", solidHeader = TRUE,
    width = 4, height = "33em",
    column(
      width = 11,
      selectInput("upsetMode", "Upset Plot Mode", choices = c(
        "Distinct" = "distinct",
        "Intersect" = "intersect",
        "Union" = "union"
      ))
    ),
    column(
      width = 1,
      # link to help
      tags$a(
        href = "https://jokergoo.github.io/ComplexHeatmap-reference/book/08-upset_files/figure-html/unnamed-chunk-7-1.png",
        target = "_blank", icon("question-circle")
      ),
      style = "margin-top:2em"
    ),

    # plot settings
    sliderInput("upSetDegree",
      label = "Intersection Sizes to display", min = 1, max = 5,
      value = c(1, 5), round = TRUE, step = 1, ticks = FALSE
    ),
    column(
      width = 5,
      selectInput("upSetOrder", label = "Order Sets by", choices = c(
        "Set Size" = "size",
        "Set Degree" = "degree"
      ))
    ),
    column(
      width = 3,
      div(checkboxInput("upSetInvert", label = "Invert Order", value = FALSE), style = "margin-top:2em")
    ),
    column(
      width = 4,
      div(checkboxInput("upSetColorDegrees", label = "Color Degrees", value = TRUE), style = "margin-top:2em")
    ),

    # download of results
    checkboxGroupInput("upSetDownloadSelection", h3("Download Genes of a specific subset"),
      choices = NULL, inline = TRUE
    ),
    downloadButton("upSetDownloadButton", label = "Download Subset Genes"),
    downloadButton("upSetPlotDownloadButton", label = "Download plot as PDF")
  )

  # Signature Refinement Boxes ----------------------------------------------
  refinementHeatmapBox <- shinydashboard::box(
    id = "tour_refinementHeatmap",
    title = "Signature", solidHeader = TRUE, width = 12, status = "info",
    column(
      width = 2,
      selectInput("refinementHeatmapScore", "Gene Score", choices = c("Entropy" = "entropy", "Gini Index" = "gini"))
    ),
    column(
      width = 2,
      selectInput("refinementHeatmapScorePlotType", "Score Plot Type", choices = c("Bars" = "bar", "Line" = "line"))
    ),
    shinydashboard::valueBoxOutput("refinementGenes", width = 2),
    shinydashboard::valueBoxOutput("refinementCellTypes", width = 2),
    shinydashboard::valueBoxOutput("refinementKappa", width = 2),
    shinydashboard::valueBoxOutput("refinementMeanEntropy", width = 2),
    column(
      width = 12,
      withSpinner(
        plotOutput("refinementHeatmapPlot")
      )
    )
  )

  refinementSettingsBox <- shinydashboard::box(
    title = span("Settings", icon("question-circle", id = "refSettingsQ")),
    solidHeader = TRUE, width = 4, status = "info",
    column(
      width = 8,
      selectInput("signatureToRefine", "Choose a signature to refine", choices = NULL)
    ),
    column(
      width = 4,
      actionButton("loadRefinementSignature", "Load", style = "margin-top: 1.7em")
    ),
    column(
      width = 8,
      textInput("refinementNewName", "New Signature Name")
    ),
    column(
      width = 4,
      actionButton("saveRefinedSignature", "Save", style = "margin-top: 1.7em")
    ),
    column(
      width = 12,
      helpText(icon("arrow-up"), "  Load a signature to run refinements and choose a new name when saving"),
      div(helpText(icon("arrow-down"), "  If required rename specific cell types"), style = "margin-top:2.5em")
    ),
    column(
      width = 8,
      selectInput("cellTypeToRename", "Cell Type to Rename", choices = NULL),
      textInput("cellTypeNewName", "New cell type name")
    ),
    column(
      width = 4,
      div(actionButton("renameCellTypeGo", "Rename"), style = "margin-top:4.5em")
    )
  )

  refSettingsPopover <-
    shinyBS::bsPopover(
      id = "refSettingsQ",
      title = "",
      content = "Load a signature for modifications. A new name is required for saving. If necessary, cell-types can be renamed."
    )

  refinementUnzeroBox <- shinydashboard::box(
    solidHeader = FALSE, width = NULL, background = "aqua",
    fluidRow(
      column(
        width = 4,
        h1("Unzero")
      ),
      column(
        width = 7,
        sliderInput("refinePercentZero", "Maximum percentage of zeroes allowed for each gene",
                    min = 0, max = 100, value = 90, step = 1, post = "%"
        )
      ),
      column(
        width = 1,
        actionButton("refinePercentZeroGo", "Run", style = "margin-top: 1.7em")
      )
    ),
    fluidRow(
      column(10, p("Remove genes that contain mostly 0 in their expression profile. All genes with more zeros than the provided percentage are removed. "))
    )
  )

  refinementRemoveUnspecificBox <- shinydashboard::box(
    solidHeader = FALSE, width = NULL, background = "yellow",
    fluidRow(
      column(
        width = 4,
        h1("Remove Unspecific")
      ),
      column(
        width = 7,
        numericInput("refineUnspecific", "Remove unspecific genes", 1)
      ),
      column(
        width = 1,
        actionButton("refineUnspecificGo", "Run", style = "margin-top: 1.7em")
      )
    ),
    fluidRow(
      column(10, p("The expression of each gene over all cell types is binned into 'high', 'medium' and 'low'. A gene is considered 'specific' if its expression is 'high' in no more than n cell types; a gene with a 'high' expression in more than n cell types is removed. The 'n' parameter can be modified in the field on the side."))
    )
  )

  refinementBestNBox <- shinydashboard::box(
    solidHeader = FALSE, width = NULL, background = "red",
    fluidRow(
      column(
        width = 4,
        h1("Best n genes")
      ),
      column(
        width = 5,
        numericInput("refineBestN", "Number of genes to select for each cell type", 20, 1)
      ),
      column(
        width = 2,
        selectInput("refineBestNScore", "How to score genes", choices = c("Entropy" = "entropy", "Gini Index" = "gini"))
      ),
      column(
        width = 1,
        actionButton("refineBestNGo", "Run", style = "margin-top: 1.7em")
      )
    ),
    fluidRow(
      column(10, p("Based on the selected scoring method, the top ranking genes are selected for each cell type. "))
    )
  )


  refinementManualBox <- shinydashboard::box(
    solidHeader = FALSE, width = NULL, background = "purple",
    fluidRow(
      column(
        width = 4,
        h1("Remove manually")
      ),
      column(
        width = 7,
        textInput("refinementManualGene", "Type in a Gene Identifier to remove")
      ),
      column(
        width = 1,
        actionButton("refinementManualGo", "Run", style = "margin-top: 1.7em")
      )
    ),
    fluidRow(
      column(10, p("Manually remove genes by providing a Gene Identifier. "))
    )
  )

  refManuallyPopover <-
    shinyBS::bsPopover(
      id = "refManuallyQ",
      title = "",
      content = 
    )

  # Info Boxes --------------------------------------------------------------
  info_overview <- shinydashboard::box(
    id = "tour_infoOverview",
    title = "Overview", solidHeader = TRUE,
    h3("DeconvExplorer is an interactive web interface to perform, evaluate and enhance cell type deconvolution
      from trancsriptome data with the omnideconv framework.", style = "margin-top:-10px;"), br(),
    p("The app contains multiple modules further explained below."), br(),
    column(
      width = 12,
      div(
        h3(icon("github"), style = "display:inline; margin-right:1em"),
        a("Omnideconv", href = "https://github.com/omnideconv", target = "_blank", style = "margin-right:1em"),
        h3(icon("link"), style = "display: inline; margin-right:1em"),
        a("omnideconv.org", href = "https://www.omnideconv.org", target = "_blank"),
        style = "display:inline; font-size: 1.4em;"
      ),
    ),
    column(
      width = 12,
      div(
        h3(icon("envelope"), style = "display: inline; margin-right:1em"),
        a("Francesca Finotello", href = "mailto:francesca.finotello@uibk.ac.at", style = "margin-right:1em"),
        a("Markus List", href = "mailto:markus.list@wzw.tum.de", style = "margin-right:1em"),
        a("Federico Marini", href = "mailto:marinif@uni-mainz.de", style = "margin-right:1em"),
        a("Constantin Zackl", href = "mailto:Constantin.Zackl@student.uibk.ac.at", style = "margin-right:1em"),
        a("Lorenzo Meretto", href = "mailto:Lorenzo.Merotto@uibk.ac.at", style = "margin-right:1em"),
        a("Alexander Dietrich", href = "mailto:alex.dietrich@tum.de", style = "margin-right:1em"),
        style = "display:block; font-size:1.4em; margin-top:0.7em; "
      )
    ), br()
  )

  info_link <- shinydashboard::box(
    title = NULL, solidHeader = TRUE,
    column(
      width = 12,
      imageOutput("logoInfo", width = "100%", height = "auto")
    ),
    column(
      width = 12,
      h2("Robust deconvolution of cell types from any tissue", style = "text-align: center; font-weight: bold; color:#003F5C;")
    ),
  )

  info_modules <- shinydashboard::box(
    title = "Information about each module", solidHeader = TRUE,
    width = 12,
    status = "primary",
    div(includeMarkdown(system.file("extdata", "app_information.md", package = "DeconvExplorer")), style = "padding:1em; padding-top:0em")
  )


  # ui definition  ----------------------------------------------------------
  deconvexplorer_ui <- dashboardPage(
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
          icon = icon(NULL, verify_fa = FALSE)
        ),
        notificationItem(text = actionButton("githubLink", "View the Code",
          onclick = "window.open('https://github.com/omnideconv', '_blank')",
          icon = icon("github")
        ), icon = icon(NULL, verify_fa = FALSE))
      ),
      dropdownMenu(
        type = "task",
        icon = icon("cog", lib = "glyphicon"),
        headerText = "Set your CIBERSORTx Credentials",
        badgeStatus = NULL,
        notificationItem(
          text = textInput("csxEmail", "Email Adress"),
          icon = icon(NULL, verify_fa = FALSE)
        ),
        notificationItem(
          text = textInput("csxToken", "Token"),
          icon = icon(NULL, verify_fa = FALSE)
        ),
        notificationItem(
          text = actionButton("setCSX", "Set CIBERSORTx Credentials"),
          icon = icon(NULL, verify_fa = FALSE)
        )
      ),
      dropdownMenu(
        type = "task",
        icon = icon("bookmark", lib = "glyphicon"),
        headerText = "Use session to proceed your work later",
        badgeStatus = NULL,
        notificationItem(
          text = downloadButton("downloadSession", "Download Session"),
          icon = icon(NULL, verify_fa = FALSE)
        ),
        notificationItem(
          text = fileInput("uploadSession", "Upload Session File", accept = c(".rds")),
          icon = icon(NULL, verify_fa = FALSE), status = "primary"
        )
      )
    ),
    dashboardSidebar(sidebarMenu(
      id = "tabs",
      shinyjs::useShinyjs(),
      rintrojs::introjsUI(),
      waiter::use_waiter(),
      waiter::waiter_show_on_load(html = tagList(waiter::spin_rotating_plane(), "Loading necessary packages ...")),
      menuItem("Data Upload", tabName = "data"),
      menuItem("Deconvolution", tabName = "deconv"),
      menuItem("Signature Exploration", tabName = "signatureExploration"),
      menuItem("Signature Refinement", tabName = "signatureRefinement"),
      menuItem("Benchmark", tabName = "benchmark"),
      menuItem("Further Information", tabName = "fInfo"),
      selectInput("globalColor", "Select Plot Color Palette",
        choices = c("Set1", "Set2", "Set3", "Paired", "Dark2", "Spectral", "Accent"),
        selected = "Dark2"
      )
    )),
    dashboardBody(
      tags$head(tags$style(
        HTML("
        .wrapper {height: auto !important;
             position:relative; overflow-x:hidden; overflow-y:hidden}
        /* logo */
        .skin-blue .main-header .logo {
                              background-color: #11415d;
        }

        /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: #11415d;
        }

        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #11415d;
        }

        /* main sidebar */
        .skin-blue .main-sidebar {
                              background-color: #11415d;
        }

        /* active selected tab in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #73bfa7;
        }

        /* other links in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #3687ba;
                              color: #000000;
        }

        /* other links in the sidebarmenu when hovered */
         .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #afdca4;
         }

        /* toggle button when hovered  */
         .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #3687ba;
         }

        /* primary box header and border */
        .box.box-solid.box-primary>.box-header {
                              color:#ffffff;
                              background:#11415d
        }
        .box.box-solid.box-primary{
                            border-bottom-color:#11415d;
                            border-left-color:#11415d;
                            border-right-color:#11415d;
                            border-top-color:#11415d;
        }

        /* info box header and border */
        .box.box-solid.box-info>.box-header {
                              color:#ffffff;
                              background:#3687ba
        }
        .box.box-solid.box-info{
                            border-bottom-color:#3687ba;
                            border-left-color:#3687ba;
                            border-right-color:#3687ba;
                            border-top-color:#3687ba;
        }

        /* warning box header and border */
        .box.box-solid.box-warning>.box-header {
                              color:#ffffff;
                              background:#ee6d3d
        }
        .box.box-solid.box-warning{
                            border-bottom-color:#ee6d3d;
                            border-left-color:#ee6d3d;
                            border-right-color:#ee6d3d;
                            border-top-color:#ee6d3d;
        }

        .popover{
                            color:#000000
        }
        ")
      )),
      tabItems(
        tabItem(tabName = "data", fluidPage(
          fluidRow(
            column(
              width = 6,
              data_load_signature, signatureUploadPopover,
              data_load_fractions, fractionsUploadPopover,
              data_load_sample, exampleDataPopover,
              data_deconvolution, deconvUploadPopover
            ),
            column(
              width = 6,
              imageOutput("logoDeconvExplorer", height = "auto")
            )
          ),
          fluidRow(
            column(
              width = 8,
              offset = 2,
              data_info
            )
          )
        )),
        tabItem(tabName = "deconv", fluidPage(
          fluidRow(
            data_upload_box, deconvSelectDataPopover,
            settings_box, deconvSettingsPopover
          ),
          fluidRow(
            deconv_all_results, deconvResultPopover
          ),
          fluidRow(
            deconv_plot_box, deconvPlotPopover,
            deconv_table_box,
            deconv_signature_box
          )
        )),
        tabItem(tabName = "signatureExploration", fluidPage(
          fluidRow(
            signature_genesPerMethod, sigGenesPopover,
            signature_kappaPerMethod, sigConditionNumberPopover,
            signature_entropyPerMethod, sigEntropyPopover
          ),
          fluidRow(
            signature_clusteredHeatmap, sigHeatmapPopover
          ),
          fluidRow(
            signature_clusteredHeatmapSubPlot
          ),
          fluidRow(
            signature_clusteredHeatmapSubTable
          ),
          fluidRow(
            signature_upsetPlot, sigUpsetPopover,
            signature_upsetPlotSettings
          )
        )),
        tabItem(tabName = "signatureRefinement", fluidPage(
          fluidRow(
            refinementHeatmapBox
          ),
          fluidRow(
            column(
              width = 8,
              refinementUnzeroBox,
              refinementRemoveUnspecificBox,
              refinementBestNBox,
              refinementManualBox
            ),
            refinementSettingsBox, refSettingsPopover
          )
        )),
        tabItem(tabName = "benchmark", fluidPage(
          fluidRow(
            benchmark_deconvolutionSelection, benchSettingsPopover
          ),
          fluidRow(
            benchmark_plot_box
          )
        )),
        tabItem(tabName = "fInfo", fluidPage(
          fluidRow(
            info_overview,
            info_link
          ),
          fluidRow(
            info_modules
          )
        ))
      )
    )
  )

  # server definition  ------------------------------------------------------

  deconvexplorer_server <- shinyServer(function(input, output, session) {
    # nocov start
    waiter::waiter_hide()

    # General Setup -----------------------------------------------------------
    internal <- shiny::reactiveValues(
      # signatures = list("dwls" = readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))),
      # deconvolutions = list("dwls_dwls" = readRDS(system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer"))),
      signatures = list(),
      deconvolutions = list(),
      bulk = list(),
      singleCell = list(),
      annotation = list(),
      batch = list(),
      markers = list()
    ) # this is new

    # options
    options(shiny.maxRequestSize = 10 * 1024^2 * 100) # 1GB

    waitress <- waiter::Waitress$new("#deconvolute", infinite = TRUE)

    # Reactives ---------------------------------------------------------------
    # collect all available signature option (Calculate new one and already available)
    allSignatureOptions <- reactive({
      # collect all available precalculated signatures
      precalcSignatures <- NULL

      # add token to make clear that these represent precalculated signatures
      for (name in names(internal$signatures)) {
        token <- stringr::str_to_title(name)
        if (grepl("precalculated_", name)) {
          token <- stringr::str_split(token, "_")[[1]][2] # remove precalc and
          precalcSignatures[token] <- paste0("precalculated_", name)
        }
        precalcSignatures[token] <- paste0("precalculated_", name)
      }

      list(
        "Calculate New Signature" = produces_signature,
        "Available Signatures" = precalcSignatures
      )
    })

    # reactiveVal of refinable signature, separated from the rest of signatures
    signatureRefined <- reactiveVal("")

    # init for later
    signatureSelectedGenesDownloadContent <- reactiveVal("") # set empty reactiveVal

    # Observers and Eventhandling ---------------------------------------------

    # start the tour
    observeEvent(input$startTour, {
      tour_steps <- read.delim(
        system.file("extdata", "tour_intro.txt",
          package = "DeconvExplorer"
        ),
        sep = ";", stringsAsFactors = FALSE,
        row.names = NULL, quote = ""
      )
      introjs(session, options = list(
        steps = tour_steps,
        "nextLabel" = ">",
        "prevLabel" = "<",
        "skipLabel" = "X"
      ), events = list(onbeforechange = rintrojs::readCallback("switchTabs")))
    })

    observe({
      updateSelectInput(session, "bulkSelection", choices = names(internal$bulk))
      updateSelectInput(session, "singleCellSelection", choices = names(internal$singleCell))
      updateSelectInput(session, "annotationSelection", choices = names(internal$annotation))
      updateSelectInput(session, "batchSelection", choices = names(internal$batch))
      updateSelectInput(session, "markerSelection", choices = names(internal$markers))
    })

    observe({
      updateSelectInput(session, "benchmark_reference", choices = names(internal$deconvolutions))
      updateSelectInput(session, "benchmark_ToPlot", choices = names(internal$deconvolutions))
    })

    observeEvent(input$selectDeconvolution, {
      updateTabItems(session, inputId = "tabs", selected = "deconv")
    })

    observeEvent(input$selectSigExploration, {
      updateTabItems(session, inputId = "tabs", selected = "signatureExploration")
    })

    observeEvent(input$selectSigRefinement, {
      updateTabItems(session, inputId = "tabs", selected = "signatureRefinement")
    })

    observeEvent(input$selectBenchmark, {
      updateTabItems(session, inputId = "tabs", selected = "benchmark")
    })

    observeEvent(input$loadSample, {
      waiter::waiter_show(html = tagList(waiter::spin_rotating_plane(), "Loading example data ..."), color = overlay_color)
      internal$bulk[["Example Bulk"]] <- omnideconv::bulk
      internal$singleCell[["Example Single-cell"]] <- omnideconv::single_cell_data_1
      internal$annotation[["Example Cell-type annotation"]] <- omnideconv::cell_type_annotations_1
      internal$batch[["Example Batch-IDs"]] <- omnideconv::batch_ids_1
      internal$deconvolutions[["Example Ground-truth"]] <- omnideconv::RefData
      internal$signatures[["Example Signature (DWLS)"]] <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))

      showNotification("Loaded Sample Data")
      waiter::waiter_hide()
    })

    # update Signature Select Options
    observeEvent(allSignatureOptions(), {
      updateSelectInput(session, "signatureMethod", choices = allSignatureOptions())
    })

    # for selecting a signature to refine, users should be able to select from all already available signatures
    observe({
      updateSelectInput(session, "signatureToRefine", choices = names(internal$signatures))
    })

    # when "load Refinement" is clicked, load signature in reactive Value
    observeEvent(input$loadRefinementSignature, {
      req(input$signatureToRefine)
      showNotification(paste0("Loading Signature for Refinement: ", input$signatureToRefine))
      signatureRefined(internal$signatures[[input$signatureToRefine]])
    })

    # run signature refinement "unzero"
    observeEvent(input$refinePercentZeroGo, {
      req(signatureRefined(), input$refinePercentZero)

      sig_pre <- signatureRefined()
      percentage <- input$refinePercentZero / 100
      shiny::showNotification(paste0("Removing genes with more than ", percentage * 100, "% zeroes in a row"))

      signatureRefined(removePercentZeros(signatureRefined(), input$refinePercentZero / 100)) # update reactive Value with result
      sig_post <- signatureRefined()

      shiny::showNotification(
        paste0("Removed a total of ", nrow(sig_pre) - nrow(sig_post), " genes")
      )
    })

    # run signature refinement "unspecific"
    observeEvent(input$refineUnspecificGo, {
      req(signatureRefined(), input$refineUnspecific)

      sig_pre <- signatureRefined()
      shiny::showNotification("removing unspecific genes from signature")

      signatureRefined(
        removeUnspecificGenes(signatureRefined(),
          number_of_bins = 3,
          max_count = input$refineUnspecific
        )
      )

      sig_post <- signatureRefined()
      shiny::showNotification(
        paste(
          "Removed", nrow(sig_pre) - nrow(sig_post),
          "unspecific genes from the signature."
        )
      )
    })

    # run signature refinement "bestN"
    observeEvent(input$refineBestNGo, {
      req(signatureRefined(), input$refineBestN)

      sig_pre <- signatureRefined()
      ref_method <- "entropy" ## hard coded for now, could be a widget itself?
      shiny::showNotification(paste0("Refining Signature by score: ", ref_method))

      signatureRefined(
        selectGenesByScore(signatureRefined(),
          scoring_method = ref_method,
          genes_per_cell_type = input$refineBestN
        )
      )

      sig_post <- signatureRefined()
      shiny::showNotification(
        paste0("Removed a total of ", nrow(sig_pre) - nrow(sig_post), " genes")
      )
    })

    # run signature refinement "manual"
    observeEvent(input$refinementManualGo, {
      if (input$refinementManualGene == "") {
        showNotification("Please provide a Gene Identifier", type = "warning")
      }

      req(input$refinementManualGene, signatureRefined())

      genes <- rownames(signatureRefined())

      if (!(input$refinementManualGene %in% genes)) {
        showNotification("Gene not in Signature!", type = "error")
      } else {
        # remove unwanted gene from gene list
        genes <- genes[!genes %in% input$refinementManualGene]

        # subselect Signature and save to reactiveVal
        signatureRefined(signatureRefined()[genes, ])

        showNotification(
          paste0("Removed Gene ", input$refinementManualGene, " from the signature"),
          type = "message"
        )
      }
    })

    # save refinedSignature
    observeEvent(input$saveRefinedSignature, {
      if (is.null(input$refinementNewName) | input$refinementNewName == "") {
        showNotification("Please provide a new signature name!", type = "error")
      }

      req(signatureRefined(), input$refinementNewName)

      internal$signatures[[input$refinementNewName]] <- isolate(signatureRefined())

      showNotification(paste0("Successfully saved signature ", input$refinementNewName), type = "message")
    })

    # load cell types of currently loaded signature
    observe({
      updateSelectInput(session, "cellTypeToRename", choices = colnames(signatureRefined()))
    })

    # rename cell type if button is clicked
    observeEvent(input$renameCellTypeGo, {
      req(input$cellTypeNewName, signatureRefined(), input$cellTypeToRename)

      if (input$cellTypeToRename %in% colnames(signatureRefined())) {
        signatureRefined(renameCellType(isolate(signatureRefined()), input$cellTypeToRename, input$cellTypeNewName))
        showNotification("Renamed cell type", type = "message")
      } else {
        showNotification("Cell Type does not exist in signature", type = "error")
      }
    })

    # delete signatures
    observeEvent(input$signatureToTableDelete, {
      req(input$signatureToTable)
      internal$signatures[[input$signatureToTable]] <- NULL
      showNotification("Deleted Signature")
    })

    # delete deconvolution results
    observeEvent(input$deconvolutionToTableDelete, {
      req(input$deconvolutionToTable)
      internal$deconvolutions[[input$deconvolutionToTable]] <- NULL

      # check if plot currently loaded, if yes, update variable
      # userData$deconvolution_result = userData$deconvolution_result[!userData$deconvolution_result %in% input$deconvolutionToTable]

      showNotification("Deleted Deconvolution Result")
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

      nDeconvolutions <- length(session_deconvolutions)
      nSignatures <- length(session_signatures)

      # works, i checked that
      for (name in names(session_deconvolutions)) {
        internal$deconvolutions[[name]] <- session_deconvolutions[[name]]
      }

      for (name in names(session_signatures)) {
        internal$signatures[[name]] <- session_signatures[[name]]
      }

      showNotification(paste0("Loaded ", nDeconvolutions, " deconvolutions and ", nSignatures, " signatures"))
    })

    # deconvolute when button is clicked
    observeEvent(input$deconvolute, {
      # reqs

      waiter::waiter_show(html = tagList(waiter::spin_rotating_plane(), "Starting deconvolution ..."), color = overlay_color)

      bulkData <- NULL
      singleCellData <- NULL
      cellTypeAnnotations <- NULL
      batchIDs <- NULL
      markers <- NULL

      # check signature method interchangeability
      # when not interchangeable set signatureMethod = DeconvMethod
      signature_Method <- input$signatureMethod
      if (!(input$deconvMethod %in% two_step_methods)) {
        signature_Method <- input$deconvMethod
      }

      # input$deconvMethod, signature_Method#
      if (is.null(input$bulkSelection) | input$bulkSelection == "") {
        waiter::waiter_hide()
        showNotification("Bulk Data Missing", type = "error")
      }
      req(input$bulkSelection)
      bulkData <- internal$bulk[[input$bulkSelection]]

      # check if Single Cell Data Necessary
      if (input$deconvMethod %in% c("autogenes", "momf", "bisque", "music", "bseqsc", "cdseq", "cpm", "scdc", "scaden") | signature_Method %in% c("cibersortx", "dwls", "momf")) {
        if (is.null(input$singleCellSelection) | input$singleCellSelection == "") {
          waiter::waiter_hide()
          showNotification("Single Cell Data Missing", type = "error")
        }
        req(input$singleCellSelection)
        singleCellData <- internal$singleCell[[input$singleCellSelection]]
      }
      # check if annotation necessary
      if (input$deconvMethod %in% c("music", "bisque", "autogenes", "bseqsc", "cdseq", "cpm", "scdc", "scaden") | signature_Method %in% c("cibersortx", "dwls", "momf")) {
        if (is.null(input$annotationSelection) | input$annotationSelection == "") {
          waiter::waiter_hide()
          showNotification("Cell type annotation missing", type = "error")
        }
        req(input$annotationSelection)
        cellTypeAnnotations <- internal$annotation[[input$annotationSelection]]
      }

      # check if batch ids necessary
      if (input$deconvMethod %in% c("music", "bisque", "bseqsc", "cdseq", "scdc")) {
        if (is.null(input$batchSelection) | input$batchSelection == "") {
          waiter::waiter_hide()
          showNotification("BatchIDs Missing", type = "error")
        }
        req(input$batchSelection)
        batchIDs <- internal$batch[[input$batchSelection]]
      }

      if (input$deconvMethod %in% c("bseqsc")) {
        if (is.null(input$markerSelection) | input$markerSelection == "") {
          waiter::waiter_hide()
          showNotification("Markers Missing", type = "error")
        }
        req(input$markerSelection)
        markers <- internal$markers[[input$markerSelection]]
      }


      # check if signature needs to be calculated or loaded
      if (grepl("precalculated", signature_Method)) {
        # load signature
        token <- stringr::str_split(signature_Method, "_")[[1]][2] # get the signature name
        signature_Method <- token
        signature <- isolate(internal$signatures[[token]])
        showNotification(paste0("Using Available Signature ", signature_Method, " for deconvolution"))
      } else {
        # calculate signature from signature method
        waiter::waiter_update(html = tagList(waiter::spin_rotating_plane(), paste0("Building Signature: ", signature_Method)))

        tryCatch(
          {
            signature <- omnideconv::build_model(
              single_cell_object = singleCellData,
              bulk_gene_expression = bulkData,
              method = signature_Method,
              batch_ids = batchIDs,
              cell_type_annotations = cellTypeAnnotations,
              markers = markers,
              verbose = TRUE
            )

            # only add signature if not null
            if (!is.null(signature) && signature_Method != "autogenes" && signature_Method != "scaden") {
              internal$signatures[[signature_Method]] <- signature
            }

            message("Finished Signature") # debug reasons

            waiter::waiter_update(html = tagList(waiter::spin_rotating_plane(), paste0("Deconvolution started: ", input$deconvMethod)))

            deconvolution_result <-
              omnideconv::deconvolute(
                bulk_gene_expression = bulkData,
                signature = signature,
                method = input$deconvMethod,
                single_cell_object = singleCellData,
                cell_type_annotations = cellTypeAnnotations,
                batch_ids = batchIDs,
                verbose = TRUE
              )

            # insert result into the internal$deconvolutions reactive ValuelogoInfo
            internal$deconvolutions[[paste0(input$deconvMethod, "_", signature_Method)]] <- deconvolution_result

            waiter::waiter_hide()
            message("Finished Deconvolution") # debug reasons
          },
          error = function(e) {
            showModal(errorModal(error_message = e$message))
            waiter_hide()
          }
        )
      }
    })

    # update avaible deconvolutions for plotting
    observe({
      selection <- input$deconvolutionToPlot
      # selection <- intersect(selection, names(internal$deconvolutions)) # remove deleted ones
      updateSelectInput(session, inputId = "deconvolutionToPlot", choices = names(internal$deconvolutions), selected = selection)
    })

    observe({
      updateSelectInput(session, "deconvolutionToDelete", "Choose a deconvolution to delete", choices = names(internal$deconvolutions))
    })

    observeEvent(input$deconvolutionToDeleteButton, {
      req(input$deconvolutionToDelete)
      internal$deconvolutions[[input$deconvolutionToDelete]] <- NULL

      showNotification("Deleted Deconvolution", type = "message")
    })

    # update Signature Tab Choices when new Deconvolution Added
    observe({
      updateSelectInput(session, inputId = "signatureToHeatmap", choices = names(internal$signatures)) # used to be allSignatures()
    })

    # update selection inputs if deconvolution gets added
    observe({
      updateSelectInput(session, inputId = "deconvolutionToTable", choices = names(internal$deconvolutions))
      updateSelectInput(session, inputId = "signatureToTable", choices = names(internal$signatures))
    })

    # Plots -------------------------------------------------------------------
    output$plotBox <- plotly::renderPlotly({
      # req(userData$deconvolution_result)
      req(input$deconvolutionToPlot)
      suppressWarnings(omnideconv::plot_deconvolution(
        returnSelectedDeconvolutions(input$deconvolutionToPlot, isolate(internal$deconvolutions)),
        input$plotMethod,
        input$facets,
        input$globalColor
      ))
    })

    # barplots
    barplotReactive <- reactive({
      req(length(internal$signatures) > 0)
      signatures <- shiny::isolate(internal$signatures)
      nGenesPlot <- plot_signatureGenesPerMethod(signatures, input$globalColor)
      conditionNumberPlot <- plot_conditionNumberPerMethod(signatures, input$globalColor)
      entropyPlot <- plot_meanEntropyPerMethod(signatures, input$globalColor)

      return(list(
        nGenesPlot = nGenesPlot,
        conditionNumberPlot = conditionNumberPlot,
        entropyPlot = entropyPlot
      ))
    })

    # Number of genes Plot
    output$signatureGenesPerMethod <- renderPlot({
      req(barplotReactive)
      barplotReactive()$nGenesPlot
    })

    output$downloadSignatureGenesPerMethod <- downloadHandler(
      filename = function() {
        "signature_genes_plot.pdf"
      },
      content = function(file) {
        req(barplotReactive)
        ggsave(file, plot = barplotReactive()$nGenesPlot, device = "pdf", width = 6, height = 6)
      }
    )

    # Condition Number Plot
    output$kappaPerMethod <- renderPlot({
      req(barplotReactive)
      barplotReactive()$conditionNumberPlot
    })

    output$downloadKappaPerMethod <- downloadHandler(
      filename = function() {
        "condition_number_plot.pdf"
      },
      content = function(file) {
        req(barplotReactive)
        ggsave(file, plot = barplotReactive()$conditionNumberPlot, device = "pdf", width = 6, height = 6)
      }
    )

    # Entropy Plot
    output$signatureEntropyPerMethod <- renderPlot({
      req(barplotReactive)
      barplotReactive()$entropyPlot
    })

    output$downloadSignatureEntropyPerMethod <- downloadHandler(
      filename = function() {
        "condition_number_plot.pdf"
      },
      content = function(file) {
        req(barplotReactive)
        ggsave(file, plot = barplotReactive()$entropyPlot, device = "pdf", width = 6, height = 6)
      }
    )

    # plot interactive heatmap
    observe({
      req(
        input$signatureToHeatmap,
        input$signatureAnnotationScore,
        input$signatureAnnotationPlotType,
        input$clusterCelltypes,
        input$clusterGenes
      )
      signature <- isolate(internal$signatures[[input$signatureToHeatmap]])
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input,
        output,
        session,
        plot_signatureClustered(signature,
          scoring_method = input$signatureAnnotationScore,
          annotation_type = input$signatureAnnotationPlotType,
          color_palette = input$globalColor,
          order_rows = input$clusterCelltypes, 
          order_columns = input$clusterGenes
        ),
        "clusteredHeatmapOneSignature",
        brush_action = brush_action
      )
    })

    upsetReactive <- reactive({
      req(length(internal$signatures) > 0, input$upSetDegree, input$upSetOrder)

      # update checkbox of setting box before rendering the plot
      # needs to be done with every plot rerendering, data could have been changed!
      updateCheckboxGroupInput(session, "upSetDownloadSelection", choices = names(isolate(internal$signatures)), inline = TRUE)

      # get upset Degree Choices from slider Input
      minDegree <- input$upSetDegree[[1]]
      maxDegree <- input$upSetDegree[[2]]

      # calculate the plot
      result <- plot_signatureUpset(shiny::isolate(internal$signatures),
        upset_mode = input$upsetMode,
        min_degree = minDegree,
        max_degree = maxDegree,
        order_sets = input$upSetOrder,
        invert_sets = input$upSetInvert,
        color_by_degrees = input$upSetColorDegrees,
        color_palette = input$globalColor
      )

      # update settings
      # probably going with a preselected range of values
      # might also be possible to update to min=1, max=numberofsamples, and if maxDegree now higher than selected: value=c(min, newmax)
      # updateSliderInput(session, inputId = "upSetDegree", max=max(ComplexHeatmap::comb_degree(result[[2]])))

      # show the plot
      return(list("plot" = result[[1]]))
    })

    # UpSet Plot
    output$signatureUpset <- renderPlot({
      req(upsetReactive)
      upsetReactive()$plot
    })

    output$upSetPlotDownloadButton <- downloadHandler(
      filename = function() {
        "upset_signatures.pdf"
      },
      content = function(file) {
        req(upsetReactive)
        pdf(file = file, width = 9, height = 6)
        print(upsetReactive()$plot)
        dev.off()
        # ggsave(file, plot = upsetReactive()$plot, device = "pdf", width = 9, height = 6)
      }
    )

    output$refinementHeatmapPlot <- renderPlot({
      req(input$refinementHeatmapScore, input$refinementHeatmapScorePlotType, signatureRefined()) # und die signatur

      plot_signatureClustered(signatureRefined(),
        scoring_method = input$refinementHeatmapScore,
        annotation_type = input$refinementHeatmapScorePlotType,
        color_palette = input$globalColor
      )
    })

    output$benchmark_scatter <- renderPlot({
      req(input$benchmark_reference, input$benchmark_ToPlot)
      reference <- internal$deconvolutions[[input$benchmark_reference]]
      estimates <- returnSelectedDeconvolutions(input$benchmark_ToPlot, isolate(internal$deconvolutions))
      print(estimates)
      tryCatch(
        {
          plot_benchmark_scatter(reference, estimates, input$globalColor)
        },
        error = function(e) {
          print(e$message)
          showModal(errorModal(e$message))
        }
      )
    })

    output$benchmark_correlation <- renderPlot({
      req(input$benchmark_reference, input$benchmark_ToPlot, input$correlationPlotType, input$correlationAnnotationType, input$correlationAnntotationColor)
      reference <- internal$deconvolutions[[input$benchmark_reference]]
      estimates <- returnSelectedDeconvolutions(input$benchmark_ToPlot, isolate(internal$deconvolutions))
      tryCatch(
        {
          plot_benchmark_correlation(
            reference,
            estimates,
            plot_method = input$correlationPlotType,
            pvalue_type = input$correlationAnnotationType,
            pvalue_color = input$correlationAnntotationColor
          )
        },
        error = function(e) {
          print(e$message)
          showModal(errorModal(e$message))
        }
      )
    })

    output$benchmark_rmse <- renderPlot({
      req(input$benchmark_reference, input$benchmark_ToPlot, input$rmsePlotType, input$globalColor)
      reference <- internal$deconvolutions[[input$benchmark_reference]]
      estimates <- returnSelectedDeconvolutions(input$benchmark_ToPlot, isolate(internal$deconvolutions))

      tryCatch(
        {
          plot_benchmark_rmse(reference,
            estimates,
            plot_type = input$rmsePlotType,
            hm_method = input$rmseHeatmapMethod,
            color_palette = input$globalColor
          )
        },
        error = function(e) {
          print(e$message)
          showModal(errorModal(e$message))
        }
      )
    })

    # ValueBoxes --------------------------------------------------------------
    output$refinementGenes <- shinydashboard::renderValueBox({
      req(signatureRefined())

      shinydashboard::valueBox(
        value = nrow(signatureRefined()),
        subtitle = "Number of Genes",
        icon = icon("dna"),
        color = "blue"
      )
    })

    output$refinementCellTypes <- shinydashboard::renderValueBox({
      req(signatureRefined())

      shinydashboard::valueBox(
        value = ncol(signatureRefined()),
        subtitle = "Number of Cell Types",
        icon = icon("disease"),
        color = "light-blue"
      )
    })

    output$refinementKappa <- shinydashboard::renderValueBox({
      req(signatureRefined())
      kappa <- round(kappa(signatureRefined(), exact = TRUE), 2)
      shinydashboard::valueBox(
        value = kappa,
        subtitle = "Condition Number",
        icon = icon("hashtag"),
        color = "blue"
      )
    })

    output$refinementMeanEntropy <- shinydashboard::renderValueBox({
      req(signatureRefined())
      meanEntropy <- round(mean(BioQC::entropySpecificity(signatureRefined())), 2)
      shinydashboard::valueBox(
        value = meanEntropy,
        subtitle = "Mean Entropy",
        icon = icon("hashtag"),
        color = "light-blue"
      )
    })


    # Uploads -----------------------------------------------------------------
    # load simbu simulation from user and extract bulk and ground truth
    observeEvent(input$data_simbu_simulation, {
      # read rds
      file <- input$data_simbu_simulation$datapath
      name <- input$data_simbu_simulation$name

      simulation <- readRDS(file)

      tryCatch(
        {
          bulk <- as.matrix(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]])
          reference <- as.matrix(simulation$cell_fractions)

          internal$deconvolutions[[paste0("simbu_reference_", name)]] <- reference
          internal$bulk[[paste0("simbu_bulk_", name)]] <- bulk
          showNotification("Successfully loaded Simulation", type = "message")
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    observeEvent(input$userBulkUpload, {
      name <- input$userBulkUpload$name
      tryCatch(
        {
          internal$bulk[[paste0("bulk_", name)]] <- loadFile(input$userBulkUpload)
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    observeEvent(input$userSingleCellUpload, {
      name <- input$userSingleCellUpload$name
      tryCatch(
        {
          internal$singleCell[[paste0("singleCell_", name)]] <- loadFile(input$userSingleCellUpload)
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    observeEvent(input$userAnnotationUpload, {
      name <- input$userAnnotationUpload$name
      tryCatch(
        {
          internal$annotation[[paste0("annotation_", name)]] <- loadFile(input$userAnnotationUpload, type = "vector")
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    observeEvent(input$userBatchUpload, {
      name <- input$userAnnotationUpload$name
      tryCatch(
        {
          internal$batch[[paste0("batchID_", name)]] <- loadFile(input$userBatchUpload, type = "vector")
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    observeEvent(input$userMarkerUpload, {
      name <- input$userMarkerUpload$name
      tryCatch(
        {
          internal$markers[[paste0("batchID_", name)]] <- loadFile(input$userMarkerUpload, type = "vector")
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    observeEvent(input$userSignatureUpload, {
      name <- input$userSignatureUpload$name
      tryCatch(
        {
          internal$signatures[[paste0("Signature.", name)]] <- loadFile(input$userSignatureUpload)
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    observeEvent(input$userFractionsUpload, {
      name <- input$userFractionsUpload$name

      tryCatch(
        {
          internal$deconvolutions[[name]] <- loadFile(input$userFractionsUpload)
        },
        error = function(e) {
          showNotification("There was an error with your upload", type = "error")
        }
      )
    })

    # Tables ------------------------------------------------------------------
    output$tableBox <- DT::renderDataTable({
      req(
        input$deconvolutionToTable != "",
        internal$deconvolutions[[input$deconvolutionToTable]]
      )

      # load deconvolution
      deconvolution <- internal$deconvolutions[[input$deconvolutionToTable]]

      # turn rownames to column to enable DT search
      deconvolution <- data.frame("Sample" = rownames(deconvolution), deconvolution, check.names = FALSE) # check.names prevents cell type names from beeing changed
      rownames(deconvolution) <- NULL

      columns <- colnames(deconvolution)[-1]

      # render table
      DT::datatable(deconvolution,
        filter = "top",
        options = list(
          dom = "tip"
        )
      ) |>
        DT::formatPercentage(columns, 2)
    })

    output$signatureBox <- DT::renderDataTable({
      # run only if variable contains correct signature
      req(
        input$signatureToTable != "",
        input$signatureToTable != "autogenes",
        input$signatureToTable != "scaden",
        internal$signatures[[input$signatureToTable]]
      )

      # load signature
      signature <- isolate(internal$signatures[[input$signatureToTable]])

      # turn rownames to column to enable DT Search
      signature <- data.frame("Gene" = rownames(signature), signature, check.names = FALSE) # check.names prevents Cell Type names to be changed
      rownames(signature) <- NULL

      columns <- colnames(signature)[-1]

      # render table
      DT::datatable(signature, filter = "top", options = list(dom = "tip")) |>
        DT::formatRound(columns, 2)
    })


    # Downloads ---------------------------------------------------------------
    output$signatureDownload <- downloadHandler(
      filename = function() {
        paste("signature_", input$signatureToTable, ".csv", sep = "")
      },
      content = function(file) {
        # data <- deconv_list[[input$signatureToTable]][[2]]
        data <- isolate(internal$signatures[[input$signatureToTable]])
        write.csv(data, file)
      }
    )

    output$deconvolutionDownload <- downloadHandler(
      filename = function() {
        paste("deconvolution_", input$deconvolutionToTable, ".csv", sep = "")
      },
      content = function(file) {
        data <- internal$deconvolutions[[input$deconvolutionToTable]] # removed [[1]]
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
        data[["deconvolutions"]] <- shiny::isolate(internal$deconvolutions)
        data[["signatures"]] <- shiny::isolate(internal$signatures)

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
        signatures <- shiny::isolate(internal$signatures)

        data <- download_signatureUpset(signatures,
          combination_to_include = input$upSetDownloadSelection,
          upset_mode = input$upsetMode
        )

        # get genes from function
        write.table(data, file)
      }
    )

    # download selected Genes from the Interactive Signature Heatmap
    output$signatureSelectedGenesDownloadButton <- downloadHandler(
      filename = function() {
        paste0("Selection_", input$signatureToHeatmap, ".txt")
      },
      content = function(file) {
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
        content <- utils::read.table(path, header = T)
      } else if (ext == "csv") {
        content <- vroom::vroom(path, delim = ",")
      } else if (ext == "tsv") {
        content <- vroom::vroom(path, delim = "\t")
      } else if (ext == "rds" | ext == "RDS") {
        content <- readRDS(path)
      } else {
        showNotification(paste("File extension ", ext, " not supported.
                               Please view documentation for further information."), type = "error")
      }

      # file is loaded, perform checks
      if (!is.null(content)) {
        # check data ... (content, structure, gene names, etc.)

        # convert to dataframe and set colnames
        content <- as.data.frame(content)
        # content <- as.matrix(content)

        # if first column is a character column use as rownames
        if (typeof(content[, 1]) == "character" && dim(content)[2] > 1) {
          rownames(content) <- content[, 1] # first column
          content[, 1] <- NULL # remove first column
        }

        # if requested a vector turn into character vector
        if (type == "vector") {
          if (dim(content)[2] != 1) {
            message("Wrong file dimensions for file ", file$name)
          }
          # turn into vector
          content <- as.vector(t(content))
        }

        showNotification(paste("Successfully Loaded File: ", file$name), type = "default")
        print(is.data.frame(content))
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

    output$logoInfo <- renderImage(
      {
        list(
          src = system.file("www", "omnideconv_logo_info.svg", package = "DeconvExplorer"),
          contentType = "image/svg+xml",
          width = "100%"
        )
      },
      deleteFile = TRUE
    )

    output$logoDeconvExplorer <- renderImage(
      {
        list(
          src = system.file("www", "deconvExplorer.png", package = "DeconvExplorer"),
          contentType = "image/png",
          width = "100%"
        )
      },
      deleteFile = TRUE
    )


    # functions ---------------------------------------------------------------
    brush_action <- function(df, input, output, session) {
      req(internal$signatures, input$signatureToHeatmap) # used to contain deconv_list

      # ClusteredHeatmapSelectedGenes(Table)

      # get index of selected columns
      column_index <- unique(unlist(df$column_index))

      # get full dataset
      # signature <- allSignatures()[[input$signatureToHeatmap]]
      signature <- isolate(internal$signatures[[input$signatureToHeatmap]])

      # get selected subset
      selected <- signature[column_index, ]

      # Output Table of selected Genes
      output$signatureHeatmap_SelectedGenesTable <- DT::renderDataTable(DT::formatRound(DT::datatable(selected), columns = 1:ncol(selected), digits = 2))

      # Output List of Gene Names for Download
      signatureSelectedGenesDownloadContent(paste(rownames(selected), sep = "\n"))
    }
    # nocov end
  })

  shiny::shinyApp(ui = deconvexplorer_ui, server = deconvexplorer_server)
}
