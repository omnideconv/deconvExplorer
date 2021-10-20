library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dplyr)

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


### the Layout
dashboardPage(
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