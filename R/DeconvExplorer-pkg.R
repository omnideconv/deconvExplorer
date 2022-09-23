#' DeconvExplorer
#' 
#' DeconvExplorer Interactive user interface for the omnideconv deconvolution toolset
#' 
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @import omnideconv
#' @importFrom shinydashboard box dashboardBody dashboardHeader dashboardPage
#' dashboardSidebar dropdownMenu menuItem notificationItem sidebarMenu valueBox valueBoxOutput renderValueBox
#' tabItem tabItems
#' @importFrom plotly ggplotly plotlyOutput renderPlotly plot_ly layout config
#' @importFrom ggplot2 aes aes_ aes_string coord_cartesian coord_flip element_text
#' facet_wrap geom_abline geom_boxplot geom_col geom_jitter geom_point
#' geom_tile ggplot guide_colorbar guides labs scale_fill_gradient theme geom_text element_blank
#' geom_hline scale_colour_brewer scale_fill_brewer ylim
#' @importFrom shinycssloaders withSpinner 
#' @importFrom waiter Waitress
#' @importFrom rlang .data
#' @importFrom rintrojs introBox introjs introjsUI readCallback
#' @importFrom DT datatable dataTableOutput renderDataTable formatRound formatPercentage
#' @importFrom shinyjs useShinyjs hide show
#' @importFrom magrittr "%>%"
#' @importFrom utils write.csv write.table
#' @importFrom stringr str_subset
#' @importFrom tidyr pivot_longer
#' @importFrom stats sd cor.test
#' @importFrom bbplot bbc_style
#' @importFrom ComplexHeatmap Heatmap make_comb_mat UpSet comb_size upset_top_annotation extract_comb
#' @importFrom grid gpar unit
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput makeInteractiveComplexHeatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom ggforce facet_grid_paginate
#' @importFrom grDevices colorRampPalette
#' @importFrom ggpubr rotate_x_text stat_cor
#' @importFrom corrplot corrplot
#' @importFrom SummarizedExperiment assays
#' 
#' @name DeconvExplorer-pkg
#' @docType package
NULL
