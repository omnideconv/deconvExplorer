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
#' geom_hline scale_colour_brewer scale_fill_brewer ylim theme_minimal geom_rect element_rect
#' @importFrom shinycssloaders withSpinner
#' @importFrom waiter Waitress waiter_hide waiter_show waiter_show_on_load
#' waiter_update
#' @importFrom rlang .data
#' @importFrom rintrojs introBox introjs introjsUI readCallback
#' @importFrom DT datatable dataTableOutput renderDataTable formatRound
#' formatPercentage
#' @importFrom shinyjs useShinyjs hide show
#' @importFrom utils write.csv write.table read.delim
#' @importFrom stringr str_to_title str_split
#' @importFrom tidyr pivot_longer
#' @importFrom stats sd cor.test complete.cases dist hclust
#' @importFrom ComplexHeatmap Heatmap make_comb_mat UpSet comb_size
#' upset_top_annotation extract_comb
#' @importFrom grid gpar unit
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#' makeInteractiveComplexHeatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom ggforce facet_grid_paginate
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom ggpubr stat_cor
#' @importFrom corrplot corrplot
#' @importFrom SummarizedExperiment assays
#' @importFrom shinyWidgets actionBttn
#' @import shinyBS
#'
#' @name DeconvExplorer-pkg
#' @docType package
"_PACKAGE"
