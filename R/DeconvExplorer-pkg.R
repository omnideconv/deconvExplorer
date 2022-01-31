#' DeconvExplorer
#' 
#' DeconvExplorer description TODO
#' 
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @import omnideconv
#' @importFrom shinydashboard box dashboardBody dashboardHeader dashboardPage
#' dashboardSidebar dropdownMenu menuItem notificationItem sidebarMenu
#' tabItem tabItems
#' @importFrom plotly ggplotly plotlyOutput renderPlotly plot_ly layout config
#' @importFrom ggplot2 aes aes_ aes_string coord_cartesian coord_flip element_text
#' facet_wrap geom_abline geom_boxplot geom_col geom_jitter geom_point
#' geom_tile ggplot guide_colorbar guides labs scale_fill_gradient theme element_blank
#' geom_hline
#' @importFrom shinycssloaders withSpinner 
#' @importFrom waiter Waitress
#' @importFrom rintrojs introBox introjs introjsUI
#' @importFrom DT datatable dataTableOutput renderDataTable formatRound formatPercentage
#' @importFrom shinyjs useShinyjs hide show
#' @importFrom magrittr "%>%"
#' @importFrom utils write.csv
#' @importFrom stringr str_subset
#' @importFrom tidyr pivot_longer
#' 
#' @name DeconvExplorer-pkg
#' @docType package
NULL