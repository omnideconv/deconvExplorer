#' Get subset of deconvolution results
#'
#' This function returns the requested deconvolution results as a list
#'
#' @param deconv_select list of deconvolution result names which should be plotted
#' @param deconv_list List containing all available deconvolution results
#'
#' @return list of deconvolution results
#'
#' @export
#'
#' @examples
#' deconv <- readRDS(system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer"))
#'
#' # list containting deconvolution results
#' deconvList <- list("bisque" = deconv, "momf" = deconv)
#'
#' returnSelectedDeconvolutions(c("momf"), deconvList)
returnSelectedDeconvolutions <- function(deconv_select, deconv_list) {
  deconvolutions <- list()

  for (deconvolution in deconv_select) {
    deconvolutions[deconvolution] <- deconv_list[deconvolution]
  }

  return(deconvolutions)
}


#' Modal window to print error messages or other warnings
#'
#' @param error_message
#'
#' @return NULL
#' @export
errorModal <- function(error_message = NULL) {
  modalDialog(
    p(error_message, style = "color:red;"),
    easyClose = T,
    modalButton("Cancel")
  )
}
