#  basic script to assemble the rintrojs tour for DeconvExplorer

getTour <- function() {
  element <- intro <- character(0)

  # basic usage: add elements to the tour based on their ID
  # element <- c(element, '#ID')
  # intro <- c(intro, "Tour content")
  # elements can be addressed by id, in case a specific id is necessary, add id = "#id" to the elements


  element <- c(element, "#logo")
  intro <- c(intro, "Welcome to DeconvExplorer, an interactive framework to perform, evaluate and enhance cell type deconvolution from transcriptome data. <br/>
             In this tour we will give an overview over DeconvExplorers functionality. DeconvExplorer performs deconvolution with the omnideconv toolset. ")

  element <- c(element, ".sidebar")
  intro <- c(intro, "DeconvExplorer is devided into six modules. You can always switch between modules but keep in mind that, while a deconvolution is running, plots are not updating.")

  element <- c(element, "#tour_upload")
  intro <- c(intro, "Upload your data here as a csv, tsv, txt or rds file. Further data requirements are listed below.")

  element <- c(element, "#tour_sample")
  intro <- c(intro, "To test the interface you can load three different sample datasets.")

  element <- c(element, "#tour_simbu")
  intro <- c(intro, "If you want to use a simulated pseudo-bulk sample from SimBu you can upload your simulation in rds format.")

  element <- c(element, "#tour_signatureUpload")
  intro <- c(intro, "If required you can upload a custom signature.")

  element <- c(element, "#tour_deconvSettings")
  intro <- c(intro, "To perform a deconvolution, select a deconvolution method. If applicable you can choose a custom signature as well.")

  element <- c(element, "#tour_deconvPlotSettings")
  intro <- c(intro, "All deconvolution results are collected here and can be selected for plotting below.")

  element <- c(element, "#tour_deconvPlot")
  intro <- c(intro, "All chosen deconvolution results are visualized here. Feel free to explore the different available plotting options!")

  element <- c(element, "#tour_genesPlot")
  intro <- c(intro, "The Signature Exploration module plots multiple metrics to analyse and compare expression signatures.")

  element <- c(element, "#tour_refinementHeatmap")
  intro <- c(intro, "In the Signature Refinement module you can further subset the genes of a signature to analyse the impact on deconvolution performance.")

  element <- c(element, "#tour_infoOverview")
  intro <- c(intro, "More information about each module can be found here.")

  # finalize tour
  data.frame(element = element, intro = intro, stringsAsFactors = FALSE)
}