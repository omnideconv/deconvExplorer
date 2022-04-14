#  basic script to assemble the rintrojs tour for DeconvExplorer

getTour <- function(){
  element <- intro <- character(0)
  
  # basic usage: add elements to the tour based on their ID
  # element <- c(element, '#ID')
  # intro <- c(intro, "Tour content")
  # elements can be addressed by id, in case a specific id is necessary, add id = "#id" to the elements
  
  
  element <- c(element, "#logo")
  intro <- c(intro, "Welcome to DeconvExplorer, an interactive Web Interface to the omnideconv deconvolution toolset. <br/>
             In this tour we will give an overview over DeconvExplorers functionality.")
  
  element <- c(element, ".sidebar")
  intro <- c(intro, "DeconvExplorer is devided into different sections. On the right side you can see the deconvolution section where your data
             can be uploaded, you can choose different different deconvolution and signature calculation options and 
             compare the results.")
  
  element <- c(element, "#tour_upload")
  intro <- c(intro, "Upload your data here")
  
  element <- c(element, "#tour_deconvSettings")
  intro <- c(intro, "Here you can select your preferred deconvolution method. If applicable, your can choose 
             a cell-type specific deconvolution signature, or choose a method to calculate a new one.")
  
  element <- c(element, "#tour_deconvPlotSettings")
  intro <- c(intro, "All deconvolution results are collected here and can be identified by selecting the deconvolution method and the used signature. 
             After the selection you can load the result to be plotted (load result) or add it to an existing visualization to simplify the comparison (add to plot):")
  
  element <- c(element, "#tour_deconvPlot")
  intro <- c(intro, "All chosen deconvolution results are visualized here. Feel free to explore the different available plotting options!")
  
  
  
  # finalize tour
  data.frame(element=element, intro=intro, stringsAsFactors = FALSE)
  
}