library("DeconvExplorer")

message("--- Loading packages...")
suppressPackageStartupMessages({
  library("omnideconv")
})
message("- Done!")



message("--- Generating objects for the testing setup...")
data("RefData", package = "omnideconv")
RefData <- as.data.frame(RefData)
deconv_example <- readRDS(system.file("extdata", "deconvolution_example.rds", package = "DeconvExplorer"))
deconv_list <- list(
  "momf" = deconv_example,
  "bisque" = deconv_example
)

signature_example <- readRDS(system.file("extdata", "signature_example.rds", package = "DeconvExplorer"))
signature_list <- list(
  "bisque" = signature_example,
  "momf" = signature_example
)
message("- Done!")

message("--- Test setup script completed!")
