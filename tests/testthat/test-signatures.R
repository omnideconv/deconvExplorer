test_that("Signature refinement operations are correct", {
  signature_renamed <- renameCellType(
    signature = signature_list$bisque,
    cellType = "B",
    newName = "B.cells"
  )

  expect_true("B" %in% colnames(signature_list$bisque))
  expect_true("B.cells" %in% colnames(signature_renamed))

  

  signature_nozeros <- removePercentZeros(baseSignature = signature_list$bisque,
                                          percentage = 0.5)
  expect_equal(nrow(signature_list$bisque) - nrow(signature_nozeros), 8)
  
  signature_nounspecific <- removeUnspecificGenes(signature = signature_list$momf,
                                                  numberOfBins = 3,
                                                  maxCount = 1)
  expect_equal(nrow(signature_list$momf) - nrow(signature_nounspecific), 199)
  
  signature_selgenes <- selectGenesByScore(signature = signature_list$bisque,
                                           method = "gini",
                                           genesPerCellType = 50)
  expect_equal(nrow(signature_list$bisque) - nrow(signature_selgenes), 126)
  
  signature_selgenes_entropy <- selectGenesByScore(signature = signature_list$bisque,
                                                   method = "entropy",
                                                   genesPerCellType = 50)
  
  
  expect_error(renameCellType())
  expect_error(renameCellType(signature_list$bisque, cellType = "", newName = "B cells"))
  expect_error(renameCellType(signature_list$bisque, cellType = "B"))
  expect_error(renameCellType(signature_list$bisque, cellType = "C cells", newName = "B cells"))
  
  expect_error(removePercentZeros())
  expect_error(removePercentZeros(signature_list$bisque, percentage = -1))
  
  expect_error(removeUnspecificGenes())
  expect_error(removeUnspecificGenes(signature_list$bisque, numberOfBins = 1))
  expect_error(removeUnspecificGenes(signature_list$bisque, maxCount = -2))
  expect_error(removeUnspecificGenes(signature_list$bisque, 
                                     maxCount = 3, 
                                     labels = c("A", "B", "C", "D")))
  
  
  entropy_calc <- scoreEntropy(signature_list$bisque[1, ])
  expect_true(is.numeric(entropy_calc))
})
