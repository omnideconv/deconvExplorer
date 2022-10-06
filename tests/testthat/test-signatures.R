test_that("Signature refinement operations are correct", {
  signature_renamed <- renameCellType(
    signature = signature_list$bisque,
    cellType = "B",
    newName = "B.cells"
  )

  expect_true("B" %in% colnames(signature_list$bisque))
  expect_true("B.cells" %in% colnames(signature_renamed))

  # COSTODO: will need to put away the notifications
  # signature_nozeros <- removePercentZeros(baseSignature = signature_list$bisque,
  #                                         percentage = 0.5)
  #
  # signature_nounspecific <- removeUnspecificGenes(signature = signature_list$momf,
  #                                                 numberOfBins = 3,
  #                                                 maxCount = 1)
  #
  # signature_selgenes <- selectGenesByScore(signature = signature_list,
  #                                          method = "gini",
  #                                          genesPerCellType = 50)

  entropy_calc <- scoreEntropy(signature_list$bisque[1, ])
  expect_true(is.numeric(entropy_calc))
})
