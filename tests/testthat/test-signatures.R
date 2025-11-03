test_that("Signature refinement operations are correct", {
  signature_renamed <- renameCellType(
    signature_mat = signature_list$bisque,
    cell_type = "B",
    new_celltype_name = "B.cells"
  )

  expect_true("B" %in% colnames(signature_list$bisque))
  expect_true("B.cells" %in% colnames(signature_renamed))


  signature_nozeros <- removePercentZeros(
    signature_mat = signature_list$bisque,
    max_percentage_zeroes = 0.5
  )
  expect_equal(nrow(signature_list$bisque) - nrow(signature_nozeros), 8)

  signature_nounspecific <- removeUnspecificGenes(
    signature_mat = signature_list$momf,
    number_of_bins = 3,
    max_count = 1
  )
  expect_equal(nrow(signature_list$momf) - nrow(signature_nounspecific), 199)

  signature_selgenes <- selectGenesByScore(
    signature_mat = signature_list$bisque,
    scoring_method = "gini",
    genes_per_cell_type = 50
  )
  expect_equal(nrow(signature_list$bisque) - nrow(signature_selgenes), 126)

  signature_selgenes_entropy <- selectGenesByScore(
    signature_mat = signature_list$bisque,
    scoring_method = "entropy",
    genes_per_cell_type = 50
  )


  expect_error(renameCellType())
  expect_error(renameCellType(signature_list$bisque, cell_type = "", new_celltype_name = "B cells"))
  expect_error(renameCellType(signature_list$bisque, cell_type = "B"))
  expect_error(renameCellType(signature_list$bisque, cell_type = "C cells", new_celltype_name = "B cells"))

  expect_error(removePercentZeros())
  expect_error(removePercentZeros(signature_list$bisque, max_percentage_zeroes = -1))

  expect_error(removeUnspecificGenes())
  expect_error(removeUnspecificGenes(signature_list$bisque, number_of_bins = 1))
  expect_error(removeUnspecificGenes(signature_list$bisque, max_count = -2))
  expect_error(removeUnspecificGenes(signature_list$bisque,
    max_count = 3,
    labels = c("A", "B", "C", "D")
  ))
})
