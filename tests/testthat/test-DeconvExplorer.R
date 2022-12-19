test_that("Shiny app is generated", {
  expect_true(
    is(DeconvExplorer::DeconvExplorer(), "shiny.appobj")
  )
})

test_that("Helper funcs", {
  internal_list <- returnSelectedDeconvolutions(c("momf"), deconv_list)
  expect_true(
    is(internal_list, "list")
  )
  expect_equal(length(internal_list), 1)
})
