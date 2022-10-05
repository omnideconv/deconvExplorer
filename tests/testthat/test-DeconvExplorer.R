test_that("Shiny app is generated", {
  expect_true(
    is(DeconvExplorer::DeconvExplorer(), "shiny.appobj")
  )
  
  expect_true(
    is(DeconvExplorer::DeconvExplorer(), "shiny.appobj")
  )
})
