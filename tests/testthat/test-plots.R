test_that("Benchmark scatter plots are correct", {
  p <- plot_benchmark_scatter(gtruth = RefData,
                              estimate = deconv_list)
  
  expect_s3_class(p, "gg")
  
  expect_error(
    plot_benchmark_scatter(gtruth = RefData,
                           estimate = deconv_example)
  )
})


test_that("Benchmark correlation plots are correct", {
  p <- plot_benchmark_correlation(gtruth = RefData,
                                  estimate = deconv_list)
  
  expect_true(is(p, "list"))
  expect_equal(names(p), c("corr", "corrPos", "arg"))
  
  expect_error(
    plot_benchmark_correlation(gtruth = RefData,
                               estimate = deconv_example)
  )
  
  expect_error(
    plot_benchmark_correlation(gtruth = RefData,
                               estimate = deconv_list, 
                               plot_method = "boxplot")
  )
  file.remove("Rplots.pdf")
})


test_that("Benchmark rmse are correct", {
  p_hm <- plot_benchmark_rmse(gtruth = RefData,
                              estimate = deconv_list,
                              plot_type = "heatmap")
  
  p_box <- plot_benchmark_rmse(gtruth = RefData,
                               estimate = deconv_list,
                               plot_type = "boxplot")
  
  expect_true(is(p_hm, "list"))
  expect_equal(names(p_hm), c("corr", "corrPos", "arg"))
  
  expect_s3_class(p_box, "gg")
  
  expect_error(
    plot_benchmark_rmse(gtruth = RefData,
                        estimate = deconv_list,
                        plot_type = "pie")
  )
})
