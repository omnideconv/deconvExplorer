test_that("Signature exploration plots are correct", {
  p_sig_met <- plot_signatureGenesPerMethod(signatures = signature_list)

  expect_s3_class(p_sig_met, "gg")

  expect_error(
    plot_signatureGenesPerMethod(signatures = signature_list[[1]])
  )

  p_cond_nr_met <- plot_conditionNumberPerMethod(signature_list)
  expect_s3_class(p_cond_nr_met, "gg")
  expect_error(
    plot_conditionNumberPerMethod(signatures = signature_list[[1]])
  )

  p_meanentropy <- plot_meanEntropyPerMethod(signature_list)
  expect_s3_class(p_meanentropy, "gg")
  expect_error(
    plot_meanEntropyPerMethod(signatures = signature_list[[1]])
  )

  p_sig_clust <- plot_signatureClustered(
    signature = signature_list$bisque,
    score = "gini", annotation_type = "bar"
  )

  expect_s4_class(p_sig_clust, "HeatmapList")
  expect_error(
    plot_signatureClustered(
      signature = signature_list,
      score = "gini", annotation_type = "bar"
    )
  )

  # generating a slghtly larger list with 3 elements
  siglist2 <- list(
    "dwls" = signature_example,
    "momf" = signature_example,
    "bisque" = signature_example
  )
  p_sig_upset <- plot_signatureUpset(siglist2, mode = "union")
  expect_true(is(p_sig_upset, "list"))

  sig_vec <- download_signatureUpset(
    signatures = siglist2,
    combination = c("dwls", "bisque"),
    mode = "intersect"
  )

  expect_true(is(sig_vec, "vector"))
})
