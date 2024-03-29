test_that("Signature exploration plots are correct", {
  p_sig_met <- plot_signatureGenesPerMethod(signature_list = signature_list)

  expect_s3_class(p_sig_met, "gg")

  expect_error(
    plot_signatureGenesPerMethod(signature_list = signature_list[[1]])
  )

  p_cond_nr_met <- plot_conditionNumberPerMethod(signature_list)
  expect_s3_class(p_cond_nr_met, "gg")
  expect_error(
    plot_conditionNumberPerMethod(signature_list = signature_list[[1]])
  )

  p_meanentropy <- plot_meanEntropyPerMethod(signature_list)
  expect_s3_class(p_meanentropy, "gg")
  expect_error(
    plot_meanEntropyPerMethod(signature_list = signature_list[[1]])
  )

  p_sig_clust <- plot_signatureClustered(
    signature_mat = signature_list$bisque,
    scoring_method = "gini", annotation_type = "bar"
  )

  expect_s4_class(p_sig_clust, "HeatmapList")
  expect_error(
    plot_signatureClustered(
      signature_mat = signature_list,
      scoring_method = "gini", annotation_type = "bar"
    )
  )
  expect_error(
    plot_signatureClustered(
      signature_mat = signature_list,
      scoring_method = "gino", annotation_type = "bar"
    )
  )
  expect_error(
    plot_signatureClustered(
      signature_mat = signature_list,
      scoring_method = "gini", annotation_type = "bart"
    )
  )
  # some variations on that
  p_sig_clust2 <- plot_signatureClustered(
    signature_mat = signature_list$bisque,
    scoring_method = "gini", annotation_type = "line"
  )
  p_sig_clust3 <- plot_signatureClustered(
    signature_mat = signature_list$bisque,
    scoring_method = "entropy", annotation_type = "bar"
  )
  p_sig_clust4 <- plot_signatureClustered(
    signature_mat = signature_list$bisque,
    scoring_method = "entropy", annotation_type = "bar"
  )

  expect_s4_class(p_sig_clust2, "HeatmapList")
  expect_s4_class(p_sig_clust3, "HeatmapList")
  expect_s4_class(p_sig_clust4, "HeatmapList")

  # generating a slghtly larger list with 3 elements
  siglist2 <- list(
    "dwls" = signature_example,
    "momf" = signature_example,
    "bisque" = signature_example
  )
  p_sig_upset <- plot_signatureUpset(siglist2, upset_mode = "union")
  expect_true(is(p_sig_upset, "list"))

  # some variations on that
  p_sig_upset2 <- plot_signatureUpset(siglist2, upset_mode = "union", order_sets = "degree")
  expect_true(is(p_sig_upset2, "list"))
  p_sig_upset3 <- plot_signatureUpset(siglist2, upset_mode = "union", color_by_degrees = FALSE)
  expect_true(is(p_sig_upset3, "list"))


  sig_vec <- download_signatureUpset(
    signature_list = siglist2,
    combination_to_include = c("dwls", "bisque"),
    upset_mode = "intersect"
  )
  expect_true(is(sig_vec, "vector"))

  sig_vec_null <- download_signatureUpset(
    signature_list = siglist2,
    combination_to_include = NULL,
    upset_mode = "intersect"
  )
  expect_null(sig_vec_null)
})
