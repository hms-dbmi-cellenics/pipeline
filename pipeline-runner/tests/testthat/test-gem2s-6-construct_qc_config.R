test_that("construct_qc_config works with bpcells", {
  scdata_list <- mock_scdata_list(use_bpcells = TRUE)
  expect_no_error(
    construct_qc_config(
      scdata_list,
      unfiltered_samples = names(scdata_list),
      technology = "10X"
    )
  )
})


test_that("cellsize filter is disabled by default and classifier is pre-filtered", {
  scdata_list <- mock_scdata_list()

  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = unfiltered_samples,
    technology = "10X"
  )

  for (sample in names(scdata_list)) {
    if (sample %in% unfiltered_samples) {
      expect_true(qc_config$classifier[[sample]]$enabled)
      expect_false(qc_config$classifier[[sample]]$prefiltered)
    } else {
      expect_false(qc_config$classifier[[sample]]$enabled)
      expect_true(qc_config$classifier[[sample]]$prefiltered)
    }

    expect_false(qc_config$cellSizeDistribution[[sample]]$enabled)
  }
})


test_that("cellsize filter is disabled by default and classifier is not pre-filtered", {
    scdata_list <- mock_scdata_list()

  unfiltered_samples <- c()
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples, technology = "10X")

  for (sample in names(scdata_list)) {
    expect_false(qc_config$cellSizeDistribution[[sample]]$enabled)
    if (sample %in% unfiltered_samples) {
      expect_true(qc_config$classifier[[sample]]$enabled)
      expect_false(qc_config$classifier[[sample]]$prefiltered)
    } else {
      expect_false(qc_config$classifier[[sample]]$enabled)
      expect_true(qc_config$classifier[[sample]]$prefiltered)
    }
  }
})

test_that("customize_doublet_config sets threshold to 0 when there are no singlets", {
    scdata_list <- mock_scdata_list()

  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples, technology = "10X")

  for (sample in names(scdata_list)) {
    scdata_list[[sample]]$doublet_class <- "doublet"
    config <- customize_doublet_config(scdata_list[[sample]], qc_config)
    expect_equal(config$filterSettings$probabilityThreshold, 0)
  }
})


test_that("classifier filter config is enabled for unfiltered samples and disabled for pre-filtered samples", {
    scdata_list <- mock_scdata_list()

  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples, technology = "10X")

  for (sample in names(scdata_list)) {
    if (sample %in% unfiltered_samples) {
      expect_true(qc_config$classifier[[sample]]$enabled)
      expect_false(qc_config$classifier[[sample]]$prefiltered)
    } else {
      expect_false(qc_config$classifier[[sample]]$enabled)
      expect_true(qc_config$classifier[[sample]]$prefiltered)
    }
  }
})

test_that("NumGenesVsUmis filter config has spline as default for Parse Datasets", {
    scdata_list <- mock_scdata_list()

  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples, "parse")

  for (sample in names(scdata_list)) {
    expect_true(qc_config$numGenesVsNumUmis[[sample]]$filterSettings$regressionType == "spline")
  }
})

test_that("NumGenesVsUmis filter config has linear as default for 10x datasets", {
    scdata_list <- mock_scdata_list()

  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples, "10X")

  for (sample in names(scdata_list)) {
    expect_true(qc_config$numGenesVsNumUmis[[sample]]$filterSettings$regressionType == "linear")
  }
})
