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

test_that("construct_qc_config sets geosketch to default", {
  scdata_list <- mock_scdata_list()
  qc_config <- expect_no_error(
    construct_qc_config(
      scdata_list,
      unfiltered_samples = names(scdata_list),
      technology = "10X"
    )
  )

  expect_equal(qc_config$dataIntegration$downsampling$method, "default")
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
  }
})

test_that("NumGenesVsUmis filter config has spline as default for Parse Datasets", {
    scdata_list <- mock_scdata_list()

  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = unfiltered_samples,
    "parse"
  )

  for (sample in names(scdata_list)) {
    expect_true(
      qc_config$numGenesVsNumUmis[[sample]]$filterSettings$regressionType == "spline"
    )
  }
})

test_that("visium_hd config builds the three spatial filters with correct defaults", {
  scdata_list <- mock_scdata_list()

  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = names(scdata_list),
    technology = "visium_hd"
  )

  spatial_steps <- list(
    spatialUmiOutlier = "lower",
    spatialNumGenesOutlier = "lower",
    spatialMitoOutlier = "upper"
  )

  for (step in names(spatial_steps)) {
    expect_true(step %in% names(qc_config))
    for (sample in names(scdata_list)) {
      cfg <- qc_config[[step]][[sample]]
      expect_true(cfg$enabled)
      expect_true(cfg$auto)
      expect_equal(cfg$filterSettings$method, "localOutlier")
      expect_equal(cfg$filterSettings$methodSettings$localOutlier$cutoff, 3)

      # absolute defaults are this sample's 3-MAD thresholds; no max for UMIs/genes
      absolute <- cfg$filterSettings$methodSettings$absoluteThreshold
      if (spatial_steps[[step]] == "lower") {
        expect_true(is.numeric(absolute$min))
        expect_true("max" %in% names(absolute))
        expect_null(absolute$max)
      } else {
        expect_false("min" %in% names(absolute))
        expect_true(is.numeric(absolute$max))
      }
    }
  }

  # integration + embedding still present
  expect_true(all(c("dataIntegration", "configureEmbedding") %in% names(qc_config)))
})

test_that("visium_hd config disables the single-cell-only steps", {
  scdata_list <- mock_scdata_list()

  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = names(scdata_list),
    technology = "visium_hd"
  )

  single_cell_steps <- c(
    "cellSizeDistribution", "mitochondrialContent", "classifier",
    "numGenesVsNumUmis", "doubletScores"
  )
  # the single-cell filter steps are not part of the spatial config at all
  expect_false(any(single_cell_steps %in% names(qc_config)))
})

test_that("xenium config builds count density, then UMI and gene filters that are off", {
  scdata_list <- mock_scdata_list()

  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = names(scdata_list),
    technology = "xenium"
  )

  expect_equal(
    names(qc_config),
    c("spatialCountDensity", "spatialUmiOutlier", "spatialNumGenesOutlier",
      "dataIntegration", "configureEmbedding")
  )

  for (sample in names(scdata_list)) {
    density <- qc_config$spatialCountDensity[[sample]]
    expect_true(density$enabled)
    expect_true(density$auto)
    expect_equal(density$filterSettings$nmads, 3)
    expect_true(density$filterSettings$filterControls)
    expect_equal(density$filterSettings$maxControlCounts, 0)

    for (step in c("spatialUmiOutlier", "spatialNumGenesOutlier")) {
      cfg <- qc_config[[step]][[sample]]
      expect_false(cfg$enabled)
      expect_true(cfg$auto)
      expect_equal(cfg$filterSettings$method, "absoluteThreshold")
      expect_null(cfg$filterSettings$methodSettings$localOutlier)
      expect_true(is.numeric(cfg$filterSettings$methodSettings$absoluteThreshold$min))
    }
  }
})

test_that("xenium config disables the single-cell-only steps", {
  scdata_list <- mock_scdata_list()

  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = names(scdata_list),
    technology = "xenium"
  )

  single_cell_steps <- c(
    "cellSizeDistribution", "mitochondrialContent", "classifier",
    "numGenesVsNumUmis", "doubletScores"
  )
  expect_false(any(single_cell_steps %in% names(qc_config)))
})

test_that("non-spatial (10X) config shape is unchanged (no spatial keys)", {
  scdata_list <- mock_scdata_list()

  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = names(scdata_list),
    technology = "10X"
  )

  expect_setequal(
    names(qc_config),
    c("cellSizeDistribution", "mitochondrialContent", "classifier",
      "numGenesVsNumUmis", "doubletScores", "dataIntegration", "configureEmbedding")
  )
  expect_false(any(startsWith(names(qc_config), "spatial")))
})

test_that("NumGenesVsUmis filter config has linear as default for 10x datasets", {
    scdata_list <- mock_scdata_list()

  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(
    scdata_list,
    unfiltered_samples = unfiltered_samples,
    "10X"
  )

  for (sample in names(scdata_list)) {
    expect_true(
      qc_config$numGenesVsNumUmis[[sample]]$filterSettings$regressionType == "linear"
    )
  }
})
