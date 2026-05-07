mock_config <- function() {
  config <- list(
    auto = TRUE,
    enabled = TRUE,
    filterSettings = list(
      regressionType = "linear",
      regressionTypeSettings = list(
        linear = list(p.level = 0.1)
      )
    )
  )

  return(config)
}


test_that("filter_gene_umi_outlier updates linear p.level in config if auto", {
  scdata_list <- mock_scdata_list()
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)

  config$filterSettings$regressionTypeSettings$linear$p.level <- 1
  out <- filter_gene_umi_outlier(scdata_list, config, sample2_id, cells_id)
  new <- out$config$filterSettings$regressionTypeSettings$linear$p.level

  expect_lt(new, 1)
})


test_that("filter_gene_umi_outlier is sample specific", {

  # one outlier in first sample
  scdata_list <- mock_scdata_list(with_outlier = TRUE)
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)

  cell_ids1 <- scdata_list[[sample1_id]]$cells_id
  cell_ids2 <- scdata_list[[sample2_id]]$cells_id
  out1 <- filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)
  out2 <- filter_gene_umi_outlier(scdata_list, config, sample2_id, cells_id)

  # filtered something
  expect_lt(length(out1$new_ids[[sample1_id]]), length(cells_id[[sample1_id]]))
  expect_lt(length(out2$new_ids[[sample2_id]]), length(cells_id[[sample2_id]]))

  # didn't filter other sample
  expect_true(all(cell_ids1 %in% out2$new_ids[[sample1_id]]))
  expect_true(all(cell_ids2 %in% out1$new_ids[[sample2_id]]))
})

test_that("filter_gene_umi_outlier can filter cells and works with bpcells", {
  for (use_bpcells in c(FALSE, TRUE)) {
     # single outlier in single sample
    scdata_list <- mock_scdata_list(with_outlier = TRUE, use_bpcells = use_bpcells)
    sample1_id <- names(scdata_list)[1]
    sample2_id <- names(scdata_list)[2]

    config <- mock_config()
    cells_id <- mock_ids(scdata_list)
    nstart <- ncol(scdata_list[[sample1_id]])

    out <-
      filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)
    expect_lt(length(unlist(out$new_ids[[sample1_id]])), nstart)
  }
})


test_that("filter_gene_umi_outlier can be disabled", {

  # single outlier in single sample
  scdata_list <- mock_scdata_list(with_outlier = TRUE)
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)
  nstart <- ncol(scdata_list[[sample1_id]])

  # disabled
  config$enabled <- FALSE

  out <- filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)
  expect_equal(ncol(out$data[[sample1_id]]), nstart)
  expect_equal(out$new_ids[[sample1_id]], 0:39)
  expect_equal(out$new_ids[[sample2_id]], 40:79)
})



test_that("filter_gene_umi_outlier return the appropriate plot data", {
  scdata_list <- mock_scdata_list()
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)

  out <- filter_gene_umi_outlier(scdata_list, config, sample2_id, cells_id)
  plot_data <- out$plotData[[1]]

  # points and lines data
  expect_equal(names(plot_data), c("pointsData", "linesData"))

  # one value per cell
  points_data <- plot_data$pointsData
  lines_data <- plot_data$linesData
  expect_equal(length(points_data), ncol(scdata_list[[sample2_id]]))

  # has the right names
  expect_equal(names(points_data[[1]]), c("log_molecules", "log_genes"))
  expect_equal(
    names(lines_data[[1]][[1]]),
    c("log_molecules", "lower_cutoff", "upper_cutoff")
  )
})



test_that("filter_gene_umi_outlier gives different results with spline fit", {
  # single outlier in single sample
  scdata_list <- mock_scdata_list(with_outlier = TRUE)
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)
  nstart <- ncol(scdata_list[[sample1_id]])

  out1 <- filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)
  expect_equal(ncol(out1$data[[sample1_id]]), nstart)
  expect_lt(length(out1$new_ids[[sample1_id]]), length(cells_id[[sample1_id]]))

  config$filterSettings$regressionType <- "spline"

  out2 <- filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)

  expect_true(!identical(out1$plotData, out2$plotData))
  expect_true(
    !identical(
      out1$new_ids[[sample1_id]],
      out2$new_ids[[sample1_id]]
    )
  )
})

test_that("filter_gene_umi_outlier spline fit works with bpcells", {
  # single outlier in single sample
  scdata_list <- mock_scdata_list(with_outlier = TRUE, use_bpcells = TRUE)
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)
  config$filterSettings$regressionType <- "spline"

  expect_no_error(
    filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)
  )
})

test_that("Gene UMI filter works if input is a a float-interpretable string", {
  # single outlier in single sample
  scdata_list <- mock_scdata_list(with_outlier = TRUE)
  sample1_id <- names(scdata_list)[1]
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)
  nstart <- ncol(scdata_list[[sample2_id]])
  config$auto <- FALSE
  config$filterSettings$predictionInterval <- "0.3"

  out_number <- filter_gene_umi_outlier(scdata_list, config, sample2_id, cells_id)

  config$filterSettings$regressionTypeSettings$linear$p.level <- "0.1"
  out_string <- filter_gene_umi_outlier(scdata_list, config, sample2_id, cells_id)

  expect_equal(ncol(out_number$data[[sample2_id]]), nstart)
  expect_equal(
    out_string$new_ids[[sample2_id]],
    out_number$new_ids[[sample2_id]]
  )
})

test_that("Gene UMI filter throws error if input is a non float-interpretable string", {
  scdata_list <- mock_scdata_list()
  sample2_id <- names(scdata_list)[2]

  config <- mock_config()
  config$auto <- FALSE
  config$filterSettings$regressionTypeSettings$linear$p.level <-
    "asd"

  expect_error(
    filter_gene_umi_outlier(scdata, config, sample2_id, cells_id)
  )
})

test_that("Gene UMI filter works with manual settings and default prediction interval value", {
  scdata_list <- mock_scdata_list()
  sample1_id <- names(scdata_list)[1]

  config <- mock_config()
  cells_id <- mock_ids(scdata_list)
  type <- "spline"
  config$filterSettings$regressionType <- type

  config$auto <- TRUE
  out_auto <- filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)

  config$auto <- FALSE
  out_manual <- filter_gene_umi_outlier(scdata_list, config, sample1_id, cells_id)

  expect_null(config$filterSettings$predictionInterval)
  expect_equal(
    out_auto$config$filterSettings$regressionTypeSettings[[type]]$p.level,
    out_manual$config$filterSettings$regressionTypeSettings[[type]]$p.level
  )
})
