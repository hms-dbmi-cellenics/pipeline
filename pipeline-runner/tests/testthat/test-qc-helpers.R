test_that("downsample_plotdata keeps at most max_number_of_cells", {
    ncol_sample <- 1001
    max_number_of_cells <- 1000
    nkeep <- downsample_plotdata(ncol_sample, max_number_of_cells)

    expect_equal(nkeep, max_number_of_cells)
})


test_that("downsample_plotdata keeps all cells if less than max_number_of_cells", {
    ncol_sample <- 999
    max_number_of_cells <- 1000
    nkeep <- downsample_plotdata(ncol_sample, max_number_of_cells)

    expect_equal(nkeep, ncol_sample)
})

test_that("get_qc_task_order keeps only single-cell steps for non-spatial technologies", {
  task_names <- names(QC_TASK_LIST)

  for (technology in list(NULL, "10x", "rhapsody", "parse")) {
    order <- get_qc_task_order(task_names, technology)
    expect_false(any(startsWith(order, "spatial")))
    expect_equal(
      order,
      c(
        "classifier", "cellSizeDistribution", "mitochondrialContent",
        "numGenesVsNumUmis", "doubletScores", "dataIntegration", "configureEmbedding"
      )
    )
  }
})


test_that("get_qc_task_order keeps only each spatial technology's steps", {
  task_names <- names(QC_TASK_LIST)

  expect_equal(
    get_qc_task_order(task_names, "visium_hd"),
    c(
      "spatialUmiOutlier", "spatialNumGenesOutlier", "spatialMitoOutlier",
      "dataIntegration", "configureEmbedding"
    )
  )
  expect_equal(
    get_qc_task_order(task_names, "xenium"),
    c(
      "spatialCountDensity", "spatialUmiOutlier", "spatialNumGenesOutlier",
      "dataIntegration", "configureEmbedding"
    )
  )
})


test_that("single-cell doubletScores hands filtered cell ids off to dataIntegration", {
  order <- get_qc_task_order(names(QC_TASK_LIST), "10x")
  expect_equal(get_next_qc_task("doubletScores", order), "dataIntegration")
})


test_that("the last spatial filter hands filtered cell ids off to dataIntegration", {
  order <- get_qc_task_order(names(QC_TASK_LIST), "visium_hd")
  expect_equal(get_next_qc_task("spatialMitoOutlier", order), "dataIntegration")

  order <- get_qc_task_order(names(QC_TASK_LIST), "xenium")
  expect_equal(get_next_qc_task("spatialNumGenesOutlier", order), "dataIntegration")
})


test_that("get_next_qc_task returns NULL for the last or an unknown task", {
  order <- c("a", "b")
  expect_null(get_next_qc_task("b", order))
  expect_null(get_next_qc_task("z", order))
})


test_that("get_scdata_technology reads the technology from a Seurat object or list", {
  scdata <- mock_spatial_scdata(technology = "visium_hd")

  expect_equal(get_scdata_technology(scdata), "visium_hd")
  expect_equal(get_scdata_technology(list(s1 = scdata)), "visium_hd")
  expect_null(get_scdata_technology(list()))
  expect_null(get_scdata_technology(NULL))
})

test_that("columnar_plot_data rounds metrics to 3 and coordinates to 1 decimal", {
  out <- columnar_plot_data(list(
    cellId = c(1L, 2L),
    value = c(1.23456, NA),
    x = c(12345.678, 2.25),
    y = c(0.04, 9.99)
  ))

  expect_named(out, "columnar")
  expect_equal(unclass(out$columnar$cellId), c(1L, 2L))
  expect_equal(unclass(out$columnar$value), c(1.235, NA))
  expect_equal(unclass(out$columnar$x), c(12345.7, 2.2))
  expect_equal(unclass(out$columnar$y), c(0, 10))
})

test_that("columnar_plot_data serialises every column as a JSON array, even length 1", {
  json <- RJSONIO::toJSON(columnar_plot_data(list(cellId = 7L, value = 1.5, x = 3, y = 4)))
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_equal(parsed$columnar$cellId, list(7L))
  expect_equal(parsed$columnar$value, list(1.5))
})
