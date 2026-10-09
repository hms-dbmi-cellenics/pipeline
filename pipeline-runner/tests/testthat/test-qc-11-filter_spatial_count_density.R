mock_density_config <- function(auto = FALSE, nmads = 3, filter_controls = TRUE,
                                 max_control_counts = 0) {
  list(
    enabled = TRUE,
    auto = auto,
    filterSettings = list(
      nmads = nmads,
      filterControls = filter_controls,
      maxControlCounts = max_control_counts,
      binStep = 0.05
    )
  )
}

# a Xenium-like sample with known counts, cell areas and control counts
mock_density_sample <- function(counts, area = rep(100, length(counts)),
                                controls = rep(0, length(counts))) {
  scdata_list <- mock_spatial_scdata_list(ncells = length(counts))
  scdata_list$spatial1$nCount_RNA <- counts
  scdata_list$spatial1$cell_area <- area
  scdata_list$spatial1$control_counts <- controls
  scdata_list
}

# 60 cells around 5 counts/um^2 plus a few low-density ones
mock_density_counts <- function() {
  set.seed(1)
  c(round(exp(rnorm(60, log(500), 0.2))), 20, 30)
}

test_that("removes 3-MAD low count-density outliers", {
  counts <- mock_density_counts()
  scdata_list <- mock_density_sample(counts)
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_count_density(scdata_list, mock_density_config(), "spatial1", cells_id)

  expected <- scuttle::isOutlier(counts / 100, log = TRUE, type = "lower", nmads = 3)
  expect_equal(sum(expected), 2)
  expect_equal(out$new_ids$spatial1, cells_id$spatial1[!expected])
})

test_that("density is per area: a small cell with few counts is kept", {
  counts <- c(mock_density_counts()[1:60], 50)
  area <- c(rep(100, 60), 10) # 5 counts/um^2, like the rest
  scdata_list <- mock_density_sample(counts, area)
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_count_density(scdata_list, mock_density_config(), "spatial1", cells_id)

  expect_true(cells_id$spatial1[61] %in% out$new_ids$spatial1)
})

test_that("removes cells with control counts above the maximum when enabled", {
  counts <- mock_density_counts()[1:60]
  controls <- c(1, 2, rep(0, 58))
  scdata_list <- mock_density_sample(counts, controls = controls)
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_count_density(scdata_list, mock_density_config(), "spatial1", cells_id)
  expect_equal(out$new_ids$spatial1, cells_id$spatial1[-(1:2)])

  out <- filter_spatial_count_density(
    scdata_list, mock_density_config(max_control_counts = 1), "spatial1", cells_id
  )
  expect_equal(out$new_ids$spatial1, cells_id$spatial1[-2])

  out <- filter_spatial_count_density(
    scdata_list, mock_density_config(filter_controls = FALSE), "spatial1", cells_id
  )
  expect_equal(out$new_ids$spatial1, cells_id$spatial1)
})

test_that("nmads widens or narrows the density cut", {
  counts <- mock_density_counts()
  scdata_list <- mock_density_sample(counts)
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_count_density(
    scdata_list, mock_density_config(nmads = 100), "spatial1", cells_id
  )
  expect_equal(out$new_ids$spatial1, cells_id$spatial1)
  expect_equal(out$config$filterSettings$nmads, 100)
})

test_that("auto uses 3 MADs and removes controls above 0", {
  counts <- mock_density_counts()[1:60]
  scdata_list <- mock_density_sample(counts, controls = c(1, rep(0, 59)))
  cells_id <- mock_ids(scdata_list)

  config <- mock_density_config(
    auto = TRUE, nmads = 100, filter_controls = FALSE, max_control_counts = 10
  )
  out <- filter_spatial_count_density(scdata_list, config, "spatial1", cells_id)

  expect_equal(out$config$filterSettings$nmads, 3)
  expect_true(out$config$filterSettings$filterControls)
  expect_equal(out$config$filterSettings$maxControlCounts, 0)
  expect_false(cells_id$spatial1[1] %in% out$new_ids$spatial1)
})

test_that("zero-count cells are outliers, cells without an area are kept", {
  counts <- c(mock_density_counts()[1:60], 0, 500)
  area <- c(rep(100, 61), NA)
  scdata_list <- mock_density_sample(counts, area)
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_count_density(scdata_list, mock_density_config(), "spatial1", cells_id)

  expect_false(cells_id$spatial1[61] %in% out$new_ids$spatial1)
  expect_true(cells_id$spatial1[62] %in% out$new_ids$spatial1)
})

test_that("disabled filter keeps every cell", {
  scdata_list <- mock_density_sample(mock_density_counts())
  config <- mock_density_config()
  config$enabled <- FALSE

  out <- filter_spatial_count_density(scdata_list, config, "spatial1", mock_ids(scdata_list))

  expect_equal(length(out$new_ids$spatial1), 62)
})

test_that("plotData has per-cell densities and the stats to recompute the threshold", {
  counts <- mock_density_counts()
  scdata_list <- mock_density_sample(counts, controls = c(3, rep(0, 61)))
  out <- filter_spatial_count_density(
    scdata_list, mock_density_config(), "spatial1", mock_ids(scdata_list),
    num_cells_to_downsample = 10
  )

  plot0 <- out$plotData[[generate_gui_uuid("spatial1", "spatialCountDensity", 0)]]
  plot1 <- out$plotData[[generate_gui_uuid("spatial1", "spatialCountDensity", 1)]]
  stats <- out$plotData[[generate_gui_uuid("spatial1", "spatialCountDensity", 2)]]

  expect_equal(length(plot0), 62)
  expect_setequal(names(plot0[[1]]), c("cellId", "value", "controls", "x", "y"))
  expect_equal(plot0[[1]]$value, log(counts[1] / 100))
  expect_equal(plot0[[1]]$controls, 3)

  expect_equal(length(plot1$values), 10)
  expect_equal(plot1$median, stats::median(log(counts / 100)))
  expect_equal(plot1$mad, stats::mad(log(counts / 100)))

  expect_setequal(names(stats), c("before", "after"))
})

test_that("missing cell areas are a no-op with a warning message", {
  scdata_list <- mock_spatial_scdata_list(ncells = 64)
  cells_id <- mock_ids(scdata_list)

  expect_message(
    out <- filter_spatial_count_density(scdata_list, mock_density_config(), "spatial1", cells_id),
    "No cell areas"
  )
  expect_identical(out$new_ids, cells_id)
  expect_equal(length(out$plotData), 0)
})

test_that("a NULL config (step not configured, e.g. single-cell data) is a no-op", {
  scdata_list <- mock_scdata_list()
  cells_id <- mock_ids(scdata_list)

  out <- filter_spatial_count_density(scdata_list, NULL, names(scdata_list)[1], cells_id)

  expect_identical(out$new_ids, cells_id)
  expect_equal(length(out$plotData), 0)
})
