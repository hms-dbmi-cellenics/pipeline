mock_ids <- function() {
  return(list("123abc" = 0:39, "123def" = 40:79))
}

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

mock_scdata <- function(with_outlier = FALSE) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  if (with_outlier) {
    pbmc_raw[, 1] <- 0
    pbmc_raw[1:10, 1] <- 100
  }

  sample_1_id <- "123abc"
  sample_2_id <- "123def"

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  # add samples
  scdata$samples <- rep(c(sample_1_id, sample_2_id), each = 40)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))
  scdata$emptyDrops_FDR <- NA

  scdata_sample1 <- subset(scdata, samples == sample_1_id)
  scdata_sample2 <- subset(scdata, samples == sample_2_id)

  scdata_list <- list(scdata_sample1, scdata_sample2)
  names(scdata_list) <- c(sample_1_id, sample_2_id)

  return(list(scdata_list, sample_1_id, sample_2_id))
}


test_that("filter_gene_umi_outlier updates linear p.level in config if auto", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config()
  cells_id <- mock_ids()

  config$filterSettings$regressionTypeSettings$linear$p.level <- 1
  out <- filter_gene_umi_outlier(scdata_list, config, sample_2_id, cells_id)
  new <- out$config$filterSettings$regressionTypeSettings$linear$p.level

  expect_lt(new, 1)
})


test_that("filter_gene_umi_outlier is sample specific", {

  # one outlier in first sample
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(with_outlier = T)
  config <- mock_config()
  cells_id <- mock_ids()

  cell_ids1 <- scdata_list[[sample_1_id]]$cells_id
  cell_ids2 <- scdata_list[[sample_2_id]]$cells_id
  out1 <- filter_gene_umi_outlier(scdata_list, config, sample_1_id, cells_id)
  out2 <- filter_gene_umi_outlier(scdata_list, config, sample_2_id, cells_id)

  # filtered something
  expect_lt(length(out1$new_ids[[sample_1_id]]), length(cells_id[[sample_1_id]]))
  expect_lt(length(out2$new_ids[[sample_2_id]]), length(cells_id[[sample_2_id]]))

  # didn't filter other sample
  expect_true(all(cell_ids1 %in% out2$new_ids[[sample_1_id]]))
  expect_true(all(cell_ids2 %in% out1$new_ids[[sample_2_id]]))
})

test_that("filter_gene_umi_outlier can filter cells", {
  # single outlier in single sample
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(with_outlier = T)
  config <- mock_config()
  cells_id <- mock_ids()
  nstart <- ncol(scdata_list[[sample_1_id]])

  out <-
    filter_gene_umi_outlier(scdata_list, config, sample_1_id, cells_id)
  expect_lt(length(unlist(out$new_ids[[sample_1_id]])), nstart)
})


test_that("filter_gene_umi_outlier can be disabled", {

  # single outlier in single sample
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(with_outlier = T)
  config <- mock_config()
  cells_id <- mock_ids()
  nstart <- ncol(scdata_list[[sample_1_id]])

  # disabled
  config$enabled <- FALSE

  out <- filter_gene_umi_outlier(scdata_list, config, sample_1_id, cells_id)
  expect_equal(ncol(out$data[[sample_1_id]]), nstart)
  expect_equal(out$new_ids[[sample_1_id]], 0:39)
  expect_equal(out$new_ids[[sample_2_id]], 40:79)
})



test_that("filter_gene_umi_outlier return the appropriate plot data", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config()
  cells_id <- mock_ids()

  out <- filter_gene_umi_outlier(scdata_list, config, sample_2_id, cells_id)
  plot_data <- out$plotData[[1]]

  # points and lines data
  expect_equal(names(plot_data), c("pointsData", "linesData"))

  # one value per cell
  points_data <- plot_data$pointsData
  lines_data <- plot_data$linesData
  expect_equal(length(points_data), ncol(scdata_list[[sample_2_id]]))

  # has the right names
  expect_equal(names(points_data[[1]]), c("log_molecules", "log_genes"))
  expect_equal(names(lines_data[[1]][[1]]), c("log_molecules", "lower_cutoff", "upper_cutoff"))
})



test_that("filter_gene_umi_outlier gives different results with spline fit", {
  # single outlier in single sample
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(with_outlier = T)
  config <- mock_config()
  cells_id <- mock_ids()
  nstart <- ncol(scdata_list[[sample_1_id]])

  out1 <- filter_gene_umi_outlier(scdata_list, config, sample_1_id, cells_id)
  expect_equal(ncol(out1$data[[sample_1_id]]), nstart)
  expect_lt(length(out1$new_ids[[sample_1_id]]), length(cells_id[[sample_1_id]]))

  config$filterSettings$regressionType <- "spline"

  out2 <- filter_gene_umi_outlier(scdata_list, config, sample_1_id, cells_id)

  expect_true(!identical(out1$plotData, out2$plotData))
  expect_true(!identical(out1$new_ids[[sample_1_id]], out2$new_ids[[sample_1_id]]))
})

test_that("Gene UMI filter works if input is a a float-interpretable string", {
  # single outlier in single sample
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata(with_outlier = T)
  config <- mock_config()
  cells_id <- mock_ids()
  nstart <- ncol(scdata_list[[sample_2_id]])
  config$auto <- FALSE

  out_number <- filter_gene_umi_outlier(scdata_list, config, sample_2_id, cells_id)

  config$filterSettings$regressionTypeSettings$linear$p.level <- "0.1"
  out_string <- filter_gene_umi_outlier(scdata_list, config, sample_2_id, cells_id)

  expect_equal(ncol(out_number$data[[sample_2_id]]), nstart)
  expect_equal(out_string$new_ids[[sample_2_id]], out_number$new_ids[[sample_2_id]])
})

test_that("Gene UMI filter throws error if input is a non float-interpretable string", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  config <- mock_config()
  config$auto <- FALSE
  config$filterSettings$regressionTypeSettings$linear$p.level <-
    "asd"

  expect_error(filter_gene_umi_outlier(scdata, config, sample_2_id, cells_id))
})
