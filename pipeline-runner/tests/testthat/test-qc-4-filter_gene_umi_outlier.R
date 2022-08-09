mock_ids <- function() {
  return(list("123abc" = 0:39, "123def" = 40:79))
}

mock_config <- function() {
    config <- list(
        auto = TRUE,
        enabled = TRUE,
        filterSettings = list(
            regressionType = 'linear',
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
        as.is = TRUE)

    if (with_outlier) {
        pbmc_raw[, 1] <- 0
        pbmc_raw[1:10, 1] <- 100
    }

    scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
    scdata$cells_id <- 0:(ncol(scdata)-1)

    # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    return(scdata)
}

test_that("filter_gene_umi_outlier updates linear p.level in config if auto", {
    scdata <- mock_scdata()
    config <- mock_config()
    cells_id <- mock_ids()
    config$filterSettings$regressionTypeSettings$linear$p.level <- 1
    out <- filter_gene_umi_outlier(scdata, config, '123def', cells_id)
    new <- out$config$filterSettings$regressionTypeSettings$linear$p.level

  expect_lt(new, 1)
})


test_that("filter_gene_umi_outlier is sample specific", {

    # one outlier in first sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    cells_id <- mock_ids()

    cell_ids1 <- scdata$cells_id[scdata$samples == '123abc']
    cell_ids2 <- scdata$cells_id[scdata$samples == '123def']

    out1 <- filter_gene_umi_outlier(scdata, config, '123abc', cells_id)
    out2 <- filter_gene_umi_outlier(scdata, config, '123def', cells_id)

    # filtered something
    expect_lt(length(out1$new_ids$`123abc`), length(cells_id$`123abc`))
    expect_lt(length(out2$new_ids$`123def`), length(cells_id$`123def`))

    # didn't filter other sample
    expect_true(all(cell_ids1 %in% out2$new_ids$`123abc`))
    expect_true(all(cell_ids2 %in% out1$new_ids$`123def`))
})

test_that("filter_gene_umi_outlier can filter cells", {

    # single outlier in single sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    cells_id <- mock_ids()
    nstart <- ncol(scdata)
    scdata$samples <- '123abc'

    out <- filter_gene_umi_outlier(scdata, config, '123abc', cells_id)
    expect_lt(length(unlist(out$new_ids)), nstart)
})

test_that("filter_gene_umi_outlier uses linear if regressionType is gam (legacy)", {

    # single outlier in single sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    cells_id <- mock_ids()
    nstart <- ncol(scdata)
    scdata$samples <- '123abc'

    out1 <- filter_gene_umi_outlier(scdata, config, '123abc', cells_id)
    expect_lt(length(unlist(out1$new_ids)), nstart)


    config$filterSettings$regressionType <- 'gam'
    out2 <- filter_gene_umi_outlier(scdata, config, '123abc', cells_id)
    expect_identical(out1$new_ids, out2$new_ids)
})

test_that("filter_gene_umi_outlier can be disabled", {

    # single outlier in single sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    nstart <- ncol(scdata)
    cells_id <- mock_ids()

    # disabled
    config$enabled <- FALSE

    out <- filter_gene_umi_outlier(scdata, config, '123abc',cells_id)
    expect_equal(ncol(out$data), nstart)
    expect_equal(out$new_ids$`123abc`,0:39)
    expect_equal(out$new_ids$`123def`,40:79)
})



test_that("filter_gene_umi_outlier return the appropriate plot data", {
  scdata <- mock_scdata()
  cells_id <- mock_ids()
  config <- mock_config()

    out <- filter_gene_umi_outlier(scdata, config, '123def',cells_id)
    plot_data <- out$plotData[[1]]

    # points and lines data
    expect_equal(names(plot_data), c('pointsData', 'linesData'))

    # one value per cell
    points_data <- plot_data$pointsData
    lines_data <- plot_data$linesData
    expect_equal(length(points_data), sum(scdata$samples == '123def'))

   # has the right names
   expect_equal(names(points_data[[1]]), c('log_molecules', 'log_genes'))
   expect_equal(names(lines_data[[1]][[1]]), c('log_molecules', 'lower_cutoff', 'upper_cutoff'))

   # has the correct number of pre-caluclated values
   expect_equal(length(lines_data), (length(c(seq(0,0.99,0.01), 0.999, 0.9999, 0.99999, 0.999999)) + 1))
})


test_that("filter_gene_umi_outlier skips if no barcodes for this sample", {
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    cells_id <- mock_ids()
    out <- filter_gene_umi_outlier(scdata, config, 'not a sample',cells_id)

    expect_equal(ncol(out$data), ncol(scdata))
    expect_equal(out$new_ids$`123abc`,0:39)
    expect_equal(out$new_ids$`123def`,40:79)
})

test_that("filter_gene_umi_outlier gives different results with spline fit", {
  # single outlier in single sample
  scdata <- mock_scdata(with_outlier = TRUE)
  cells_id <- mock_ids()
  config <- mock_config()
  nstart <- ncol(scdata)
  scdata$samples <- '123abc'

  out1 <- filter_gene_umi_outlier(scdata, config, '123abc',cells_id)
  expect_equal(ncol(out1$data), nstart)
  expect_lt(length(out1$new_ids$`123abc`),length(cells_id$`123abc`))

  config$filterSettings$regressionType <- 'spline'

  out2 <- filter_gene_umi_outlier(scdata, config, '123abc',cells_id)

  expect_true(!identical(out1$plotData, out2$plotData))
  expect_true(!identical(out1$new_ids$`123abc`, out2$new_ids$`123abc`))
})

test_that("Gene UMI filter works if input is a a float-interpretable string", {
  scdata <- mock_scdata()
  cells_id <- mock_ids()
  config <- mock_config()
  nstart <- ncol(scdata)
  config$auto <- FALSE

  out_number <- filter_gene_umi_outlier(scdata, config, "123def", cells_id)

  config$filterSettings$regressionTypeSettings$linear$p.level <- "0.1"
  out_string <- filter_gene_umi_outlier(scdata, config, '123def',cells_id)

  expect_equal(ncol(out_number$data), nstart)
  expect_equal(out_string$new_ids$`123def`, out_number$new_ids$`123def`)
})

test_that("Gene UMI filter throws error if input is a non float-interpretable string", {
    scdata <- mock_scdata()
    config <- mock_config()
    config$auto <- FALSE
    config$filterSettings$regressionTypeSettings$linear$p.level <- "asd"

  expect_error(filter_gene_umi_outlier(scdata, config, "123def", cells_id))
})

