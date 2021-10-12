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
        pbmc_raw[1:10, 1] <- 20
    }

    scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

    # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    return(scdata)
}


test_that("filter_gene_umi_outlier updates linear p.level in config if auto", {
    scdata <- mock_scdata()
    config <- mock_config()
    config$filterSettings$regressionTypeSettings$linear$p.level <- 1
    out <- filter_gene_umi_outlier(scdata, config, '123def')
    new <- out$config$filterSettings$regressionTypeSettings$linear$p.level

    expect_lt(new, 1)
})


test_that("filter_gene_umi_outlier is sample specific", {

    # one outlier in first sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()

    barcodes1 <- colnames(scdata)[scdata$samples == '123abc']
    barcodes2 <- colnames(scdata)[scdata$samples == '123def']

    out1 <- filter_gene_umi_outlier(scdata, config, '123abc')
    out2 <- filter_gene_umi_outlier(scdata, config, '123def')

    # filtered something
    expect_lt(ncol(out1$data), ncol(scdata))
    expect_lt(ncol(out2$data), ncol(scdata))

    # didn't filter other sample
    expect_true(all(barcodes1 %in% colnames(out2$data)))
    expect_true(all(barcodes2 %in% colnames(out1$data)))

})

test_that("filter_gene_umi_outlier can filter cells", {

    # single outlier in single sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    nstart <- ncol(scdata)
    scdata$samples <- '123abc'

    out <- filter_gene_umi_outlier(scdata, config, '123abc')
    expect_lt(ncol(out$data), nstart)
})

test_that("filter_gene_umi_outlier uses linear if regressionType is gam (legacy)", {

    # single outlier in single sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    nstart <- ncol(scdata)
    scdata$samples <- '123abc'

    out1 <- filter_gene_umi_outlier(scdata, config, '123abc')
    expect_lt(ncol(out1$data), nstart)


    config$filterSettings$regressionType <- 'gam'
    out2 <- filter_gene_umi_outlier(scdata, config, '123abc')
    expect_identical(out1$data, out2$data)
})

test_that("filter_gene_umi_outlier can be disabled", {

    # single outlier in single sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    nstart <- ncol(scdata)
    scdata$samples <- '123abc'

    # disabled
    config$enabled <- FALSE

    out <- filter_gene_umi_outlier(scdata, config, '123abc')
    expect_equal(ncol(out$data), nstart)
})



test_that("filter_gene_umi_outlier return the appropriate plot data", {
    scdata <- mock_scdata()
    config <- mock_config()

    out <- filter_gene_umi_outlier(scdata, config, '123def')
    plot_data <- out$plotData[[1]]

    # points and lines data
    expect_equal(names(plot_data), c('pointsData', 'linesData'))

    # one value per cell
    points_data <- plot_data$pointsData
    lines_data <- plot_data$linesData
    expect_equal(length(points_data), sum(scdata$samples == '123def'))

    # has the right names
    expect_equal(names(points_data[[1]]), c('log_molecules', 'log_genes'))
    expect_equal(names(lines_data[[1]]), c('log_molecules', 'lower_cutoff', 'upper_cutoff'))
})


test_that("filter_gene_umi_outlier skips if no barcodes for this sample", {
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    out <- filter_gene_umi_outlier(scdata, config, 'not a sample')

    expect_equal(ncol(out$data), ncol(scdata))
})


test_that("filter_gene_umi_outlier gives different results with spline fit", {
    # single outlier in single sample
    scdata <- mock_scdata(with_outlier = TRUE)
    config <- mock_config()
    nstart <- ncol(scdata)
    scdata$samples <- '123abc'

    out1 <- filter_gene_umi_outlier(scdata, config, '123abc')
    expect_lt(ncol(out1$data), nstart)

    config$filterSettings$regressionType <- 'spline'

    out2 <- filter_gene_umi_outlier(scdata, config, '123abc')

    expect_true(!identical(out1$plotData, out2$plotData))
    expect_true(!identical(out1$data, out2$data))
})

