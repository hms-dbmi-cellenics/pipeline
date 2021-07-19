mock_config <- function() {
    config <- list(
        auto = TRUE,
        enabled = TRUE,
        filterSettings = list(
            regressionType = 'gam',
            regressionTypeSettings = list(
                gam = list(p.level = 0.1)
            )
        )
    )

    return(config)
}


mock_scdata <- function() {
    pbmc_raw <- read.table(
        file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
        as.is = TRUE)

    scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

    # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ''))
    return(scdata)
}


test_that("filter_gene_umi_outlier updates gam p.level in config if auto", {
    scdata <- mock_scdata()
    config <- mock_config()
    config$filterSettings$regressionTypeSettings$gam$p.level <- 1
    out <- filter_gene_umi_outlier(scdata, config, '123def')
    new <- out$config$filterSettings$regressionTypeSettings$gam$p.level

    expect_lt(new, 1)
})

test_that("filter_gene_umi_outlier can filter cells", {
    scdata <- mock_scdata()
    config <- mock_config()
    nstart <- ncol(scdata)

    out <- filter_gene_umi_outlier(scdata, config, '123def')
    expect_lt(ncol(out$data), nstart)
})

test_that("filter_gene_umi_outlier is sample specific", {
    scdata <- mock_scdata()
    config <- mock_config()

    out1 <- filter_gene_umi_outlier(scdata, config, '123abc')
    out2 <- filter_gene_umi_outlier(scdata, config, '123def')
    expect_lt(ncol(out2$data), ncol(out1$data))
})

test_that("filter_gene_umi_outlier can be disabled", {
    scdata <- mock_scdata()
    config <- mock_config()
    nstart <- ncol(scdata)
    config$enabled <- FALSE

    out <- filter_gene_umi_outlier(scdata, config, '123def')
    expect_equal(ncol(out$data), nstart)
})



test_that("filter_gene_umi_outlier return the appropriate plot data", {
    scdata <- mock_scdata()
    config <- mock_config()

    out <- filter_gene_umi_outlier(scdata, config, '123def')
    pdat <- out$plotData[[1]]

    # one value per cell
    expect_equal(length(pdat), sum(scdata$samples == '123def'))

    # has the right names
    expected_names <- c('log_molecules', 'log_genes', 'lower_cutoff', 'upper_cutoff' )
    expect_equal(names(pdat[[1]]), expected_names)

    # is numeric
    expect_equal(class(pdat[[1]]), 'numeric')
})
