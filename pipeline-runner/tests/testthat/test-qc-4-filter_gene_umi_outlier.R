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
    scdata$cells_id <- 0:(ncol(scdata) - 1)

    # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ''))
    return(scdata)
}


test_that("filter_gene_umi_outlier updates gam p.level in config if auto", {
    scdata <- mock_scdata()
    cells_id <- mock_ids()
    config <- mock_config()

    config$filterSettings$regressionTypeSettings$gam$p.level <- 1
    out <- filter_gene_umi_outlier(scdata, config, '123def',cells_id)
    new <- out$config$filterSettings$regressionTypeSettings$gam$p.level

    expect_lt(new, 1)
})

test_that("filter_gene_umi_outlier can filter cells", {
    scdata <- mock_scdata()
    cells_id <- mock_ids()
    config <- mock_config()
    nstart <- ncol(scdata)

    out <- filter_gene_umi_outlier(scdata, config, '123def',cells_id)
    expect_equal(ncol(out$data), nstart)
    expect_lt(length(out$new_ids$`123def`),length(cells_id$`123def`))
})

test_that("filter_gene_umi_outlier is sample specific", {
    scdata <- mock_scdata()
    cells_id <- mock_ids()
    config <- mock_config()

    out1 <- filter_gene_umi_outlier(scdata, config, '123abc',cells_id)

    expect_equal(length(out1$new_ids$`123def`),length(cells_id$`123def`))

    out2 <- filter_gene_umi_outlier(scdata, config, '123def',cells_id)

    expect_equal(length(out2$new_ids$`123abc`),length(cells_id$`123abc`))
    expect_lt(length(out2$new_ids$`123def`),length(cells_id$`123def`))

    expect_lt(length(out2$new_ids$`123def`), length(out1$new_ids$`123abc`))
})

test_that("filter_gene_umi_outlier can be disabled", {
    scdata <- mock_scdata()
    cells_id <- mock_ids()
    config <- mock_config()
    nstart <- ncol(scdata)
    config$enabled <- FALSE

    out <- filter_gene_umi_outlier(scdata, config, '123def',cells_id)
    expect_equal(out$new_ids$`123abc`,0:39)
    expect_equal(out$new_ids$`123def`,40:79)
    expect_equal(ncol(out$data), nstart)
})



test_that("filter_gene_umi_outlier return the appropriate plot data", {
    scdata <- mock_scdata()
    cells_id <- mock_ids()
    config <- mock_config()

    out <- filter_gene_umi_outlier(scdata, config, '123def',cells_id)
    pdat <- out$plotData[[1]]

    # one value per cell
    expect_equal(length(pdat), sum(scdata$samples == '123def'))

    # has the right names
    expected_names <- c('log_molecules', 'log_genes', 'lower_cutoff', 'upper_cutoff' )
    expect_equal(names(pdat[[1]]), expected_names)

    # is numeric
    expect_equal(class(pdat[[1]]), 'numeric')
})

test_that("Gene UMI filter works if input is a a float-interpretable string", {
    scdata <- mock_scdata()
    config <- mock_config()
    nstart <- ncol(scdata)
    config$auto <- FALSE
    out_number <- filter_gene_umi_outlier(scdata, config, '123def')

    config$filterSettings$regressionTypeSettings$gam$p.level <- "0.1"
    out_string <- filter_gene_umi_outlier(scdata, config, '123def')

    expect_lt(ncol(out_number$data),nstart)
    expect_equal(ncol(out_number$data),ncol(out_string$data))
})

test_that("Gene UMI filter throws error if input is a non float-interpretable string", {
    scdata <- mock_scdata()
    config <- mock_config()
    config$auto <- FALSE
    config$filterSettings$regressionTypeSettings$gam$p.level <- "asd"

    expect_error(filter_gene_umi_outlier(scdata, config, '123def'))
})
