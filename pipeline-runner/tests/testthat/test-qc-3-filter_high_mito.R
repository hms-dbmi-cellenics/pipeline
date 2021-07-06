mock_config <- function(max_fraction = 0.1) {
    config <- list(
        auto = TRUE,
        enabled = TRUE,
        filterSettings = list(
            method = 'absolute_threshold',
            methodSettings = list(
                absolute_threshold = list(maxFraction = max_fraction)
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

    # add mitochondrial percent
    scdata@meta.data$percent.mt <- 5
    scdata@meta.data$percent.mt[1:10] <- 11

    # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ''))
    return(scdata)
}



test_that("filter_high_mito filters based on threshold", {
    scdata <- mock_scdata()

    # should filter first 10 cells
    config <- mock_config(0.1)
    config$auto <- FALSE
    out <- filter_high_mito(scdata, config, '123abc')
    expect_equal(ncol(out$data), 70)

    # should filter all cells in first sample
    config <- mock_config(0.01)
    config$auto <- FALSE

    out <- filter_high_mito(scdata, config, '123abc')
    expect_equal(ncol(out$data), 40)
})

test_that("filter_high_mito can be disabled", {
    scdata <- mock_scdata()
    nstart <- ncol(scdata)

    # should filter first 10 cells
    config <- mock_config(0.1)
    config$enabled <- FALSE
    out <- filter_high_mito(scdata, config, '123abc')
    expect_equal(ncol(out$data), nstart)
})

test_that("filter_high_mito is sample specific", {
    scdata <- mock_scdata()
    nstart <- ncol(scdata)

    # should filter first 10 cells of 123abc only
    config <- mock_config(0.1)
    config$auto <- FALSE

    out1 <- filter_high_mito(scdata, config, '123abc')
    out2 <- filter_high_mito(scdata, config, '123def')

    expect_equal(ncol(out1$data), 70)
    expect_equal(ncol(out2$data), 80)
})

test_that("filter_high_mito can be set to auto", {
    scdata <- mock_scdata()
    nstart <- ncol(scdata)

    # would filter all 40 cells in first sample unless auto
    config <- mock_config(0.01)
    out <- filter_high_mito(scdata, config, '123abc')
    expect_equal(ncol(out$data), 70)
})


