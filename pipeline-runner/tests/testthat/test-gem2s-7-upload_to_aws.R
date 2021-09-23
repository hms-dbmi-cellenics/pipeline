mock_scdata <- function() {
    pbmc_raw <- read.table(
        file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
        as.is = TRUE)

    scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
    scdata$cells_id <- seq(0, ncol(scdata)-1)

    # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    return(scdata)
}


test_that("samples_sets adds cells to correct sample", {
    scdata <- mock_scdata()

    ids.123abc <- scdata$cells_id[scdata$samples == '123abc']
    ids.123def <- scdata$cells_id[scdata$samples == '123def']
    color_pool <- get_color_pool()

    sample_ids <- c('123abc', '123def')
    input <- list(sampleIds = sample_ids, sampleNames = sample_ids)

    res <- samples_sets(input, scdata, color_pool)

    # one child for each sample
    expect_length(res$children, 2)

    keys <- sapply(res$children, `[[`, 'key')

    key.123abc <- which(keys == sample_ids[1])
    key.123def <- which(keys == sample_ids[2])

    expect_equal(res$children[[key.123abc]]$cellIds, unname(ids.123abc))
    expect_equal(res$children[[key.123def]]$cellIds, unname(ids.123def))
})
