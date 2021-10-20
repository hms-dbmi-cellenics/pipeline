mock_scdata <- function() {
    pbmc_raw <- read.table(
        file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
        as.is = TRUE)

    scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

        # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    return(scdata)
}


test_that("harmony integration works", {
    scdata <- mock_scdata()
    config <- list(dimensionalityReduction = list(numPCs = 2),
                   dataIntegration = list(method = 'harmony', methodSettings = list(harmony = list(numGenes = 10, normalisation = 'logNormalize'))))

    scdata <- run_dataIntegration(scdata, config)
    expect_s4_class(scdata, 'Seurat')
})
