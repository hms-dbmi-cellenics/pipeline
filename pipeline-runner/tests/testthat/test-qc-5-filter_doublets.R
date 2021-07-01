mock_config <- function(thr = 0.1,auto=FALSE,enabled=TRUE) {
    config <- list(
        auto = auto,
        enabled = enabled,
        filterSettings = list(
            probabilityThreshold = thr
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
    scdata@meta.data$doublet_scores <- 0.01
    scdata@meta.data$doublet_scores[1:10] <- 0.9

    scdata@meta.data$doublet_class <- rep("singlet",80)
    scdata@meta.data$doublet_class[1:10] <- rep("doublet",10)

    # add samples
    scdata$samples <- rep(c('123abc', '123def'), each = 40)
    scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ''))
    return(scdata)
}



test_that("filter_doublets filters based on threshold", {
    scdata <- mock_scdata()
    # should filter first 10 cells
    config <- mock_config(0.5)
    out <- filter_doublets(scdata, config, '123abc')
    expect_equal(ncol(out$data),70)
})

test_that("filter_doublets is sample aware", {
    scdata <- mock_scdata()
    scdata@meta.data$doublet_scores[71:80] <- 0.9
    config <- mock_config(0.5)

    out <- filter_doublets(scdata, config, '123abc')
    expect_equal(ncol(out$data),70)

    out<-filter_doublets(out$data,config,"123def")
    expect_equal(ncol(out$data),60)
})

test_that("filter_doublets filters works with auto", {
    scdata <- mock_scdata()
    # should filter first 10 cells
    config <- mock_config(0.001,auto=TRUE)
    out <- filter_doublets(scdata, config, '123abc')
    expect_equal(ncol(out$data),70)

    config <- mock_config(0.001,auto=FALSE)
    out <- filter_doublets(scdata, config, '123abc')
    expect_equal(ncol(out$data),40)

})

test_that("filter_doublets can be disabled", {
    scdata <- mock_scdata()

    config <- mock_config(0.5,enabled=FALSE)
    out <- filter_doublets(scdata, config, '123abc')
    expect_equal(ncol(out$data),80)

    config <- mock_config(0.5,enabled=TRUE)
    out <- filter_doublets(scdata, config, '123abc')
    expect_equal(ncol(out$data),70)
})
