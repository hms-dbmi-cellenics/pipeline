mock_config <- function(mcs=100,auto_settings=FALSE,enabled_on=TRUE) {
    config <- list(
        auto = auto_settings,
        enabled = enabled_on,
        filterSettings = list(
            minCellSize = mcs
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


test_that("filter_low_cellsize removes cells", {
    scdata <- mock_scdata()
    config <- mock_config(mcs=10000)

    out<-filter_low_cellsize(scdata,config,"123abc")
    expect_lt(ncol(out$data),ncol(scdata))
})

test_that("filter_low_cellsize filters only appropiate cells", {
    mcs<-100
    scdata <- mock_scdata()
    config <- mock_config(mcs)
    out<-filter_low_cellsize(scdata,config,"123abc")

    data<- out$data
    barcode_names_this_sample <- rownames(data@meta.data[grep("123abc", rownames(data@meta.data)), ])
    sample_subset <- subset(data, cells = barcode_names_this_sample)
    expect_false(all(sample_subset@meta.data$nCount_RNA <= mcs))
})

test_that("filter_low_cellsize is sample aware", {
    mcs<-10000
    scdata <- mock_scdata()
    config <- mock_config(mcs)
    out<-filter_low_cellsize(scdata,config,"123abc")

    data<- out$data
    barcode_names_this_sample <- rownames(data@meta.data[grep("123abc", rownames(data@meta.data)), ])

    expect_lt(length(barcode_names_this_sample),40)

    barcode_names_this_sample <- rownames(data@meta.data[grep("123def", rownames(data@meta.data)), ])
    expect_equal(length(barcode_names_this_sample),40)

    out<-filter_low_cellsize(scdata,config,"123def")
    data<- out$data

    barcode_names_this_sample <- rownames(data@meta.data[grep("123def", rownames(data@meta.data)), ])
    expect_lt(length(barcode_names_this_sample),40)
})

test_that("filter_low_cellsize works on empty data/wrong sample", {
    scdata <- mock_scdata()
    config <- mock_config(100000)

    out<-filter_low_cellsize(scdata,config,"123abc")
    out<-filter_low_cellsize(out$data,config,"123def")

    expect_equal(1,ncol(out$data))
    out<-filter_low_cellsize(out$data,config,"123abc")
    expect_equal(1,ncol(out$data))
})

test_that("filter_low_cellsize works with auto", {
    scdata <- mock_scdata()
    config <- mock_config(auto_settings =TRUE)
    out<-filter_low_cellsize(scdata,config,"123def")

    expect_false(config$filterSettings$minCellSize == out$config$filterSettings$minCellSize)

    data<- out$data
    barcode_names_this_sample <- rownames(data@meta.data[grep("123def", rownames(data@meta.data)), ])
    sample_subset <- subset(data, cells = barcode_names_this_sample)

    expect_true(all(sample_subset@meta.data$nCount_RNA >= out$config$filterSettings$minCellSize))
})


test_that("filter_low_cellsize can be disabled", {
    scdata <- mock_scdata()
    config <- mock_config(mcs=10000,enabled=FALSE)
    out<-filter_low_cellsize(scdata,config,"123def")

    expect_equal(ncol(out$data),80)
})

