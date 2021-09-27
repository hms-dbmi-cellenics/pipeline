mock_config <- function() {
    config <- list(
        sampleNames = list('WT1', 'WT2'),
        sampleIds = list('123abc', '123def')
    )

    return(config)
}


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




test_that("get_cell_sets creates scratchpad and sample sets if no metadata", {
    scdata <- mock_scdata()
    config <- mock_config()

    cell_sets <- get_cell_sets(scdata, config)
    keys <- sapply(cell_sets$cellSets, `[[`, 'key')

    expect_setequal(keys, c('scratchpad', 'sample'))
})


test_that("get_cell_sets adds correct cell ids for each sample", {
    scdata <- mock_scdata()
    config <- mock_config()

    cell_sets <- get_cell_sets(scdata, config)
    sets_key <- sapply(cell_sets$cellSets, `[[`, 'key')

    sample_sets <- cell_sets$cellSets[[which(sets_key == 'sample')]]
    samples_key <- sapply(sample_sets$children, `[[`, 'key')

    for (sample_id in config$sampleIds) {
        sample_cells <- sample_sets$children[[which(samples_key == sample_id)]]$cellIds
        expected_cells <- unname(scdata$cells_id)[scdata$samples == sample_id]

        expect_equal(sample_cells, expected_cells)
    }
})


test_that("get_cell_sets adds correct cell ids for each sample", {
    scdata <- mock_scdata()
    config <- mock_config()

    cell_sets <- get_cell_sets(scdata, config)
    sets_key <- sapply(cell_sets$cellSets, `[[`, 'key')

    sample_sets <- cell_sets$cellSets[[which(sets_key == 'sample')]]
    samples_key <- sapply(sample_sets$children, `[[`, 'key')

    # ids are correct for each child
    for (sample_id in config$sampleIds) {
        sample_cells <- sample_sets$children[[which(samples_key == sample_id)]]$cellIds
        expected_cells <- unname(scdata$cells_id)[scdata$samples == sample_id]

        expect_equal(sample_cells, expected_cells)
    }
})


test_that("get_cell_sets without metadata matches initial snapshot", {
    scdata <- mock_scdata()
    config <- mock_config()

    cell_sets <- get_cell_sets(scdata, config)$cellSets[[1]]
    expect_equal(cell_sets$type, 'cellSets')
    expect_true(cell_sets$rootNode)

    expect_snapshot(cell_sets)
})


test_that("get_cell_sets adds a single metadata column", {
    scdata <- mock_scdata()
    config <- mock_config()

    # mockup metadata and lookups
    config$metadata <- list(Group = list('Hello', 'WT2'))

    scdata$Group <- NA
    scdata$Group[scdata$samples == '123abc'] <- 'Hello'
    scdata$Group[scdata$samples == '123def'] <- 'WT2'

    scdata@misc <- list(metadata_lookups = c(Group = 'Group'))

    cell_sets <- get_cell_sets(scdata, config)

    # have it as a key
    keys <- sapply(cell_sets$cellSets, `[[`, 'key')
    expect_setequal(keys, c('scratchpad', 'sample', 'Group'))

    group_set <- cell_sets$cellSets[[which(keys == 'Group')]]
    group_names <- sapply(group_set$children, `[[`, 'name')


    # cell ids are correct for each child
    for (group_name in group_names) {
        group_cells <- group_set$children[[which(group_names == group_name)]]$cellIds
        expected_cells <- unname(scdata$cells_id)[scdata$Group == group_name]

        expect_equal(group_cells, expected_cells)
    }
})

test_that("get_cell_sets adds two metadata columns", {
    scdata <- mock_scdata()
    config <- mock_config()

    # mockup metadata and lookups
    config$metadata <- list(Group1 = list('Hello', 'WT2'),
                            Group2 = list('WT', 'WT'))

    scdata$Group1 <- NA
    scdata$Group1[scdata$samples == '123abc'] <- 'Hello'
    scdata$Group1[scdata$samples == '123def'] <- 'WT2'
    scdata$Group2 <- 'WT'

    scdata@misc <- list(metadata_lookups = c(Group1 = 'Group1', Group2 = 'Group2'))
    cell_sets <- get_cell_sets(scdata, config)

    # have as keys
    keys <- sapply(cell_sets$cellSets, `[[`, 'key')
    expect_setequal(keys, c('scratchpad', 'sample', 'Group1', 'Group2'))

    # check that Group2 has all cells
    group2_set <- cell_sets$cellSets[[which(keys == 'Group2')]]
    group2_cells <- group2_set$children[[1]]$cellIds
    expect_equal(group2_cells, unname(scdata$cells_id))
})


test_that("get_cell_sets uses unique colors for each cell set", {
    scdata <- mock_scdata()
    config <- mock_config()

    # mockup metadata and lookups
    config$metadata <- list(Group1 = list('Hello', 'WT2'),
                            Group2 = list('WT', 'WT'))

    scdata$Group1 <- NA
    scdata$Group1[scdata$samples == '123abc'] <- 'Hello'
    scdata$Group1[scdata$samples == '123def'] <- 'WT2'
    scdata$Group2 <- 'WT'

    scdata@misc <- list(metadata_lookups = c(Group1 = 'Group1', Group2 = 'Group2'))
    cell_sets <- get_cell_sets(scdata, config)

    flat_cell_sets <- unlist(cell_sets)
    colors <- flat_cell_sets[grepl('[.]color', names(flat_cell_sets))]
    colors <- unname(colors)

    expect_equal(unique(colors), colors)
})
