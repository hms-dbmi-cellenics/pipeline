mock_counts <- function() {
    read.table(
        file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
        as.is = TRUE
    )
}


test_that("construct_metadata works with syntactically invalid column names", {

    sample <- 'hello'
    counts <- mock_counts()

    config <- list(
        samples = 'hello',
        metadata = list('TRUE' = list('a'),
                        '1first' = list('b'),
                        'with space' = list('c'),
                        'with-dash' = list('d'))
        )

    metadata <- construct_metadata(counts, sample, config)

    # doesn't change column names of metadata
    expect_equal(colnames(metadata), c('samples', 'TRUE', '1first', 'with space', 'with-dash'))

    # stores correct values in metadata
    expect_true(all(metadata[['TRUE']] == 'a'))
    expect_true(all(metadata[['1first']] == 'b'))
    expect_true(all(metadata[['with space']] == 'c'))
    expect_true(all(metadata[['with-dash']] == 'd'))
})


test_that("result of get_metadata_lookups can be used to access syntactically invalid slots", {

    sample <- 'hello'
    counts <- mock_counts()

    config <- list(
        samples = 'hello',
        metadata = list('TRUE' = list('a'),
                        '1first' = list('b'),
                        'with space' = list('c'),
                        'with-dash' = list('d'))
    )

    metadata <- construct_metadata(counts, sample, config)
    lookups <- get_metadata_lookups(metadata)

    # seurat stores as syntactically valid names
    scdata <- SeuratObject::CreateSeuratObject(counts, meta.data = metadata)
    expect_true(all(lookups %in% colnames(scdata@meta.data)))

    # use lookups to access
    keys <- c('TRUE', '1first', 'with space', 'with-dash')
    values <- c('a', 'b', 'c', 'd')

    for (i in seq_along(keys)) {
        key <- keys[i]
        val <- values[i]
        seurat_key <- lookups[[key]]

        expect_true(seurat_key %in% colnames(scdata@meta.data))
        expect_true(all(scdata[[seurat_key]] == val))
    }
})
