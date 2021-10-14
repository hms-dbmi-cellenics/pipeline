mock_counts <- function() {
    read.table(
        file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
        as.is = TRUE
    )
}

mock_doublet_scores <- function(counts) {

    doublet_scores <- runif(ncol(counts))
    doublet_class <- ifelse(doublet_scores < 0.8, 'singlet', 'doublet')

    data.frame(
        row.names = colnames(counts),
        barcodes = colnames(counts),
        doublet_class = doublet_class,
        doublet_scores = doublet_scores
    )
}

mock_prev_out <- function(samples = 'sample_a', counts = NULL) {

    if (is.null(counts)) {
        counts <- DropletUtils:::simCounts()
        colnames(counts) <- paste0('cell', seq_len(ncol(counts)))
    }

    eout <- DropletUtils::emptyDrops(counts)

    counts_list <- list()
    edrops <- list()
    doublet_scores <- list()

    for (sample in samples) {
        counts_list[[sample]] <- counts
        edrops[[sample]] <- eout
        doublet_scores[[sample]] <- mock_doublet_scores(counts)
    }

    list(
        counts_list = counts_list,
        edrops = edrops,
        doublet_scores = doublet_scores,
        annot = data.frame(name = row.names(counts), input = row.names(counts)),
        config = list(name = 'project name')
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
    valid_names <- make.names(colnames(metadata), unique = TRUE)

    # changes column names of metadata
    expect_equal(colnames(metadata), valid_names)

    # stores correct values in metadata
    expect_true(all(metadata$TRUE. == 'a'))
    expect_true(all(metadata$X1first == 'b'))
    expect_true(all(metadata$with.space == 'c'))
    expect_true(all(metadata$with.dash == 'd'))
})



test_that("create_seurat works without emptyDrops result", {

    prev_out <- mock_prev_out()
    prev_out$edrops$sample_a <- NULL

    expect_message(create_seurat(NULL, NULL, prev_out), '^emptyDrops .+ skipping')
})

test_that("create_seurat fails with missing items in prev_out", {

    prev_out <- mock_prev_out()

    for (item in names(prev_out)) {
        obj <- prev_out[[item]]

        prev_out[[item]] <- NULL
        expect_error(create_seurat(NULL, NULL, prev_out))

        prev_out[[item]] <- obj
    }
})

test_that("create_seurat adds mitochondrial percentage", {

    # without mitcondrial genes - gets set to 0
    prev_out <- mock_prev_out()
    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_true(all(scdata$percent.mt == 0))

    # with mitochondrial genes - some not zero
    counts <- DropletUtils:::simCounts()
    row.names(counts)[1:10] <- paste0('mt-', 1:10)
    colnames(counts) <- paste0('cell', seq_len(ncol(counts)))
    prev_out <- mock_prev_out(counts = counts)

    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_false(all(scdata$percent.mt == 0))

    # case insensitive
    row.names(counts)[1:10] <- paste0('MT-', 1:10)
    prev_out <- mock_prev_out(counts = counts)

    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_false(all(scdata$percent.mt == 0))
})

test_that("create_seurat adds emptyDrops_FDR to SeuratObject", {

    prev_out <- mock_prev_out()
    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_true('emptyDrops_FDR' %in% colnames(scdata@meta.data))
})

test_that("create_seurat adds doublet_scores and doublet_class to SeuratObject", {

    prev_out <- mock_prev_out()
    out <- create_seurat(NULL, NULL, prev_out)$output
    scdata <- out$scdata_list[[1]]

    expect_true(all(c('doublet_scores', 'doublet_class') %in% colnames(scdata@meta.data)))
})

test_that("create_seurat works with multiple samples", {

    prev_out <- mock_prev_out(samples = c('a', 'b'))
    scdata_list <- create_seurat(NULL, NULL, prev_out)$output$scdata_list

    # scdata_list names are right
    expect_true(all(c('a', 'b') %in% names(scdata_list)))

    # samples column added to SeuratObjects
    expect_true(all(scdata_list[['a']]$samples == 'a'))
    expect_true(all(scdata_list[['b']]$samples == 'b'))
})
