
mock_cellranger_files <- function(counts, features, sample_dir) {

    # save features
    features_path <- file.path(sample_dir, 'features.tsv')
    write.table(features, features_path, col.names = FALSE, quote = FALSE, sep = '\t', row.names = FALSE)
    R.utils::gzip(features_path)

    # save barcodes
    barcodes <- colnames(counts)
    barcodes_path <- file.path(sample_dir, 'barcodes.tsv')
    writeLines(barcodes, barcodes_path)
    R.utils::gzip(barcodes_path)

    # save Matrix
    if (is(counts, 'data.frame')) counts <- as.matrix(counts)
    sparse.mat <- Matrix::Matrix(counts, sparse = TRUE)
    matrix_path <- file.path(sample_dir,  'matrix.mtx')
    Matrix::writeMM(sparse.mat, matrix_path)
    R.utils::gzip(matrix_path)
}


mock_counts <- function() {
    pbmc_raw <- read.table(
        file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
        as.is = TRUE)
}


test_that("format_annot keeps unique rows", {

    annot_list <- list(
        sample1 = data.frame(ENSID = 1:5, SYMBOL = paste0('gene', 1:5)),
        sample2 = data.frame(ENSID = 1:5, SYMBOL = paste0('gene', 1:5))
    )

    annot <- pipeline:::format_annot(annot_list)

    expect_s3_class(annot, 'data.frame')
    expect_true(nrow(annot) == nrow(annot_list$sample1))
})

test_that("format_annot deduplicates name column", {

    annot_list <- list(
        sample1 = data.frame(ENSID = 1:6, SYMBOL = paste0('gene', c(1, 1:5)))
    )

    annot <- pipeline:::format_annot(annot_list)

    expect_true(length(annot$name) == length(unique(annot$name)))
})


test_that("load_cellranger loads a count matrix", {

    counts <- mock_counts()
    features <- data.frame(ensid = paste0('ENSFAKE', seq_len(nrow(counts))),
                           symbol = row.names(counts))

    outdir <- tempdir()
    sample <- 'sample_a'
    sample_dir <- file.path(outdir, sample)
    dir.create(sample_dir)

    mock_cellranger_files(counts, features, sample_dir)

    prev_out <- list(config = list(samples = sample))
    out <- load_cellranger_files(NULL, NULL, prev_out, outdir)$output

    expect_true('counts_list' %in% names(out))
    expect_true(sample %in% names(out$counts_list))

    expect_s4_class(out$counts_list[[1]], 'dgCMatrix')
    unlink(sample_dir, recursive = TRUE)
})


test_that("load_cellranger generates feature annotation", {

    counts <- mock_counts()
    features <- data.frame(ensid = paste0('ENSFAKE', seq_len(nrow(counts))),
                           symbol = row.names(counts))

    outdir <- tempdir()
    sample <- 'sample_a'
    sample_dir <- file.path(outdir, sample)
    dir.create(sample_dir)

    mock_cellranger_files(counts, features, sample_dir)

    prev_out <- list(config = list(samples = sample))
    out <- load_cellranger_files(NULL, NULL, prev_out, outdir)$output

    expect_true('annot' %in% names(out))
    expect_true(
        all(c('input', 'name', 'original_name') %in% colnames(out$annot))
    )

    unlink(sample_dir, recursive = TRUE)
})


test_that("load_cellranger deduplicates gene symbols", {

    counts <- mock_counts()

    symbols <- row.names(counts)
    symbols[1:5] <- 'DUPLICATED'

    features <- data.frame(ensid = paste0('ENSFAKE', seq_len(nrow(counts))),
                           symbol = symbols)

    outdir <- tempdir()
    sample <- 'sample_a'
    sample_dir <- file.path(outdir, sample)
    dir.create(sample_dir)

    mock_cellranger_files(counts, features, sample_dir)

    prev_out <- list(config = list(samples = sample))
    annot <- load_cellranger_files(NULL, NULL, prev_out, outdir)$output$annot

    # unique gene names is same as number of gene names
    expect_length(unique(annot$name), length(symbols))

    # unique original names is same as unique gene names
    expect_length(unique(annot$original_name), length(unique(symbols)))
    unlink(sample_dir, recursive = TRUE)
})

test_that("load_cellranger uses appropriate feature columns", {

    counts <- mock_counts()

    symbols <- row.names(counts)
    features <- data.frame(ensid = paste0('ENSFAKE', seq_len(nrow(counts))),
                           symbol = symbols)

    outdir <- tempdir()
    sample <- 'sample_a'
    sample_dir <- file.path(outdir, sample)
    dir.create(sample_dir)

    mock_cellranger_files(counts, features, sample_dir)

    prev_out <- list(config = list(samples = sample))
    out <- load_cellranger_files(NULL, NULL, prev_out, outdir)$output

    # ensembl ids are counts row names
    expect_equal(
        features$ensid,
        row.names(out$counts_list[[1]])
    )

    # ensembl ids are in column 'input' of annot
    expect_equal(out$annot$input, features$ensid)

    # symbols are in column 'name' of annot
    expect_equal(out$annot$name, symbols)

    unlink(sample_dir, recursive = TRUE)
})


test_that("load_cellranger loads multisample experiments", {

    counts <- mock_counts()

    symbols <- row.names(counts)
    ensids <- paste0('ENSFAKE', seq_len(nrow(counts)))

    features <- data.frame(ensid = ensids, symbol = symbols)

    # most overlapping, some unique
    ensids2 <- ensids
    symbols2 <- symbols
    ensids2[1:20] <- paste0(ensids2[1:20], '_2')
    symbols2[1:20] <- paste0(symbols2[1:20], '_2')

    features2 <- data.frame(ensid = ensids2, symbols = symbols2)

    outdir <- tempdir()
    samples <- c('sample_a', 'samble_b')
    sample_dirs <- file.path(outdir, samples)
    dir.create(sample_dirs[1])
    dir.create(sample_dirs[2])

    mock_cellranger_files(counts, features, sample_dirs[1])
    mock_cellranger_files(counts, features2, sample_dirs[2])

    prev_out <- list(config = list(samples = samples))
    out <- load_cellranger_files(NULL, NULL, prev_out, outdir)$output

    # loaded both
    expect_equal(names(out$counts_list), samples)

    # kept only unique ensembl ids
    unique_ensids <- unique(c(ensids, ensids2))
    expect_length(out$annot$input, length(unique_ensids))

    # removed duplicated ensembl ids
    expect_lt(
        length(unique_ensids),
        length(c(ensids, ensids2))
    )

    unlink(sample_dirs, recursive = TRUE)
})
