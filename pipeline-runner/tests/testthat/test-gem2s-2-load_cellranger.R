
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
    sparse.mat <- as(counts, 'dgCMatrix')
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

    annot <- format_annot(annot_list)

    expect_s3_class(annot, 'data.frame')
    expect_true(nrow(annot) == nrow(annot_list$sample1))
})

test_that("format_annot deduplicates name column", {

    annot_list <- list(
        sample1 = data.frame(ENSID = 1:6, SYMBOL = paste0('gene', c(1, 1:5)))
    )

    annot <- format_annot(annot_list)

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
    expect_true(all(c('input', 'name', 'original_name') %in% colnames(out$annot)))
})
