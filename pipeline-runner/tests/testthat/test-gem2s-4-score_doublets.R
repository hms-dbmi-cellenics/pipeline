mock_counts <- function(...) {
    set.seed(0)
    sce <- scDblFinder::mockDoubletSCE(...)
    sce@assays@data$counts
}


test_that("score_doublets returns expected columns", {
    counts <- mock_counts()

    prev_out <- list(counts_list = list(sample1=counts),
                     config = list(),
                     annot = list())
    out <- score_doublets(NULL, NULL, prev_out)$output

    expect_setequal(colnames(out$doublet_scores$sample1),
                    c('barcodes', 'doublet_class', 'doublet_scores'))
})


test_that("score_doublets filters cells to avoid warning of extremely low read counts", {
    counts <- mock_counts(ncells = c(200, 300, 400, 200, 500, 300))
    counts <- round(counts/2)

    expect_warning(scDblFinder:::.checkSCE(counts), 'extremely low read counts')

    prev_out <- list(counts_list = list(sample1=counts),
                     config = list(),
                     annot = list())
    expect_warning(score_doublets(NULL, NULL, prev_out)$output, NA)
})
