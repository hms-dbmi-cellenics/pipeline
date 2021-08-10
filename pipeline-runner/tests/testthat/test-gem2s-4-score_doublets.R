mock_counts <- function() {
    set.seed(0)
    sce <- scDblFinder::mockDoubletSCE()
    sce@assays@data$counts
}


test_that("score_doublets returns expected columns", {
    counts <- mock_counts()

    prev_out <- list(counts_list = list(sample1=counts))
    out <- score_doublets(NULL, NULL, prev_out)$output

    expect_setequal(colnames(out$doublet_scores$sample1),
                    c('barcodes', 'doublet_class', 'doublet_scores'))
})
