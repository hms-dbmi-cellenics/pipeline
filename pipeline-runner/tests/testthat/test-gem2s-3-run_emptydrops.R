mock_counts <- function() {
    set.seed(0)
    DropletUtils:::simCounts()
}


test_that("run_emptydrops is skipped when dataset is pre-filtered", {
    counts <- mock_counts()
    ntot <- Matrix::colSums(counts)
    counts <- counts[, ntot > 100]

    prev_out <- list(counts_list = list(sample1=counts))
    expect_message(
        out <- run_emptydrops(NULL, NULL, prev_out),
        'filtered'
    )

    expect_equal(out$output$edrops, list())
})

test_that("run_emptydrops runs when not pre-filtered", {
    counts <- mock_counts()
    prev_out <- list(counts_list = list(sample1=counts, sample2=counts))
    out <- run_emptydrops(NULL, NULL, prev_out)$output

    # names of stored things are right
    expect_named(out, c('counts_list', 'edrops'))
    expect_named(out$edrops, c('sample1', 'sample2'))

    # emptyDrops stores something
    expect_s4_class(out$edrops$sample1, 'DFrame')

    # preserves counts list
    expect_equal(prev_out$counts_list, out$counts_list)
})
