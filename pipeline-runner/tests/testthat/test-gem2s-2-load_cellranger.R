


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

    expect_true(length(annot$name) == length(unique(annot$name)))
})
