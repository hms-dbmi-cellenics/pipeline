test_that("downsample_plotdata keeps at most max_number_of_cells", {
    ncol_sample <- 1001
    max_number_of_cells <- 1000
    nkeep <- downsample_plotdata(ncol_sample, max_number_of_cells)

    expect_equal(nkeep, max_number_of_cells)
})
