mock_scdata_list <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )

  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)
  # add samples
  scdata$samples <- rep("123abc", 80)
  scdata <- Seurat::RenameCells(scdata, paste(scdata$samples, colnames(scdata), sep = ""))

  # add doublet scores
  scdata$doublet_scores <- rep(c(0.01, 0.9), each = 40)
  scdata$doublet_class <- rep(c("singlet", "doublet"), each = 40)

  # add mitochondrial percent
  scdata$percent.mt <- rnorm(ncol(scdata), mean = 6)

  # create an scdata_list with duplicated samples
  scdata_list <- list()
  for (sample_id in scdata$samples) {
      scdata_list[[sample_id]] <- scdata
  }
  return(scdata_list)
}

test_that("cellsize filter is disabled by default and classifier is pre-filtered", {
  scdata_list <- mock_scdata_list()
  qc_config <- construct_qc_config(scdata_list, any_filtered = TRUE, subset_experiment = FALSE)

  for (sample in names(scdata_list)) {
    expect_false(qc_config$classifier[[sample]]$enabled)
    expect_true(qc_config$classifier[[sample]]$prefiltered)
    expect_false(qc_config$cellSizeDistribution[[sample]]$enabled)
  }
})


test_that("cellsize filter is disabled by default and classifier is not pre-filtered", {
  scdata_list <- mock_scdata_list()
  qc_config <- construct_qc_config(scdata_list, any_filtered = FALSE, subset_experiment = FALSE)

  for (sample in names(scdata_list)) {
    expect_false(qc_config$cellSizeDistribution[[sample]]$enabled)
    expect_true(qc_config$classifier[[sample]]$enabled)
    expect_false(qc_config$classifier[[sample]]$prefiltered)
  }
})

