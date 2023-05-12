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
  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples)

  for (sample in names(scdata_list)) {
    if (sample %in% unfiltered_samples) {
      expect_true(qc_config$classifier[[sample]]$enabled)
      expect_false(qc_config$classifier[[sample]]$prefiltered)
    } else {
      expect_false(qc_config$classifier[[sample]]$enabled)
      expect_true(qc_config$classifier[[sample]]$prefiltered)
    }

    expect_false(qc_config$cellSizeDistribution[[sample]]$enabled)
  }
})


test_that("cellsize filter is disabled by default and classifier is not pre-filtered", {
  scdata_list <- mock_scdata_list()
  unfiltered_samples <- c()
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples)

  for (sample in names(scdata_list)) {
    expect_false(qc_config$cellSizeDistribution[[sample]]$enabled)
    if (sample %in% unfiltered_samples) {
      expect_true(qc_config$classifier[[sample]]$enabled)
      expect_false(qc_config$classifier[[sample]]$prefiltered)
    } else {
      expect_false(qc_config$classifier[[sample]]$enabled)
      expect_true(qc_config$classifier[[sample]]$prefiltered)
    }
  }
})

test_that("customize_doublet_config sets threshold to 0 when there are no singlets", {
  scdata_list <- mock_scdata_list()
  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples)

  for (sample in names(scdata_list)) {
    scdata_list[[sample]]$doublet_class <- "doublet"
    config <- customize_doublet_config(scdata_list[[sample]], qc_config)
    expect_equal(config$filterSettings$probabilityThreshold, 0)
  }
})


test_that("classifier filter config is enabled for unfiltered samples and disabled for pre-filtered samples", {
  scdata_list <- mock_scdata_list()
  unfiltered_samples <- c("123abc")
  qc_config <- construct_qc_config(scdata_list, unfiltered_samples = unfiltered_samples)

  for (sample in names(scdata_list)) {
    if (sample %in% unfiltered_samples) {
      expect_true(qc_config$classifier[[sample]]$enabled)
      expect_false(qc_config$classifier[[sample]]$prefiltered)
    } else {
      expect_false(qc_config$classifier[[sample]]$enabled)
      expect_true(qc_config$classifier[[sample]]$prefiltered)
    }
  }
})
