library(Seurat)

mock_scdata <- function(rename_genes = c(), n_rep = 1) {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  # replicate matrix columns n times to create a bigger mock dataset
  pbmc_raw <- do.call("cbind", replicate(n_rep, pbmc_raw, simplify = FALSE))

  if (length(rename_genes) > 0) {
    # rename some genes to match cell cycle genes
    some_genes <- sample(1:nrow(pbmc_raw), length(rename_genes))
    rownames(pbmc_raw)[some_genes] <- rename_genes
  }

  sample_1_id <- "123abc"
  sample_2_id <- "123def"

  pbmc_raw <- as(as.matrix(pbmc_raw), 'dgCMatrix')
  scdata <- Seurat::CreateSeuratObject(counts = pbmc_raw)

  # add samples
  scdata$samples <- rep(c(sample_1_id, sample_2_id), each = ncol(scdata) / 2)
  scdata$cells_id <- 0:(ncol(scdata) - 1)
  scdata@misc$gene_annotations <- data.frame(input = rownames(scdata), name = paste("SYMBOL -", rownames(scdata)))
  rownames(scdata@misc$gene_annotations) <- rownames(scdata)

  scdata_sample1 <- subset(scdata, samples == sample_1_id)
  scdata_sample2 <- subset(scdata, samples == sample_2_id)


  scdata_list <- list(scdata_sample1, scdata_sample2)
  names(scdata_list) <- c(sample_1_id, sample_2_id)

  return(list(scdata_list, sample_1_id, sample_2_id))
}

mock_ids <- function() {
  return(list("123abc" = 0:39, "123def" = 40:79))
}

test_that("FastMNN works", {
  c(scdata_list, sample_1_id, sample_2_id) %<-% mock_scdata()
  cells_id <- mock_ids()

  config <- list(
    dimensionalityReduction = list(numPCs = 5),
    dataIntegration = list(method = "fastmnn", methodSettings = list(fastmnn = list(numGenes = 1000, normalisation = "logNormalize")))
  )

  integrated_scdata <- suppressWarnings(run_fastmnn(scdata_list, config = config, cells_id))
  integrated_scdata <- clean_timestamp(integrated_scdata)
  integrated_scdata <- remove_commands_functions(integrated_scdata)

  expect_s4_class(integrated_scdata, "Seurat")
  expect_snapshot(str(integrated_scdata))
})
