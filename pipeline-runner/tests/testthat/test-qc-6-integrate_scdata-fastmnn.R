library(Seurat)



test_that("FastMNN works", {
  scdata_list <- mock_scdata_list()
  cells_id <- mock_ids(scdata_list)

  config <- list(
    dimensionalityReduction = list(numPCs = 5),
    dataIntegration = list(
      method = "fastmnn",
      methodSettings = list(
        fastmnn = list(numGenes = 1000, normalisation = "logNormalize")
      )
    )
  )

  integrated_scdata <- suppressWarnings(
    run_fastmnn(scdata_list, config = config, cells_id)
  )
  integrated_scdata <- clean_timestamp(integrated_scdata)
  integrated_scdata <- remove_commands_functions(integrated_scdata)

  expect_s4_class(integrated_scdata, "Seurat")
  skip_if(is_bpcells())
  expect_snapshot(str(integrated_scdata))
})
