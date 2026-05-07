library(Seurat)



test_that("FastMNN works with or without bpcells", {

  # use_bpcells = FALSE last for snapshot
  for (use_bpcells in c(TRUE, FALSE)) {
     scdata_list <- mock_scdata_list(use_bpcells = use_bpcells)
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

   }
  
  # snapshot is on result without bpcells
  skip_on_ci()
  expect_snapshot(str(integrated_scdata))
})
