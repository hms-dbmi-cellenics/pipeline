mock_scdata <- function(){
  processed_path <- "/Users/german/bm/cellenics/data/8ecc9d20-30e4-49eb-b536-a0d1f0ba420d/processed_r.rds"
  readRDS(processed_path)
}

mock_cellsets <- function(){
  cellsets_path <- "/Users/german/bm/cellenics/data/8ecc9d20-30e4-49eb-b536-a0d1f0ba420d/cellsets.json"
  jsonlite::fromJSON(cellsets_path, flatten = TRUE)
}

mock_input <- function() {
  input <- list(
    name = "mock_subset_experiment_name",
    parentExperimentId = "mock_parent_experiment_id",
    subsetExperimentId = "mock_subset_experiment_id",
    cellSetKeys =  c("louvain-0", "louvain-1")
  )

  return(input)
}

parent_scdata <- mock_scdata()
parent_cellsets <- parse_cellsets(mock_cellsets())
sample_mapping <- mock_sample_id_mapping()
input <- mock_input()
