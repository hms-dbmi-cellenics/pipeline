get_mock_cell_sets <- function () {
  cell_sets_path <- file.path(setup_test_paths()$mock_data, "cell_sets")
  path <- file.path(cell_sets_path, "cell_sets_2_samples.json")
  cell_sets <- list(cellSets = jsonlite::fromJSON(path, flatten = T))

  return(cell_sets)
}