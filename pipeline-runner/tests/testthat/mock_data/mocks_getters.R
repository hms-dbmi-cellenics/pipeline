source("../helper-stub_functions")

get_mock_cell_sets <- function(flatten = TRUE) {
  cell_sets_path <- file.path(setup_test_paths()$mock_data, "cell_sets")
  path <- file.path(cell_sets_path, "cell_sets_2_samples.json")
  cell_sets <- jsonlite::fromJSON(path, flatten = flatten)

  return(cell_sets)
}