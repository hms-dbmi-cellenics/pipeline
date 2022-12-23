# mock_scdata <- function(){
#   paths <- path_setup()
#   source("tests/testthat/_snaps/qc/mock_experiment_id-integrated_scdata.R")
#   scdata <- snap_list$data
#   rm(snap_list, envir = parent.frame())
#   return(scdata)
# }
#
# mock_cellsets <- function(){
#
#   jsonlite::fromJSON("tests/testthat/_snaps/gem2s/gem2s-7-mock_experiment_id-cellsets.json", flatten = TRUE)
#
# }
#
# mock_cluster_cellsets <- function(cellsets) {
#
# }
#
# mock_input <- function() {
#   input <- list(
#     name = "mock_subset_experiment_name",
#     parentExperimentId = "mock_parent_experiment_id",
#     subsetExperimentId = "mock_subset_experiment_id",
#     cellSetKeys =  c("louvain-0", "louvain-1")
#   )
#
#   return(input)
# }
#
# parent_scdata <- mock_scdata()
# cellsets <- mock_cellsets()
# parent_cellsets <- parse_cellsets(mock_cellsets())
#
# parent <- list(scdata = parent_scdata, cellsets = parent_cellsets)
#
# input <- mock_input()
