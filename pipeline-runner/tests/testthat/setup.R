library(httptest)
library(mockery)
library(withr)
library(zeallot)
library(Seurat)
options(Seurat.object.assay.version = "v5")

source("mock_data/send_output_to_api_mock_data.R")
withr::defer(options(testthat.summary.max_reports = 30))

# prevent deletion of snapshot files
announce_snapshot_file(name = "_snaps/gem2s-7-upload_to_aws/cellsets.json")
announce_snapshot_file(name = "_snaps/gem2s/gem2s-6-mock_experiment_id-out.R")
announce_snapshot_file(name = "_snaps/gem2s/gem2s-6-mock_experiment_id-qc_config.R")
announce_snapshot_file(name = "_snaps/gem2s/gem2s-7-mock_experiment_id-cellsets.json")
announce_snapshot_file(name = "_snaps/qc/cluster_cell_sets.json")
announce_snapshot_file(name = "_snaps/qc/mock_experiment_id-integrated_scdata.R")