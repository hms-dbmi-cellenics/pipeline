source("data-raw/processing_config_template.R")

bucket_list <- list(
  originals_bucket = "biomage-originals",
  source_bucket = "biomage-source",
  processed_bucket = "processed-matrix",
  results_bucket = "worker-results",
  cells_id_bucket = "biomage-filtered-cells",
  plot_data_bucket = "plots-tables",
  cell_sets_bucket = "cell-sets",
  debug_bucket = "biomage-pipeline-debug"
)

# list of task functions named by task name
GEM2S_TASK_LIST <- list(
  "downloadGem" = "download_user_files",
  "preproc" = "load_user_files",
  "emptyDrops" = "run_emptydrops",
  "doubletScores" = "score_doublets",
  "createSeurat" = "create_seurat",
  "prepareExperiment" = "prepare_experiment",
  "uploadToAWS" = "upload_to_aws"
)

SUBSET_SEURAT_TASK_LIST <- list(
  "subsetSeurat" = "subset_seurat",
  "prepareExperiment" = "prepare_experiment",
  "uploadToAWS" = "upload_to_aws"
)

# list of task functions named by task name
QC_TASK_LIST <- list(
  "classifier" = "filter_emptydrops",
  "cellSizeDistribution" = "filter_low_cellsize",
  "mitochondrialContent" = "filter_high_mito",
  "numGenesVsNumUmis" = "filter_gene_umi_outlier",
  "doubletScores" = "filter_doublets",
  "dataIntegration" = "integrate_scdata",
  "configureEmbedding" = "embed_and_cluster"
)

# directory where download_user_files downloads user files
INPUT_DIR <- "/input"

# constants used in GEM2S
gem2s <- list(
  max.edrops.fdr = 0.001,
  max.empty.counts = 100,
  max.empty.drops = 50
)


# minimum number of cells required in a sample to have the pipeline not break.
MIN_CELLS_IN_SAMPLE <- 15


RANDOM_SEED <- 42

# path where dump/log files are saved
# mounted as a volume outside container to local-runner/debug
DEBUG_PATH <- "/debug"

# File management, it needs to match the sample_file_type enum in sql
# (they are originally defined in 20220304184711_schema.js in the api)
file_types_by_technology <- list(
  "10x" = list("barcodes10x", "features10x", "matrix10x"),
  "rhapsody" = list("rhapsody")
)

file_names <- list(
  barcodes10x = "barcodes.tsv.gz",
  features10x = "features.tsv.gz",
  matrix10x = "matrix.mtx.gz",
  rhapsody = "expression_data.st.gz"
)

source("data-raw/cell_cycle_genes.R")

cc_genes <- list(
  "human" = human_cc_genes,
  "mouse" = mouse_cc_genes
)

# annotation type constants
SYM_IDS <- "sym_ids"
SYM_SYM <- "sym_sym"
IDS_SYM <- "ids_sym"
IDS_IDS <- "ids_ids"

pipeline_version <- 2

UNISAMPLE <- "unisample"

usethis::use_data(
  debug_timestamp,
  processing_config_template,
  bucket_list,
  gem2s,
  SUBSET_SEURAT_TASK_LIST,
  GEM2S_TASK_LIST,
  QC_TASK_LIST,
  INPUT_DIR,
  RANDOM_SEED,
  MIN_CELLS_IN_SAMPLE,
  DEBUG_PATH,
  file_names,
  file_types_by_technology,
  SYM_IDS,
  SYM_SYM,
  IDS_SYM,
  IDS_IDS,
  cc_genes,
  pipeline_version,
  UNISAMPLE,
  internal = TRUE,
  overwrite = TRUE
)
