source("data-raw/processing_config_template.R")

bucket_list <- list(
  originals_bucket = "biomage-originals",
  source_bucket = "biomage-source",
  processed_bucket = "processed-matrix",
  results_bucket = "worker-results",
  cells_id_bucket = "biomage-filtered-cells",
  plot_data_bucket = "plots-tables",
  cell_sets_bucket = "cell-sets",
  debug_bucket = "biomage-pipeline-debug",
  cl_metadata_bucket = "cellenics-cell-level-metadata"
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

COPY_TASK_LIST <- list(
  "copyS3Objects" = "copy_s3_objects"
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

SEURAT_TASK_LIST <- list(
  "downloadSeurat" = "download_user_files",
  "processSeurat" = "load_obj2s_file",
  "uploadSeuratToAWS" = "upload_obj2s_to_aws"
)

# directory where download_user_files downloads user files
INPUT_DIR <- "/input"

# constants used in GEM2S
gem2s <- list(
  max.edrops.fdr = 0.001,
  max.empty.counts = 100,
  max.empty.drops = 100
)


# minimum number of cells required in a sample to have the pipeline not break.
MIN_CELLS_IN_SAMPLE <- 15


RANDOM_SEED <- 42

# path where dump/log files are saved
# mounted as a volume outside container to local-runner/debug
DEBUG_PATH <- "/debug"

# File management, it needs to match the sample_file_type enum in sql
# (they are originally defined in the structure of the SQL database)
file_types_by_technology <- list(
  "10x" = list("barcodes10x", "features10x", "matrix10x"),
  "seurat" = list("seurat"),
  "single_cell_experiment" = list("singleCellExperiment"),
  "rhapsody" = list("rhapsody"),
  "10x_h5" = list("10XH5"),
  "parse" = list("barcodesParse", "featuresParse", "matrixParse")
)

file_names <- list(
  barcodes10x = "barcodes.tsv.gz",
  features10x = "features.tsv.gz",
  matrix10x = "matrix.mtx.gz",
  seurat = "r.rds",
  singleCellExperiment = "r.rds",
  rhapsody = "expression_data.st.gz",
  "10XH5" = "matrix.h5.gz",
  barcodesParse = "cell_metadata.csv.gz",
  featuresParse = "all_genes.csv.gz",
  matrixParse = "DGE.mtx.gz"
)

MITOCHONDRIAL_REGEX <- "^mt[-:]"
RIBOSOMAL_REGEX <- "^M?RP[LS]|FAU|UBA52|DAP3"

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

# Parse kits empirical doublet rates
DOUBLET_RATE_PARSE <- list(mini = 0.046, WT = 0.034, mega = 0.064)

# pipeline error constants
errors <- list(
  ERROR_SEURAT_RDS = 'ERROR_SEURAT_RDS',
  ERROR_SEURAT_COUNTS = 'ERROR_SEURAT_COUNTS',
  ERROR_SEURAT_HVFINFO = 'ERROR_SEURAT_HVFINFO',
  ERROR_SEURAT_METADATA = 'ERROR_SEURAT_METADATA',
  ERROR_SEURAT_CLUSTERS = 'ERROR_SEURAT_CLUSTERS',
  ERROR_SEURAT_REDUCTION = 'ERROR_SEURAT_REDUCTION',
  ERROR_SEURAT_LOGCOUNTS = 'ERROR_SEURAT_LOGCOUNTS'
)

usethis::use_data(
  processing_config_template,
  bucket_list,
  gem2s,
  SUBSET_SEURAT_TASK_LIST,
  COPY_TASK_LIST,
  GEM2S_TASK_LIST,
  QC_TASK_LIST,
  SEURAT_TASK_LIST,
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
  MITOCHONDRIAL_REGEX,
  RIBOSOMAL_REGEX,
  pipeline_version,
  UNISAMPLE,
  DOUBLET_RATE_PARSE,
  errors,
  internal = TRUE,
  overwrite = TRUE
)
