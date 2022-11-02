# time stamp used for directory to store log/dump files in event of error
debug_timestamp <- format(Sys.time(), format = "%Y-%m-%d_at_%H-%M-%OS3")


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

# constants used in GEM2S
gem2s <- list(
  max.edrops.fdr = 0.001,
  max.empty.counts = 100,
  max.empty.drops = 50
)

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

usethis::use_data(
  debug_timestamp,
  bucket_list,
  gem2s,
  RANDOM_SEED,
  DEBUG_PATH,
  file_names,
  file_types_by_technology,
  SYM_IDS,
  SYM_SYM,
  IDS_SYM,
  IDS_IDS,
  cc_genes,
  pipeline_version,
  internal = TRUE,
  overwrite = TRUE
)

# add comment to a file that's never run
