#  - Read input folder {10x data}
#  - Prepare rds with a list of raw counts matrix per sample to compute emptyDrops

load_cellranger <- function(input, pipeline_config) {
  message("Loading configuration...")
  config <- RJSONIO::fromJSON("/input/meta.json")

  print_config(2,"Load Cellranger",input,pipeline_config,config)

  # We include in variable scdata_list all the sparse matrix per sample
  message("Creating raw dataframe...")
  scdata_list <- call_read10x(config)

  # We store the raw scdata_list for the emptyDrops, since to compute the background we cannot remove any cells.
  message("Exporting raw scdata for emptyDrops...")
  saveRDS(scdata_list, file = "/output/pre-doublet-scdata_list.rds", compress = FALSE)

  message("Step 2 completed.")

  return(list())
}

# checking_10x_structure
#' @description The design of the input data needs to be in a particular way. With this function we are going to check if
#' the input folder is designed correctly:
#' intput/
#' ------ sample_name/
#' ------------------- features.tsv.gz
#' ------------------- barcodes.tsv.gz
#' ------------------- matrix.mtx.gz
#'
#' cell ranger output
#' https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
#' STAR solo conventions (drop in replacement for cell ranger):
#' https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#' @param samples samples to check
#'
#' @return TRUE if the design is correct FALSE otherwise
check_10x_input <- function(samples) {
  #The UI will upload files with these names, regardless of actual Cellranger version. 
  fnames <- c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")
  fpaths <- file.path("/input", samples, fnames)
  all(file.exists(fpaths))
}


#' Calls Read10X
#'
#' Cellranger outputs from V2 and V3 kits were renamed to look like V3 (features.tsv.gz).
#'
#' @param config experiment settings.
#' @return list with an element per sample with dgCMatrix with genes as rows and cells as column
#'
call_read10x <- function(config) {
  data_type <- config$input["type"]
  scdata <- list()

  # Check config format. It requires to have a 'samples' key. For the multisample experiment, it needs to be the same
  # as the name of the folders that are inside the input folder. In the case of unisample it should be [].
  if (!"samples" %in% names(config)) {
    stop("The format of the config is wrong. There should be a category for the sample.
     Current content of config:", names(config))
  }

  # If no samples have been added, we will take all the samples of the project.
  if (length(config$samples) == 0) {
    warning("All samples in input folder will be included in the analysis!")
    samples <- list.dirs("/input", full.names = F)[-1]
  } else {
    samples <- config$samples

    if (!all(samples %in% list.dirs("/input", full.names = FALSE))) {
      stop(
        "Check samples to be used in the analysis,
      since there are some of them that hasn't got a folder with the files: ",
        samples[!samples %in% list.dirs("/input", full.names = FALSE)]
      )
    }
  }

  message("Samples to include in the analysis: ", paste(samples, collapse = " - "))

  if (data_type == "10x") {
    message("Loading 10x data set from input folder.")

    # Checking design
    if (!check_10x_input(samples)) {
      stop(
        "Please! Check in files inside the input folder.",
        "There should be the files {features.tsv.gz, barcodes.tsv.gz and matrix.mtx.gz}"
      )
    }

    annotation_features <- list()
    # overall feature annotation is derived from input data saved in genes.tsv features.tsv.gz
    # since each sample only carries a subset of annotation for its expressed genes, the annotation for all samples is merged.
    # More information about genes.tsv features.tsv.gz: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices

    for (sample in samples) {
      sample_dir <- file.path("/input", sample)
      sample_fpaths <- list.files(sample_dir, full.names = TRUE)
      message("Reading files --> ", sample_fpaths)

      scdata[[sample]] <- Seurat::Read10X(sample_dir, gene.column = 1)
      annot_fpath <- sample_fpaths[grepl("genes.tsv$|features.tsv.gz$", sample_fpaths)]
      annotation_features[[sample]] <- read.delim(annot_fpath, header = FALSE)

      message(
        "Found ", nrow(scdata[[sample]]), " genes and ",
        ncol(scdata[[sample]]), " cells in sample ", sample, "."
      )
    }
    annotation_features_df <- unique(do.call("rbind", annotation_features))
    annotation_features_df <- annotation_features_df[, c(1, 2)]
    colnames(annotation_features_df) <- c("input", "name")
    write.table(annotation_features_df, "/output/features_annotations.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }

  if (data_type == "table") {
    stop("data_type of table not yet implemented")

    message("Loading table-type data set from ", path)
    for (sample in samples) {
      scdata[[sample]] <- as.matrix(read.table(paste("/input", sample, sep = "/")))
      message(
        "Found ", nrow(scdata$raw[[sample]]), " genes and ",
        ncol(scdata$raw[[sample]]), " cells in sample ", sample, "."
      )
    }
  }

  return(scdata)
}
