#  - Create a seurat object per sample
#  - Adding emtpyDrops,  scrublet and MT-content
#  - Create a file with flag_filtered

create_seurat <- function(input, pipeline_config, prev_out) {
  # destructure previous output
  counts_list <- prev_out$counts_list
  config <- prev_out$config
  annot <- prev_out$annot

  samples <- names(counts_list)

  message("Creating Seurat Object...")
  scdata_list <- list()
  for (sample in samples) {
    message("\tConverting into seurat object sample --> ", sample)
    construct_scdata(counts_list[[sample]], sample, annot, config)
  }

  message("Step 5 completed.")

  return(list())
}

#' construct metadata for each SeuratObject
#' @param counts matrix with barcodes as columns to assign metadata information
#' @param sample name of the sample
#' @param config configuration for experiment
#'
#' @export
#' @return dataframe with metadata information of the sample

construct_metadata <- function(counts, sample, config) {
  metadata <- data.frame(row.names = colnames(counts), samples = rep(sample, ncol(counts)))

  # Add "metadata" if exists in config
  if ("metadata" %in% names(config)) {
    rest_metadata <- as.data.frame(config$metadata)
    rest_metadata <- rest_metadata[match(metadata$samples, config$samples),, drop = FALSE]
    metadata <- cbind(metadata, rest_metadata)
  }

  return(metadata)
}

construct_scdata <- function(counts, sample, annot, config) {
  message("Converting into seurat object sample --> ", sample)

  metadata <- construct_metadata(counts, sample, config)
  scdata <- Seurat::CreateSeuratObject(counts, meta.data = metadata, project = config$name)

  if (any(grepl("^mt-", annot$name, ignore.case = TRUE))) {
    message("[", sample, "] \t Adding MT information...")
    mt.features <- annot$input[grep("^mt-", annot$name, ignore.case = TRUE)]
    mt.features <- mt.features[mt.features %in% rownames(scdata)]
    if (length(mt.features)) {
      scdata <- PercentageFeatureSet(scdata, features = mt.features, col.name = "percent.mt")
    }
  }

  if (is.null(scdata@meta.data$percent.mt)) scdata$percent.mt <- 0

  message("[", sample, "] \t Getting scrublet results...")
  scores <- get_doublet_score(sample)
  rownames(scores) <- scores$barcodes

  idt <- scores$barcodes[scores$barcodes %in% rownames(scdata@meta.data)]
  scdata@meta.data[idt, "doublet_scores"] <- scores[idt, "doublet_scores"]
  # Doublet class is the classifications that scDblFinder does to set the threshold of doublet_scores
  # (https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/2_scDblFinder.html#thresholding-and-local-calibration)
  scdata@meta.data[idt, "doublet_class"] <- scores[idt, "doublet_class"]

  message("[", sample, "] \t Adding emptyDrops...")
  file_ed <- paste("/output/pre-emptydrops-", sample, ".rds", sep = "")


  if (file.exists(file_ed)) {
    scdata@tools$flag_filtered <- FALSE
    message("\t \t getting emptyDrops results...")
    emptydrops_out <- readRDS(file = file_ed)

    emptydrops_out_df <- emptydrops_out %>%
      as.data.frame() %>%
      rlang::set_names(~ paste0("emptyDrops_", .)) %>%
      tibble::rownames_to_column("barcode")

    # adding emptydrops data to meta.data for later filtering (using left join)
    meta.data <- scdata@meta.data %>%
      tibble::rownames_to_column("barcode") %>%
      dplyr::left_join(emptydrops_out_df)
    rownames(meta.data) <- meta.data$barcode

    message("\t \t Adding emptyDrops scores information...")
    scdata@meta.data <- meta.data
    # previously (before joining into meta.data), results were just dumped as a additional slot
    # leaving the code here in case bugs arise from above solution
    # scdata@tools$CalculateEmptyDrops <- emptydrops_out
  } else {
    # Later on, when creating the config file, the enable will look the value of flag_filtered to   deactivate the classifier filter
    message("\t \t emptyDrops results not present, skipping...")
    scdata@meta.data$emptyDrops_FDR <- NA
    scdata@tools$flag_filtered <- TRUE
  }

  if (!dir.exists("/output/rds_samples")) {
    dir.create("/output/rds_samples")
  }

  message("[", sample, "] Saving R object...")
  saveRDS(scdata, file = paste("/output/rds_samples/", sample, ".rds", sep = ""), compress = FALSE)

  return(scdata@tools$flag_filtered)
}
