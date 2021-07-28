#' Create a SeuratObject per sample
#'
#' @inheritParams download_cellranger
#' @param prev_out  'output' slot from call to \code{score_doublets}
create_seurat <- function(input, pipeline_config, prev_out) {
  message("Creating Seurat Objects...")
  message("items in prev_out: ", paste(names(prev_out), collapse = ' - '))

  # destructure previous output
  list2env(prev_out, envir = environment())

  samples <- names(counts_list)
  scdata_list <- list()
  for (sample in samples) {
    message("Creating SeuratObject for sample --> ", sample)

    scdata_list[[sample]] <- construct_scdata(
      counts = counts_list[[sample]],
      doublet_score = doublet_scores[[sample]],
      edout = edrops[[sample]],
      sample = sample,
      annot = annot,
      config = config)
  }


  prev_out$scdata_list <- scdata_list
  res <- list(
    data = list(),
    output = prev_out)

  message("Step 5 completed.")
  return(res)
}

# construct SeuratObject
construct_scdata <- function(counts, doublet_score, edout, sample, annot, config) {

  metadata <- construct_metadata(counts, sample, config)
  scdata <- Seurat::CreateSeuratObject(counts, meta.data = metadata, project = config$name)
  scdata <- add_mito(scdata, annot)
  scdata <- add_dblscore(scdata, doublet_score)
  scdata <- add_edrops(scdata, edout)

  return(scdata)
}

# construct metadata for each SeuratObject
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

# add mitochondrial percent to SeuratObject
add_mito <- function(scdata, annot) {

  if (any(grepl("^mt-", annot$name, ignore.case = TRUE))) {
    message("\tAdding MT information...")
    mt.features <- annot$input[grep("^mt-", annot$name, ignore.case = TRUE)]
    mt.features <- mt.features[mt.features %in% rownames(scdata)]
    if (length(mt.features)) {
      scdata <- PercentageFeatureSet(scdata, features = mt.features, col.name = "percent.mt")
    }
  }

  if (is.null(scdata@meta.data$percent.mt)) scdata$percent.mt <- 0
  return(scdata)
}

# add emptyDrops result to SeuratObject
add_edrops <- function(scdata, edout) {

  scdata@tools$flag_filtered <- is.null(edout)

  if (!scdata@tools$flag_filtered) {
    message("\tAdding emptyDrops scores...")

    edout <- emptydrops_out %>%
      as.data.frame() %>%
      rlang::set_names(~ paste0("emptyDrops_", .)) %>%
      tibble::rownames_to_column("barcode")

    # adding emptydrops data to meta.data for later filtering (using left join)
    meta.data <- scdata@meta.data %>%
      tibble::rownames_to_column("barcode") %>%
      dplyr::left_join(edout)
    rownames(meta.data) <- meta.data$barcode

    scdata@meta.data <- meta.data

  } else {
    message("\temptyDrops results not present, skipping...")
    scdata@meta.data$emptyDrops_FDR <- NA
  }

  return(scdata)
}

# add scDblFinder result to SeuratObject
add_dblscore <- function(scdata, score) {
  message("\tAdding doublet scores...")

  idt <- scores$barcodes[scores$barcodes %in% rownames(scdata@meta.data)]
  scdata@meta.data[idt, "doublet_scores"] <- scores[idt, "doublet_scores"]
  scdata@meta.data[idt, "doublet_class"] <- scores[idt, "doublet_class"]
  return(scdata)
}
