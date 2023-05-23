#' Create a SeuratObject per sample
#'
#' @inheritParams download_user_files
#' @param prev_out  'output' slot from call to \code{score_doublets}
#'
#' @return \code{prev_out} with added slot 'scdata_list' containing \code{SeuratObject}'s for each sample.
#'
#' @export
#'
create_seurat <- function(input, pipeline_config, prev_out) {
  message("Creating Seurat Objects...")

  # NOTE: edrops can be empty list
  check_names <- c("config", "counts_list", "annot", "doublet_scores", "edrops")
  check_prev_out(prev_out, check_names)

  # destructure previous output: config, counts_list, annot, and doublet_scores
  config <- prev_out$config
  counts_list <- prev_out$counts_list
  annot <- prev_out$annot
  doublet_scores <- prev_out$doublet_scores
  edrops <- prev_out$edrops

  scdata_list <- list()
  samples <- names(counts_list)
  for (sample in samples) {
    message("\nCreating SeuratObject for sample --> ", sample)

    scdata_list[[sample]] <- construct_scdata(
      counts = counts_list[[sample]],
      doublet_score = doublet_scores[[sample]],
      edrops_out = edrops[[sample]],
      sample = sample,
      annot = annot,
      config = config
    )
  }


  prev_out$scdata_list <- scdata_list
  prev_out$disable_qc_filters <- FALSE
  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nCreation of Seurat objects step complete.")
  return(res)
}

# construct SeuratObject
construct_scdata <- function(counts, doublet_score, edrops_out, sample, annot, config, min.cells = 0, min.features = 10) {
  metadata <- construct_metadata(counts, sample, config)

  scdata <- Seurat::CreateSeuratObject(
    counts,
    meta.data = metadata,
    project = config$name,
    min.cells = min.cells,
    min.features = min.features
  )

  scdata <- scdata %>%
    add_mito(annot) %>%
    add_dblscore(doublet_score) %>%
    add_edrops(edrops_out)

  return(scdata)
}


#' Construct metadata for each SeuratObject
#'
#' This function creates a `data.frame` with the barcodes, sampleIDs and user supplied
#' metadata that corresponds to each one.
#'
#' @param counts count matrix
#' @param sample character sample ID
#' @param config list containing experiment config
#'
#' @return data.frame of sample metadata
#'
construct_metadata <- function(counts, sample, config) {
  message("Constructing metadata data.frame...")
  metadata_df <- data.frame(row.names = colnames(counts), samples = rep(sample, ncol(counts)))

  # Add "metadata" if exists in config
  user_metadata <- config$metadata
  if (!is.null(user_metadata)) {
    user_metadata <- lapply(user_metadata, unlist)
    user_metadata <- data.frame(user_metadata, row.names = config$samples, check.names = FALSE)
    metadata_df[names(user_metadata)] <- user_metadata[sample, ]
  }

  # make syntactically valid column names
  colnames(metadata_df) <- make.names(colnames(metadata_df), unique = TRUE)

  return(metadata_df)
}

# add mitochondrial percent to SeuratObject
add_mito <- function(scdata, annot) {
  if (any(grepl(MITOCHONDRIAL_REGEX, annot$name, ignore.case = TRUE))) {
    message("Adding MT information...")
    mt.features <- annot$input[grep(MITOCHONDRIAL_REGEX, annot$name, ignore.case = TRUE)]
    mt.features <- mt.features[mt.features %in% rownames(scdata)]
    if (length(mt.features)) {
      scdata <- Seurat::PercentageFeatureSet(scdata, features = mt.features, col.name = "percent.mt")
    }
  }

  if (is.null(scdata@meta.data$percent.mt)) scdata$percent.mt <- 0
  return(scdata)
}

# add emptyDrops result to SeuratObject
add_edrops <- function(scdata, edout) {
  scdata@tools$flag_filtered <- is.null(edout)

  if (!scdata@tools$flag_filtered) {
    message("Adding emptyDrops scores...")

    edout <- edout %>%
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
    message("emptyDrops results not present, skipping...")
    scdata@meta.data$emptyDrops_FDR <- NA
  }

  return(scdata)
}

# add scDblFinder result to SeuratObject
add_dblscore <- function(scdata, score) {
  message("Adding doublet scores...")

  idt <- score$barcodes[score$barcodes %in% rownames(scdata@meta.data)]
  scdata@meta.data[idt, "doublet_scores"] <- score[idt, "doublet_scores"]
  scdata@meta.data[idt, "doublet_class"] <- score[idt, "doublet_class"]
  return(scdata)
}
