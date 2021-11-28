#' Create a SeuratObject per sample
#'
#' @inheritParams download_cellranger_files
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
  list2env(prev_out, envir = environment())

  samples <- names(counts_list)
  scdata_list <- list()
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
  res <- list(
    data = list(),
    output = prev_out
  )

  message("\nCreation of Seurat objects step complete.")
  return(res)
}

# construct SeuratObject
construct_scdata <- function(counts, doublet_score, edrops_out, sample, annot, config, min.cells = 3, min.features = 10) {
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



# NOTE: any changes here must be reflected in meta_sets

# construct metadata for each SeuratObject
construct_metadata <- function(counts, sample, config) {
  message("Constructing metadata df...")
  metadata <- data.frame(row.names = colnames(counts), samples = rep(sample, ncol(counts)))

  # Add "metadata" if exists in config
  rest <- config$metadata
  if (!is.null(rest)) {
    rest <- lapply(rest, unlist)
    rest <- data.frame(rest, row.names = config$samples, check.names = FALSE)
    metadata[names(rest)] <- rest[sample, ]
  }

  # make syntactically valid column names
  colnames(metadata) <- make.names(colnames(metadata), unique = TRUE)

  return(metadata)
}

# add mitochondrial percent to SeuratObject
add_mito <- function(scdata, annot) {
  mt_regex <- "^mt[-:]"
  if (any(grepl(mt_regex, annot$name, ignore.case = TRUE))) {
    message("Adding MT information...")
    mt.features <- annot$input[grep(mt_regex, annot$name, ignore.case = TRUE)]
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