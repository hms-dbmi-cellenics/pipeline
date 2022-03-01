#' Read input folder of 10x data
#'
#' @inheritParams download_user_files
#' @param prev_out list with experiment configuration settings
#'
#' @return list with 'output' slot containing \itemize{
#'   \item{"counts_list"}{named list of dgCMatrix per sample}
#'   \item{"annot"}{data.frame from reading features.tsv.gz for each sample}
#' }
#' @export
#'
load_cellranger_files <- function(input, pipeline_config, prev_out, input_dir = "/input") {
  message("Loading cellranger output ...")
  check_prev_out(prev_out, "config")

  # destructure previous output
  config <- prev_out$config

  output <- c(prev_out, call_read10x(config, input_dir))

  res <- list(
    data = list(),
    output = output
  )

  message("\nLoading of cellranger files step complete.")
  return(res)
}

#' Calls Read10X
#'
#' Cellranger outputs from V2 and V3 kits were renamed to look like V3 (features.tsv.gz).
#'
#' @param config experiment settings.
#'
call_read10x <- function(config, input_dir) {
  counts_list <- list()
  annot_list <- list()

  samples <- config$samples
  message("Samples to include in the analysis:\n- ", paste(samples, collapse = "\n- "))
  message("Loading 10x data set from input folder.")


  for (sample in samples) {
    sample_dir <- file.path(input_dir, sample)
    sample_fpaths <- list.files(sample_dir)
    annot_fpath <- file.path(sample_dir, "features.tsv.gz")

    message("\nSample --> ", sample)
    message("Reading files from ", sample_dir, " --> ", paste(sample_fpaths, collapse = " - "))

    counts <- Seurat::Read10X(sample_dir, gene.column = 1)

    if (is(counts, "list")) {
      slot <- "Gene Expression"
      # questionable: grab first slot if no gene expression
      if (!slot %in% names(counts)) slot <- names(counts)[1]
      counts <- counts[[slot]]
    }

    annot <- read.delim(annot_fpath, header = FALSE)

    # Equalizing number of columns in case theres no Gene Expression column
    annot <- annot[, c(1, 2)]

    message(
      sprintf("Sample %s has %s genes and %s droplets.", sample, nrow(counts), ncol(counts))
    )

    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annot
  }

  annot <- format_annot(annot_list)

  return(list(counts_list = counts_list, annot = annot))
}

format_annot <- function(annot_list) {
  annot <- unique(do.call("rbind", annot_list))
  colnames(annot) <- c("input", "name")

  message("Deduplicating gene annotations...")

  # add ENSEMBL ID for genes that are duplicated (geneNameDuplicated-ENSEMBL)
  # original name kept in 'original_name' column
  gname <- annot$name
  annot$original_name <- gname
  is.dup <- duplicated(gname) | duplicated(gname, fromLast = TRUE)

  # We need to convert the gene inputs from _ to - bc when we create the Seurat
  # object we do this, and the match would return NA values if any of the inputs still has _.
  annot$input <- gsub("_", "-", annot$input)
  annot$name[is.dup] <- paste(gname[is.dup], annot$input[is.dup], sep = " - ")

  annot <- annot[!duplicated(annot$input), ]

  rownames(annot) <- annot$input
  return(annot)
}


read_rhapsody_matrix <- function(config, input_dir) {
  counts_list <- list()
  # annot_list <- list()

  samples <- config$samples
  message("Samples to include in the analysis:\n- ", paste(samples, collapse = "\n- "))
  message("Loading rhapsody data set from input folder.")


  for (sample in samples) {
    sample_dir <- file.path(input_dir, sample)
    sample_fpaths <- list.files(sample_dir)

    message("\nSample --> ", sample)
    message("Reading files from ", sample_dir, " --> ", paste(sample_fpaths, collapse = " - "))

    counts <- data.table::fread(sample_fpaths)

    # catch absent DBEC column
    if ("DBEC_Adjusted_Molecules" %in% names(counts)) {
      counts <- counts[, c("Cell_Index", "Gene", "DBEC_Adjusted_Molecules")]
    } else {
      counts <- counts[, c("Cell_Index", "Gene", "RSEC_Adjusted_Molecules")]
    }

    # order by cell indices and gene, to ensure correct cell_index_j to
    # column name association when using frank
    setorder(counts, Cell_Index, Gene)

    counts[, Gene := factor(Gene)]
    counts[, gene_i := as.integer(Gene)]

    # to create small sparse matrix, and retain original cell indices ("barcodes")
    counts[, cell_index_j := frank(counts[, Cell_Index], ties.method = "dense")]

    # create dimnames
    dimnames_i <- levels(counts[, Gene])
    dimnames_j <- unique(counts[, Cell_Index])

    counts <- Matrix::sparseMatrix(
      i = counts[, gene_i],
      j = counts[, cell_index_j],
      x = counts[, DBEC_Adjusted_Molecules],
      dimnames = list(dimnames_i, dimnames_j)
    )

    message(
      sprintf("Sample %s has %s genes and %s wells", sample, nrow(counts), ncol(counts))
    )

    counts_list[[sample]] <- counts
    # annot_list[[sample]] <- annot
  }

  # annot <- format_annot(annot_list)

  return(list(counts_list = counts_list)) # , annot = annot))
}
