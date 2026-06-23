#' Read user input files
#'
#' Checks technology used and dispatches call to correct reading function.
#' In case of 10x, `annot` contains a data.frame from reading the
#' features.tsv.gz file, while in rhapsody it contains the gene symbols as read
#' from the count matrix.
#'
#' @inheritParams download_user_files
#' @param prev_out list with experiment configuration settings
#'
#' @return list with 'output' slot containing \itemize{
#'   \item{"counts_list"}{named list of dgCMatrix per sample}
#'   \item{"annot"}{data.frame with gene ids and/or symbols}
#'   }
#' @export
load_user_files <- function(
  input, pipeline_config, prev_out, input_dir = INPUT_DIR
) {
  message("Loading user files...")
  check_prev_out(prev_out, "config")

  # destructure previous output
  config <- prev_out$config

  technology <- ifelse(
    config$input$type %in% c("rhapsody", "10x_h5", "parse", "visium_hd", "xenium"),
    config$input$type,
    "10x"
  )

  read_fun <- switch(technology,
    "10x" = read_10x_files,
    "rhapsody" = read_rhapsody_files,
    "10x_h5" = read_10x_h5_file,
    "parse" = read_parse_files,
    "visium_hd" = read_visium_hd_files,
    "xenium" = read_xenium_files
  )

  message(
    "Samples to include in the analysis:\n- ",
    paste(config$samples, collapse = "\n- ")
  )
  message("Loading ", technology, " data set from input folder.")

  user_matrices <- c(prev_out, read_fun(config, input_dir))

  res <- list(
    data = list(),
    output = user_matrices
  )

  message("\nLoading of ", technology, " files step complete.")
  return(res)
}

#' Combine and name results from parallel sample processing
#'
#' @param results list of results from bplapply, each with slots: counts,
#'   annotations, matrix_dir
#' @param samples character vector of sample names
#' @param extract_nested logical. If TRUE, extracts nested annotations
#'   (assumes annotations$annot and annotations$feature_types structure)
#'
#' @return list with named counts_list, annot_list, matrix_dir_list, and
#'   optionally feature_types_list
#'
combine_parallel_results <- function(results, samples, extract_nested = FALSE) {
  counts_list <- lapply(results, function(x) x$counts)
  matrix_dir_list <- lapply(results, function(x) x$matrix_dir)

  if (extract_nested) {
    annot_list <- lapply(results, function(x) x$annotations$annot)
    feature_types_list <- lapply(
      results,
      function(x) x$annotations$feature_types
    )
  } else {
    annot_list <- lapply(results, function(x) x$annotations)
    feature_types_list <- NULL
  }

  names(counts_list) <- samples
  names(annot_list) <- samples
  names(matrix_dir_list) <- samples

  if (!is.null(feature_types_list)) {
    names(feature_types_list) <- samples
  }

  result <- list(
    counts_list = counts_list,
    annot_list = annot_list,
    matrix_dir_list = matrix_dir_list
  )

  if (!is.null(feature_types_list)) {
    result$feature_types_list <- feature_types_list
  }

  return(result)
}

#' Read Parse data for a single sample
#'
#' @param sample character sample ID
#' @param input_dir character path to input directory
#'
#' @return list with slots:
#'   \item{counts}{count matrix}
#'   \item{annotations}{data.frame with annotations}
#'   \item{matrix_dir}{directory where matrix is stored}
#'
read_parse_sample <- function(sample, input_dir) {
  sample_dir <- file.path(input_dir, sample)
  sample_fpaths <- list.files(sample_dir)
  annot_fpath <- file.path(sample_dir, "all_genes.csv.gz")

  message("\nSample --> ", sample)
  message(
    "Reading files from ",
    sample_dir,
    " --> ",
    paste(sample_fpaths, collapse = " - ")
  )

  annotations <- read_parse_annotations(annot_fpath, sample)

  mtx_path <- file.path(sample_dir, "DGE.mtx.gz")
  barcodes_path <- file.path(sample_dir, "cell_metadata.csv.gz")
  features_path <- file.path(sample_dir, "all_genes.csv.gz")

  # read barcodes (cells) from CSV
  cell_barcodes <- as.data.frame(
    data.table::fread(barcodes_path, header = FALSE, skip = 1)
  )
  cell_names <- cell_barcodes[[1]]

  # read features (genes) from CSV, using column 1 (parse format)
  feature_names <- as.data.frame(
    data.table::fread(features_path, header = FALSE, skip = 1)
  )
  gene_names <- feature_names[[1]]
  gene_names <- make.unique(gene_names)

  # import matrix from parse format (cell x gene orientation)
  # then transpose to get standard gene x cell format
  counts_mat <- BPCells::import_matrix_market(
    mtx_path,
    row_names = cell_names,
    col_names = gene_names
  )
  counts_mat <- Matrix::t(counts_mat)

  # write the matrix to a directory
  matrix_dir <- file.path(
    tempdir(),
    paste0(sample, "_matrix_dir")
  )
  unlink(matrix_dir, recursive = TRUE)

  # hangs indefinitely in BPCells not attached (see init.R)
  counts <- BPCells::write_matrix_dir(counts_mat, dir = matrix_dir)

  message(
    sprintf(
      "Sample %s has %s genes and %s barcodes",
      sample, nrow(counts), ncol(counts)
    )
  )

  c(counts, annotations) %<-%
    filter_unnamed_features(counts, annotations, sample)

  list(
    counts = counts,
    annotations = annotations,
    matrix_dir = matrix_dir
  )
}

#' Read Parse files. Calls readMTX
#'
#' @param config experiment settings
#' @param input_dir
#'
#' @return
#' @export
#'
read_parse_files <- function(config, input_dir) {
  samples <- config$samples

  # set up parallel processing parameters
  nworkers <- min(length(samples), BATCH_POD_CPUS)
  bpparam <- BiocParallel::MulticoreParam(workers = nworkers)

  # use bplapply for parallel processing
  results <- BiocParallel::bplapply(
    samples,
    read_parse_sample,
    input_dir = input_dir,
    BPPARAM = bpparam
  )

  # combine and name results from parallel processing
  combined <- combine_parallel_results(
    results,
    samples,
    extract_nested = FALSE
  )

  annot <- format_annot(combined$annot_list)

  list(
    counts_list = combined$counts_list,
    annot = annot,
    matrix_dir_list = combined$matrix_dir_list
  )
}

#' Read 10X h5 data for a single sample
#'
#' @param sample character sample ID
#' @param input_dir character path to input directory
#'
#' @return list with slots:
#'   \item{counts}{count matrix}
#'   \item{annotations}{data.frame with annotations}
#'   \item{matrix_dir}{directory where matrix is stored}
#'
read_10x_h5_sample <- function(sample, input_dir, feature_type = NULL) {
  sample_dir <- file.path(input_dir, sample)
  sample_fpaths <- list.files(sample_dir, pattern = "[.]h5(.gz)?$")
  sample_counts_path <- file.path(
    sample_dir,
    sample_fpaths[[1]]
  )

  message("\nSample --> ", sample)
  message(
    "Reading files from ",
    sample_dir,
    " --> ",
    paste(sample_fpaths, collapse = " - ")
  )

  if (length(sample_fpaths) > 1) {
    stop("Only one h5 expected. More files detected.")
  }

  if (R.utils::isGzipped(sample_counts_path)) {
    message("Detected gzipped h5 file. Unzipping...")
    sample_counts_path <- R.utils::gunzip(sample_counts_path)
  }

  # feature_type (e.g. "Gene Expression") subsets multi-modal h5s to a single
  # feature type, dropping the others; NULL keeps every feature (default).
  gene_names <- read_10x_h5_feature_names(
    sample_counts_path,
    feature_type = feature_type
  )
  counts_mat <- BPCells::open_matrix_10x_hdf5(
    sample_counts_path,
    feature_type = feature_type
  )

  # write the matrix to a directory
  matrix_dir <- file.path(
    tempdir(),
    paste0(sample, "_matrix_dir")
  )
  unlink(matrix_dir, recursive = TRUE)

  # hangs indefinitely if BPCells not attached (see init.R)
  counts <- BPCells::write_matrix_dir(counts_mat, dir = matrix_dir)

  # use Gene Expression modality if multiple
  if (methods::is(gene_names, "list")) {
    gene_names <- gene_names$`Gene Expression`
  }

  annotations <- data.frame(
    input = rownames(counts),
    symbol = gene_names
  )

  list(
    counts = counts,
    annotations = annotations,
    matrix_dir = matrix_dir
  )
}

read_10x_h5_file <- function(config, input_dir) {
  samples <- config$samples

  # set up parallel processing parameters
  nworkers <- min(length(samples), BATCH_POD_CPUS)
  bpparam <- BiocParallel::MulticoreParam(workers = nworkers)

  # use bplapply for parallel processing
  results <- BiocParallel::bplapply(
    samples,
    read_10x_h5_sample,
    input_dir = input_dir,
    BPPARAM = bpparam
  )

  # combine and name results from parallel processing
  combined <- combine_parallel_results(
    results,
    samples,
    extract_nested = FALSE
  )

  annot <- format_annot(combined$annot_list)

  list(
    counts_list = combined$counts_list,
    annot = annot,
    matrix_dir_list = combined$matrix_dir_list
  )
}

read_10x_h5_feature_names <- function(
  filename,
  use.names = TRUE,
  unique.features = TRUE,
  feature_type = NULL
) {
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(infile)

  # Determine feature slot based on file structure and use.names argument.
  # feature_type filtering relies on the CellRanger >= 3 "features/feature_type"
  # dataset, which the old-style ("gene_names"/"genes") format lacks.
  if (hdf5r::existsGroup(infile, "matrix")) {
    feature_slot <- ifelse(use.names, "features/name", "features/id")
    type_slot <- "features/feature_type"
  } else {
    feature_slot <- ifelse(use.names, "gene_names", "genes")
    type_slot <- NULL
  }

  output <- list()
  for (genome in genomes) {
    features <- infile[[paste0(genome, "/", feature_slot)]][]

    # keep only the requested feature type(s), matching the row subset
    # BPCells::open_matrix_10x_hdf5 applies for the same feature_type
    if (!is.null(feature_type) && !is.null(type_slot)) {
      feature_types <- infile[[paste0(genome, "/", type_slot)]][]
      features <- features[feature_types %in% feature_type]
    }

    if (unique.features) {
      features <- make.unique(names = features)
    }

    output[[genome]] <- features
  }

  infile$close_all()

  if (length(output) == 1) {
    return(output[[genome]])
  } else {
    return(output)
  }
}

#' Calls Read10X
#'
#' Cellranger outputs from V2 and V3 kits were renamed to look like V3
#' (features.tsv.gz).
#'
#' @param config experiment settings.
#'
#' Read 10X data for a single sample
#'
#' @param sample character sample ID
#' @param input_dir character path to input directory
#'
#' @return list with slots:
#'   \item{counts}{count matrix}
#'   \item{annotations}{list with annot and feature_types}
#'   \item{matrix_dir}{directory where matrix is stored}
#'
read_10x_sample <- function(sample, input_dir) {
  sample_dir <- file.path(input_dir, sample)
  sample_fpaths <- list.files(sample_dir)
  annot_fpath <- file.path(sample_dir, "features.tsv.gz")

  message("\nSample --> ", sample)
  message(
    "Reading files from ",
    sample_dir,
    " --> ",
    paste(sample_fpaths, collapse = " - ")
  )

  annotations <- read_10x_annotations(annot_fpath, sample)

  # previously used Seurat::Read10X
  counts_mat <- read_10x_bpcells(
    sample_dir,
    gene.column = annotations[["gene_column"]],
    unique.features = TRUE
  )

  # handle case where read_10x_bpcells returns a list (multiple feature types)
  if (is.list(counts_mat)) {
    slot <- "Gene Expression"
    # questionable: grab first slot if no slot named gene expression
    if (!(slot %in% names(counts_mat))) slot <- names(counts_mat)[1]
    counts_mat <- counts_mat[[slot]]
  }

  # write the matrix to a directory
  matrix_dir <- file.path(
    tempdir(),
    paste0(sample, "_matrix_dir")
  )
  unlink(matrix_dir, recursive = TRUE)
  counts <- BPCells::write_matrix_dir(counts_mat, dir = matrix_dir)

  if (methods::is(counts, "list")) {
    slot <- "Gene Expression"
    # questionable: grab first slot if no slot named gene expression
    if (!(slot %in% names(counts))) slot <- names(counts)[1]
    counts <- counts[[slot]]
  }

  message(
    sprintf(
      "Sample %s has %s genes and %s droplets.",
      sample, nrow(counts), ncol(counts)
    )
  )

  c(counts, annotations) %<-%
    filter_unnamed_features(counts, annotations, sample)

  list(
    counts = counts,
    annotations = annotations,
    matrix_dir = matrix_dir
  )
}

read_10x_files <- function(config, input_dir) {
  samples <- config$samples

  # set up parallel processing parameters
  nworkers <- min(length(samples), BATCH_POD_CPUS)
  bpparam <- BiocParallel::MulticoreParam(workers = nworkers)

  # use bplapply for parallel processing
  results <- BiocParallel::bplapply(
    samples,
    read_10x_sample,
    input_dir = input_dir,
    BPPARAM = bpparam
  )

  # combine and name results from parallel processing
  combined <- combine_parallel_results(
    results,
    samples,
    extract_nested = TRUE
  )

  c(combined$counts_list, combined$annot_list) %<-%
    normalize_annotation_types(
      combined$annot_list,
      combined$counts_list,
      combined$feature_types_list,
      samples
    )

  annot <- format_annot(combined$annot_list)
  list(
    counts_list = combined$counts_list,
    annot = annot,
    matrix_dir_list = combined$matrix_dir_list
  )
}

# mimics Seurat::Read10X
read_10x_bpcells <- function(
  data.dir,
  gene.column = 2,
  unique.features = TRUE
) {
  if (!dir.exists(data.dir)) {
    stop("Directory provided does not exist: ", data.dir)
  }

  # Define file paths (only newer format with .gz files)
  barcode.loc <- file.path(data.dir, "barcodes.tsv.gz")
  features.loc <- file.path(data.dir, "features.tsv.gz")
  matrix.loc <- file.path(data.dir, "matrix.mtx.gz")

  # Check files exist
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing. Expecting ", basename(barcode.loc))
  }
  if (!file.exists(features.loc)) {
    stop("Features file missing. Expecting ", basename(features.loc))
  }
  if (!file.exists(matrix.loc)) {
    stop("Matrix file missing. Expecting ", basename(matrix.loc))
  }

  # Read barcodes and features using data.table::fread
  cell.barcodes <- as.data.frame(
    data.table::fread(barcode.loc, header = FALSE)
  )
  cell.names <- cell.barcodes[[1]]

  feature.names <- as.data.frame(
    data.table::fread(features.loc, header = FALSE)
  )

  # Validate gene.column (3)
  fcols <- ncol(feature.names)
  if (fcols < gene.column) {
    stop(paste0(
      "gene.column was set to ", gene.column,
      " but features.tsv.gz only has ", fcols, " columns.",
      " Try setting the gene.column argument to a value <= to ", fcols, "."
    ))
  }

  row.names <- feature.names[[gene.column]]

  # Handle NA feature names (2)
  if (any(is.na(row.names))) {
    warning(
      "Some features names are NA. ",
      "Replacing NA names with ID from the opposite column requested",
      call. = FALSE, immediate. = TRUE
    )
    na.features <- which(is.na(row.names))
    replacement.column <- ifelse(gene.column == 2, 1, 2)
    row.names[na.features] <- feature.names[na.features, replacement.column]
  }

  # Read matrix using BPCells
  data <- BPCells::import_matrix_market(
    matrix.loc,
    row_names = row.names,
    col_names = cell.names
  )

  # Handle unique features
  if (unique.features) {
    rownames(data) <- make.unique(rownames(data))
  }

  # Handle multiple feature types (e.g., Gene Expression, Protein Expression)
  if (ncol(feature.names) > 2) {
    data_types <- factor(feature.names[[3]])
    lvls <- levels(data_types)

    # Protein Expression messaging
    if ("Protein Expression" %in% lvls) {
      message(
        "Xenium protein expression detected, but no scaling factor",
        " is supplied with the MEX matrices (vs HDF5). The loaded",
        " matrix is scaled by a constant from the original values."
      )
    }

    if (length(lvls) > 1) {
      message(
        "10X data contains more than one type and is being ",
        "returned as a list containing matrices of each type."
      )
    }

    # Reorder to prioritize Gene Expression
    expr_name <- "Gene Expression"
    if (expr_name %in% lvls) {
      lvls <- c(expr_name, lvls[-which(lvls == expr_name)])
    }

    data <- lapply(lvls, function(l) {
      data[data_types == l, , drop = FALSE]
    })
    names(data) <- lvls
    return(data)
  } else {
    return(data)
  }
}

#' Calls BD rhapsody data parsing functions
#'
#' Currently we only have support for sparse expression matrices.
#'
#' @inheritParams download_user_files
#'
#' @return list containing \itemize{
#'   \item{"counts_list"}{named list of dgCMatrix per sample}
#'   \item{"annot"}{data.frame with gene symbols}
#'   }
#'
read_rhapsody_files <- function(config, input_dir) {
  # if we add support for other rhapsody file types (csv matrices) we should
  # check filetypes here and dispatch accordingly.

  out <- parse_rhapsody_matrix(config, input_dir)
  return(out)
}

#' Read BD Rhapsody expression_matrix.st files
#'
#' Parses sparse count matrices generated by BD Rhapsody primary analysis.
#'
#' @inheritParams download_user_files
#' @return list containing \itemize{
#'   \item{"counts_list"}{named list of dgCMatrix per sample}
#'   \item{"annot"}{data.frame with gene symbols}
#'   }
#'
parse_rhapsody_matrix <- function(config, input_dir) {
  counts_list <- list()
  annot_list <- list()

  samples <- config$samples
  sample_options <- config$sampleOptions


  for (sample in samples) {
    sample_dir <- file.path(input_dir, sample)
    sample_fpaths <- file.path(sample_dir, file_names[["rhapsody"]])
    include_abseq <- sample_options[[sample]]$includeAbSeq

    message("\nSample --> ", sample)
    message(
      "Reading files from ",
      sample_dir,
      " --> ",
      paste(sample_fpaths, collapse = " - ")
    )

    counts <- data.table::fread(sample_fpaths)

    # catch absent DBEC column
    adjusted_col <- ifelse(
      "DBEC_Adjusted_Molecules" %in% colnames(counts),
      "DBEC_Adjusted_Molecules",
      "RSEC_Adjusted_Molecules"
    )

    # AbSeq has a "Bioproduct" col instead of "Gene" to account for proteins
    if ("Bioproduct" %in% colnames(counts)) {
      data.table::setnames(counts, "Bioproduct", "Gene")
    }

    # The ..keep is data.table syntax to grab the keep columns
    keep <- c("Cell_Index", "Gene", adjusted_col)
    counts <- counts[, ..keep]

    # convert Cell_Index to string! we parse strings from jsons a lot, and
    # having ints alone breaks things, as they are coerced to numbers
    counts[, Cell_Index := paste0("cell_", Cell_Index)]

    # clean AbSeq names, removing symbols
    counts[, Gene := gsub("[\\|:]", "_", Gene)]

    if (!include_abseq) {
      message("Remove abseq genes from sample ", sample)
      counts <- counts[!grepl("(p_?ab_?o)$", Gene, ignore.case = TRUE), ]
    }

    # we need the genes as ints to create the sparse matrix
    counts[, Gene := factor(Gene)]
    counts[, gene_i := as.integer(Gene)]

    features <- levels(counts$Gene)
    barcodes <- unique(counts$Cell_Index)

    # to create small sparse matrix, and keep original cell indices ("barcodes")
    counts[, cell_index_j := match(Cell_Index, barcodes)]

    counts <- Matrix::sparseMatrix(
      i = counts$gene_i,
      j = counts$cell_index_j,
      x = counts[[adjusted_col]],
      dimnames = list(features, barcodes)
    )

    message(
      sprintf(
        "Sample %s has %s genes and %s wells",
        sample, nrow(counts), ncol(counts)
      )
    )

    # Rhapsody data does not have ensemblIDs, but format_annot needs 2 cols
    annot <- data.frame(input = features, name = features)

    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annot
  }

  annot <- format_annot(annot_list)

  list(counts_list = counts_list, annot = annot)
}

#' Loads annotations and formats the annotation file with the column names required by format_annot
#'
#' @param annot_fpath Path to annotations file
#' @param sample Sample name
#'
#' @return
#' @export
#'
read_parse_annotations <- function(annot_fpath, sample) {
  annot <- read.delim(annot_fpath, header = TRUE, sep = ",")

  # Make the names in annot the same as the ones in the Read10x generated count matrix
  # Since Seurat uses makes.unique, we need to as well.
  # Only for the first column, at this stage column 1 are the counts matrix rownames.
  annot[, 1] <- make.unique(annot[, 1])

  # Equalizing number of columns
  # We are removing the genome column, so far it's not used.
  annot <- annot[, c(1, 2)]
  colnames(annot) <- c("input", "name")

  return(annot)
}

#' Load and parse annotations from the feature files
#'
#' This function reads the features file into a data.frame, and takes steps
#' towards its standardization, converting it to a two-column data.frame, the first column
#' being the same as the rownames of the count matrix (ideally ensemblIDs), and
#' the second either the gene symbols if they are available or a copy of the gene IDs.
#'
#' @param annot_fpath chraracter path to features file
#' @param sample character sample id
#'
#' @return list of annotations data.frame
#' @export
#'
read_10x_annotations <- function(annot_fpath, sample) {
  gene_column <- 1

  annot <- read.delim(annot_fpath, header = FALSE)

  # Remove features that are not "Gene Expression"
  if (ncol(annot) > 2 && length(grep("Gene Expression", annot$V3)) > 0) {
    annot <- annot |> dplyr::filter(V3 == "Gene Expression")
  }

  # Some feature files have less columns than expected.
  # Duplicate first column if there is only one column with gene names/ids, or
  # if there is a "Gene Expression" column
  if (ncol(annot) == 1 || annot[1, 2] == "Gene Expression") {
    annot[, 2] <- annot[, 1]
  }

  feature_types <- get_feature_types(annot)

  message("Feature types are ", feature_types, " for sample ", sample)

  # reverse annot cols if symbols are first
  if (feature_types == SYM_IDS) {
    annot[, c(1, 2)] <- annot[, c(2, 1)]
    # Read10x reads from the feature file which is not modified by us, but contains IDS
    # in the second column. Use gene_column to indicate to Read10x from where to read the ids.
    gene_column <- 2
    feature_types <- IDS_SYM
  }

  # Make the names in annot the same as the ones in the Read10x generated count matrix
  # Since Seurat uses makes.unique, we need to as well.
  # Only for the first column, at this stage column 1 are the counts matrix rownames.
  annot[, 1] <- make.unique(annot[, 1])

  # Equalizing number of columns in case there's no Gene Expression column
  annot <- annot[, c(1, 2)]
  colnames(annot) <- c("input", "name")

  list(
    "annot" = annot,
    "feature_types" = feature_types,
    "gene_column" = gene_column
  )
}

format_annot <- function(annot_list) {
  annot <- unique(do.call("rbind", annot_list))
  colnames(annot) <- c("input", "name")

  message("Deduplicating gene annotations...")

  # store original name before dedup
  gname <- annot$name
  annot$original_name <- gname

  # if empty gene name, use input as name (e.g. ENSG123)
  is_empty <- gname == ""
  gname[is_empty] <- annot$input[is_empty]
  annot$name[is_empty] <- gname[is_empty]

  # add ENSEMBL ID for genes that are duplicated (geneNameDuplicated-ENSEMBL)
  # original name kept in 'original_name' column
  is.dup <- duplicated(gname) | duplicated(gname, fromLast = TRUE)

  # We need to convert the gene inputs from _ to - bc when we create the Seurat
  # object we do this, and the match would return NA values if any
  # of the inputs still has _.
  annot$input <- gsub("_", "-", annot$input)
  annot$name[is.dup] <- paste(gname[is.dup], annot$input[is.dup], sep = " - ")

  annot <- annot[!duplicated(annot$input), ]

  rownames(annot) <- annot$input
  return(annot)
}


#' normalize annotation types
#'
#' This function makes annotations compatible between samples with different types.
#' There are three possible options for feature_types at this stage:
#' 0. SYM_SYM
#' 1. IDS_SYM
#' 2. IDS_IDS
#'
#' In the case of combination of homogeneous samples (SYM_SYM and IDS_IDS) we need
#' a "key" sample (IDS_SYM), in which case it'll try to infer symbols from ids or
#' viceversa. In case of absence of a key sample, this will throw an error, due
#' to incompatible features.
#'
#' @param annot_list list of annotation data.frames
#' @param counts_list list of count matrices
#' @param feature_types_list list of feature types
#' @param samples list of samples
#'
#' @return list with normalized annotations and count matrices
#' @export
#'
normalize_annotation_types <- function(annot_list, counts_list, feature_types_list, samples) {
  if (any(feature_types_list == IDS_IDS) &&
    any(feature_types_list == SYM_SYM) &&
    !any(feature_types_list == IDS_SYM)) {
    stop("Incompatible features detected.")
  }

  if (any(feature_types_list == IDS_SYM) &&
    any(feature_types_list == IDS_IDS) ||
    any(feature_types_list == SYM_SYM)) {
    annot_with_ids <-
      make_annot_with_ids(annot_list, feature_types_list)

    for (sample in samples) {
      sample_annot <- annot_list[[sample]]

      # Try to replace input column (currently symbols) in sample_annot
      # with ids from annot_with_ids
      if (feature_types_list[[sample]] == SYM_SYM) {
        sample_annot <- sym_to_ids(sample_annot, annot_with_ids)

        # counts were loaded with symbols, we need to replace rownames with ids
        rownames(counts_list[[sample]]) <- sample_annot$input
      }

      # Try to replace names column (currently ids) in sample_annot
      # with symbols from annot_with_ids
      else if (feature_types_list[[sample]] == IDS_IDS) {
        sample_annot <- ids_to_sym(sample_annot, annot_with_ids)
      }
      annot_list[[sample]] <- sample_annot
    }
  }
  list(counts_list = counts_list, annot_list = annot_list)
}

#' Determine the type of features in the annot data frame
#'
#' Classifies the features file columns into either ensemblIds or symbols
#'
#' @param annot data.frame read from features file
#' @return character vector indicating feature types
#'
#' @export
get_feature_types <- function(annot) {
  is_ens_col1 <- startsWith(annot[[1]], "ENS")
  pct_ens_col1 <- sum(is_ens_col1) / nrow(annot)

  is_ens_col2 <- startsWith(annot[[2]], "ENS")
  pct_ens_col2 <- sum(is_ens_col2) / nrow(annot)

  is_ens <- c(pct_ens_col1, pct_ens_col2) >= 0.5

  # reverse case, sym in first and id in second column
  if (!is_ens[1] && is_ens[2]) {
    return(SYM_IDS)
  }

  # regular cases. sum of booleans returns ints. convert to char to string match
  feature_type <- switch(as.character(sum(is_ens)),
    "0" = SYM_SYM,
    "1" = IDS_SYM,
    "2" = IDS_IDS
  )

  return(feature_type)
}


#' Make key annotation data.frame
#'
#' This function takes the annotation data.frames which contain both symbols and
#' ids and creates a key data.frame, to be used for conversion of symbols to ids
#' or viceversa.
#'
#' @param annot_list list of annotation data.frames
#' @param feature_types_list list of feature types
#'
#' @return data.frame of ids and symbols
#' @export
#'
make_annot_with_ids <- function(annot_list, feature_types_list) {
  annot_with_ids <-
    unique(do.call("rbind", annot_list[feature_types_list == IDS_SYM]))

  annot_with_ids <-
    annot_with_ids[!duplicated(annot_with_ids$input), ]

  return(annot_with_ids)
}

#' Convert symbols to ids using key data.frame
#'
#' This function tries to convert symbols to ids using the key data.frame created
#' with `make_annot_with_ids`. In case there's no ID available, it keeps the
#' existing symbol.
#'
#' @param sample_annot data.frame of annotations
#' @param annot_with_ids data.frame of annotations with IDs. Key data.frame
#'
#' @return data.frame of annotations
#' @export
#'
sym_to_ids <- function(sample_annot, annot_with_ids) {
  symbol_idx_in_annot <-
    match(sample_annot$input, annot_with_ids$name)
  symbols_with_ids_in_annot <-
    which(sample_annot$input %in% annot_with_ids$name)

  sample_annot$input[symbols_with_ids_in_annot] <-
    annot_with_ids$input[na.omit(symbol_idx_in_annot)]

  # This avoids duplicates after combining with the annotated df.
  # Leads to a mismatch in genes between samples but it seems like the best solution
  sample_annot$input <- make.unique(sample_annot$input)

  return(sample_annot)
}

#' Convert IDS to symbols using key data.frame
#'
#' This function tries to convert IDs to symbols using the key data.frame created
#' with `make_annot_with_ids`. In case there's no symbol available, it keeps the
#' existing ID
#'
#' @param sample_annot data.frame of annotations
#' @param annot_with_ids data.frame of annotations with IDs. Key data.frame
#'
#' @return data.frame of annotations
#' @export
#'
ids_to_sym <- function(sample_annot, annot_with_ids) {
  ids_index_in_annot <- match(sample_annot$input, annot_with_ids$input)
  ids_with_symbols_in_annot <-
    which(sample_annot$input %in% annot_with_ids$input)

  sample_annot$name[ids_with_symbols_in_annot] <-
    annot_with_ids$name[na.omit(ids_index_in_annot)]

  return(sample_annot)
}


#' Try to fix genes without name or remove them
#'
#' This function checks if genes whose annotations in the count matrix (`rownames(counts)`)
#' are empty, can be annotated with a different column from the features file,
#' read into the `annotations[["annot"]]` table. If it can, they are replaced everywhere
#' (both in the count matrix and in all columns of the annotations table). In case
#' the other columns do not provide a better annotation, they are removed from
#' the count matrix and annotations table.
#'
#' @param counts count matrix
#' @param annotations list of annotations data.frame, feature types and gene_column
#' @param sample character specifying current sample
#'
#' @return list of counts and annotations
#' @export
#'
filter_unnamed_features <- function(counts, annotations, sample) {
  # Check existence of empty gene symbols in count rownames
  # first will be empty; then ".1" because make.unique.
  unnamed_pat <- "^\\.[0-9]+$|^$"
  unnamed_genes_idx <- grep(unnamed_pat, rownames(counts))

  # check if unnamed genes have a gene symbol
  keep_ids <-
    !grepl(unnamed_pat, annotations$annot$name[unnamed_genes_idx])

  # replace count rownames and gene annotations
  # Using two masks to avoid storing large boolean vectors
  if (any(keep_ids)) {
    available_gene_symbols <- annotations$annot$name[unnamed_genes_idx][keep_ids]

    rownames(counts)[unnamed_genes_idx][keep_ids] <- available_gene_symbols

    rownames(annotations$annot)[unnamed_genes_idx][keep_ids] <- available_gene_symbols

    annotations$annot$input[unnamed_genes_idx][keep_ids] <- available_gene_symbols

    message(
      "Replaced ", length(which(keep_ids)),
      " empty gene names with available not empty annotation."
    )
  }

  # remove remaining rows if there are any
  if (any(!keep_ids)) {
    genes_to_remove <- unnamed_genes_idx[!keep_ids]
    counts <- counts[-genes_to_remove, ]
    annotations$annot <- annotations$annot[-genes_to_remove, ]
    message(
      "Removed ",
      length(unnamed_genes_idx[!keep_ids]),
      " genes without annotations"
    )
  }

  list("counts" = counts, "annotations" = annotations)
}

#' Read a single Visium HD sample (SpaceRanger >= 4.0 with cell segmentations)
#'
#' Expects the following files in the sample directory:
#'   - filtered_feature_cell_matrix.h5
#'   - cell_segmentations.geojson
#'   - spatial/tissue_hires_image.png
#'   - spatial/scalefactors_json.json
#'
#' @param sample character sample ID
#' @param input_dir character path to input directory
#'
#' @return list with:
#'   \item{counts}{BPCells matrix}
#'   \item{annotations}{data.frame with input and symbol columns}
#'   \item{matrix_dir}{directory where BPCells matrix is stored}
#'   \item{scdata}{Seurat spatial object created by Load10X_Spatial}
#'
read_visium_hd_sample <- function(sample, input_dir) {
  sample_dir <- file.path(input_dir, sample)
  h5_path <- file.path(sample_dir, "filtered_feature_cell_matrix.h5")

  message("\nSample --> ", sample)
  message("Reading Visium HD files from ", sample_dir)

  # Read10X_Segmentations expects segmentations
  # to be in "segmented_outputs/" subdirectory
  move_segmentations_to_subdir(sample_dir)

  results <- read_10x_h5_sample(sample, input_dir)
  segmentations <- Seurat::Read10X_Segmentations(
    image.dir = sample_dir,
    data.dir = sample_dir,
    image.name = "tissue_hires_image.png",
    assay = "RNA",
    slice = paste0(sample, ".polygons"),
    compact = FALSE
  )

  # simplify polygons and rotate if width is less than height
  results$segmentations <- segmentations |>
    simplify_segmentations() |>
    pivot_image_wide()

  # edrops and doublet scores not calculated for visium HD
  # return empty lists for compatibility with downstream functions
  results$edrops <- results$doublet_scores <- list()
  return(results)
}

move_segmentations_to_subdir <- function(sample_dir) {
  seg_dir <- file.path(sample_dir, "segmented_outputs")
  dir.create(seg_dir)
  file.copy(
    from = file.path(sample_dir, "cell_segmentations.geojson"),
    to = file.path(seg_dir, "cell_segmentations.geojson"),
    overwrite = TRUE
  )
}

simplify_segmentations <- function(segmentations) {
  # tol determined empircally for visium hd
  segmentations[["simplified.segmentations"]] <- SeuratObject::Simplify(
    coords = segmentations[["segmentation"]],
    tol = 5
  )

  # obtain compact coords
  simplified_coords <- get_simplified_coords(segmentations)

  # remove polygons to save memory
  segmentations[["simplified.segmentations"]] <- CreateSegmentation(
    coords = simplified_coords,
    compact = TRUE
  )
  segmentations[["segmentations"]] <- CreateSegmentation(
    coords = segmentations[["segmentation"]]@sf.data,
    compact = TRUE
  )

  return(segmentations)
}

pivot_image_wide <- function(segmentations) {
  image <- segmentations@image

  # retrieve image dimensions
  image_height <- dim(image)[1]
  image_width <- dim(image)[2]

  if (image_width < image_height) {
    message(
      "Image width is smaller than height.",
      "Rotating image and coordinates to match expected orientation."
    )
    segmentations <- rotate_visiumv2(segmentations)
  }
  return(segmentations)
}

# adapted from satijalab/seurat#9344
rotate_visiumv2 <- function(segmentations) {

  # extract cell coordinates and image RGB array
  bounds <- segmentations@boundaries
  coords_list <- list(
    centroids = bounds[["centroids"]]@coords,
    polygons = bounds[["segmentations"]]@sf.data,
    polygons_simple = bounds[["simplified.segmentations"]]@sf.data
  )

  image <- segmentations@image
  scale_factor <- segmentations@scale.factors$hires

  # retrieve image dimensions
  image_height <- dim(image)[1]
  image_width <- dim(image)[2]

  # retrieve cell coordinates dimensions
  coord_width <- image_width / scale_factor
  coord_height <- image_height / scale_factor

  # define center in coordinate space
  center_x <- coord_width / 2
  center_y <- coord_height / 2

  # Rotate image: transpose then flip vertically
  rotated_image <- aperm(image, c(2, 1, 3))
  rotated_image <- rotated_image[dim(rotated_image)[1]:1, , ]

  # new coordinate width/height
  new_coord_width <- coord_height
  new_coord_height <- coord_width

  # define center in new coordinate space
  new_center_x <- new_coord_width / 2
  new_center_y <- new_coord_height / 2

  coords_rot_list <- lapply(
    coords_list,
    function(coords) {
      # center cell coordinates at the origin
      coords_centered <- coords
      coords_centered[, "x"] <- coords[, "x"] - center_x
      coords_centered[, "y"] <- coords[, "y"] - center_y

      # 90-degree counterclockwise: (x, y) -> (-y, x)
      coords_rot <- coords_centered
      coords_rot[, "x"] <- coords_centered[, "y"]
      coords_rot[, "y"] <- -coords_centered[, "x"]

      # shift according to coordinates origin
      coords_rot[, "x"] <- coords_rot[, "x"] + new_center_x
      coords_rot[, "y"] <- coords_rot[, "y"] + new_center_y
      return(coords_rot)
    }
  )

  # update segmentations object
  bounds[["centroids"]]@coords <- coords_rot_list$centroids
  bounds[["segmentations"]]@sf.data <- coords_rot_list$polygons
  bounds[["simplified.segmentations"]]@sf.data <- coords_rot_list$polygons_simple
  segmentations@boundaries <- bounds
  segmentations@image <- rotated_image

  return(segmentations)
}

get_simplified_coords <- function(segmentations) {
  polygons <- segmentations[["simplified.segmentations"]]@polygons
  coords_list <- lapply(
    segmentations[["simplified.segmentations"]]@polygons,
    function(x) {
      x@Polygons[[1]]@coords
    }
  )

  cells <- Seurat::Cells(segmentations)

  ncoords <- sapply(coords_list, nrow)
  coords <- do.call(rbind, coords_list)
  colnames(coords) <- c("x", "y")
  coords <- as.data.frame(coords)
  coords$cell <- rep(cells, ncoords)
  coords
}


#' Read Visium HD files for all samples
#'
#' @param config experiment settings
#' @param input_dir character path to input directory
#'
#' @return list with counts_list, annot, matrix_dir_list, scdata_list
#' @export
#'
read_visium_hd_files <- function(config, input_dir) {
  samples <- config$samples

  nworkers <- min(length(samples), BATCH_POD_CPUS)
  bpparam <- BiocParallel::MulticoreParam(workers = nworkers)

  results <- BiocParallel::bplapply(
    samples,
    read_visium_hd_sample,
    input_dir = input_dir,
    BPPARAM = bpparam
  )
  names(results) <- samples

  counts_list <- lapply(results, function(x) x$counts)
  annot_list <- lapply(results, function(x) x$annotations)
  matrix_dir_list <- lapply(results, function(x) x$matrix_dir)
  segmentations_list <- lapply(results, function(x) x$segmentations)

  annot <- format_annot(annot_list)

  list(
    counts_list = counts_list,
    annot = annot,
    matrix_dir_list = matrix_dir_list,
    segmentations_list = segmentations_list,
    edrops = list(),
    doublet_scores = list()
  )
}


#' Read a single Xenium sample
#'
#' Expects the following files in the sample directory:
#'   - cell_feature_matrix.h5
#'   - cells.parquet
#'   - cell_boundaries.parquet
#'
#' Counts are read directly from \code{cell_feature_matrix.h5} into a disk-backed
#' BPCells matrix (reusing the 10X H5 reader), keeping only the Gene Expression
#' feature type. The count matrix is therefore never fully held in memory, unlike
#' the \code{LoadXenium} shortcut this replaces. Cell centroids and boundaries are
#' read as plain data frames from the parquet files; FOV assembly is deferred to
#' \code{create_seurat} (gem2s-5), matching the load-vs-assemble separation the
#' other technologies follow.
#'
#' @param sample character sample ID
#' @param input_dir character path to input directory
#'
#' @return list with:
#'   \item{counts}{BPCells matrix (Gene Expression feature type only)}
#'   \item{annotations}{data.frame with input and symbol columns}
#'   \item{matrix_dir}{directory where BPCells matrix is stored}
#'   \item{segmentations}{list of raw centroid/boundary/segmentation_method frames}
#'
read_xenium_sample <- function(sample, input_dir) {
  sample_dir <- file.path(input_dir, sample)

  message("\nSample --> ", sample)
  message("Reading Xenium files from ", sample_dir)

  # Counts: disk-backed BPCells matrix, Gene Expression feature type only.
  # Control/codeword feature types (Negative Control Probe/Codeword, Unassigned
  # Codeword, Genomic Control, ...) are intentionally dropped for v1 — including
  # them would corrupt QC/normalization.
  results <- read_10x_h5_sample(
    sample,
    input_dir,
    feature_type = "Gene Expression"
  )

  # Centroids + cell boundaries as raw data frames; assembled into a FOV later.
  results$segmentations <- read_xenium_segmentations(sample_dir)

  # edrops and doublet scores are not calculated for spatial technologies;
  # return empty lists for compatibility with downstream functions
  results$edrops <- results$doublet_scores <- list()

  message(
    sprintf(
      "Sample %s has %s genes and %s cells",
      sample, nrow(results$counts), ncol(results$counts)
    )
  )

  return(results)
}

#' Read Xenium centroids and cell boundaries from the raw parquet files
#'
#' Replicates the cell-level parts of \code{Seurat::ReadXenium} we support,
#' reading the parquet files directly via \code{arrow} and returning plain data
#' frames (no Seurat object). \code{cells.parquet} is read once for both the
#' centroids and the per-cell segmentation method.
#'
#' @param sample_dir character path to the sample directory
#'
#' @return list with:
#'   \item{centroids}{data.frame with columns x, y, cell}
#'   \item{segmentations}{data.frame of polygon vertices, columns cell, x, y}
#'   \item{segmentation_method}{data.frame keyed by cell, or NULL if absent}
#'   \item{transcripts}{data.frame of molecules (x, y, gene, and qv if present),
#'     or NULL when the optional transcripts.parquet is absent}
#'
read_xenium_segmentations <- function(sample_dir) {
  cells_path <- file.path(sample_dir, "cells.parquet")
  boundaries_path <- file.path(sample_dir, "cell_boundaries.parquet")

  # cells.parquet: one row per cell (x_centroid/y_centroid/cell_id, and
  # optionally segmentation_method)
  cells <- as.data.frame(arrow::read_parquet(cells_path))
  cells$cell_id <- binary_to_string(cells$cell_id)

  centroids <- data.frame(
    x = cells$x_centroid,
    y = cells$y_centroid,
    cell = cells$cell_id
  )

  # per-cell segmentation method, indexed by cell id (attached as metadata in
  # create_seurat). ReadXenium defaults this to "cell" when the column is absent.
  segmentation_method <- NULL
  if ("segmentation_method" %in% colnames(cells)) {
    segmentation_method <- data.frame(
      segmentation_method = cells$segmentation_method,
      row.names = cells$cell_id
    )
  }

  # cell_boundaries.parquet: long-format polygon vertices (cell, vertex_x,
  # vertex_y) -> columns cell, x, y
  boundaries <- as.data.frame(arrow::read_parquet(boundaries_path))
  colnames(boundaries) <- c("cell", "x", "y")
  boundaries$cell <- binary_to_string(boundaries$cell)

  # transcripts.parquet is an OPTIONAL Xenium input: one row per molecule
  # (x_location/y_location/feature_name, and optionally qv). Read it as a plain
  # frame here; the molecule pyramid is built straight from it in gem2s-7. Don't
  # build a Seurat object / CreateMolecules here (loader stays assembly-free).
  transcripts <- read_xenium_transcripts(sample_dir)

  list(
    centroids = centroids,
    segmentations = boundaries,
    segmentation_method = segmentation_method,
    transcripts = transcripts
  )
}

#' Read the optional Xenium transcripts (molecules) parquet
#'
#' \code{transcripts.parquet} is not a required input; returns NULL when absent.
#' Columns \code{x_location}/\code{y_location}/\code{feature_name} (+ \code{qv}
#' if present) -> \code{x}/\code{y}/\code{gene} (+ \code{qv}). \code{feature_name}
#' is decoded the same way as cell ids. The installed \code{ReadXenium} does not
#' apply a QV threshold; the filter is applied at pyramid-build time (gem2s-7).
#'
#' @param sample_dir character path to the sample directory
#'
#' @return data.frame with columns x, y, gene (+ qv), or NULL if the file absent
read_xenium_transcripts <- function(sample_dir) {
  transcripts_path <- file.path(sample_dir, "transcripts.parquet")
  if (!file.exists(transcripts_path)) {
    message("No transcripts.parquet found; skipping molecule pyramid.")
    return(NULL)
  }

  message("Reading Xenium transcripts from ", transcripts_path)
  tx <- as.data.frame(arrow::read_parquet(transcripts_path))

  transcripts <- data.frame(
    x = tx$x_location,
    y = tx$y_location,
    gene = binary_to_string(tx$feature_name)
  )
  if ("qv" %in% colnames(tx)) {
    transcripts$qv <- tx$qv
  }

  transcripts
}

#' Decode Xenium hex-encoded binary cell ids
#'
#' Xenium cell ids may be stored as a list of hex-encoded bytes; pass them
#' through unchanged when they are already character (the common parquet case).
#' Mirrors the helper inside \code{Seurat::ReadXenium}.
binary_to_string <- function(arrow_binary) {
  if (typeof(arrow_binary) == "list") {
    unlist(lapply(
      arrow_binary,
      function(x) rawToChar(as.raw(strtoi(x, 16L)))
    ))
  } else {
    arrow_binary
  }
}


#' Read Xenium files for all samples
#'
#' @param config experiment settings
#' @param input_dir character path to input directory
#'
#' @return list with counts_list, annot, matrix_dir_list, segmentations_list
#' @export
#'
read_xenium_files <- function(config, input_dir) {
  samples <- config$samples

  nworkers <- min(length(samples), BATCH_POD_CPUS)
  bpparam <- BiocParallel::MulticoreParam(workers = nworkers)

  results <- BiocParallel::bplapply(
    samples,
    read_xenium_sample,
    input_dir = input_dir,
    BPPARAM = bpparam
  )
  names(results) <- samples

  counts_list <- lapply(results, function(x) x$counts)
  annot_list <- lapply(results, function(x) x$annotations)
  matrix_dir_list <- lapply(results, function(x) x$matrix_dir)
  segmentations_list <- lapply(results, function(x) x$segmentations)

  annot <- format_annot(annot_list)

  list(
    counts_list = counts_list,
    annot = annot,
    matrix_dir_list = matrix_dir_list,
    segmentations_list = segmentations_list,
    edrops = list(),
    doublet_scores = list()
  )
}
