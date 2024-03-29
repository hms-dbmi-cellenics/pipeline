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
load_user_files <- function(input, pipeline_config, prev_out, input_dir = INPUT_DIR) {
  message("Loading user files...")
  check_prev_out(prev_out, "config")

  # destructure previous output
  config <- prev_out$config

  technology <- ifelse(config$input$type %in% c("rhapsody", "10x_h5", "parse"), config$input$type, "10x")

  read_fun <- switch(technology,
    "10x" = read_10x_files,
    "rhapsody" = read_rhapsody_files,
    "10x_h5" = read_10x_h5_file,
    "parse" = read_parse_files
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

#' Read Parse files. Calls readMTX
#'
#' @param config experiment settings
#' @param input_dir
#'
#' @return
#' @export
#'
read_parse_files <- function(config, input_dir) {
  counts_list <- list()
  annot_list <- list()
  feature_types_list <- list()

  samples <- config$samples

  for (sample in samples) {
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


    # We use readMtx instead of Seurat::ReadParse because the feature.column needs to be 1 instead of 2.
    counts <- Seurat::ReadMtx(
      mtx = mtx_path, cells = barcodes_path, features = features_path,
      cell.column = 1, feature.column = 1, cell.sep = ",",
      feature.sep = ",", skip.cell = 1, skip.feature = 1, mtx.transpose = TRUE
    )

    message(
      sprintf(
        "Sample %s has %s genes and %s barcodes",
        sample, nrow(counts), ncol(counts)
      )
    )

    c(counts, annotations) %<-% filter_unnamed_features(counts, annotations, sample)

    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annotations
  }

  annot <- format_annot(annot_list)

  return(list(counts_list = counts_list, annot = annot))
}

#' Read h5 file
#'
#' Calls read10x_h5
#'
#' @param config list of configuration parameters
#' @param input_dir character input dir
#'
#' @return list with counts and annotations
#' @export
#'
read_10x_h5_file <- function(config, input_dir) {
  counts_list <- list()
  annot_list <- list()

  samples <- config$samples

  for (sample in samples) {
    sample_dir <- file.path(input_dir, sample)
    sample_fpaths <- list.files(sample_dir)
    sample_counts_path <- file.path(sample_dir, sample_fpaths[[1]])

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

    ungzipped_counts_path <- R.utils::gunzip(sample_counts_path)

    counts_names <- Seurat::Read10X_h5(ungzipped_counts_path)
    counts <- Seurat::Read10X_h5(ungzipped_counts_path, use.names = FALSE)

    # use Gene Expression modality if multiple
    if (methods::is(counts, "list")) {
      counts_names <- counts_names$`Gene Expression`
      counts <- counts$`Gene Expression`
    }

    gene_names <- row.names(counts_names)

    annotations <-
      data.frame(input = rownames(counts), symbol = gene_names)
    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annotations
  }

  annot <- format_annot(annot_list)

  return(list(counts_list = counts_list, annot = annot))
}

#' Calls Read10X
#'
#' Cellranger outputs from V2 and V3 kits were renamed to look like V3
#' (features.tsv.gz).
#'
#' @param config experiment settings.
#'
read_10x_files <- function(config, input_dir) {
  counts_list <- list()
  annot_list <- list()
  feature_types_list <- list()

  samples <- config$samples

  for (sample in samples) {
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
    counts <- Seurat::Read10X(sample_dir, gene.column = annotations[["gene_column"]], unique.features = TRUE)

    if (is(counts, "list")) {
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

    c(counts, annotations) %<-% filter_unnamed_features(counts, annotations, sample)

    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annotations[["annot"]]
    feature_types_list[[sample]] <- annotations[["feature_types"]]
  }

  c(counts_list, annot_list) %<-% normalize_annotation_types(annot_list, counts_list, feature_types_list, samples)
  annot <- format_annot(annot_list)

  return(list(counts_list = counts_list, annot = annot))
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

  return(list(counts_list = counts_list, annot = annot))
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
    annot <- annot %>% dplyr::filter(V3 == "Gene Expression")
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

  return(list("annot" = annot, "feature_types" = feature_types, "gene_column" = gene_column))
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
    any(feature_types_list == IDS_IDS) || any(feature_types_list == SYM_SYM)) {
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
  return(list(counts_list = counts_list, annot_list = annot_list))
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

  return(list("counts" = counts, "annotations" = annotations))
}
