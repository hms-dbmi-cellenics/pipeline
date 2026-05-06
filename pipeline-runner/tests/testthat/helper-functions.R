#' Fix time stamps in seurat object
#'
#' Seurat adds logs to certain command runs, with time stamps that make snapshot
#' tests fail. This function replaces them with a fixed datetime object.
#'
#' @param scdata seurat object
#'
#' @return seurat object with fixed time stamps
#'
clean_timestamp <- function(scdata) {
  fixed_datetime <- as.POSIXct("1991-12-19 05:23:00", tz = "UTC")

  for (slot in names(scdata@commands)) {
    scdata@commands[[slot]]@time.stamp <- fixed_datetime
  }

  if ("ingestionDate" %in% names(scdata@misc)) {
    scdata@misc$ingestionDate <- fixed_datetime
  }

  return(scdata)
}


# Seurat object commands are sometimes functions which
# will have a new environment for every call which breaks snapshots
remove_commands_functions <- function(data) {
  command_names <-  names(data@commands)
  for (command_name in command_names) {
    param_names <- names(data@commands[[command_name]]@params)

    for (param_name in param_names) {
      param <- data@commands[[command_name]]@params[[param_name]]

      if (methods::is(param, "function"))
        data@commands[[command_name]]@params[[param_name]] <- NULL
    }
  }
  return(data)
}

#' Check if BPCELLS testing mode is enabled
#'
#' Returns TRUE if the BPCELLS environment variable is set to "true".
#' Used to conditionally enable BPCells matrix conversion in tests.
#'
#' @return logical TRUE if BPCELLS="true", FALSE otherwise
#'
is_bpcells <- function() {
  Sys.getenv("BPCELLS") == "true"
}

skip_if_bpcells <- function() {
  if (is_bpcells()) {
    skip("Skipping test because BPCELLS=true")
  }
}

#' Convert a matrix to disk-backed BPCells format if BPCELLS env var is set
#'
#' Simple wrapper that converts a single matrix to BPCells IterableMatrix
#' for testing with disk-backed matrices. Only converts if BPCELLS="true".
#'
#' @param mat A matrix (dgCMatrix, matrix, etc.)
#'
#' @return BPCells IterableMatrix if BPCELLS=true, otherwise original matrix
#'
#' @examples
#' \dontrun{
#'   counts <- matrix(1:10, nrow = 2)
#'   counts <- maybe_bpcells(counts)
#' }
#'
maybe_bpcells <- function(mat, matrix_dir) {
  if (!is_bpcells()) {
    return(mat)
  }

  # Convert to BPCells
  BPCells::write_matrix_dir(mat, dir = matrix_dir)
}

mock_scdata_list <- function(
  rename_genes = c(),
  n_rep = 1,
  with_outlier = FALSE
) {

  counts <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
  rownames(counts)[grep("^RP[LS]", rownames(counts))] <- "SOX1"

  if (with_outlier) {
    counts[, 1] <- 0
    counts[1:10, 1] <- 100
  }

  # replicate matrix columns n times to create a bigger mock dataset
  counts <- do.call("cbind", replicate(n_rep, counts, simplify = FALSE))
  colnames(counts) <- make.unique(colnames(counts))
  ncells_each <- ncol(counts) / 2

  if (length(rename_genes) > 0) {
    # rename some genes to match cell cycle genes
    some_genes <- sample(nrow(counts), length(rename_genes))
    rownames(counts)[some_genes] <- rename_genes
  }

  # setup samples
  sample_ids <- c("123abc", "123def")
  samples <- rep(sample_ids, each = ncells_each)
  colnames(counts) <- paste0(samples, colnames(counts))

  gene_annotations <- data.frame(
    input = rownames(counts),
    name = rownames(counts)
  )

  counts <- as(as.matrix(counts), "dgCMatrix")

  # Set seed once before generating all samples
  set.seed(123)

  # create seurat objects and add metadata using loop
  scdata_list <- list()

  for (i in seq_along(sample_ids)) {

    # split counts for this sample and convert to bpcells if testing
    start_col <- (i - 1) * ncells_each + 1
    end_col <- i * ncells_each
    counts_split <- maybe_bpcells(
      counts[, seq(start_col, end_col)],
      withr::local_tempfile(.local_envir = parent.frame())
    )

    scdata <- Seurat::CreateSeuratObject(counts = counts_split)

    scdata$cells_id <- seq((i - 1) * ncells_each, i * ncells_each - 1)
    scdata$samples <- sample_ids[i]
    scdata@misc$gene_annotations <- gene_annotations

    # add doublet scores and class
    scdata$doublet_scores <- 0.01
    scdata$doublet_class <- "singlet"

    # add mitochondrial percent
    scdata$percent.mt <- rnorm(ncol(scdata), mean = 6)

    # add empty drops FDR
    scdata$emptyDrops_FDR <- 0.009

    scdata_list[[sample_ids[i]]] <- scdata
  }

  return(scdata_list)
}

mock_ids <- function(scdata_list) {
  # Generate IDs based on scdata_list structure
  sample_names <- names(scdata_list)
  cells_id <- list()
  total_cells <- 0

  for (sample_name in sample_names) {
    ncells <- ncol(scdata_list[[sample_name]])
    cells_id[[sample_name]] <- total_cells:(total_cells + ncells - 1)
    total_cells <- total_cells + ncells
  }

  return(cells_id)
}

materialize_scdata_list <- function(res) {
  # materialize bpcells matrices for transmission to test-qc.R
  res$output$scdata_list <- lapply(
    res$output$scdata_list,
    function(scdata) {
      scdata[["RNA"]]$counts <- as(scdata[["RNA"]]$counts, "dgCMatrix")
      return(scdata)
    }
  )
  return(res)
}

materialize_counts_list <- function(res) {
  res$output$counts_list <- lapply(
    res$output$counts_list,
    function(x) {
      as(x, "dgCMatrix")
    }
  )
  return(res)
}

materialize_res <- function(res) {
  if (!is.null(res$output$counts_list)) {
    res <- materialize_counts_list(res)
  }
  if (!is.null(res$output$scdata_list)) {
    res <- materialize_scdata_list(res)
  }

  if (!is.null(res$output$matrix_dir_list)) {
    res$output$matrix_dir_list <- list()
  }
  return(res)
}