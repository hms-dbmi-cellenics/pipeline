#' Read input folder of 10x data
#'
#' @inheritParams download_cellranger_files
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
  experiment_id <- input$experimentId

  output <- c(prev_out, call_read10x(config, input_dir, experiment_id))

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
call_read10x <- function(config, input_dir, experiment_id) {
  counts_list <- list()
  annot_list <- list()

  samples <- config$samples
  message("Samples to include in the analysis:\n- ", paste(samples, collapse = "\n- "))
  message("Loading 10x data set from input folder.")


  for (sample in samples) {
    sample_dir <- file.path(input_dir, sample)
    sample_fpaths <- list.files(sample_dir)

    spatial <- NULL
    is.spatial <- 'aligned_fiducials.jpg.gz' %in% sample_fpaths
    if (is.spatial) {
      # process spatial files
      setup_visium_files(sample_dir)
      spatial <- process_spaceranger_files(sample_dir)
      spatial <- prepare_spatial_experiment(spatial, sample, experiment_id)

      # mock RNA-seq data to run GEM2S on
      unlink(sample_dir, recursive = TRUE)
      mock_10x_counts(sample_dir)
    }

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

    message(
      sprintf("Sample %s has %s genes and %s droplets.", sample, nrow(counts), ncol(counts))
    )

    counts_list[[sample]] <- counts
    annot_list[[sample]] <- annot
  }

  annot <- format_annot(annot_list)

  return(list(counts_list = counts_list, annot = annot, spatial = spatial))
}

mock_10x_counts <- function(sample_dir) {
  sce <- scDblFinder::mockDoubletSCE()
  my.counts <- SingleCellExperiment::counts(sce)
  my.counts <- as(my.counts, "dgCMatrix")

  ngenes <- nrow(my.counts)
  gene.ids <- paste0("ENSG0000", seq_len(ngenes))
  gene.symb <- paste0("GENE", seq_len(ngenes))
  cell.ids <- colnames(my.counts)

  # Writing this to file:
  DropletUtils::write10xCounts(
    sample_dir,
    my.counts,
    gene.id = gene.ids,
    gene.symbol = gene.symb,
    barcodes = cell.ids,
    version = '3')
}

setup_visium_files <- function(sample_dir) {
  sample_fpaths <- list.files(sample_dir)
  for (sample_fpath in sample_fpaths) {
    R.utils::gunzip(file.path(sample_dir, sample_fpath))
  }

  sample_fpaths <- list.files(sample_dir)
  spatial_files <- grep('[.]h5$', sample_fpaths, value = TRUE, invert = TRUE)
  spatial_dir <- file.path(sample_dir, 'spatial')
  dir.create(spatial_dir)
  file.copy(file.path(sample_dir, spatial_files),
            file.path(spatial_dir, spatial_files))

  unlink(file.path(sample_dir, spatial_files))
}


prepare_spatial_experiment <- function(scdata, sample_name, experiment_id) {

  # add cell ids
  scdata$cells_id <- 0:(ncol(scdata) - 1)

  # move SCT assay to RNA for expression plots
  scdata[['RNA']] <- scdata[['SCT']]

  # add samples/other things for platform
  scdata$samples <- sample_name
  scdata@misc$experimentId <- experiment_id
  scdata@misc$color_pool <- get_color_pool()
  scdata@misc$numPCs <- 30

  scdata <- add_dispersions(scdata)

  return(scdata)
}

process_spaceranger_files <- function(data.dir) {

  # load and normalize
  scdata <- Seurat::Load10X_Spatial(data.dir, use.names = FALSE)
  scdata <- add_annotation(scdata, data.dir)

  scdata <- Seurat::SCTransform(scdata, assay = "Spatial", verbose = TRUE)

  # cluster and run umap
  scdata <- Seurat::RunPCA(scdata, assay = "SCT", verbose = FALSE)
  scdata <- Seurat::FindNeighbors(scdata, reduction = "pca", dims = 1:30)
  scdata <- Seurat::FindClusters(scdata, verbose = FALSE)
  scdata <- Seurat::RunUMAP(scdata, reduction = "pca", dims = 1:30)

  return(scdata)
}

add_annotation <- function(scdata, data.dir) {
  # get gene symbols
  h5file <- file.path(data.dir, 'filtered_feature_bc_matrix.h5')
  infile <- hdf5r::H5File$new(filename = h5file, mode = "r")
  features <- infile[[paste0('matrix/features/name')]][]

  annot <- data.frame(
    input = row.names(scdata),
    name = features,
    row.names = row.names(scdata)
  )

  scdata@misc[["gene_annotations"]] <- annot
  return(scdata)
}

add_dispersions <- function(scdata) {
  scdata <- Seurat::FindVariableFeatures(scdata)
  vars <- Seurat::HVFInfo(object = scdata, selection.method = 'vst')

  # fake for sctransform
  colnames(vars) <- c('mean', 'variance', 'variance.standardized')

  annotations <- scdata@misc[["gene_annotations"]]
  vars$SYMBOL <- annotations$name[match(rownames(vars), annotations$input)]
  vars$ENSEMBL <- rownames(vars)
  scdata@misc[["gene_dispersion"]] <- vars
  return(scdata)
}

format_annot <- function(annot_list) {
  annot <- unique(do.call("rbind", annot_list))
  annot <- annot[, c(1, 2)]
  colnames(annot) <- c("input", "name")

  message("Deduplicating gene annotations...")

  # add ENSEMBL ID for genes that are duplicated (geneNameDuplicated-ENSEMBL)
  # original name kept in 'original_name' column
  gname <- annot$name
  annot$original_name <- gname
  is.dup <- duplicated(gname) | duplicated(gname, fromLast = TRUE)

  # We need to convert the gene inputs from _ to - bc when we create the Seurat object we do this, and the match would return NA values if any of the inputs still has _.
  annot$input <- gsub("_", "-", annot$input)
  annot$name[is.dup] <- paste(gname[is.dup], annot$input[is.dup], sep = " - ")

  annot <- annot[!duplicated(annot$input), ]

  rownames(annot) <- annot$input
  return(annot)
}
