mock_scdata <- function(data_dir, obj2s_type, remove_reductions = NULL) {
  data("pbmc_small", package = "SeuratObject")
  pbmc_small$seurat_clusters <- pbmc_small$RNA_snn_res.1

  for (red in remove_reductions) {
    pbmc_small@reductions[[red]] <- NULL
  }

  if (obj2s_type == "sce_object") {
    pbmc_small <- Seurat::as.SingleCellExperiment(pbmc_small)
  } else if (obj2s_type == "anndata_object") {
    pbmc_small <- Seurat::as.SingleCellExperiment(pbmc_small)

    # setup X with counts under '/X'
    SingleCellExperiment::logcounts(pbmc_small) <- NULL
    SummarizedExperiment::assayNames(pbmc_small, "counts") <- "X"

    # change reduction names to match AnnData convention
    red_names <- SingleCellExperiment::reducedDimNames(pbmc_small)
    if (length(red_names)) {
      SingleCellExperiment::reducedDimNames(pbmc_small) <- paste0(
        "X_", red_names
      )
    }

    # write to h5ad
    anndata <- reticulate::import("anndata")
    adata <- zellkonverter::SCE2AnnData(pbmc_small)

    anndata$AnnData$write_h5ad(adata, file.path(data_dir, "adata.h5ad"))
    return(pbmc_small)
  }

  saveRDS(pbmc_small, file.path(data_dir, "r.rds"))
  return(pbmc_small)
}

obj2s_types <- c("seurat_object", "sce_object", "anndata_object")


test_that("load_obj2s_file works", {
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    mock_scdata(data_dir, obj2s_type)
    prev_out <- list(
      config = list(samples = "pbmc_small",
      input = list(type = obj2s_type))
    )
    res <- load_obj2s_file(
      input = NULL,
      pipeline_config = NULL,
      prev_out = prev_out,
      input_dir = input_dir
    )
    expect_type(res, "list")
    expect_s4_class(res$output$scdata, "Seurat")
  }
})

test_that("load_obj2s_file has the correct structure", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    orig_scdata <- mock_scdata(data_dir, obj2s_type)
    mock_scdata(data_dir, obj2s_type)
    prev_out <- list(
      config = list(samples = "pbmc_small",
      input = list(type = obj2s_type))
    )
    res <- load_obj2s_file(
      input = NULL,
      pipeline_config = NULL,
      prev_out = prev_out,
      input_dir = input_dir
    )
    scdata <- res$output$scdata

    # misc has dispersions and annotations
    expect_setequal(
      names(scdata@misc),
      c("gene_dispersion", "gene_annotations")
    )

    # dispersions have expected columns
    dispersions <- scdata@misc$gene_dispersion
    expect_setequal(
      colnames(dispersions),
      c("mean", "variance", "variance.expected", "variance.standardized", "SYMBOL", "ENSEMBL")
    )

    # annotations have expected columns
    annot <- scdata@misc$gene_annotations
    expect_setequal(colnames(annot), c("input", "name", "original_name"))

    # only the default and 'pca' reduction is present
    expect_equal(c("tsne", "pca"), names(scdata@reductions))

    if (obj2s_type == "seurat_object") {
      # all metadata is preserved and active.ident is added
      expect_equal(
        colnames(scdata@meta.data),
        c(colnames(orig_scdata@meta.data), "active.ident")
      )
    } else if (obj2s_type %in% c("sce_object", "anndata_object")) {
      # all metadata is preserved
      expect_equal(colnames(scdata@meta.data), colnames(orig_scdata@colData))
    }
  }
})

test_that("load_obj2s_file fails if the file isn't a .rds file", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)


  for (obj2s_type in obj2s_types) {
    dataset_fname <- switch(obj2s_type,
      "seurat_object" = "r.rds",
      "sce_object" = "r.rds",
      "anndata_object" = "adata.h5ad"
    )

    mock_scdata(data_dir, obj2s_type)
    # overwrite with mtcars csv
    write.csv(mtcars, file.path(data_dir, dataset_fname))

    prev_out <- list(
      config = list(
        samples = "pbmc_small",
        input = list(type = obj2s_type)
      )
    )
    expect_error(
      load_obj2s_file(
        input = NULL,
        pipeline_config = NULL,
        prev_out = prev_out,
        input_dir = input_dir
      ),
      regexp = "ERROR_OBJ2S_READ"
    )
  }
})

test_that("load_obj2s_file fails if there is no RNA assay for seurat_object", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, obj2s_type = "seurat_object")

  # rename RNA assay
  names(orig_scdata@assays) <- "blah"
  saveRDS(orig_scdata, file.path(data_dir, "r.rds"))

  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "seurat_object")
    )
  )
  expect_error(
    load_obj2s_file(
      input = NULL,
      pipeline_config = NULL,
      prev_out = prev_out,
      input_dir = input_dir
    ),
    regexp = "ERROR_OBJ2S_COUNTS"
  )
})

test_that("load_obj2s_file does not fail if there is no seurat_clusters column in metadata", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, "seurat_object")

  # remove seurat_clusters
  orig_scdata$seurat_clusters <- NULL
  saveRDS(orig_scdata, file.path(data_dir, "r.rds"))

  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "seurat_object")
    )
  )
  expect_error(
    load_obj2s_file(
      input = NULL,
      pipeline_config = NULL,
      prev_out = prev_out,
      input_dir = input_dir
    ),
    NA
  )
})

test_that("load_obj2s_file fails if there is no reduction", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    scdata <- mock_scdata(
      data_dir, obj2s_type, remove_reductions = c("pca", "tsne")
    )
    prev_out <- list(
      config = list(
        samples = "pbmc_small",
        input = list(type = obj2s_type)
      )
    )

    expect_error(
      load_obj2s_file(
        input = NULL,
        pipeline_config = NULL,
        prev_out = prev_out,
        input_dir = input_dir
      ),
      regexp = "ERROR_OBJ2S_REDUCTION"
    )
  }
})

test_that("load_obj2s_file does not fail if there is no pca reduction", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    scdata <- mock_scdata(data_dir, obj2s_type, remove_reductions = "pca")
    prev_out <- list(
      config = list(
        samples = "pbmc_small",
        input = list(type = obj2s_type)
      )
    )

    expect_error(
      load_obj2s_file(
        input = NULL,
        pipeline_config = NULL,
        prev_out = prev_out,
        input_dir = input_dir
      ),
      NA
    )
  }
})


test_that("load_obj2s_file generates HVFInfo if it is not present", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, "seurat_object")


  # mess up meta.features (used by HVFInfo)
  colnames(orig_scdata@assays$RNA@meta.features) <- paste0(
    "blah.",
    colnames(orig_scdata@assays$RNA@meta.features)
  )
  saveRDS(orig_scdata, file.path(data_dir, "r.rds"))

  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "seurat_object")
    )
  )
  expect_error(
    res <- load_obj2s_file(
      input = NULL,
      pipeline_config = NULL,
      prev_out = prev_out,
      input_dir = input_dir
    ),
    NA
  )

  scdata <- res$output$scdata
  disp_colnames <- c(
    "mean",
    "variance",
    "variance.expected",
    "variance.standardized",
    "SYMBOL",
    "ENSEMBL"
  )
  expect_true(methods::is(scdata@misc$gene_dispersion, "data.frame"))
  expect_setequal(colnames(scdata@misc$gene_dispersion), disp_colnames)
})

test_that("load_obj2s_file normalizes logcounts from counts instead of using the uploaded logcounts", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, "seurat_object")

  # replace logcounts with counts
  orig_counts <- orig_scdata[["RNA"]]$counts
  orig_scdata[["RNA"]]$data <- orig_counts

  saveRDS(orig_scdata, file.path(data_dir, "r.rds"))
  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "seurat_object")
    )
  )
  res <- load_obj2s_file(
    input = NULL,
    pipeline_config = NULL,
    prev_out = prev_out,
    input_dir = input_dir
  )
  scdata <- res$output$scdata

  # logcounts are log-normalized counts, as for count matrix uploads
  counts <- as(scdata[["RNA"]]$counts, "dgCMatrix")
  expect_equal(counts, orig_counts)
  expect_equal(
    as(scdata[["RNA"]]$data, "dgCMatrix"),
    Seurat::NormalizeData(orig_counts),
    tolerance = 1e-6
  )
})

test_that("load_obj2s_file calculates logcounts for sce_object if not present", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, "sce_object")

  # remove logcounts
  SingleCellExperiment::logcounts(orig_scdata) <- NULL

  saveRDS(orig_scdata, file.path(data_dir, "r.rds"))
  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "sce_object")
    )
  )
  res <- load_obj2s_file(
    input = NULL,
    pipeline_config = NULL,
    prev_out = prev_out,
    input_dir = input_dir
  )
  scdata <- res$output$scdata

  # logcounts are log-transformed counts
  orig_counts <- SingleCellExperiment::counts(orig_scdata)
  counts <- as(scdata[["RNA"]]$counts, "dgCMatrix")
  expect_equal(counts, orig_counts)
  expect_equal(
    as(scdata[["RNA"]]$data, "dgCMatrix"),
    Seurat::NormalizeData(orig_counts),
    tolerance = 1e-6
  )
})

test_that("load_obj2s_file works when default reducion for seurat_object is different than umap or tsne", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, "seurat_object")
  withr::defer(unlink(data_dir, recursive = TRUE))

  # simulate the condition where renaming is needed
  names(orig_scdata@reductions)[names(orig_scdata@reductions) == "tsne"] <- "tsne.ref"
  orig_scdata@reductions$tsne <- orig_scdata@reductions$tsne.ref
  saveRDS(orig_scdata, file.path(data_dir, "r.rds"))

  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "seurat_object")
    )
  )
  res <- load_obj2s_file(
    input = NULL,
    pipeline_config = NULL,
    prev_out = prev_out,
    input_dir = input_dir
  )
  scdata <- res$output$scdata

  expect_true("tsne" %in% names(scdata@reductions))
  expect_equal("tsne", SeuratObject::DefaultDimReduc(scdata))
})

test_that("reconstruct_sce loads an sce_object serialized with a newer Bioconductor (rowRanges Seqinfo)", {
  # sce_seqinfo.rds was saved with Bioconductor >= 3.21, where the Seqinfo S4
  # class lives in a standalone 'Seqinfo' package not installed in this image.
  # Its rowRanges reference that package, so touching them via updateObject()
  # errored with "unable to find required package 'Seqinfo'". reconstruct_sce
  # must drop rowRanges before any accessor runs. See mock_data/obj2s_newer_bioc.
  # Skipped on CI: CI runs under covr::package_coverage(), which itself executes
  # the suite in a subprocess. Nesting the callr process below inside it leaves a
  # child that covr can't reap at shutdown ("unable to terminate some child
  # processes"). This still runs in local `make test`, like the snapshot tests.
  skip_on_ci()

  fixture <- testthat::test_path(
    "mock_data", "obj2s_newer_bioc", "sce_seqinfo.rds"
  )

  # Run in a fresh R process: the failure only reproduces with a cold S4 method
  # dispatch cache, as in production. Within the shared test session earlier
  # tests warm dispatch for the Seqinfo/GRanges hierarchy and mask the error.
  reductions <- callr::r(
    function(pkg, fpath) {
      pkgload::load_all(pkg, quiet = TRUE)
      scdata <- reconstruct_sce(fpath)
      names(scdata@reductions)
    },
    args = list(pkg = pkgload::pkg_path(), fpath = normalizePath(fixture))
  )

  expect_true("tsne" %in% reductions)
})

test_that("load_obj2s_file works when reduction for sce_object is a string match for umap or tsne", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  orig_scdata <- mock_scdata(data_dir, "sce_object")
  withr::defer(unlink(data_dir, recursive = TRUE))

  # simulate the condition where renaming is needed
  SingleCellExperiment::reducedDimNames(orig_scdata) <- c("PCA", "TSNE.REF")
  saveRDS(orig_scdata, file.path(data_dir, "r.rds"))

  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "sce_object")
    )
  )
  res <- load_obj2s_file(
    input = NULL,
    pipeline_config = NULL,
    prev_out = prev_out,
    input_dir = input_dir
  )
  scdata <- res$output$scdata

  expect_true("tsne" %in% names(scdata@reductions))
  expect_equal("tsne", SeuratObject::DefaultDimReduc(scdata))
})

test_that("update_reduction_name correctly renames dimensionality reductions", {
  # setup
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  scdata <- mock_scdata(data_dir, "seurat_object")
  withr::defer(unlink(data_dir, recursive = TRUE))

  # simulate the condition where renaming is needed
  names(scdata@reductions)[names(scdata@reductions) == "tsne"] <- "tsne.ref"
  scdata@reductions$tsne <- scdata@reductions$tsne.ref

  red_name <- SeuratObject::DefaultDimReduc(scdata)
  red_match <- grep("umap|tsne", red_name, value = TRUE)
  is_umap <- grepl("umap", red_match)
  is_tsne <- grepl("tsne", red_match)
  new_red_name <- ifelse(
    is_umap,
    "umap",
    ifelse(is_tsne, "tsne", NA)
  )

  updated_scdata <- update_reduction_name(scdata, red_name, new_red_name)
  expect_equal(SeuratObject::DefaultDimReduc(updated_scdata), new_red_name)
  expect_true("tsne.ori" %in% names(updated_scdata@reductions))
})

test_that("load_obj2s_file works with an Assay5 RNA for seurat_object", {
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  scdata <- mock_scdata(data_dir, "seurat_object")

  # convert to Assay5
  scdata[["RNA"]] <- as(scdata[["RNA"]], "Assay5")
  saveRDS(scdata, file.path(data_dir, "r.rds"))

  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "seurat_object")
    )
  )
  res <- load_obj2s_file(
    input = NULL,
    pipeline_config = NULL,
    prev_out = prev_out,
    input_dir = input_dir
  )
  expect_type(res, "list")
  expect_s4_class(res$output$scdata, "Seurat")
  expect_s4_class(res$output$scdata[["RNA"]], "Assay5")
})

test_that("load_obj2s_file works with split Assay5 RNA for seurat_object", {
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  scdata <- mock_scdata(data_dir, "seurat_object")

  # convert to Assay5 and split
  scdata[["RNA"]] <- as(scdata[["RNA"]], "Assay5")
  scdata$batch <- sample(
    c("batchA", "batchB", "batchC"),
    ncol(scdata),
    replace = TRUE
  )
  scdata[["RNA"]] <- split(scdata[["RNA"]], f = scdata$batch)

  # was split
  expect_setequal(
    names(scdata[["RNA"]]@layers),
    c(
      "data.batchB",
      "data.batchC",
      "data.batchA",
      "scale.data",
      "counts.batchB",
      "counts.batchC",
      "counts.batchA"
    )
  )

  saveRDS(scdata, file.path(data_dir, "r.rds"))

  prev_out <- list(
    config = list(
      samples = "pbmc_small",
      input = list(type = "seurat_object")
    )
  )
  res <- load_obj2s_file(
    input = NULL,
    pipeline_config = NULL,
    prev_out = prev_out,
    input_dir = input_dir
  )
  expect_type(res, "list")
  expect_s4_class(res$output$scdata, "Seurat")
  expect_s4_class(res$output$scdata[["RNA"]], "Assay5")

  # is joined back
  expect_setequal(
    names(res$output$scdata[["RNA"]]@layers),
    c("counts", "data")
  )
})

load_mock_obj2s <- function(input_dir, obj2s_type) {
  prev_out <- list(
    config = list(samples = "pbmc_small", input = list(type = obj2s_type))
  )
  res <- load_obj2s_file(
    input = NULL,
    pipeline_config = NULL,
    prev_out = prev_out,
    input_dir = input_dir
  )
  res$output$scdata
}

test_that("load_obj2s_file stores counts in a single BPCells matrix dir, as count matrix uploads", {
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)

  for (obj2s_type in obj2s_types) {
    mock_scdata(data_dir, obj2s_type)
    scdata <- load_mock_obj2s(input_dir, obj2s_type)

    expect_s4_class(scdata[["RNA"]]$counts, "IterableMatrix")
    expect_s4_class(scdata[["RNA"]]$data, "IterableMatrix")

    # data is normalized lazily from counts: one dir for the processed object
    expect_equal(
      get_matrix_dirs(scdata),
      normalizePath(file.path(tempdir(), "matrix_dir"))
    )
  }
})

test_that("load_obj2s_file fails with disk-backed (BPCells) counts in the uploaded object", {
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)

  scdata <- mock_scdata(data_dir, "seurat_object")
  scdata[["RNA"]] <- as(scdata[["RNA"]], "Assay5")
  scdata[["RNA"]]$counts <- BPCells::write_matrix_dir(
    scdata[["RNA"]]$counts,
    dir = withr::local_tempfile()
  )
  saveRDS(scdata, file.path(data_dir, "r.rds"))

  expect_error(
    load_mock_obj2s(input_dir, "seurat_object"),
    errors$ERROR_OBJ2S_COUNTS
  )
})


# Spatial objects shaped as the Seurat loaders create them: Xenium a FOV with
# centroids, segmentation and molecules (LoadXenium), Visium HD a VisiumV2 with
# cell segmentations (Load10X_Spatial with segmentations) and Visium a VisiumV2
# with spot centroids. One image per sample.
mock_obj2s_square <- function(cell, cx, cy, r = 2) {
  data.frame(
    x = c(cx - r, cx + r, cx + r, cx - r, cx - r),
    y = c(cy - r, cy - r, cy + r, cy + r, cy - r),
    cell = cell
  )
}

# segmentation: "cell" or "nucleus" polygons (named as LoadXenium / Seurat's
# Visium HD loaders may name them), or "none" for centroids only (Xenium loaded
# with LoadXenium's default, or Visium HD bins)
mock_obj2s_image <- function(cells, technology, assay, key, segmentation = "cell") {
  side <- ceiling(sqrt(length(cells)))
  grid <- expand.grid(x = seq_len(side), y = seq_len(side))[seq_along(cells), ]
  centroids <- data.frame(x = 20 + grid$x * 10, y = 20 + grid$y * 10, cell = cells)
  segmentations <- do.call(rbind, lapply(seq_along(cells), function(i) {
    r <- if (segmentation == "nucleus") 1 else 2
    mock_obj2s_square(cells[i], centroids$x[i], centroids$y[i], r = r)
  }))
  segmentation_boundary <- if (segmentation == "nucleus") "nucleus" else "segmentation"

  if (technology == "xenium") {
    molecules <- data.frame(
      x = centroids$x,
      y = centroids$y,
      gene = rep(c("MS4A1", "CD79B"), length.out = length(cells))
    )
    coords <- list(centroids = SeuratObject::CreateCentroids(centroids))
    if (segmentation != "none") {
      coords[[segmentation_boundary]] <- SeuratObject::CreateSegmentation(segmentations)
    }
    return(SeuratObject::CreateFOV(
      coords = coords,
      molecules = molecules,
      assay = assay,
      key = key
    ))
  }

  if (technology == "visium_v1") {
    # older Seurat (Load10X_Spatial before v5.1): spot coordinates, no boundaries
    scale_factors <- Seurat::scalefactors(spot = 1, fiducial = 1, hires = 0.1, lowres = 0.01)
    image <- methods::new(
      "VisiumV1",
      image = array(0.5, dim = c(30, 40, 3)),
      scale.factors = scale_factors,
      coordinates = data.frame(
        tissue = 1L,
        row = as.integer(grid$y),
        col = as.integer(grid$x),
        imagerow = centroids$y,
        imagecol = centroids$x,
        row.names = cells
      ),
      spot.radius = 0.01,
      assay = assay,
      key = key
    )
    return(image)
  }

  fov <- SeuratObject::CreateFOV(
    coords = list(centroids = SeuratObject::CreateCentroids(centroids)),
    type = "centroids",
    assay = assay,
    key = key
  )

  scale_factors <- list(spot = 1, fiducial = 1, hires = 0.1, lowres = 0.01)
  class(scale_factors) <- "scalefactors"

  image <- methods::new(
    "VisiumV2",
    molecules = fov@molecules,
    boundaries = fov@boundaries,
    assay = fov@assay,
    key = fov@key,
    image = array(0.5, dim = c(30, 40, 3)),
    scale.factors = scale_factors
  )

  if (technology == "visium_hd" && segmentation != "none") {
    image[[segmentation_boundary]] <- SeuratObject::CreateSegmentation(segmentations)
  }

  image
}

mock_spatial_obj2s <- function(technology, assay = "Spatial", nimages = 2, segmentation = "cell") {
  data("pbmc_small", package = "SeuratObject")
  scdata <- pbmc_small
  scdata$seurat_clusters <- scdata$RNA_snn_res.1

  if (assay != "RNA") {
    # key collisions of the renamed assay are irrelevant here
    scdata <- suppressWarnings(SeuratObject::RenameAssays(scdata, RNA = assay))
  }

  image_cells <- split(colnames(scdata), rep(seq_len(nimages), length.out = ncol(scdata)))
  for (i in seq_len(nimages)) {
    scdata[[paste0("slice", i)]] <- mock_obj2s_image(
      image_cells[[i]], technology, assay, paste0("slice", i, "_"), segmentation
    )
  }

  scdata
}

save_mock_spatial_obj2s <- function(scdata) {
  input_dir <- withr::local_tempdir(.local_envir = parent.frame())
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)
  saveRDS(scdata, file.path(data_dir, "r.rds"))
  input_dir
}

test_that("load_obj2s_file detects the spatial technology of seurat_spatial_object", {
  for (technology in c("visium", "visium_hd", "xenium")) {
    input_dir <- save_mock_spatial_obj2s(mock_spatial_obj2s(technology))
    scdata <- load_mock_obj2s(input_dir, "seurat_spatial_object")

    expect_equal(scdata@misc$technology, technology)
    expect_equal(SeuratObject::Assays(scdata), "RNA")
    expect_s4_class(scdata[["RNA"]]$counts, "IterableMatrix")

    # one image per sample, named by sample and tied to the RNA assay
    expect_setequal(Seurat::Images(scdata), c("slice1", "slice2"))
    expect_setequal(unique(scdata$samples), c("slice1", "slice2"))
    for (image_name in Seurat::Images(scdata)) {
      image <- scdata[[image_name]]
      expect_equal(SeuratObject::DefaultAssay(image), "RNA")
      expect_setequal(
        SeuratObject::Cells(image),
        colnames(scdata)[scdata$samples == image_name]
      )
    }
  }
})

test_that("load_obj2s_file sets up Xenium objects as gem2s does", {
  input_dir <- save_mock_spatial_obj2s(mock_spatial_obj2s("xenium", assay = "Xenium"))
  scdata <- load_mock_obj2s(input_dir, "seurat_spatial_object")

  fov <- scdata[["slice1"]]
  expect_s4_class(fov, "FOV")
  expect_setequal(SeuratObject::Boundaries(fov), c("centroids", "segmentations"))
  expect_equal(SeuratObject::DefaultBoundary(fov), "centroids")

  # polygons accessor used to upload segmentations works
  coords <- Seurat::GetTissueCoordinates(fov, which = "segmentations")
  expect_setequal(unique(coords$cell), SeuratObject::Cells(fov))

  # molecules are kept apart for the molecule artifact
  expect_null(names(fov@molecules))
  expect_setequal(names(scdata@misc$transcripts), c("slice1", "slice2"))
  expect_setequal(colnames(scdata@misc$transcripts$slice1), c("x", "y", "gene"))
})

test_that("load_obj2s_file sets up Visium HD objects as gem2s does", {
  input_dir <- save_mock_spatial_obj2s(mock_spatial_obj2s("visium_hd"))
  scdata <- load_mock_obj2s(input_dir, "seurat_spatial_object")

  image <- scdata[["slice1"]]
  expect_s4_class(image, "VisiumV2")
  expect_setequal(SeuratObject::Boundaries(image), c("centroids", "segmentations"))

  # compact segmentations, as read by crop_to_capture_area and gem2s-7
  expect_true(nrow(image[["segmentations"]]@sf.data) > 0)

  # cropped to the cells' bounding box
  expect_lt(prod(dim(image@image)[1:2]), 30 * 40)

  # no micron scale on objects from the Seurat loaders
  expect_null(scdata@misc$microns_per_pixel)
})

test_that("load_obj2s_file loads Visium objects with VisiumV1 images (older Seurat)", {
  input_dir <- save_mock_spatial_obj2s(mock_spatial_obj2s("visium_v1"))
  scdata <- load_mock_obj2s(input_dir, "seurat_spatial_object")

  expect_equal(scdata@misc$technology, "visium")
  expect_setequal(Seurat::Images(scdata), c("slice1", "slice2"))
  expect_s4_class(scdata[["slice1"]], "VisiumV1")
  expect_setequal(unique(scdata$samples), c("slice1", "slice2"))
})

test_that("load_obj2s_file loads Visium HD bins (no segmentations) as Visium spots", {
  # Load10X_Spatial(bin.size = 8): VisiumV2 images with centroids only
  scdata <- mock_spatial_obj2s("visium_hd", assay = "Spatial.008um", segmentation = "none")
  input_dir <- save_mock_spatial_obj2s(scdata)
  res <- load_mock_obj2s(input_dir, "seurat_spatial_object")

  expect_equal(res@misc$technology, "visium")
  expect_equal(SeuratObject::Boundaries(res[["slice1"]]), "centroids")
  expect_equal(ncol(res), ncol(scdata))
})

test_that("load_obj2s_file uses the default assay for Visium HD with bins and segmentations", {
  # bins and cell segmentations loaded as separate assays, each with its image
  scdata <- mock_spatial_obj2s("visium_hd", assay = "Spatial.008um", nimages = 1, segmentation = "none")
  cells <- colnames(scdata)
  suppressWarnings(scdata[["Spatial.Polygons"]] <- scdata[["Spatial.008um"]])
  scdata[["slice1.polygons"]] <- mock_obj2s_image(cells, "visium_hd", "Spatial.Polygons", "polygons_")

  SeuratObject::DefaultAssay(scdata) <- "Spatial.Polygons"
  res <- load_mock_obj2s(save_mock_spatial_obj2s(scdata), "seurat_spatial_object")
  expect_equal(res@misc$technology, "visium_hd")
  expect_equal(Seurat::Images(res), "slice1.polygons")

  SeuratObject::DefaultAssay(scdata) <- "Spatial.008um"
  res <- load_mock_obj2s(save_mock_spatial_obj2s(scdata), "seurat_spatial_object")
  expect_equal(res@misc$technology, "visium")
  expect_equal(Seurat::Images(res), "slice1")
})

test_that("load_obj2s_file loads Visium HD and Xenium nucleus segmentations", {
  for (technology in c("visium_hd", "xenium")) {
    scdata <- mock_spatial_obj2s(technology, segmentation = "nucleus")
    res <- load_mock_obj2s(save_mock_spatial_obj2s(scdata), "seurat_spatial_object")

    expect_equal(res@misc$technology, technology)
    # found by class whatever the boundary is named, stored as gem2s names it
    expect_setequal(SeuratObject::Boundaries(res[["slice1"]]), c("centroids", "segmentations"))
  }
})

test_that("load_obj2s_file loads Xenium objects with centroids only", {
  # LoadXenium's default (segmentations = NULL)
  scdata <- mock_spatial_obj2s("xenium", assay = "Xenium", segmentation = "none")
  res <- load_mock_obj2s(save_mock_spatial_obj2s(scdata), "seurat_spatial_object")

  expect_equal(res@misc$technology, "xenium")
  fov <- res[["slice1"]]
  expect_equal(SeuratObject::Boundaries(fov), "centroids")
  expect_setequal(SeuratObject::Cells(fov), colnames(res)[res$samples == "slice1"])

  # molecules are still kept for the molecule artifact
  expect_setequal(names(res@misc$transcripts), c("slice1", "slice2"))
})

test_that("load_obj2s_file loads Xenium objects with segmentations only", {
  scdata <- mock_spatial_obj2s("xenium", assay = "Xenium", nimages = 1)
  fov <- scdata[["slice1"]]
  fov@boundaries$centroids <- NULL
  SeuratObject::DefaultBoundary(fov) <- "segmentation"
  scdata[["slice1"]] <- NULL
  scdata[["slice1"]] <- fov

  res <- load_mock_obj2s(save_mock_spatial_obj2s(scdata), "seurat_spatial_object")

  # centroids are derived from the polygons
  expect_setequal(SeuratObject::Boundaries(res[["slice1"]]), c("centroids", "segmentations"))
  centroids <- Seurat::GetTissueCoordinates(res[["slice1"]], which = "centroids")
  expect_setequal(centroids$cell, colnames(res))
})

test_that("get_obj2s_microns_per_pixel keeps the micron scale persisted by Cellenics", {
  scdata <- mock_spatial_obj2s("visium_hd")
  scdata@misc$microns_per_pixel <- 0.27

  expect_equal(get_obj2s_microns_per_pixel(scdata, scdata[["slice1"]], "Spatial", "visium_hd"), 0.27)
})

test_that("get_obj2s_microns_per_pixel derives the micron scale of Visium HD bins from the assay", {
  # spot is the bin diameter in full-resolution pixels (here 1px = 8um)
  scdata <- mock_spatial_obj2s("visium_hd", assay = "Spatial.008um", segmentation = "none")
  image <- scdata[["slice1"]]

  expect_equal(get_obj2s_microns_per_pixel(scdata, image, "Spatial.008um", "visium"), 8)
})

test_that("get_obj2s_microns_per_pixel derives the micron scale of Visium spots", {
  # mock spots are 1px across and 10px apart, as spots (not touching bins)
  scdata <- mock_spatial_obj2s("visium")
  image <- scdata[["slice1"]]

  expect_equal(
    get_obj2s_microns_per_pixel(scdata, image, "Spatial", "visium"),
    VISIUM_SPOT_DIAMETER_MICRONS
  )

  # VisiumV1 images too
  scdata <- mock_spatial_obj2s("visium_v1")
  expect_equal(
    get_obj2s_microns_per_pixel(scdata, scdata[["slice1"]], "Spatial", "visium"),
    VISIUM_SPOT_DIAMETER_MICRONS
  )
})

test_that("get_obj2s_microns_per_pixel returns NULL when the micron scale is unknown", {
  scdata <- mock_spatial_obj2s("visium")
  image <- scdata[["slice1"]]

  # touching bins without a bin size in the assay name
  image@scale.factors$spot <- 10
  expect_null(get_obj2s_microns_per_pixel(scdata, image, "Spatial", "visium"))

  # spot size isn't a cell size for cell-segmented Visium HD
  scdata <- mock_spatial_obj2s("visium_hd")
  expect_null(get_obj2s_microns_per_pixel(scdata, scdata[["slice1"]], "Spatial", "visium_hd"))
})

test_that("load_obj2s_file persists the micron scale of Visium objects for Spaco", {
  scdata <- mock_spatial_obj2s("visium_hd", assay = "Spatial.008um", segmentation = "none")
  res <- load_mock_obj2s(save_mock_spatial_obj2s(scdata), "seurat_spatial_object")

  expect_equal(res@misc$technology, "visium")
  expect_equal(res@misc$microns_per_pixel, 8)
})

test_that("load_obj2s_file uses the assay of the images for counts", {
  scdata <- mock_spatial_obj2s("visium", assay = "Spatial.016um")
  suppressWarnings(scdata[["SCT"]] <- scdata[["Spatial.016um"]])
  SeuratObject::DefaultAssay(scdata) <- "SCT"

  input_dir <- save_mock_spatial_obj2s(scdata)
  res <- load_mock_obj2s(input_dir, "seurat_spatial_object")

  expect_equal(
    as(res[["RNA"]]$counts, "dgCMatrix"),
    scdata[["Spatial.016um"]]$counts
  )
})

test_that("load_obj2s_file fails for images of different spatial technologies", {
  scdata <- mock_spatial_obj2s("visium")
  scdata[["slice2"]] <- NULL
  cells <- colnames(scdata)[seq(2, ncol(scdata), by = 2)]
  scdata[["slice2"]] <- mock_obj2s_image(cells, "xenium", "Spatial", "xen_")

  input_dir <- save_mock_spatial_obj2s(scdata)
  expect_error(
    load_mock_obj2s(input_dir, "seurat_spatial_object"),
    errors$ERROR_OBJ2S_IMAGES
  )
})

test_that("load_obj2s_file fails for seurat_spatial_object without images", {
  data("pbmc_small", package = "SeuratObject")
  input_dir <- save_mock_spatial_obj2s(pbmc_small)

  expect_error(
    load_mock_obj2s(input_dir, "seurat_spatial_object"),
    errors$ERROR_OBJ2S_IMAGES
  )
})

# what a Cellenics download of the processed object looks like: counts
# materialized in memory (materialize_bpcells_matrix in the worker), sample
# names in 'samples' (add_clusters in the worker)
as_cellenics_download <- function(scdata, input) {
  for (layer in SeuratObject::Layers(scdata[["RNA"]])) {
    SeuratObject::LayerData(scdata[["RNA"]], layer) <-
      as(SeuratObject::LayerData(scdata[["RNA"]], layer), "dgCMatrix")
  }

  sample_names <- setNames(input$sampleNames, input$sampleIds)
  scdata$samples <- unname(sample_names[scdata$samples])
  scdata
}

test_that("seurat_spatial_object downloaded from Cellenics can be re-uploaded", {
  for (technology in c("visium", "visium_hd", "xenium")) {
    input_dir <- save_mock_spatial_obj2s(mock_spatial_obj2s(technology))
    scdata <- load_mock_obj2s(input_dir, "seurat_spatial_object")

    # as uploadObj2sToAWS sets up the processed object
    input <- add_samples_to_input(format_obj2s(scdata, "1234"), list())
    processed <- format_obj2s(scdata, "1234")
    processed <- change_sample_names_to_ids(processed, input)
    processed <- rename_obj2s_images(processed)
    expect_setequal(Seurat::Images(processed), paste0(input$sampleIds, ".polygons"))
    processed@misc$transcripts <- NULL
    if (technology == "visium_hd") processed@misc$microns_per_pixel <- 0.27

    download <- as_cellenics_download(processed, input)
    reupload_dir <- save_mock_spatial_obj2s(download)
    reuploaded <- load_mock_obj2s(reupload_dir, "seurat_spatial_object")

    expect_equal(reuploaded@misc$technology, technology)

    # samples keep their names, images are named by them again
    expect_setequal(unique(reuploaded$samples), c("slice1", "slice2"))
    expect_setequal(Seurat::Images(reuploaded), c("slice1", "slice2"))

    # same cells, counts and image setup
    expect_setequal(colnames(reuploaded), colnames(scdata))
    expect_equal(
      as(reuploaded[["RNA"]]$counts, "dgCMatrix")[, colnames(scdata)],
      as(scdata[["RNA"]]$counts, "dgCMatrix")
    )
    expect_equal(
      SeuratObject::Boundaries(reuploaded[["slice1"]]),
      SeuratObject::Boundaries(scdata[["slice1"]])
    )
    expect_equal(
      Seurat::GetTissueCoordinates(reuploaded[["slice1"]]),
      Seurat::GetTissueCoordinates(scdata[["slice1"]])
    )

    if (technology == "visium_hd") {
      expect_equal(reuploaded@misc$microns_per_pixel, 0.27)
      expect_equal(dim(reuploaded[["slice1"]]@image), dim(scdata[["slice1"]]@image))
    }
  }
})

test_that("seurat_object downloaded from Cellenics can be re-uploaded", {
  input_dir <- withr::local_tempdir()
  data_dir <- file.path(input_dir, "pbmc_small")
  dir.create(data_dir)

  mock_scdata(data_dir, "seurat_object")
  scdata <- load_mock_obj2s(input_dir, "seurat_object")

  input <- add_samples_to_input(format_obj2s(scdata, "1234"), list())
  processed <- format_obj2s(scdata, "1234")
  processed <- change_sample_names_to_ids(processed, input)

  saveRDS(as_cellenics_download(processed, input), file.path(data_dir, "r.rds"))
  reuploaded <- load_mock_obj2s(input_dir, "seurat_object")

  expect_equal(
    as(reuploaded[["RNA"]]$counts, "dgCMatrix"),
    as(scdata[["RNA"]]$counts, "dgCMatrix")
  )
  expect_equal(names(reuploaded@reductions), names(scdata@reductions))
  expect_equal(
    SeuratObject::Embeddings(reuploaded, "tsne"),
    SeuratObject::Embeddings(scdata, "tsne")
  )
})

test_that("get_obj2s_gene_annotations reuses gene symbols from a Cellenics download", {
  data("pbmc_small", package = "SeuratObject")
  rns <- c("ENSG1", "ENSG2", "ENSG3", "ENSG4")
  # as gem2s stores them: duplicated names disambiguated, original kept
  pbmc_small@misc$gene_annotations <- data.frame(
    input = rns,
    name = c("CD3E", "MS4A1 - ENSG2", "MS4A1 - ENSG3", "ENSG4"),
    original_name = c("CD3E", "MS4A1", "MS4A1", ""),
    row.names = rns
  )

  annot <- get_obj2s_gene_annotations(pbmc_small, rev(rns))

  expect_equal(annot$input, rev(rns))
  expect_equal(annot$name, c("ENSG4", "MS4A1 - ENSG3", "MS4A1 - ENSG2", "CD3E"))
  expect_equal(annot$original_name, c("", "MS4A1", "MS4A1", "CD3E"))
  expect_equal(rownames(annot), rev(rns))
})

test_that("get_obj2s_gene_annotations uses row names without matching annotations", {
  data("pbmc_small", package = "SeuratObject")
  rns <- rownames(pbmc_small)

  # no annotations
  annot <- get_obj2s_gene_annotations(pbmc_small, rns)
  expect_equal(annot$name, rns)

  # annotations missing some genes
  pbmc_small@misc$gene_annotations <- data.frame(input = rns[-1], name = tolower(rns[-1]))
  annot <- get_obj2s_gene_annotations(pbmc_small, rns)
  expect_equal(annot$name, rns)
  expect_equal(annot$original_name, rns)
})
