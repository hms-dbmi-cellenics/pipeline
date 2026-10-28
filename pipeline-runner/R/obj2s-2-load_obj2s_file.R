#' Loads a Seurat/SingleCellExperiment/AnnData object provided by the user
#'
#' This functions is a thin wrapper around `reconstruct_*` functions, which are used to
#' safely read, check, and reconstruct a Seurat object for downstream use.
#'
#'
#' @param input The input params received from the api.
#' @param pipeline_config List with S3 bucket name and SNS topics.
#' @param prev_out 'output' slot from call to \code{download_user_files}
#' @param input_dir A string, where the folder containing the downloaded
#' Seurat object is stored
#'
#' @export
load_obj2s_file <- function(input, pipeline_config, prev_out, input_dir = "/input") {

  config <- prev_out$config
  dataset_dir <- config$samples[1]
  obj2s_type <- config$input$type

  dataset_fname <- switch(
    obj2s_type,
    "seurat_spatial_object" = "r.rds",
    "seurat_object" = "r.rds",
    "sce_object" = "r.rds",
    "anndata_object" = "adata.h5ad"
  )

  dataset_fpath <- file.path(input_dir, dataset_dir, dataset_fname)

  reconstruct_fun <- switch(
    obj2s_type,
    "seurat_object" = reconstruct_seurat,
    "seurat_spatial_object" = reconstruct_seurat_spatial,
    "sce_object" = reconstruct_sce,
    "anndata_object" = reconstruct_anndata,
  )

  scdata <- reconstruct_fun(dataset_fpath)

  prev_out$scdata <- scdata

  res <- list(
    data = list(),
    output = prev_out
  )
  return(res)
}

#' Reconstructs a spatial Seurat object provided by the user
#'
#' Supports Visium (VisiumV1/VisiumV2 spots), Visium HD (VisiumV2 with cell
#' segmentations) and Xenium (FOV without a tissue image) objects, as created by
#' the Seurat loaders or downloaded from Cellenics. The technology is detected
#' from the images and persisted to \code{@misc$technology}, as gem2s does.
#'
#' Each image is a sample. Images are named by their sample here and renamed to
#' \code{<sample id>.polygons} (as gem2s names them) once sample ids are known.
#'
#' @param dataset_fpath The path to the r.rds file containing the Seurat object.
#'
#' @return reconstructed SeuratObject
reconstruct_seurat_spatial <- function(dataset_fpath) {

  # read it
  tryCatch({
    user_scdata <- readRDS(dataset_fpath)
    stopifnot(methods::is(user_scdata, "Seurat"))
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_READ, call. = FALSE)
  })

  # get images, the assay they measure, and the technology they come from
  tryCatch({
    images <- get_obj2s_images(user_scdata)
    assay <- get_spatial_assay(user_scdata, images)
    images <- images[vapply(images, SeuratObject::DefaultAssay, "") == assay]
    technology <- detect_spatial_technology(images)
    message("Detected spatial technology: ", technology, " (assay ", assay, ")")
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_IMAGES, call. = FALSE)
  })

  # get counts
  tryCatch({

    # if V5 object, ensure layers are rejoined
    if (methods::is(user_scdata[[assay]], "Assay5"))
      user_scdata[[assay]] <- SeuratObject::JoinLayers(user_scdata[[assay]])

    counts <- SeuratObject::LayerData(user_scdata, assay = assay, layer = "counts")
    test_user_sparse_mat(counts)
    rns <- row.names(counts)
    check_type_is_safe(rns)

  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_COUNTS, call. = FALSE)
  })

  # get meta data
  tryCatch({
    metadata <- user_scdata@meta.data
    metadata$active.ident <- user_scdata@active.ident
    test_user_df(metadata)
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_METADATA, call. = FALSE)
  })

  # each image is a sample: keep only cells that are in an image
  image_samples <- get_image_samples(images, metadata)
  image_cells <- lapply(images, function(image) {
    intersect(SeuratObject::Cells(image), colnames(counts))
  })

  keep_cells <- colnames(counts)[colnames(counts) %in% unlist(image_cells)]
  if (length(keep_cells) < ncol(counts)) {
    message("Removing ", ncol(counts) - length(keep_cells), " cells not in any image")
  }

  counts <- counts[, keep_cells]
  metadata <- metadata[keep_cells, ]

  metadata$samples <- NA
  for (image_name in names(images)) {
    metadata[image_cells[[image_name]], "samples"] <- image_samples[[image_name]]
  }

  # reconstruct Seurat object
  scdata <- create_obj2s_scdata(counts, metadata)

  # add gene annotations
  scdata@misc$gene_annotations <- get_obj2s_gene_annotations(user_scdata, rns)

  # add dispersions (need for gene list)
  scdata <- add_obj2s_dispersions(scdata)

  # add default dimensionality reduction
  red_names <- tolower(names(user_scdata@reductions))
  names(user_scdata@reductions) <- red_names

  tryCatch({
    red_name <- SeuratObject::DefaultDimReduc(user_scdata)
    check_type_is_safe(red_name)
    red_match <- grep("umap|tsne", red_name, value = TRUE)

    if (length(red_match) && !(red_match %in% c("umap", "tsne"))) {
      is_umap <- grepl("umap", red_match)
      new_red_name <- ifelse(is_umap, "umap", "tsne")

      message("Found reduction name ", red_match," containing ", new_red_name)
      user_scdata <- update_reduction_name(user_scdata, red_name, new_red_name)
      red_name <- SeuratObject::DefaultDimReduc(user_scdata)
      message("Updated default reduction: ", red_name)
    }

    stopifnot(red_name %in% c("umap", "tsne"))

    red <- user_scdata@reductions[[red_name]]
    test_user_df(red@cell.embeddings)
    test_user_df(red@feature.loadings)
    test_user_df(red@feature.loadings.projected)

    scdata@reductions[[red_name]] <- SeuratObject::CreateDimReducObject(
      embeddings = red@cell.embeddings[keep_cells, , drop = FALSE],
      assay = "RNA",
      key = paste0(red_name, "_")
    )
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  # add pca dimensionality reduction
  tryCatch({
    if ("pca" %in% red_names) {

      pca <- user_scdata@reductions[["pca"]]@cell.embeddings
      test_user_df(pca)
      red <- SeuratObject::CreateDimReducObject(
        embeddings = pca[keep_cells, , drop = FALSE],
        assay = "RNA",
        key = "pca_"
      )
      scdata@reductions[["pca"]] <- red
    }

  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  # add images, one per sample
  tryCatch({
    build_image <- switch(
      technology,
      visium = build_obj2s_visium_image,
      visium_hd = build_obj2s_visium_hd_image,
      xenium = build_obj2s_xenium_fov
    )

    transcripts <- list()
    for (i in seq_along(images)) {
      image_name <- names(images)[i]
      sample <- image_samples[[image_name]]
      image <- images[[image_name]]

      # Xenium molecules are uploaded as their own artifact, not kept on the FOV
      if (technology == "xenium") {
        transcripts[[sample]] <- get_obj2s_transcripts(image)
      }

      image <- build_image(image)
      image <- subset(image, cells = image_cells[[image_name]])
      SeuratObject::Key(image) <- paste0("image", i, "_")
      scdata[[sample]] <- image
    }

    if (length(transcripts)) scdata@misc$transcripts <- transcripts
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_IMAGES, call. = FALSE)
  })

  # persist the technology and micron scale as gem2s does, so the worker can
  # dispatch on them (the micron scale sizes the Spaco colour neighbourhood)
  scdata@misc$technology <- technology
  if (technology %in% c("visium", "visium_hd")) {
    scdata@misc$microns_per_pixel <- get_obj2s_microns_per_pixel(
      user_scdata, images[[1]], assay, technology
    )
  }

  return(scdata)
}

# images of supported classes: VisiumV2 and Xenium FOVs are both FOVs
get_obj2s_images <- function(user_scdata) {
  image_names <- Seurat::Images(user_scdata)
  if (!length(image_names)) stop("Object has no images.")

  images <- user_scdata@images[image_names]
  is_supported <- vapply(images, function(image) {
    methods::is(image, "VisiumV1") || methods::is(image, "FOV")
  }, logical(1))

  if (!all(is_supported)) {
    stop(
      "Unsupported image class: ",
      paste(unique(sapply(images[!is_supported], class)), collapse = ", ")
    )
  }

  return(images)
}

# the assay the images measure: the default assay if it has images, otherwise
# the first assay with images (e.g. 'Spatial.016um' for binned Visium HD)
get_spatial_assay <- function(user_scdata, images) {
  image_assays <- unique(vapply(images, SeuratObject::DefaultAssay, ""))
  image_assays <- intersect(image_assays, SeuratObject::Assays(user_scdata))
  if (!length(image_assays)) stop("No assay found for the images.")

  default_assay <- SeuratObject::DefaultAssay(user_scdata)
  if (default_assay %in% image_assays) return(default_assay)

  if (length(image_assays) > 1) {
    message(
      "Images measure several assays (", paste(image_assays, collapse = ", "),
      "), using ", image_assays[1]
    )
  }

  return(image_assays[1])
}

# Xenium is a FOV without a tissue image, Visium HD a VisiumV2 with cell
# segmentations. Spots (VisiumV1, binned Visium HD) are Visium.
get_image_technology <- function(image) {
  if (methods::is(image, "VisiumV1")) return("visium")
  if (!methods::is(image, "VisiumV2")) return("xenium")
  if (!is.null(get_segmentation_boundary(image))) return("visium_hd")
  return("visium")
}

detect_spatial_technology <- function(images) {
  technologies <- unique(vapply(images, get_image_technology, ""))
  if (length(technologies) > 1) {
    stop(
      "Images are from different spatial technologies: ",
      paste(technologies, collapse = ", ")
    )
  }

  return(technologies)
}

# name of the cell polygons boundary: 'segmentations' in Cellenics objects,
# 'segmentation' in objects from the Seurat loaders
get_segmentation_boundary <- function(image) {
  boundaries <- SeuratObject::Boundaries(image)
  is_segmentation <- vapply(boundaries, function(boundary) {
    methods::is(image[[boundary]], "Segmentation")
  }, logical(1))

  segmentation_boundaries <- boundaries[is_segmentation]
  if (!length(segmentation_boundaries)) return(NULL)
  if ("segmentations" %in% segmentation_boundaries) return("segmentations")
  return(segmentation_boundaries[1])
}

get_centroids_boundary <- function(image) {
  boundaries <- SeuratObject::Boundaries(image)
  is_centroids <- vapply(boundaries, function(boundary) {
    methods::is(image[[boundary]], "Centroids")
  }, logical(1))

  centroids_boundaries <- boundaries[is_centroids]
  if (!length(centroids_boundaries)) return(NULL)
  if ("centroids" %in% centroids_boundaries) return("centroids")
  return(centroids_boundaries[1])
}

# x/y/cell coordinates of a boundary, type checked
get_boundary_coords <- function(image, boundary) {
  coords <- SeuratObject::GetTissueCoordinates(image[[boundary]])
  coords <- coords[, c("x", "y", "cell")]
  test_user_df(coords)
  return(coords)
}

# Sample name of each image: the user's 'samples'/'sample' column (e.g. objects
# downloaded from Cellenics) when each image has a single distinct value of it,
# otherwise the image name.
get_image_samples <- function(images, metadata) {
  image_names <- names(images)
  sample_col <- intersect(c("samples", "sample"), colnames(metadata))[1]
  if (is.na(sample_col)) return(as.list(setNames(image_names, image_names)))

  samples <- lapply(images, function(image) {
    cells <- intersect(SeuratObject::Cells(image), row.names(metadata))
    unique(as.character(metadata[cells, sample_col]))
  })

  is_single <- vapply(samples, function(x) length(x) == 1 && !is.na(x), logical(1))
  samples <- unlist(samples[is_single])

  if (!all(is_single) || anyDuplicated(samples)) {
    return(as.list(setNames(image_names, image_names)))
  }

  return(as.list(samples))
}

# Visium spots: keep the image, tied to the RNA assay of the new object
build_obj2s_visium_image <- function(image) {
  check_type_is_safe(image)
  SeuratObject::DefaultAssay(image) <- "RNA"
  return(image)
}

# Visium HD cells: keep the tissue image with compact centroids and
# segmentations, then crop and rotate it as gem2s does
build_obj2s_visium_hd_image <- function(image) {
  check_type_is_safe(image)

  centroids_boundary <- get_centroids_boundary(image)
  segmentation_boundary <- get_segmentation_boundary(image)
  segmentations <- get_boundary_coords(image, segmentation_boundary)

  centroids <- if (is.null(centroids_boundary)) {
    get_polygon_centroids(segmentations)
  } else {
    get_boundary_coords(image, centroids_boundary)
  }

  image@boundaries <- list(
    centroids = SeuratObject::CreateCentroids(centroids),
    segmentations = SeuratObject::CreateSegmentation(segmentations, compact = TRUE)
  )
  image@molecules <- list()
  SeuratObject::DefaultBoundary(image) <- "centroids"
  SeuratObject::DefaultAssay(image) <- "RNA"

  image |>
    crop_to_capture_area() |>
    pivot_image_wide()
}

# Xenium cells: rebuild the FOV from its coordinates as gem2s does. Cell (or
# nucleus) segmentations are optional: LoadXenium loads centroids only by default
build_obj2s_xenium_fov <- function(image) {
  centroids_boundary <- get_centroids_boundary(image)
  segmentation_boundary <- get_segmentation_boundary(image)
  if (is.null(centroids_boundary) && is.null(segmentation_boundary)) {
    stop("Xenium FOV has no cell centroids or segmentations.")
  }

  segmentations <- if (!is.null(segmentation_boundary)) {
    get_boundary_coords(image, segmentation_boundary)
  }
  centroids <- if (is.null(centroids_boundary)) {
    get_polygon_centroids(segmentations)
  } else {
    get_boundary_coords(image, centroids_boundary)
  }

  fov <- build_xenium_fov(list(centroids = centroids, segmentations = segmentations))
  SeuratObject::DefaultBoundary(fov) <- "centroids"
  return(fov)
}

get_polygon_centroids <- function(segmentations) {
  stats::aggregate(cbind(x, y) ~ cell, data = segmentations, FUN = mean)[, c("x", "y", "cell")]
}

# molecules frame (x, y, gene) of a Xenium FOV, NULL if it has none
get_obj2s_transcripts <- function(image) {
  molecules <- names(image@molecules)
  if (!length(molecules)) return(NULL)

  transcripts <- SeuratObject::GetTissueCoordinates(image[[molecules[1]]])
  transcripts <- data.frame(
    x = transcripts$x,
    y = transcripts$y,
    gene = as.character(transcripts$molecule)
  )
  test_user_df(transcripts)
  return(transcripts)
}

# Space Ranger's spot_diameter_fullres spans a theoretical 65um spot (not the
# 55um capture area)
VISIUM_SPOT_DIAMETER_MICRONS <- 65
# Visium spot centres are 100um apart (~1.54x that diameter) while Visium HD
# bins touch (1x), so their spacing tells spots from bins
MIN_SPOT_SPACING_RATIO <- 1.25

is_valid_scale <- function(x) {
  is.numeric(x) && length(x) == 1 && isTRUE(x > 0)
}

#' Micron scale of a spatial object's full-resolution pixels
#'
#' Uses the scale persisted by Cellenics (\code{@misc$microns_per_pixel}). Objects
#' from the Seurat loaders don't keep it, so for Visium it is derived from the
#' spot (or bin) diameter in full-resolution pixels (\code{scale.factors$spot}):
#' bins are sized by the assay name (e.g. \code{Spatial.008um}), spots are 65um.
#'
#' @param user_scdata uploaded SeuratObject
#' @param image the (uncropped) image of a sample
#' @param assay the assay the images measure
#' @param technology detected spatial technology
#' @return microns per full-resolution pixel, or NULL if unknown
get_obj2s_microns_per_pixel <- function(user_scdata, image = NULL, assay = NULL,
                                        technology = NULL) {
  microns_per_pixel <- user_scdata@misc$microns_per_pixel
  if (is_valid_scale(microns_per_pixel)) return(microns_per_pixel)

  # spot size isn't a cell size for cell-segmented Visium HD: only spots/bins
  if (!identical(technology, "visium") || is.null(image)) return(NULL)

  spot_diameter <- tryCatch(image@scale.factors$spot, error = function(e) NULL)
  if (!is_valid_scale(spot_diameter)) return(NULL)

  bin_microns <- get_bin_size_microns(assay)
  if (!is.null(bin_microns)) return(bin_microns / spot_diameter)

  spacing <- get_median_spot_spacing(image)
  if (is.null(spacing) || spacing / spot_diameter < MIN_SPOT_SPACING_RATIO) return(NULL)

  return(VISIUM_SPOT_DIAMETER_MICRONS / spot_diameter)
}

# bin size of a Visium HD bin assay named by Load10X_Spatial (e.g. Spatial.008um)
get_bin_size_microns <- function(assay) {
  if (is.null(assay) || !grepl("\\.[0-9]+um$", assay)) return(NULL)
  as.numeric(sub(".*\\.([0-9]+)um$", "\\1", assay))
}

# median distance between neighbouring spot centres, in full-resolution pixels
get_median_spot_spacing <- function(image, max_spots = 2000) {
  coords <- Seurat::GetTissueCoordinates(image, scale = NULL)
  coords <- as.matrix(coords[, vapply(coords, is.numeric, logical(1)), drop = FALSE][, 1:2])
  if (nrow(coords) < 2) return(NULL)

  sampled <- unique(round(seq(1, nrow(coords), length.out = min(max_spots, nrow(coords)))))
  query <- coords[sampled, , drop = FALSE]
  nearest <- RANN::nn2(coords, query, k = 2)$nn.dists[, 2]
  nearest <- nearest[nearest > 0]
  if (!length(nearest)) return(NULL)

  stats::median(nearest)
}

# zellkonverter seems to miss setting this
get_h5ad_raw_rownames <- function(dataset_fpath) {
  h5_list <- rhdf5::h5ls(dataset_fpath)
  raw_name <- h5_list |>
    dplyr::filter(group == "/raw/var") |>
    dplyr::filter(name %in% c("_index", "index")) |>
    dplyr::pull(name) |>
    dplyr::first()

  rns <- rhdf5::h5read(dataset_fpath, file.path("/raw/var", raw_name))
  rns <- as.character(rns)
}

#' Gene annotations for an uploaded Seurat object
#'
#' Seurat objects downloaded from Cellenics have Ensembl ids as row names, with
#' the gene symbols in \code{@misc$gene_annotations} (as gem2s stores them).
#' Reuse those symbols when they cover every gene, so the re-upload shows gene
#' names rather than Ensembl ids. Otherwise the row names are used as names.
#'
#' @param user_scdata uploaded SeuratObject
#' @param rns gene row names of the counts
#' @return gene annotations data.frame with input, name and original_name
get_obj2s_gene_annotations <- function(user_scdata, rns) {
  annot <- user_scdata@misc$gene_annotations

  has_symbols <- is.data.frame(annot) &&
    all(c("input", "name") %in% colnames(annot)) &&
    all(rns %in% annot$input)

  if (!has_symbols) {
    return(data.frame(
      input = rns,
      name = rns,
      original_name = rns,
      row.names = rns
    ))
  }

  annot <- annot[match(rns, annot$input), ]
  symbols <- if ("original_name" %in% colnames(annot)) annot$original_name else annot$name
  symbols <- as.character(symbols)
  symbols[is.na(symbols)] <- ""

  # dedup and fill empty names as gem2s does
  format_annot(list(data.frame(input = rns, name = symbols)))
}


reconstruct_anndata <- function(dataset_fpath) {

  # get counts
  tryCatch({
    # read using anndata via reticulate
    filename <- normalizePath(dataset_fpath, mustWork = FALSE)
    anndata <- reticulate::import("anndata")
    user_adata <- anndata$read_h5ad(filename)

    # convert to SingleCellExperiment
    user_scdata <- zellkonverter::AnnData2SCE(user_adata, raw = TRUE)
    stopifnot(methods::is(user_scdata, "SingleCellExperiment"))
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_READ, call. = FALSE)
  })

  has.raw <- "raw" %in% SingleCellExperiment::altExpNames(user_scdata)

  # get X and X_raw matrices
  tryCatch({

    # counts are in "raw/X" if present, otherwise "/X"
    if (has.raw) {
      raw <- SingleCellExperiment::altExp(user_scdata, "raw")
      counts <- SummarizedExperiment::assay(raw, "X")
      row.names(counts) <- get_h5ad_raw_rownames(dataset_fpath)

    } else {
      counts <- SummarizedExperiment::assay(user_scdata, "X")

    }

    test_user_sparse_mat(counts)
    rns <- row.names(counts)
    check_type_is_safe(rns)

  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_COUNTS, call. = FALSE)
  })

  # get meta data
  tryCatch({
    metadata <- user_scdata@colData
    metadata <- as.data.frame(metadata)
    test_user_df(metadata)
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_METADATA, call. = FALSE)
  })

  # construct Seurat object from SingleCellExperiment items
  scdata <- create_obj2s_scdata(counts, metadata)

  # add gene annotations
  gene_annotations <- data.frame(
    input = rns,
    name = rns,
    original_name = rns,
    row.names = rns
  )
  scdata@misc$gene_annotations <- gene_annotations

  scdata <- add_obj2s_dispersions(scdata)

  red_names <- tolower(SingleCellExperiment::reducedDimNames(user_scdata))

  # add dimensionality reduction
  tryCatch({
    red_match <- grep("x_umap|x_tsne", red_names)

    stopifnot(length(red_match) > 0)

    # use last reduction as default (most recent call)
    red.idx <- tail(red_match, 1)
    red_name <- ifelse(grepl("umap", red_names[red.idx]), "umap", "tsne")

    red <- SingleCellExperiment::reducedDims(user_scdata)[[red.idx]]
    class(red) <- "matrix"
    test_user_df(red)

    red <- Seurat::CreateDimReducObject(
      embeddings = red,
      assay = "RNA",
      key = paste0(red_name, "_"))

    scdata@reductions[[red_name]] <- red
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  # add pca dimensionality reduction
  tryCatch({
    pca.idx <- which(red_names == "x_pca")
    if (length(pca.idx) > 0) {
      pca <- SingleCellExperiment::reducedDims(user_scdata)[[pca.idx]]
      class(pca) <- "matrix"
      test_user_df(pca)

      pca <- SeuratObject::CreateDimReducObject(
        embeddings = pca,
        assay = "RNA",
        key = "pca_"
      )
      scdata@reductions[["pca"]] <- pca
    }
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  return(scdata)

}

reconstruct_sce <- function(dataset_fpath) {
  # read it
  tryCatch({
    user_scdata <- readRDS(dataset_fpath)
    stopifnot(methods::is(user_scdata, "SingleCellExperiment"))
    colnames(user_scdata) <- make.unique(colnames(user_scdata))

    # Objects serialized with a newer Bioconductor reference the standalone
    # `Seqinfo` package (the class was moved out of GenomeInfoDb), which isn't
    # installed here. `rowRanges` is unused downstream, but accessing it via
    # updateObject() errors with "unable to find required package 'Seqinfo'".
    # Rebuild the object without rowRanges to avoid touching that class.
    user_scdata <- SingleCellExperiment::SingleCellExperiment(
      assays = as.list(SummarizedExperiment::assays(user_scdata)),
      colData = SummarizedExperiment::colData(user_scdata),
      reducedDims = SingleCellExperiment::reducedDims(user_scdata)
    )
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_READ, call. = FALSE)
  })

  # get counts
  tryCatch({
    counts <- SingleCellExperiment::counts(user_scdata)
    test_user_sparse_mat(counts)
    rns <- row.names(counts)
    check_type_is_safe(rns)

  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_COUNTS, call. = FALSE)
  })

  # get meta data
  tryCatch({
    metadata <- user_scdata@colData
    metadata <- as.data.frame(metadata)
    test_user_df(metadata)
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_METADATA, call. = FALSE)
  })

  # construct Seurat object from SingleCellExperiment items
  scdata <- create_obj2s_scdata(counts, metadata)

  # add gene annotations
  gene_annotations <- data.frame(
    input = rns,
    name = rns,
    original_name = rns,
    row.names = rns
  )
  scdata@misc$gene_annotations <- gene_annotations

  scdata <- add_obj2s_dispersions(scdata)

  red_names <- tolower(SingleCellExperiment::reducedDimNames(user_scdata))

  # add dimensionality reduction
  tryCatch({
    red_match <- grep("umap|tsne", red_names)

    stopifnot(length(red_match) > 0)

    # use last reduction as default (most recent call)
    red.idx <- tail(red_match, 1)
    red_name <- ifelse(grepl("umap", red_names[red.idx]), "umap", "tsne")

    red <- SingleCellExperiment::reducedDims(user_scdata)[[red.idx]]
    class(red) <- "matrix"
    test_user_df(red)

    # red needs to have rownames
    row.names(red) <- colnames(scdata)

    red <- Seurat::CreateDimReducObject(
      embeddings = red,
      assay = "RNA",
      key = paste0(red_name, "_")
    )

    scdata@reductions[[red_name]] <- red
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  # add pca dimensionality reduction
  tryCatch({
    pca.idx <- which(red_names == "pca")
    if (length(pca.idx) > 0) {
      pca <- SingleCellExperiment::reducedDims(user_scdata)[[pca.idx]]
      class(pca) <- "matrix"
      test_user_df(pca)

      # pca needs to have rownames
      row.names(pca) <- colnames(scdata)

      pca <- SeuratObject::CreateDimReducObject(
        embeddings = pca,
        assay = "RNA",
        key = "pca_"
      )
      scdata@reductions[["pca"]] <- pca
    }
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  return(scdata)
}

#' Reconstructs the Seurat object provided by the user
#'
#' The slots needed from the Seurat object are extracted and then their in-built types
#' are strictly checked before creating a new Seurat object. This reconstruction and type
#' checking prevents potentially malicious executable code  (e.g. type `closure` or `environment`)
#' from being stored in the final Seurat object.
#'
#' @param dataset_fpath The path to the r.rds file containing the Seurat object.
#'
#' @return reconstructed SeuratObject
#' @export
reconstruct_seurat <- function(dataset_fpath) {

  # read it
  tryCatch({
    user_scdata <- readRDS(dataset_fpath)
    stopifnot(methods::is(user_scdata, "Seurat"))
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_READ, call. = FALSE)
  })

  # get counts
  tryCatch({
    SeuratObject::DefaultAssay(user_scdata) <- "RNA"

    # if V5 object, ensure layers are rejoined
    if (methods::is(user_scdata[["RNA"]], "Assay5"))
      user_scdata[["RNA"]] <- SeuratObject::JoinLayers(user_scdata[["RNA"]])


    counts <- user_scdata[["RNA"]]$counts
    test_user_sparse_mat(counts)
    rns <- row.names(counts)
    check_type_is_safe(rns)

  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_COUNTS, call. = FALSE)
  })

  # get meta data
  tryCatch({
    metadata <- user_scdata@meta.data
    metadata$active.ident <- user_scdata@active.ident
    test_user_df(metadata)
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_METADATA, call. = FALSE)
  })

  # reconstruct Seurat object
  scdata <- create_obj2s_scdata(counts, metadata)

  # add gene annotations
  scdata@misc$gene_annotations <- get_obj2s_gene_annotations(user_scdata, rns)

  # add dispersions (need for gene list)
  scdata <- add_obj2s_dispersions(scdata)

  # add default dimensionality reduction
  red_names <- tolower(names(user_scdata@reductions))
  names(user_scdata@reductions) <- red_names

  tryCatch({
    red_name <- SeuratObject::DefaultDimReduc(user_scdata)
    check_type_is_safe(red_name)
    red_match <- grep("umap|tsne", red_name, value = TRUE)

    if (length(red_match) && !(red_match %in% c("umap", "tsne"))) {
      is_umap <- grepl("umap", red_match)
      new_red_name <- ifelse(is_umap, "umap", "tsne")

      message("Found reduction name ", red_match, " containing ", new_red_name)
      user_scdata <- update_reduction_name(user_scdata, red_name, new_red_name)
      red_name <- SeuratObject::DefaultDimReduc(user_scdata)
      message("Updated default reduction: ", red_name)
    }

    stopifnot(red_name %in% c("umap", "tsne"))

    red <- user_scdata@reductions[[red_name]]
    test_user_df(red@cell.embeddings)
    test_user_df(red@feature.loadings)
    test_user_df(red@feature.loadings.projected)
    red@assay.used <- "RNA"

    scdata@reductions[[red_name]] <- red
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  # add pca dimensionality reduction
  tryCatch({
    if ("pca" %in% red_names) {

      pca <- user_scdata@reductions[["pca"]]@cell.embeddings
      test_user_df(pca)
      red <- SeuratObject::CreateDimReducObject(
        embeddings = pca,
        assay = "RNA"
      )
      scdata@reductions[["pca"]] <- red
    }

  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_REDUCTION, call. = FALSE)
  })

  return(scdata)
}


#' Create the Cellenics object from the user's counts and metadata
#'
#' Mirrors a processed count matrix upload (gem2s and QC): counts are stored in
#' a disk-backed BPCells matrix and logcounts are normalized lazily from them, so
#' the processed object has a single matrix dir (see
#' \code{write_merged_matrix_dir}). Logcounts in the user's object are not used.
#'
#' @param counts sparse counts matrix, type checked
#' @param metadata data.frame of cell metadata, type checked
#'
#' @return SeuratObject with an RNA assay
create_obj2s_scdata <- function(counts, metadata) {
  tryCatch({
    counts <- write_obj2s_matrix_dir(counts)
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_COUNTS, call. = FALSE)
  })

  scdata <- SeuratObject::CreateSeuratObject(
    counts,
    meta.data = metadata,
  )

  tryCatch({
    scdata <- Seurat::NormalizeData(
      scdata,
      normalization.method = "LogNormalize",
      assay = "RNA",
      verbose = FALSE
    )
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_LOGCOUNTS, call. = FALSE)
  })

  return(scdata)
}

# write counts to the matrix dir uploaded with the processed object
write_obj2s_matrix_dir <- function(counts) {
  counts <- methods::as(counts, "dgCMatrix")

  # integer counts are stored as integers, as when gem2s imports them
  is_integer <- all(counts@x >= 0 & counts@x == round(counts@x))
  counts <- BPCells::convert_matrix_type(
    counts,
    ifelse(is_integer, "uint32_t", "double")
  )

  matrix_dir <- file.path(tempdir(), "matrix_dir")
  unlink(matrix_dir, recursive = TRUE)
  message("Writing counts matrix dir to: ", matrix_dir)

  BPCells::write_matrix_dir(counts, dir = matrix_dir)
}

# dispersions (needed for gene list), computed as in QC for LogNormalize
add_obj2s_dispersions <- function(scdata) {
  tryCatch({
    scdata <- Seurat::FindVariableFeatures(scdata, assay = "RNA", verbose = FALSE)
    scdata <- add_dispersions(scdata)
  },
  error = function(e) {
    message(e$message)
    stop(errors$ERROR_OBJ2S_HVFINFO, call. = FALSE)
  })

  return(scdata)
}


test_user_df <- function(df) {
  for (col in colnames(df)) {
    check_type_is_safe(df[, col])
  }
}

test_user_sparse_mat <- function(counts) {
  stopifnot(isS4(counts))

  if (methods::is(counts, "IterableMatrix")) {
    stop("Disk-backed (BPCells) counts are not supported, convert them to a dgCMatrix.")
  }

  for (slot in methods::slotNames(counts)) {
    check_type_is_safe(methods::slot(counts, slot))
  }
}

check_slot <- function(object, slot_name) {
  slot <- methods::slot(object, slot_name)
  check_type_is_safe(slot)
}

check_type_is_safe <- function(x) {

  # typeof determines the R internal type or storage mode of any object
  # whereas methods::is can be tricked by setting class 
  # inappropriately (e.g. disguising a function as a numeric)
  safe_types <- c("character", "double", "integer", "logical", "NULL", "array")


  # recurse into lists until reach node
  if (typeof(x) == "list") {
    lapply(x, check_type_is_safe)

  } else if (isS4(x)) {
    # check slots of S4 objects
    slot_names <- methods::slotNames(x)
    lapply(slot_names, check_slot, object = x)

  } else if (!typeof(x) %in% safe_types) {
    stop(sprintf("Unsafe data type: '%s'.", typeof(x)))
  }

  invisible(NULL)
}


update_reduction_name <- function(scdata, red_name, new_name) {
  current_names <- names(scdata@reductions)
  if (new_name %in% current_names) {
    message(
      "Renaming existing reduction name ",
      names(scdata@reductions)[current_names == new_name],
      " to ", paste0(new_name, ".ori")
    )
    names(scdata@reductions)[current_names == new_name] <- paste0(new_name, ".ori")
  }
  names(scdata@reductions)[current_names == red_name] <- new_name

  return(scdata)
}
