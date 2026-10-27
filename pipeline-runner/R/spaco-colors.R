# Spatially-coherent cluster/category colors via the Spaco algorithm.
#
# IMPORTANT: the Spaco colouring here mirrors the worker package
# (rworker, R/cluster.R). If you update the Spaco logic, change both.
#
# For spatial datasets (e.g. Visium HD), this assigns colors so that
# spatially interlaced categories (clusters, cell-level metadata values) get
# perceptually distinct colors. We only use it for categoricals that vary WITHIN a slice;
# per-sample categoricals (samples, sample-level metadata) keep the color pool.
#
# This is a native R reimplementation of the Spaco algorithm
# (https://github.com/BrainStOrmics/Spaco), avoiding the Python package and its
# slow numba start-up.

#' Resolve the colors to use for a categorical track
#'
#' The single entry point for colouring any cellset track. It decides
#' internally whether to colour by spatial coherence: for spatial data it
#' returns a \code{\link{spaco_palette}} aligned to \code{categories}; for
#' non-spatial data, or on any Spaco failure, it returns the supplied
#' \code{color_pool} unchanged. Spaco is only invoked when \code{scdata} carries
#' images, so non-spatial runs never touch Python.
#'
#' @param scdata a Seurat object (or NULL)
#' @param cell_ids vector of cells_id, one per cell
#' @param labels vector of category labels, one per cell (parallel to
#'  \code{cell_ids})
#' @param categories character vector of category labels in the exact order the
#'  formatter will consume colors
#' @param color_pool character vector of fallback colors in hex
#'
#' @return character vector of colors to consume in \code{categories} order
#' @export
#'
resolve_color_pool <- function(scdata, cell_ids, labels, categories, color_pool) {
  # only spatial data (a Seurat object with images) gets Spaco colouring
  if (is.null(scdata) || !length(Seurat::Images(scdata))) {
    return(color_pool)
  }

  cells_id_labels <- stats::setNames(
    as.character(labels), as.character(cell_ids)
  )
  spaco_colors <- spaco_palette(scdata, cells_id_labels, categories)
  if (is.null(spaco_colors)) color_pool else spaco_colors
}


#' Resolve cluster colors from a cell_sets data.frame
#'
#' Clustering-specific adapter around \code{\link{resolve_color_pool}} that pulls
#' the per-cell labels and the (naturally-sorted) category order straight out of
#' a cluster \code{cell_sets} data.frame, matching how
#' \code{format_cluster_cellsets} consumes colors.
#'
#' @param scdata a Seurat object (or NULL)
#' @param cell_sets data.frame with \code{cluster} and \code{cell_ids} columns
#' @param color_pool character vector of fallback colors in hex
#'
#' @return character vector of colors, one per cluster in sorted order
#' @export
#'
cluster_color_pool <- function(scdata, cell_sets, color_pool) {
  resolve_color_pool(
    scdata,
    cell_sets$cell_ids,
    cell_sets$cluster,
    sort_cluster_names(unique(cell_sets$cluster)),
    color_pool
  )
}


#' Spaco-derived palette for a categorical track, ordered to match a formatter
#'
#' Convenience wrapper around \code{\link{get_spaco_color_map}} that returns a
#' plain color vector aligned to \code{categories} - i.e. ready to drop in as
#' the \code{color_pool} argument of a cellset formatter that consumes colors
#' sequentially in that same category order.
#'
#' @param scdata_list a Seurat object, or a list of them (one per slice)
#' @param cells_id_labels named character vector mapping cells_id (as character)
#'  to the category label for that cell
#' @param categories character vector of category labels in the exact order the
#'  formatter will consume colors
#'
#' @return character vector of hex colors aligned to \code{categories}, or
#'  \code{NULL} if colouring failed or did not produce every category.
#' @export
#'
spaco_palette <- function(scdata_list, cells_id_labels, categories) {
  color_map <- get_spaco_color_map(scdata_list, cells_id_labels)
  if (is.null(color_map)) {
    return(NULL)
  }

  palette <- unname(color_map[as.character(categories)])
  # if a color was not produced for every category, fall back entirely
  if (anyNA(palette)) {
    return(NULL)
  }

  palette
}


# Colors Spaco assigns to categories,
# most distinct first: Kelly's colors of maximum contrast (Kelly 1965; without
# white, black and grey), then colors each chosen as the one most different
# (CIEDE2000) from all before it (as in Glasbey et al. 2007), up to 100 colors.
# Candidates are an sRGB grid of 18 levels per channel excluding dark (Lab L <= 30)
# and grey (chroma <= 15) colors. Colors hard to see against what clusters are
# drawn on are excluded throughout: within CIEDE2000 25 of the UI's light grey
# background (#F0F2F5), or within 15 of H&E stain (pinks to purples, sampled
# from a mouse brain section).
SPACO_COLOR_POOL <- c(
  "#f3c300", "#f38400", "#be0032", "#008856", "#0067a5", "#f99379",
  "#f6a600", "#b3446c", "#dcd300", "#882d17", "#8db600", "#654522",
  "#e25822", "#2b3d26", "#00b4d2", "#8700ff", "#00ffc3", "#9696ff",
  "#968769", "#5a6900", "#877896", "#87b496", "#1e6969", "#00ff00",
  "#0f96f0", "#870078", "#96695a", "#e1ffa5", "#ff3c5a", "#a56900",
  "#694b5a", "#1ef0f0", "#00c35a", "#968700", "#ff1ea5", "#c3b478",
  "#5a8796", "#5a69f0", "#005a0f", "#69785a", "#3c9600", "#00a596",
  "#87965a", "#ffb478", "#c3ff00", "#c39678", "#78693c", "#871e3c",
  "#a54b00", "#7887a5", "#d20f00", "#96e196", "#78a5a5", "#4b4b1e",
  "#c3964b", "#966978", "#0fc396", "#2d694b", "#874b4b", "#c36969",
  "#ffe187", "#5a8778", "#2d7800", "#c3783c", "#4bc3ff", "#2d5a69",
  "#4bc3c3", "#c3d269", "#d20096", "#5aa569", "#873c69", "#69d22d",
  "#0000ff", "#ffff78", "#f00069", "#c3a51e", "#0096c3", "#5a87f0",
  "#004b78", "#ff001e", "#008787", "#b45a3c", "#00ff87", "#78b4ff",
  "#78e1c3", "#789600", "#693c2d", "#007896", "#5a874b", "#876900",
  "#d2785a", "#a5b487", "#ff695a", "#ff783c", "#a5a53c", "#875a2d",
  "#4b5a3c", "#a5ff78", "#96c369", "#789678"
)


# Spaco neighbourhood radius, in microns (the Spaco paper's default). Applied
# directly to Xenium (micron coords); converted to image pixels for Visium HD.
SPACO_RADIUS_MICRONS <- 50

#' Neighbourhood radius for \code{spatial_distance_r}, in the coordinate units
#' returned by \code{GetTissueCoordinates(scdata, image_name, scale = scale)}.
#'
#' The Spaco radius is defined in microns. Xenium coords are already microns
#' (\code{scale} is NULL) so the radius is used as-is. Visium HD coords are image
#' pixels at the hires/lowres scale, so convert microns -> fullres pixels (via
#' \code{microns_per_pixel}, persisted onto \code{@misc} at load) -> scaled
#' pixels (via the image's hires/lowres scale factor). Returns \code{NULL} when
#' the conversion factor is unavailable so the caller skips spatial colouring
#' rather than running the neighbourhood graph at the wrong scale.
#'
#' @param scdata a Seurat object
#' @param image_name name of the image/FOV slice
#' @param scale the scale passed to \code{GetTissueCoordinates} (NULL, "hires"
#'  or "lowres")
#'
#' @return radius in coordinate units, or \code{NULL} if it can't be derived
get_spaco_radius <- function(scdata, image_name, scale) {
  if (is.null(scale)) {
    return(SPACO_RADIUS_MICRONS)
  }

  microns_per_pixel <- scdata@misc$microns_per_pixel
  scale_factor <- switch(scale,
    hires = scdata[[image_name]]@scale.factors$hires,
    lowres = scdata[[image_name]]@scale.factors$lowres,
    NULL
  )
  if (is.null(microns_per_pixel) || is.null(scale_factor)) {
    return(NULL)
  }

  SPACO_RADIUS_MICRONS / microns_per_pixel * scale_factor
}

#' Compute spatially-coherent category colors (Spaco algorithm, in R)
#'
#' Native R reimplementation of the Spaco colouring, treating each image in the
#' Seurat object(s) as a slice. Per slice it builds a cluster interlacement
#' distance graph (\code{\link{spatial_distance_r}}), sums the graphs across
#' slices, then assigns \code{SPACO_COLOR_POOL} so the most interlaced
#' categories get the most distinct colors (\code{\link{assign_spaco_palette}}).
#'
#' When a slice is sketched, colors are computed from the sketch cells only (far
#' fewer cells, so the neighbourhood graph is much cheaper); category colors
#' apply to the whole dataset regardless.
#'
#' @param scdata_list a Seurat object, or a list of them (one per slice)
#' @param cells_id_labels named character vector mapping cells_id (as character)
#'  to the category label for that cell
#'
#' @return named character vector mapping category label to a hex color, or
#'  \code{NULL} if no slice has labelled cells or anything fails.
#' @export
#'
get_spaco_color_map <- function(scdata_list, cells_id_labels) {
  if (inherits(scdata_list, "Seurat")) {
    scdata_list <- list(scdata_list)
  }

  tstart <- Sys.time()
  tryCatch(
    {
      categories <- sort(unique(unname(cells_id_labels)))

      # per-slice cluster interlacement distance graphs
      per_slice <- list()
      for (scdata in scdata_list) {
        img_names <- Seurat::Images(scdata)
        if (length(img_names) == 0) next

        # match the worker (cluster.R): select the slice's image and apply the
        # same scale rule via get_image_scale (NULL for scaleless techs like
        # Xenium, "hires"/"lowres" for Visium HD), passing it to
        # GetTissueCoordinates so the coords frame matches the upload/worker paths.
        img <- img_names[[1]]
        scale <- get_image_scale(img, scdata)
        coords <- Seurat::GetTissueCoordinates(scdata, img, scale = scale)
        if (is.null(coords) || nrow(coords) == 0) next

        # 50um neighbourhood in this slice's coordinate units; skip the slice if
        # the micron->pixel factor is missing rather than colour at a wrong scale
        radius <- get_spaco_radius(scdata, img, scale)
        if (is.null(radius)) next

        # map slice cells (barcodes) -> cells_id -> category label
        cell_ids <- scdata@meta.data[coords$cell, "cells_id"]
        labels <- cells_id_labels[as.character(cell_ids)]
        keep <- !is.na(labels)

        # if sketched, compute from the sketch cells only (far fewer)
        if ("sketch" %in% names(scdata@assays)) {
          keep <- keep & coords$cell %in% colnames(scdata[["sketch"]])
        }
        if (sum(keep) < 4) next

        per_slice <- append(per_slice, list(spatial_distance_r(
          as.matrix(coords[keep, c("x", "y")]), as.character(labels[keep]),
          radius = radius
        )))
      }

      if (length(per_slice) == 0) {
        return(NULL)
      }

      # sum interlacement graphs across slices, aligned to the full category set
      cluster_distance <- merge_cluster_distances(per_slice, categories)

      # most interlaced categories get the most distinct colors
      color_map <- assign_spaco_palette(cluster_distance, SPACO_COLOR_POOL)

      message(
        "Spaco(R): coloring time ",
        round(difftime(Sys.time(), tstart, units = "secs"), 2),
        " seconds for ", length(categories), " categories"
      )
      color_map
    },
    error = function(e) {
      message(
        "Spaco(R) coloring failed (", conditionMessage(e),
        "), using default colors"
      )
      NULL
    }
  )
}


#' Cluster spatial interlacement distance graph (Spaco \code{spatial_distance})
#'
#' Native R port of \code{spaco.distance.spatial_distance}. For each cell, finds
#' its \code{n_neighbors} nearest neighbours within \code{radius}; cells with
#' fewer than \code{n_cells} same-cluster neighbours are "banished" (contribute
#' nothing). The interlacement score between two clusters accumulates
#' inverse-distance weights over neighbouring cells of different clusters, then
#' is symmetrised (max) with a zero diagonal.
#'
#' @param coords numeric matrix of spatial coordinates (n cells x 2)
#' @param labels character vector of cluster labels (length n)
#' @param radius,n_neighbors,n_cells Spaco neighbourhood parameters. \code{radius}
#'  is in the same units as \code{coords}; production callers always pass a value
#'  derived from \code{\link{get_spaco_radius}} (50um converted to the coords'
#'  scale). The default is only a fallback for direct/unit-test calls.
#'
#' @return symmetric cluster x cluster matrix with cluster labels as dimnames
#'
spatial_distance_r <- function(
  coords, labels, radius = 90, n_neighbors = 16L, n_cells = 3L
) {
  clusters <- sort(unique(labels))
  k_clusters <- length(clusters)
  n <- nrow(coords)

  knn <- RANN::nn2(coords, coords, k = min(n_neighbors, n))
  idx <- knn$nn.idx     # n x k, self in column 1
  dst <- knn$nn.dists   # n x k, self distance 0 in column 1

  label_idx <- match(labels, clusters)
  neighbor_label <- matrix(label_idx[idx], nrow = n)
  within <- dst <= radius

  # banish cells with too few same-cluster neighbours within radius
  same_count <- rowSums((neighbor_label == label_idx) & within)
  banished <- same_count < n_cells
  # neighbourhood size (including self) used to normalise scores
  size_n <- rowSums(within)

  # accumulate scores over non-self, within-radius neighbours of kept cells
  mask <- within
  mask[, 1] <- FALSE
  mask[banished, ] <- FALSE
  pos <- which(mask, arr.ind = TRUE)

  weights <- (1 / dst[pos]) / size_n[pos[, 1]]
  cluster_i <- label_idx[pos[, 1]]
  cluster_j <- neighbor_label[pos]
  agg <- tapply(
    weights,
    list(factor(cluster_i, seq_len(k_clusters)),
         factor(cluster_j, seq_len(k_clusters))),
    sum
  )
  agg[is.na(agg)] <- 0
  score <- matrix(as.numeric(agg), k_clusters, k_clusters)

  # keep max of the two directions, zero the diagonal
  score <- pmax(score, t(score))
  diag(score) <- 0
  dimnames(score) <- list(clusters, clusters)
  score
}


#' Sum per-slice cluster interlacement graphs, aligned to a cluster set
#'
#' @param per_slice list of cluster x cluster matrices (one per slice)
#' @param clusters character vector of all clusters to align to
#'
#' @return cluster x cluster matrix summed across slices (missing entries 0)
#'
merge_cluster_distances <- function(per_slice, clusters) {
  merged <- matrix(
    0, length(clusters), length(clusters), dimnames = list(clusters, clusters)
  )
  for (m in per_slice) {
    rn <- rownames(m)
    merged[rn, rn] <- merged[rn, rn] + m
  }
  merged
}


#' Assign palette colors to clusters by interlacement (Spaco palette mode)
#'
#' Spaco's mode for a given palette: choose which cluster gets which color so
#' that the most interlaced clusters get the most distinct colors, maximising
#' the sum over cluster pairs of interlacement x CIEDE2000 color difference.
#' Clusters are assigned greedily, most interlaced first, then pairs of clusters
#' swap colors while that improves the sum.
#'
#' Used instead of Spaco's auto-generated palette (a UMAP embedding of the graph
#' into CIE-Lab), which crowds many clusters into similar hues and changes from
#' run to run.
#'
#' @param cluster_distance symmetric cluster x cluster interlacement matrix
#' @param palette hex colors, most distinct first. The first one per cluster are
#'  used, recycled if there are more clusters than colors.
#' @param max_passes maximum passes over the clusters when swapping. Almost all
#'  of the improvement comes in the first few, and each costs ~k^3 operations
#'
#' @return named character vector mapping cluster to hex color
#'
assign_spaco_palette <- function(cluster_distance, palette, max_passes = 3) {
  clusters <- rownames(cluster_distance)
  k <- length(clusters)
  colors <- rep_len(palette, k)

  lab <- farver::convert_colour(t(grDevices::col2rgb(colors)), "rgb", "lab")
  color_diff <- farver::compare_colour(lab, lab, "lab", method = "cie2000")

  # assigned[i] is the index in colors of cluster i's color. Ties go to the
  # earlier color, so unconnected clusters take the palette in order.
  assigned <- rep(NA_integer_, k)
  for (i in order(-rowSums(cluster_distance))) {
    done <- which(!is.na(assigned))
    free <- setdiff(seq_len(k), assigned[done])
    gain <- vapply(free, function(color) {
      sum(cluster_distance[i, done] * color_diff[color, assigned[done]])
    }, numeric(1))
    assigned[i] <- free[which.max(gain)]
  }

  # swap cluster i's color with the cluster j where that improves the sum most.
  # Swapping only changes their pairs with the other clusters m:
  # gain_j = sum_m (d[i, m] - d[j, m]) * (c[j, m] - c[i, m]), m != i, j
  for (pass in seq_len(max_passes)) {
    improved <- FALSE
    for (i in seq_len(k)) {
      pair_diff <- color_diff[assigned, assigned]
      gain <- rowSums(
        sweep(-cluster_distance, 2, cluster_distance[i, ], `+`) *
          sweep(pair_diff, 2, pair_diff[i, ])
      ) + 2 * cluster_distance[i, ] * pair_diff[i, ] # remove the m = i, j terms
      gain[i] <- 0

      j <- which.max(gain)
      if (gain[j] > 1e-9) {
        assigned[c(i, j)] <- assigned[c(j, i)]
        improved <- TRUE
      }
    }
    if (!improved) break
  }

  stats::setNames(colors[assigned], clusters)
}
