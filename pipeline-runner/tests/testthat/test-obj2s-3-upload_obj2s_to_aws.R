mock_scdata <- function() {
  data("pbmc_small", package = "SeuratObject")
  rns <- row.names(pbmc_small)
  pbmc_small@misc$gene_annotations <- data.frame(
    input = rns,
    name = rns,
    original_name = rns,
    row.names = rns
  )

  return(pbmc_small)
}

mock_put_object_in_s3 <- function(pipeline_config, bucket, object, key) {
  return(NULL)
}

mock_upload_matrix_to_s3 <- function(pipeline_config, experiment_id, data) {
  return("object_key")
}

test_that("upload_obj2s_to_aws completes successfully", {
  mockery::stub(upload_obj2s_to_aws, "put_object_in_s3", mock_put_object_in_s3)
  mockery::stub(
    upload_obj2s_to_aws, "upload_matrix_to_s3", mock_upload_matrix_to_s3
  )

  scdata <- mock_scdata()
  samples <- rep(c("A", "B", "C", "D"), each = ncol(scdata) / 4)
  scdata$samples <- samples
  scdata$seurat_clusters <- rep(letters[1:8], length.out = ncol(scdata))

  input <- list(experimentId = "1234")
  prev_out <- list(
    scdata = scdata,
    config = list(input = list(type = "seurat_object"))
  )

  expect_error(upload_obj2s_to_aws(input, NULL, prev_out), NA)
})


test_that("mock_metadata can be used to add columns that are missing", {
  scdata <- mock_scdata()
  metadata_cols <- list(
    "percent.mt" = 0,
    "doublet_scores" = 0,
    "doublet_class" = "singlet"
  )

  expect_true(all(!names(metadata_cols) %in% colnames(scdata@meta.data)))

  scdata <- mock_metadata(scdata, metadata_cols)
  expect_true(all(names(metadata_cols) %in% colnames(scdata@meta.data)))
  expect_true(all(scdata$percent.mt == 0))
  expect_true(all(scdata$doublet_scores == 0))
  expect_true(all(scdata$doublet_class == "singlet"))
})

test_that("mock_metadata skips over columns that are present", {
  scdata <- mock_scdata()
  metadata_cols <- list(
    "percent.mt" = 0,
    "doublet_scores" = 0,
    "doublet_class" = "singlet"
  )

  expect_true(all(!names(metadata_cols) %in% colnames(scdata@meta.data)))
  scdata$percent.mt <- 20

  scdata <- mock_metadata(scdata, metadata_cols)
  expect_true(all(names(metadata_cols) %in% colnames(scdata@meta.data)))
  expect_true(all(scdata$percent.mt == 20))
  expect_true(all(scdata$doublet_scores == 0))
  expect_true(all(scdata$doublet_class == "singlet"))
})


test_that("add_samples_col labels all cells as one sample if no sample or samples column", {
  scdata <- mock_scdata()

  expect_null(scdata@meta.data$sample)
  expect_null(scdata@meta.data$samples)

  scdata <- add_samples_col(scdata)
  expect_true(all(scdata$samples == "NA"))
})

test_that("add_samples_col uses existing 'samples' or 'sample' metadata column", {
  scdata <- mock_scdata()
  samples1 <- rep(c("A", "B"), length.out = ncol(scdata))

  # will add 'sample' if present
  scdata$sample <- samples1
  scdata <- add_samples_col(scdata)
  expect_equal(scdata@meta.data$samples, samples1)

  # will leave 'samples' if present
  scdata <- mock_scdata()
  scdata$samples <- samples1
  scdata <- add_samples_col(scdata)
  expect_equal(scdata@meta.data$samples, samples1)

  # will prefer 'samples' or 'sample'
  scdata <- mock_scdata()
  scdata$sample <- samples1
  samples2 <- rep(c("C", "D"), length.out = ncol(scdata))
  scdata$samples <- samples2
  scdata <- add_samples_col(scdata)
  expect_equal(scdata@meta.data$samples, samples2)
})


test_that("format_obj2s adds required metadata", {
  scdata <- mock_scdata()
  scdata <- format_obj2s(scdata, "1234")

  # added samples
  expect_true(all(scdata$samples == "NA"))

  # added 0-indexed cells_id
  expect_true(all(scdata$cells_id == seq(0, ncol(scdata) - 1)))

  # added misc
  expect_equal(scdata@misc$experimentId, "1234")
  expect_setequal(
    names(scdata@misc),
    c("experimentId", "color_pool", "ingestionDate", "gene_annotations")
  )

  # added required metadata columns
  expect_true(all(
    c("percent.mt", "doublet_scores", "doublet_class") %in%
      colnames(scdata@meta.data)
  ))
})


test_that("find_group_columns finds all group columns that are superset of sample identities", {
  scdata <- mock_scdata()
  samples <- rep(c("A", "B", "C", "D"), each = ncol(scdata) / 4)
  scdata$samples <- samples

  # ignores column that is just duplicate of samples
  scdata$samples_copy <- as.numeric(as.factor(samples))
  expect_length(find_group_columns(scdata@meta.data), 0)

  # finds column that has two samples per group
  scdata$group1 <- ifelse(samples %in% c("A", "B"), "AB", "CD")
  expect_equal(find_group_columns(scdata@meta.data), "group1")

  # ignores column that has same grouping as another group column
  scdata$group1_copy <- ifelse(samples %in% c("A", "B"), "group_ab", "group_cd")
  expect_equal(find_group_columns(scdata@meta.data), "group1")

  # keeps column with same grouping as another if remove.dups is FALSE
  expect_equal(
    find_group_columns(scdata@meta.data, remove.dups = FALSE),
    c("group1", "group1_copy")
  )

  # finds column that has three samples in one group
  scdata$group2 <- ifelse(samples %in% c("A", "B", "C"), "ABC", "D")
  expect_equal(find_group_columns(scdata@meta.data), c("group1", "group2"))

  # ignores column that has all samples together (trivial)
  scdata$group3 <- "ABCD"
  expect_equal(find_group_columns(scdata@meta.data), c("group1", "group2"))

  # ignores columns that divide a sample between two+ groups
  scdata$group4 <- rep(LETTERS[1:8], each = ncol(scdata) / 8)
  expect_equal(find_group_columns(scdata@meta.data), c("group1", "group2"))
})

test_that("make_vals_numeric turns equivalent groups into identical vectors", {
  # works when differently ordered
  vals1 <- c("a", "a", "a", "a", "b", "b", "b")
  vals2 <- c("y", "y", "y", "y", "x", "x", "x")
  expect_identical(make_vals_numeric(vals1), make_vals_numeric(vals2))

  bad_make_vals_numeric <- function(vals) as.numeric(factor(vals))
  expect_failure(
    expect_identical(make_vals_numeric(vals1), bad_make_vals_numeric(vals2))
  )
})


test_that("add_metadata_to_input adds group metadata to input list", {
  scdata <- mock_scdata()
  samples <- rep(c("A", "B", "C", "D"), each = ncol(scdata) / 4)
  scdata$samples <- samples

  # add column that has two samples per group
  scdata$group1 <- ifelse(samples %in% c("A", "B"), "AB", "CD")
  # add column that has three samples in one group
  scdata$group2 <- ifelse(samples %in% c("A", "B", "C"), "ABC", "D")

  input <- add_metadata_to_input(scdata, input = list())

  expect_equal(input$metadata$group1, c("AB", "AB", "CD", "CD"))
  expect_equal(input$metadata$group2, c("ABC", "ABC", "ABC", "D"))
})


test_that("add_samples_to_input adds sampleNames and sampleIds to input", {
  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  scdata$samples <- rep(sample_names, each = ncol(scdata) / 4)

  input <- add_samples_to_input(scdata, input = list())

  expect_equal(input$sampleNames, sample_names)
  expect_length(input$sampleIds, 4)
})


test_that("change_sample_names_to_ids substitutes samples column for values in sampleIds input using order based on sampleNames input", {
  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  input <- list(sampleNames = sample_names, sampleIds = letters[1:4])
  scdata <- change_sample_names_to_ids(scdata, input)
  expect_equal(scdata@meta.data$samples, tolower(samples))

  # uses reversed order from input$sampleNames
  scdata$samples <- samples
  input$sampleNames <- rev(input$sampleNames)
  scdata <- change_sample_names_to_ids(scdata, input)
  expect_equal(scdata@meta.data$samples, rev(tolower(samples)))
})

test_that("find_cluster_columns finds cluster columns", {
  expected_cols <- c(
    "seurat_clusters",
    "RNA_snn_res.0.8",
    "letter.idents",
    "groups",
    "RNA_snn_res.1"
  )

  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  scdata$seurat_clusters <- scdata$RNA_snn_res.0.8

  cluster_cols <- find_cluster_columns(scdata)
  expect_setequal(cluster_cols, expected_cols)
})


test_that("find_cluster_columns skips the per-assay QC count columns", {
  scdata <- mock_scdata()
  scdata$samples <- rep(c("A", "B", "C", "D"), each = ncol(scdata) / 4)
  scdata$seurat_clusters <- scdata$RNA_snn_res.0.8

  # integer-valued with few values, as for a small Visium HD or Xenium object
  for (col in c("nCount_Spatial.008um", "nFeature_Spatial.008um", "nCount_Xenium", "nFeature_sketch")) {
    scdata[[col]] <- rep(1:5, length.out = ncol(scdata))
  }

  cluster_cols <- find_cluster_columns(scdata)
  expect_false(any(grepl("^n(Count|Feature)_", cluster_cols)))
  expect_true("seurat_clusters" %in% cluster_cols)
})

test_that("find_cluster_columns skips column with single value", {
  expected_cols <- c(
    "seurat_clusters",
    "RNA_snn_res.0.8",
    "letter.idents",
    "groups",
    "RNA_snn_res.1"
  )

  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  scdata$seurat_clusters <- scdata$RNA_snn_res.0.8
  scdata$blah <- "blah"

  cluster_cols <- find_cluster_columns(scdata)
  expect_setequal(cluster_cols, expected_cols)
})

test_that("find_cluster_columns skips columns with non-integer numeric values", {
  expected_cols <- c(
    "seurat_clusters",
    "RNA_snn_res.0.8",
    "letter.idents",
    "groups",
    "RNA_snn_res.1"
  )

  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  scdata$seurat_clusters <- scdata$RNA_snn_res.0.8
  scdata$blah <- rnorm(ncol(scdata))

  cluster_cols <- find_cluster_columns(scdata)
  expect_setequal(cluster_cols, expected_cols)
})

test_that("find_cluster_columns skips columns where repeated values are too infrequent", {
  expected_cols <- c(
    "seurat_clusters",
    "RNA_snn_res.0.8",
    "letter.idents",
    "groups",
    "RNA_snn_res.1"
  )

  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  scdata$seurat_clusters <- scdata$RNA_snn_res.0.8
  scdata$blah <- rep(c(letters[1:20], LETTERS[1:20]), each = 2)

  cluster_cols <- find_cluster_columns(scdata)
  expect_setequal(cluster_cols, expected_cols)
})

test_that("find_cluster_columns puts 'louvain' column first if exists", {
  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  scdata$louvain <- scdata$RNA_snn_res.0.8

  cluster_cols <- find_cluster_columns(scdata)
  expect_equal(cluster_cols[1], "louvain")
})


test_that("find_cluster_columns omits columns that start with scDblFinder", {
  expected_cols <- c(
    "RNA_snn_res.0.8",
    "letter.idents",
    "groups",
    "RNA_snn_res.1"
  )

  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  scdata$scDblFinder.cluster <- scdata$RNA_snn_res.0.8

  cluster_cols <- find_cluster_columns(scdata)
  expect_setequal(cluster_cols, expected_cols)
})

test_that("find_cluster_columns works when active.ident is a group column that is duplicate of other column", {
  expected_cols <- c(
    "RNA_snn_res.0.8",
    "letter.idents",
    "groups",
    "RNA_snn_res.1"
  )

  scdata <- mock_scdata()
  sample_names <- c("A", "B", "C", "D")
  samples <- rep(sample_names, each = ncol(scdata) / 4)
  scdata$samples <- samples

  batch <- ifelse(samples %in% c("A", "B"), "Batch 1", "Batch 2")
  scdata$active.ident <- scdata$batch <- batch


  cluster_cols <- find_cluster_columns(scdata)
  expect_setequal(cluster_cols, expected_cols)
})


test_that("rename_active_ident renames active.ident correctly when there is a duplicate column", {
  meta <- data.frame(
    active.ident = c("cluster1", "cluster2", "cluster3"),
    cell_type = c("cluster1", "cluster2", "cluster3"),
    other_clustering = c("clusterA", "clusterA", "clusterB")
  )

  expected_name <- "cell_type"

  meta_renamed <- rename_active_ident(meta)

  expect_equal(colnames(meta_renamed)[1], expected_name)
  expect_false("active.ident" %in% colnames(meta_renamed))
  expect_equal(colnames(meta_renamed), c("cell_type", "other_clustering"))
})

test_that("rename_active_ident does not rename active.ident when there isn't a duplicate column", {
  meta <- data.frame(
    active.ident = c("cluster1", "cluster2", "cluster3"),
    cell_type = c("cluster4", "cluster5", "cluster6"),
    other_clustering = c("clusterA", "clusterA", "clusterB")
  )

  meta_renamed <- rename_active_ident(meta)
  expect_identical(meta, meta_renamed)
})

test_that("rename_active_ident works when one column is factor and the other is character", {
  meta <- data.frame(
    active.ident = c("cluster1", "cluster2", "cluster3"),
    cell_type = factor(c("cluster1", "cluster2", "cluster3")),
    other_clustering = c("clusterA", "clusterA", "clusterB")
  )

  expected_name <- "cell_type"

  meta_renamed <- rename_active_ident(meta)

  expect_equal(colnames(meta_renamed)[1], expected_name)
  expect_false("active.ident" %in% colnames(meta_renamed))
  expect_equal(colnames(meta_renamed), c("cell_type", "other_clustering"))
})

# spatial object as reconstruct_seurat_spatial leaves it: a Xenium-shaped FOV
# per sample, named by sample, with molecules kept in misc
mock_spatial_obj2s_scdata <- function(technology = "xenium", segmentations = TRUE) {
  scdata <- mock_scdata()
  scdata$seurat_clusters <- rep(letters[1:4], length.out = ncol(scdata))

  image_cells <- split(colnames(scdata), rep(c("A", "B"), length.out = ncol(scdata)))
  scdata$samples <- NA
  for (sample in names(image_cells)) {
    cells <- image_cells[[sample]]
    scdata@meta.data[cells, "samples"] <- sample

    coords <- data.frame(x = seq_along(cells) * 10, y = seq_along(cells) * 10, cell = cells)
    squares <- do.call(rbind, lapply(seq_along(cells), function(i) {
      data.frame(
        x = coords$x[i] + c(-2, 2, 2, -2, -2),
        y = coords$y[i] + c(-2, -2, 2, 2, -2),
        cell = cells[i]
      )
    }))
    fov <- build_xenium_fov(list(centroids = coords, segmentations = if (segmentations) squares))
    SeuratObject::Key(fov) <- paste0("image", sample, "_")
    scdata[[sample]] <- fov
  }

  scdata@misc$technology <- technology
  scdata@misc$transcripts <- list(
    A = data.frame(x = 1, y = 1, gene = "MS4A1"),
    B = data.frame(x = 2, y = 2, gene = "CD79B")
  )
  scdata
}

test_that("rename_obj2s_images names images by sample id as gem2s does", {
  scdata <- mock_spatial_obj2s_scdata()
  input <- list(sampleNames = c("A", "B"), sampleIds = c("id-a", "id-b"))

  scdata <- change_sample_names_to_ids(scdata, input)
  scdata <- rename_obj2s_images(scdata)
  expect_setequal(Seurat::Images(scdata), c("id-a.polygons", "id-b.polygons"))
})

test_that("get_obj2s_sample_transcripts keys molecules by sample id", {
  scdata <- mock_spatial_obj2s_scdata()
  input <- list(sampleNames = c("B", "A"), sampleIds = c("id-b", "id-a"))

  transcripts <- get_obj2s_sample_transcripts(scdata, input)
  expect_equal(transcripts[["id-a"]]$gene, "MS4A1")
  expect_equal(transcripts[["id-b"]]$gene, "CD79B")
})

test_that("upload_obj2s_to_aws uploads the matrix dir and Xenium files linked to their sample", {
  scdata <- mock_spatial_obj2s_scdata()
  scdata[["RNA"]] <- suppressWarnings(as(scdata[["RNA"]], "Assay5"))
  counts <- scdata[["RNA"]]$counts
  scdata[["RNA"]]$counts <- BPCells::write_matrix_dir(counts, withr::local_tempfile())

  uploaded_rds <- NULL
  mockery::stub(upload_obj2s_to_aws, "put_object_in_s3", mock_put_object_in_s3)
  mockery::stub(upload_obj2s_to_aws, "upload_matrix_to_s3", function(pipeline_config, experiment_id, data) {
    uploaded_rds <<- data
    "object_key"
  })

  save_matrix_dir <- mockery::mock("tarfile")
  upload_matrix_dir <- mockery::mock(NULL)
  mockery::stub(upload_obj2s_to_aws, "save_matrix_dir_to_tarfile", save_matrix_dir)
  mockery::stub(upload_obj2s_to_aws, "upload_matrix_dir_to_s3", upload_matrix_dir)

  upload_spatial <- mockery::mock(NULL)
  mockery::stub(upload_obj2s_to_aws, "upload_obj2s_spatial_to_s3", upload_spatial)

  input <- list(experimentId = "1234", sampleIds = "obj2s-sample")
  prev_out <- list(
    scdata = scdata,
    config = list(input = list(type = "seurat_spatial_object"))
  )

  res <- upload_obj2s_to_aws(input, NULL, prev_out)

  # detected technology reported to be saved on the sample
  expect_equal(res$data$item$meta$type, "seurat_spatial_object")
  expect_equal(res$data$item$meta$spatialTechnology, "xenium")

  # molecules are not kept on the processed object
  expect_null(uploaded_rds@misc$transcripts)
  expect_true(all(grepl("[.]polygons$", Seurat::Images(uploaded_rds))))

  mockery::expect_called(upload_matrix_dir, 1)
  expect_equal(mockery::mock_args(upload_matrix_dir)[[1]][[3]], "tarfile")

  # spatial files uploaded for the single obj2s sample, molecules by sample id
  mockery::expect_called(upload_spatial, 1)
  args <- mockery::mock_args(upload_spatial)[[1]]
  spatial_input <- args[[2]]
  expect_equal(spatial_input$obj2sSampleId, "obj2s-sample")
  expect_setequal(names(args[[5]]), spatial_input$sampleIds)
})

test_that("upload_obj2s_to_aws does not report a spatial technology for single cell objects", {
  mockery::stub(upload_obj2s_to_aws, "put_object_in_s3", mock_put_object_in_s3)
  mockery::stub(upload_obj2s_to_aws, "upload_matrix_to_s3", mock_upload_matrix_to_s3)
  upload_spatial <- mockery::mock(NULL)
  mockery::stub(upload_obj2s_to_aws, "upload_obj2s_spatial_to_s3", upload_spatial)

  scdata <- mock_scdata()
  scdata$seurat_clusters <- rep(letters[1:8], length.out = ncol(scdata))
  prev_out <- list(scdata = scdata, config = list(input = list(type = "seurat_object")))

  res <- upload_obj2s_to_aws(list(experimentId = "1234"), NULL, prev_out)

  expect_null(res$data$item$meta$spatialTechnology)
  mockery::expect_called(upload_spatial, 0)
})

test_that("upload_obj2s_spatial_to_s3 uploads Xenium segmentations and molecules linked to their sample", {
  scdata <- mock_spatial_obj2s_scdata()
  scdata$cells_id <- seq_len(ncol(scdata)) - 1
  input <- list(
    obj2sSampleId = "obj2s-sample",
    sampleNames = c("A", "B"),
    sampleIds = c("id-a", "id-b")
  )
  transcripts <- get_obj2s_sample_transcripts(scdata, input)
  scdata <- change_sample_names_to_ids(scdata, input)
  scdata <- rename_obj2s_images(scdata)

  upload_image <- mockery::mock(NULL, cycle = TRUE)
  upload_polygons <- mockery::mock(NULL, cycle = TRUE)
  upload_molecules <- mockery::mock(NULL, cycle = TRUE)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_image_to_s3", upload_image)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_polygons_to_s3", upload_polygons)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_molecules_to_s3", upload_molecules)

  upload_obj2s_spatial_to_s3(NULL, input, "1234", scdata, transcripts)

  # Xenium has no tissue image
  mockery::expect_called(upload_image, 0)
  mockery::expect_called(upload_polygons, 2)
  mockery::expect_called(upload_molecules, 2)

  polygons_args <- mockery::mock_args(upload_polygons)
  expect_setequal(
    unname(sapply(polygons_args, `[[`, "obj2s_sample_id")),
    c("id-a", "id-b")
  )
  # attached to the single obj2s sample: the first upload replaces files from
  # previous runs, the rest are added alongside it
  expect_true(all(sapply(polygons_args, `[[`, 9) == "obj2s-sample"))
  expect_equal(unname(sapply(polygons_args, `[[`, "overwrite_existing")), c(TRUE, FALSE))

  molecules_args <- mockery::mock_args(upload_molecules)
  expect_equal(unname(sapply(molecules_args, `[[`, "overwrite_existing")), c(TRUE, FALSE))
  genes <- setNames(
    sapply(molecules_args, function(args) args[[4]]$gene),
    sapply(molecules_args, `[[`, "obj2s_sample_id")
  )
  expect_equal(genes[["id-a"]], "MS4A1")
  expect_equal(genes[["id-b"]], "CD79B")
})

test_that("upload_obj2s_spatial_to_s3 uploads only molecules for Xenium without segmentations", {
  scdata <- mock_spatial_obj2s_scdata(segmentations = FALSE)
  scdata$cells_id <- seq_len(ncol(scdata)) - 1
  input <- list(
    obj2sSampleId = "obj2s-sample",
    sampleNames = c("A", "B"),
    sampleIds = c("id-a", "id-b")
  )
  transcripts <- get_obj2s_sample_transcripts(scdata, input)
  scdata <- change_sample_names_to_ids(scdata, input)
  scdata <- rename_obj2s_images(scdata)

  upload_image <- mockery::mock(NULL, cycle = TRUE)
  upload_polygons <- mockery::mock(NULL, cycle = TRUE)
  upload_molecules <- mockery::mock(NULL, cycle = TRUE)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_image_to_s3", upload_image)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_polygons_to_s3", upload_polygons)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_molecules_to_s3", upload_molecules)

  upload_obj2s_spatial_to_s3(NULL, input, "1234", scdata, transcripts)

  # cells are shown at their centroids, so there are no polygons to upload
  mockery::expect_called(upload_image, 0)
  mockery::expect_called(upload_polygons, 0)
  mockery::expect_called(upload_molecules, 2)
})

test_that("upload_obj2s_spatial_to_s3 uploads padded Visium images keyed by sample", {
  scdata <- mock_scdata()
  scdata$cells_id <- seq_len(ncol(scdata)) - 1
  image_cells <- split(colnames(scdata), rep(c("A", "B"), length.out = ncol(scdata)))

  scale_factors <- list(spot = 1, fiducial = 1, hires = 1, lowres = 0.1)
  class(scale_factors) <- "scalefactors"
  dims <- list(A = c(700, 800), B = c(900, 650))
  for (sample in names(image_cells)) {
    cells <- image_cells[[sample]]
    scdata@meta.data[cells, "samples"] <- paste0("id-", tolower(sample))
    coords <- data.frame(x = seq_along(cells), y = seq_along(cells), cell = cells)
    fov <- SeuratObject::CreateFOV(
      coords = list(centroids = SeuratObject::CreateCentroids(coords)),
      type = "centroids", assay = "RNA", key = paste0("image", sample, "_")
    )
    scdata[[paste0("id-", tolower(sample), ".polygons")]] <- methods::new(
      "VisiumV2",
      molecules = fov@molecules, boundaries = fov@boundaries,
      assay = fov@assay, key = fov@key,
      image = array(0.5, dim = c(dims[[sample]], 3)),
      scale.factors = scale_factors
    )
  }
  scdata@misc$technology <- "visium"

  upload_image <- mockery::mock(NULL, cycle = TRUE)
  upload_polygons <- mockery::mock(NULL, cycle = TRUE)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_image_to_s3", upload_image)
  mockery::stub(upload_obj2s_spatial_to_s3, "upload_polygons_to_s3", upload_polygons)

  input <- list(obj2sSampleId = "obj2s-sample")
  upload_obj2s_spatial_to_s3(NULL, input, "1234", scdata)

  mockery::expect_called(upload_image, 2)
  mockery::expect_called(upload_polygons, 0)

  image_args <- mockery::mock_args(upload_image)
  # image file id is the sample cell set id, as before files were linked
  expect_setequal(unname(sapply(image_args, `[[`, 6)), c("id-a", "id-b"))
  expect_setequal(unname(sapply(image_args, `[[`, "obj2s_sample_id")), c("id-a", "id-b"))
  # padded to the largest height and width
  expect_true(all(sapply(image_args, `[[`, 8) == 900))
  expect_true(all(sapply(image_args, `[[`, 9) == 800))
})
