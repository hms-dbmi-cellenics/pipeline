remove_bucket_folder <- function(pipeline_config, bucket, folder) {
  keys_to_remove <- list()

  s3 <- paws::s3(config = pipeline_config$aws_config)
  object_list <- s3$list_objects(bucket, Prefix = paste0(folder, "/"))
  for (object in object_list$Contents) {
    keys_to_remove <- append(keys_to_remove, object$Key)
    s3$delete_object(
      Bucket = bucket,
      Key = object$Key
    )
  }

  if (length(keys_to_remove)) {
    message(
      "\nRemoved files from ", bucket, ":",
      "\n- ", paste0(keys_to_remove, sep = "\n- ")
    )
  }
}

remove_cell_ids <- function(pipeline_config, experiment_id) {
  remove_bucket_folder(pipeline_config, pipeline_config$cells_id_bucket, experiment_id)
}

upload_cells_id <- function(pipeline_config, object_key, cells_id) {
  object_data <- tempfile()
  saveRDS(cells_id, file = object_data)
  message("\nUploading filtered cell ids to S3:")
  put_object_in_s3(
    pipeline_config,
    pipeline_config$cells_id_bucket,
    object_data,
    object_key
  )
  return(object_key)
}

load_processed_scdata <- function(s3, pipeline_config, experiment_id) {
  bucket <- pipeline_config$processed_bucket
  message("Loading processed scdata")

  # List the objects under the experiment_id folder
  objects_response <- s3$list_objects_v2(
    Bucket = bucket,
    Prefix = paste0(experiment_id, "/")
  )

  # Extract object keys
  object_keys <- vapply(objects_response$Contents, `[[`, "", "Key")

  # Identify the appropriate file
  rds_file <- grep("r.rds$", object_keys, value = TRUE)
  qs_file <- grep("r.qs$", object_keys, value = TRUE)

  if (length(qs_file) > 0) {
    message("Found .qs file: ", qs_file)
    key_to_load <- qs_file
    load_rds <- FALSE
  } else if (length(rds_file) > 0) {
    message("Found .rds file: ", rds_file)
    key_to_load <- rds_file
    load_rds <- TRUE
  } else {
    stop("No .qs or .rds files found")
  }

  # Get the object
  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = key_to_load
  )

  if (load_rds) {
    # Load with readRDS
    conn <- gzcon(rawConnection(body))
    obj <- readRDS(conn)
    close(conn)
  } else {
    # Load with qs::qread
    tmp_file <- tempfile(fileext = ".qs")
    tryCatch({
      writeBin(body, tmp_file)
      obj <- qs::qread(tmp_file)
    }, finally = {
      # Clean up the temporary file
      unlink(tmp_file)
    })
  }

  return(obj)
}

# get_nnzero will return how many non-zero counts the count matrix has
# it is used to order samples according to their size
get_nnzero <- function(x) {
  nfeat_colname <- paste0(
    "nFeature_",
    Seurat::DefaultAssay(x)
  )

  return(sum(x@meta.data[[nfeat_colname]]))
}

order_by_size <- function(scdata_list) {
  return(scdata_list[order(sapply(scdata_list, get_nnzero))])
}


load_source_scdata_list <- function(s3, pipeline_config, experiment_id) {
  bucket <- pipeline_config$source_bucket
  objects <- s3$list_objects(
    Bucket = bucket,
    Prefix = experiment_id
  )
  experiment_files <- objects$Contents
  message(
    "\nLoading source scdata for experiment: ", experiment_id,
    "\n- bucket: ", bucket
  )

  sample_files <- group_files_by_sample(experiment_files)

  # Load each sample
  scdata_list <- list()
  for (sample_id in names(sample_files)) {
    file_keys <- sample_files[[sample_id]]
    rds_key <- file_keys[["r.rds"]]
    matrix_dir_key <- file_keys[["matrix_dir.tar.zst"]]

    scdata <- load_source_rds(s3, rds_key, bucket)

    if (!is.null(matrix_dir_key)) {
      # untar and update bpcells matrix dir
      new_matrix_dir <-
        load_source_matrix_dir(s3, matrix_dir_key, bucket, sample_id)

      scdata@assays$RNA@layers <- lapply(
        scdata@assays$RNA@layers,
        replace_matrix_dir_paths,
        new_matrix_dir
      )
    }

    scdata_list[[sample_id]] <- scdata
  }

  # order samples according to their size
  # to make the merge independent of samples order in the UI
  return(order_by_size(scdata_list))
}

load_source_rds <- function(s3, rds_key, bucket) {
  message("... loading key ", rds_key)

  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = rds_key
  )

  conn <- gzcon(rawConnection(body))
  obj <- readRDS(conn)
  # close connection explicitly or
  # R will run out of available connections and fail
  close(conn)
  return(obj)
}

load_source_matrix_dir <- function(s3, matrix_dir_key, bucket, sample_id) {
  message("... loading key ", matrix_dir_key)

  # download directly to tar file
  new_matrix_dir <- file.path(tempdir(), paste0(sample_id, "_matrix_dir"))
  tarfile <- paste0(new_matrix_dir, ".tar.zst")

  s3$download_file(
    Bucket = bucket,
    Key = matrix_dir_key,
    Filename = tarfile
  )

  # extract contents ({sample}_matrix_dir/*) to tempdir
  untar_zstd(tarfile, exdir = tempdir())
  unlink(tarfile)
  return(new_matrix_dir)
}

download_processed_matrix_dir <- function(s3, pipeline_config, experiment_id) {
  bucket <- pipeline_config$processed_bucket
  object_key <- paste0(experiment_id, "/matrix_dir.tar.zst")
  message("\nDownloading processed matrix dir: ", object_key)

  # download directly to tar file
  tarfile <- file.path(tempdir(), paste0(experiment_id, "_matrix_dir.tar.zst"))

  s3$download_file(
    Bucket = bucket,
    Key = object_key,
    Filename = tarfile
  )

  return(tarfile)
}

# Update BPCells matrix paths after extraction
replace_matrix_dir_paths <- function(obj, new_dir) {
  new_dir <- normalizePath(new_dir)
  if (!dir.exists(new_dir)) {
    stop(new_dir, " doesn't exist. Please move BPcells folder first.")
  }

  if (inherits(obj, "MatrixDir")) {
    obj@dir <- new_dir
    return(obj)
  }
  if (is.list(obj)) {
    return(lapply(obj, replace_matrix_dir_paths, new_dir))
  }
  for (sn in slotNames(obj)) {
    slot(obj, sn) <- replace_matrix_dir_paths(slot(obj, sn), new_dir)
  }
  return(obj)
}

group_files_by_sample <- function(experiment_files) {
  # Group files by sample
  sample_files <- list()
  for (experiment_file in experiment_files) {
    key <- experiment_file$Key
    parts <- strsplit(key, "/")[[1]]
    sample_id <- parts[2]
    file_name <- parts[3]
    sample_files[[sample_id]][[file_name]] <- key
  }
  return(sample_files)
}

# reload_data_from_s3 will reload:
# * scdata_list for all steps before integration (included)
# * scdata file for all steps after data integration
reload_data_from_s3 <- function(pipeline_config, experiment_id, task_name, tasks) {
  task_names <- names(tasks)
  integration_index <- match("dataIntegration", task_names)
  s3 <- paws::s3(config = pipeline_config$aws_config)

  # TODO: remove if block
  # this never runs, because embed and cluster runs in the worker if modified.
  # If the task is after data integration, we need scdata from processed_matrix
  if (match(task_name, task_names) > integration_index) {
    return(load_processed_scdata(s3, pipeline_config, experiment_id))
  }

  # Otherwise, return scdata_list
  return(load_source_scdata_list(s3, pipeline_config, experiment_id))
}

load_cells_id_from_s3 <- function(pipeline_config, experiment_id, task_name, tasks, samples) {
  s3 <- paws::s3(config = pipeline_config$aws_config)
  object_list <- s3$list_objects(
    Bucket = pipeline_config$cells_id_bucket,
    Prefix = paste0(experiment_id, "/", task_name, "/")
  )
  message("- bucket: ", pipeline_config$cells_id_bucket)
  message("- task name: ", task_name)

  cells_id <- list()
  message("\nTotal of ", length(object_list$Contents), " samples:")
  task_names <- names(tasks)
  integration_index <- match("dataIntegration", task_names)

  # after data integration the cell ids
  # are no longer used because there is no filtering
  # so they are not loaded
  if (match(task_name, task_names) <= integration_index) {
    for (object in object_list$Contents) {
      key <- object$Key
      sample_id <- tools::file_path_sans_ext(basename(key))

      if (!sample_id %in% samples) {
        message(
          "Unexpected filtered ids object for sample ", sample_id,
          ". Filtered cell ids' objects for removed samples ",
          "should be removed in the first step of the QC pipeline ",
          "after gem2s is triggered (due to sample removal). ",
          "Removing it now."
        )
        s3$delete_object(
          Bucket = pipeline_config$cells_id_bucket,
          Key = key
        )
        next
      }

      c(body, ...rest) %<-% s3$get_object(
        Bucket = pipeline_config$cells_id_bucket,
        Key = key
      )
      id_file <- tempfile()
      writeBin(body, con = id_file)
      sample <- readRDS(id_file)
      cells_id[[sample_id]] <- sample[[sample_id]]
      message(
        "- Sample ", sample_id,
        " with ", length(cells_id[[sample_id]]), " cells"
      )
    }
  }
  return(cells_id)
}

build_qc_response <- function(id, input, error, pipeline_config) {
  msg <- list(
    experimentId = input$experimentId,
    taskName = input$taskName,
    input = input,
    response = list(
      error = error
    ),
    pipelineVersion = pipeline_version,
    apiUrl = pipeline_config$api_url
  )

  if (!is.null(id)) {
    msg$output <- list(
      bucket = pipeline_config$results_bucket,
      key = id
    )
  }

  return(msg)
}

send_output_to_api <- function(pipeline_config, input, plot_data_keys, output) {
  c(config, plot_data = plotData) %<-% output

  config <- config[!names(config) %in% exclude_from_config]

  # upload output
  id <- ids::uuid()
  output <- RJSONIO::toJSON(
    list(
      config = config,
      plotDataKeys = plot_data_keys,
      plotData = plot_data
    )
  )

  message("\nUploading results to S3:")
  put_object_in_s3(
    pipeline_config,
    pipeline_config$results_bucket,
    charToRaw(output),
    id
  )

  message("\nSending to SNS topic: ", pipeline_config$sns_topic)
  sns <- paws::sns(config = pipeline_config$aws_config)

  message("... building the message")
  msg <- build_qc_response(id, input, FALSE, pipeline_config)

  message("... publishing the message")
  result <- sns$publish(
    Message = RJSONIO::toJSON(msg),
    TopicArn = pipeline_config$sns_topic,
    MessageAttributes = list(
      type = list(
        DataType = "String",
        StringValue = "PipelineResponse",
        BinaryValue = NULL
      )
    )
  )
  message("... done publishing")

  return(result$MessageId)
}

send_pipeline_update_to_api <- function(pipeline_config, experiment_id, task_name, data, input, string_value) {
  message("\nSending to SNS topic ", pipeline_config$sns_topic)
  sns <- paws::sns(config = pipeline_config$aws_config)
  job_id <- Sys.getenv("AWS_BATCH_JOB_ID", unset = "")

  # TODO -REMOVE DUPLICATE AUTHJWT IN RESPONSE
  msg <- c(
    data,
    taskName = list(task_name),
    experimentId = list(experiment_id),
    jobId = list(job_id),
    authJWT = list(input$auth_JWT),
    input = list(input),
    apiUrl = pipeline_config$api_url
  )

  result <- sns$publish(
    Message = RJSONIO::toJSON(msg),
    TopicArn = pipeline_config$sns_topic,
    MessageAttributes = list(
      type = list(
        DataType = "String",
        StringValue = string_value,
        BinaryValue = NULL
      )
    )
  )

  return(result$MessageId)
}

send_pipeline_fail_update <- function(pipeline_config, input, error_message) {
  sns <- paws::sns(config = pipeline_config$aws_config)

  string_value <- ""
  response <- ""

  process_name <- input$processName
  if (process_name == "qc") {
    string_value <- "PipelineResponse"

    global_env_config <- NULL
    # We don't have dynamic config set for non-sample steps so we can ignore it in that cases
    if (!is.null(input$sampleUuid)) {
      config_key <- paste0("config-", input$taskName, "-", input$sampleUuid)
      global_env_config <- get0(config_key, envir = globalenv(), ifnotfound = NULL)
    }

    # If the step didn't backup any config, don't upload anything
    if (is.null(global_env_config)) {
      response <- build_qc_response(NULL, input, process_name, pipeline_config)
    } else {
      # upload output
      id <- ids::uuid()
      output <- RJSONIO::toJSON(
        list(
          config = global_env_config
        )
      )

      message("\nUploading config to S3:")
      put_object_in_s3(
        pipeline_config,
        pipeline_config$results_bucket,
        charToRaw(output), id
      )

      response <- build_qc_response(id, input, process_name, pipeline_config)
    }
  } else if (process_name %in% c("gem2s", "obj2s", "subset", "copy")) {
    string_value <- ifelse(process_name == "obj2s", "OBJ2SResponse", "GEM2SResponse")

    # TODO - REMOVE THE DUPLICATE EXPERIMENT ID FROM INPUT RESPONSE
    response <- list(
      experimentId = input$experimentId,
      taskName = input$taskName,
      input = input,
      apiUrl = pipeline_config$api_url,
      response = list(
        error = process_name
      )
    )
  } else {
    message(paste("Invalid process_name given: ", process_name))
    return()
  }

  result <- sns$publish(
    Message = RJSONIO::toJSON(response),
    TopicArn = pipeline_config$sns_topic,
    MessageAttributes = list(
      type = list(
        DataType = "String",
        StringValue = string_value,
        BinaryValue = NULL
      )
    )
  )
}

send_plot_data_to_s3 <- function(pipeline_config, experiment_id, output) {
  plot_data <- output$plotData

  s3 <- paws::s3(config = pipeline_config$aws_config)

  plot_names <- names(plot_data)

  plot_data_keys <- list()

  for (plot_data_name in names(plot_data)) {
    id <- ids::uuid()

    plot_data_keys[[plot_data_name]] <- id

    output <- RJSONIO::toJSON(
      list(
        plotData = plot_data[[plot_data_name]]
      )
    )

    tags <- paste0(
      "experimentId=", experiment_id,
      "&plotUuid=", plot_data_name
    )

    message("\nUploading plot data to S3:")
    put_object_in_s3(
      pipeline_config,
      pipeline_config$plot_data_bucket,
      charToRaw(output),
      id,
      tags
    )
  }

  return(plot_data_keys)
}

upload_matrix_to_s3 <- function(pipeline_config, experiment_id, data) {
  object_key <- paste0(experiment_id, "/r.qs")

  message("\nSaving count matrix.")
  count_matrix <- tempfile()
  qs::qsave(data, file = count_matrix)

  message("Count matrix file size: ", fs::file_size(count_matrix))

  message("\nUploading count matrix to S3:")
  put_object_in_s3_multipart(
    pipeline_config,
    pipeline_config$processed_bucket,
    count_matrix,
    object_key
  )

  return(object_key)
}

save_matrix_dir_to_tarfile <- function(experiment_id, data) {
  matrix_dir <- get_matrix_dirs(data)

  message("\nTarring experiment matrix dir: ", matrix_dir)
  tarfile <- tar_matrix_dir(experiment_id, matrix_dir)
  message("Tar file size: ", fs::file_size(tarfile))
  return(tarfile)
}

upload_matrix_dir_to_s3 <- function(
  pipeline_config, experiment_id, matrix_dir_tarfile
) {
  object_key <- paste0(experiment_id, "/matrix_dir.tar.zst")
  message("\nSaving matrix dir to: ", object_key)

  message("\nUploading matrix dir to S3:")
  put_object_in_s3_multipart(
    pipeline_config,
    pipeline_config$processed_bucket,
    matrix_dir_tarfile,
    object_key
  )
}



find_matrix_dir_paths <- function(obj) {
  matrix_dir_paths <- c()
  if (inherits(obj, "MatrixDir")) {
    return(obj@dir)
  }
  if (is.list(obj)) {
    res <- lapply(obj, find_matrix_dir_paths)
    return(unlist(res))
  }
  for (sn in slotNames(obj)) {
    matrix_dir_paths <- c(
      matrix_dir_paths,
      find_matrix_dir_paths(slot(obj, sn))
    )
  }
  return(matrix_dir_paths)
}

get_matrix_dirs <- function(scdata) {
  matrix_dirs <- lapply(
    scdata@assays$RNA@layers,
    find_matrix_dir_paths
  )

  matrix_dirs <- unique(unlist(matrix_dirs))
  return(matrix_dirs)
}


upload_image_to_s3 <- function(
  pipeline_config, input, experiment_id, image_arr, image_name, image_id, sample_id,
  image_max_height = dim(image_arr)[1], overwrite_existing = TRUE
) {
  # things for api requests
  api_url <- pipeline_config$api_url
  auth_jwt <- input$authJWT

  # pad images to have same fixed height
  image_arr <- pad_image_height(image_arr, image_max_height)

  # where to save zarr folder locally
  zarr_name <- paste0(image_name, ".ome.zarr")
  output_path <- file.path(tempdir(), zarr_name)

  message("Saving image data to: ", output_path, "...")

  # save as ome zarr folder
  rgb_image_to_ome_zarr(image_arr, output_path, image_name)

  # zip all files in zarr folder
  zip_name <- paste0(zarr_name, ".zip")
  zip_path <- file.path(tempdir(), zip_name)

  workdir <- getwd()
  setwd(output_path)
  utils::zip(zip_path, files = ".", flags = "-rq0")
  setwd(workdir)

  # log file size before upload
  message("Image zip file size: ", fs::file_size(zip_path))

  # upload ome.zarr.zip to s3
  # use multipart upload for files larger than 100MB to handle large images
  message("\nUploading image data to S3:")
  put_object_in_s3_multipart(
    pipeline_config,
    bucket = pipeline_config$spatial_image_bucket,
    object = zip_path,
    key = image_id
  )

  # create sql entry in sample_file,
  # (also creates entry in sample_to_sample_file_map)
  create_sample_file(
    api_url,
    experiment_id,
    sample_id,
    "ome_zarr_zip",
    file.size(zip_path),
    # gets used as s3_path by API
    image_id,
    # FALSE for obj2s so that can have multiple ome_zarr_zip for single "sample"
    overwrite_existing,
    auth_jwt
  )

  invisible()
}

# based off of vitessce-python demo
rgb_image_to_ome_zarr <- function(image_arr, output_path, image_name) {
  np <- reticulate::import("numpy")
  zarr <- reticulate::import("zarr")
  ome_zarr <- reticulate::import("ome_zarr")

  # Convert [0, 1] → [0, 255] uint8, then transpose to (C, H, W) "cyx"
  image_arr <- np$array(image_arr * 255.0, dtype = "uint8")
  image_arr <- np$transpose(image_arr, c(2L, 0L, 1L))

  default_window <- list(
    start = 0,
    min = 0,
    max = 255,
    end = 255
  )

  z_root <- zarr$open_group(output_path, mode = "w")

  # scale_factors [2, 4, 8, 16] → 5 levels (0 = full res, 1–4 = coarse)
  # method NEAREST matches skimage order=0 resize — correct for uint8 RGB
  # (no fractional pixel values can appear so any method is safe here,
  # but NEAREST is the fastest and avoids any anti-aliasing artefacts)
  scale_factors <- 2L ^ seq_len(4L)

  ome_zarr$writer$write_image(
    image = image_arr,
    group = z_root,
    axes = "cyx",
    scale_factors = scale_factors,
    storage_options = list(
      chunks = reticulate::tuple(1L, 512L, 512L),
      dimension_separator = "/"
    )
  )

  omero <- list(
    "name" = image_name,
    "version" = "0.3",
    "rdefs" = list(),
    "channels" = list(
      list("label" = "R", "color" = "FF0000", "window" = default_window),
      list("label" = "G", "color" = "00FF00", "window" = default_window),
      list("label" = "B", "color" = "0000FF", "window" = default_window)
    )
  )
  reticulate::py_set_attr(z_root, "omero", omero)

  invisible()
}

upload_polygons_to_s3 <- function(
  pipeline_config, input, experiment_id, polygons, height, width,
  image_name, polygons_id, sample_id, overwrite_existing
) {
  # things for api requests
  api_url <- pipeline_config$api_url
  auth_jwt <- input$authJWT

  # where to save zarr folder locally
  zarr_name <- paste0(image_name, ".segmentations.ome.zarr")
  output_path <- file.path(tempdir(), zarr_name)

  message("\nSaving polygons data: ", output_path)

  # save as ome zarr folder
  polygons_to_bitmask_ome_zarr(
    polygons,
    height,
    width,
    output_path
  )

  # zip all files in zarr folder
  zip_name <- paste0(zarr_name, ".zip")
  zip_path <- file.path(tempdir(), zip_name)

  workdir <- getwd()
  setwd(output_path)
  utils::zip(zip_path, files = ".", flags = "-rq0")
  setwd(workdir)

  # log file size before upload
  message("Segmentations zip file size: ", fs::file_size(zip_path))

  # upload ome.zarr.zip to s3
  message("\nUploading polygons data to S3:")
  put_object_in_s3_multipart(
    pipeline_config,
    bucket = pipeline_config$spatial_segmentations_bucket,
    object = zip_path,
    key = polygons_id
  )

  # create sql entry in sample_file,
  # (also creates entry in sample_to_sample_file_map)
  create_sample_file(
    api_url,
    experiment_id,
    sample_id,
    "segmentations_ome_zarr_zip",
    file.size(zip_path),
    # gets used as s3_path by API
    polygons_id,
    overwrite_existing,
    auth_jwt
  )

  invisible()
}

#' Default Xenium transcript quality-value threshold (Q20).
#'
#' The installed \code{ReadXenium} does not apply \code{mols.qv.threshold}; we
#' filter explicitly here (per the format contract) when a qv column is present.
XENIUM_QV_THRESHOLD <- 20

#' Build and upload the gene-partitioned molecule artifact for a Xenium sample
#'
#' Modeled on \code{upload_polygons_to_s3}. Reads the raw transcripts frame
#' (x/y/gene, optionally qv), applies the QV filter, maps gene -> dense 0-based
#' feature_code, then writes ONE Feather file PER GENE (\code{{code}.feather},
#' columns x/y) plus our meta.json (the feature dictionary + build metadata, per
#' the format contract; the UI derives gene colours from the code). Zipped stored
#' (range-readable), uploaded to the spatial_molecules_bucket and registered as a
#' molecules_by_gene sample file. Built straight from the frame (no Seurat object).
#'
#' Gene-partitioned (not a spatial quadtree): the UI only ever renders a few genes
#' at a time and deck.gl renders that point count directly, so a gene's molecules
#' are range-read from its own entry with no level-of-detail / spatial tiling.
#'
#' @param transcripts data.frame with columns x, y, gene (+ qv if present)
upload_molecules_to_s3 <- function(
  pipeline_config, input, experiment_id, transcripts, sample_id,
  overwrite_existing
) {
  api_url <- pipeline_config$api_url
  auth_jwt <- input$authJWT

  # QV filter (contract risk #6): drop sub-threshold molecules when qv present.
  qv_threshold <- NULL
  if ("qv" %in% colnames(transcripts)) {
    qv_threshold <- XENIUM_QV_THRESHOLD
    before <- nrow(transcripts)
    transcripts <- transcripts[transcripts$qv >= qv_threshold, , drop = FALSE]
    message(
      "QV filter (Q", qv_threshold, "): kept ", nrow(transcripts),
      " of ", before, " molecules"
    )
  }

  # Dense 0-based feature dictionary: gene -> feature_code. The code is a stable,
  # panel-wide index (genes sorted, so it doesn't shift with the selection); the UI
  # derives a colour from it (utils/spatial/moleculeColors) so
  # colour assignment isn't baked here.
  genes <- sort(unique(transcripts$gene))
  feature_dict <- data.frame(
    code = seq_along(genes) - 1L,
    gene = genes,
    stringsAsFactors = FALSE
  )

  x <- as.numeric(transcripts$x)
  y <- as.numeric(transcripts$y)
  # full micron bbox over the kept molecules (drives the UI's zoomed-out fit).
  x_range <- range(x)
  y_range <- range(y)

  molecules_dir <- file.path(tempdir(), paste0(sample_id, ".molecules.bygene"))
  unlink(molecules_dir, recursive = TRUE)
  dir.create(molecules_dir, recursive = TRUE)

  # One Feather entry per gene: x/y as float32. The entry IS the gene, so no
  # feature_code column is stored inside; the code<->gene map is meta.json.
  # float32 (not float64): halves the bytes and the UI downcasts to Float32Array on
  # read anyway; at Xenium's micron scale (~1e4) float32 ULP is ~nm, far below any
  # visible scale. compression = "zstd": the Arrow IPC body is ZSTD-compressed, which
  # apache-arrow JS decompresses on read (>= v21.1.0). So the zip below stays stored.
  message("\nWriting per-gene molecule Feather entries: ", molecules_dir)
  rows_by_gene <- split(seq_along(x), transcripts$gene)
  gene_entries <- vapply(feature_dict$gene, function(gene) {
    rows <- rows_by_gene[[gene]]
    entry <- paste0(feature_dict$code[feature_dict$gene == gene], ".feather")
    gene_table <- arrow::arrow_table(
      x = arrow::Array$create(x[rows], type = arrow::float32()),
      y = arrow::Array$create(y[rows], type = arrow::float32())
    )
    arrow::write_feather(
      gene_table, file.path(molecules_dir, entry),
      compression = "zstd"
    )
    entry
  }, character(1))
  n_points <- vapply(feature_dict$gene, function(gene) {
    length(rows_by_gene[[gene]])
  }, integer(1))

  meta <- list(
    version = 2L,
    qvThreshold = qv_threshold,
    rootExtent = list(x = x_range, y = y_range),
    genes = lapply(seq_len(nrow(feature_dict)), function(i) {
      list(
        code = feature_dict$code[i],
        gene = feature_dict$gene[i],
        entry = unname(gene_entries[i]),
        nPoints = unname(n_points[i])
      )
    })
  )
  jsonlite::write_json(
    meta,
    file.path(molecules_dir, "meta.json"),
    auto_unbox = TRUE,
    null = "null"
  )

  # zip stored (-rq0): each Feather entry is already ZSTD-compressed, so deflating
  # the zip on top would be wasted work. Stored entries let the UI's ZipFileStore
  # range-read exactly the selected genes' bytes (same pattern as the OME-Zarr zips).
  zip_name <- paste0(sample_id, ".molecules.bygene.zip")
  zip_path <- file.path(tempdir(), zip_name)

  workdir <- getwd()
  setwd(molecules_dir)
  utils::zip(zip_path, files = ".", flags = "-rq0")
  setwd(workdir)

  message("Molecule artifact zip file size: ", fs::file_size(zip_path))

  molecules_id <- uuid::UUIDgenerate()

  message("\nUploading molecule artifact to S3:")
  put_object_in_s3_multipart(
    pipeline_config,
    bucket = pipeline_config$spatial_molecules_bucket,
    object = zip_path,
    key = molecules_id
  )

  create_sample_file(
    api_url,
    experiment_id,
    sample_id,
    "molecules_by_gene",
    file.size(zip_path),
    molecules_id,
    overwrite_existing,
    auth_jwt
  )

  invisible()
}

polygons_to_bitmask_ome_zarr <- function(
  polygons,
  height,
  width,
  output_path,
  max_pyramid_levels = 4L
) {
  np <- reticulate::import("numpy")
  zarr <- reticulate::import("zarr")
  ome_zarr <- reticulate::import("ome_zarr")
  pil <- reticulate::import("PIL")

  # 1-indexed: pixel value 0 reserved for background (transparent in shader)
  polygons$cells_id <- as.integer(polygons$cells_id) + 1L
  cell_ids <- unique(polygons$cells_id)

  # ── Rasterise polygons into a 32-bit label image ──────────────────────────
  # PIL mode "I" = 32-bit signed int; fill values are the 1-indexed cell IDs.
  image <- pil$Image$new(
    "I",
    reticulate::tuple(as.integer(width), as.integer(height)),
    0L
  )
  draw <- pil$ImageDraw$Draw(image)

  data.table::setDT(polygons)
  polygons_list <- vector("list", max(cell_ids))
  pairs <- polygons[, .(v = list(c(rbind(x, y)))), by = cells_id]
  polygons_list[pairs$cells_id] <- pairs$v
  for (cid in cell_ids) {
    draw$polygon(polygons_list[[cid]], fill = as.integer(cid))
  }

  # float32: integer-exact up to 2^24 ≈ 16.7 M cells.
  # XRLayer uploads this as R32F so the BitmaskLayer shader reads raw cell IDs.
  arr <- np$array(image, dtype = np$float32)
  arr <- np$expand_dims(arr, axis = 0L)   # (H, W) → (1, H, W)  "cyx"

  z_root <- zarr$open_group(output_path, mode = "w")

  # ── Pyramid via write_image ───────────────────────────────────────────────
  # write_image generates the pyramid automatically from scale_factors.
  # write_multiscale does NOT have a scale_factors parameter and only stores
  # whatever arrays you pass; using it with a single array was the bug.
  #
  # scale_factors as cumulative integers [2, 4, 8, 16]:
  #   level 0 = original full resolution
  #   level 1 = 1/2  (2×  downsample in y, x)
  #   level 2 = 1/4  (4×  downsample in y, x)
  #   level 3 = 1/8  (8×  downsample in y, x)
  #   level 4 = 1/16 (16× downsample in y, x)
  #   → 5 levels total, matching the RGB image writer
  #
  # method = Methods.NEAREST (nearest-neighbour, skimage order=0):
  #   Selects an existing pixel value at each downsampled position.
  #   Cell IDs (integers stored as float32) are preserved exactly at every level.
  #   Methods.RESIZE (the default) uses anti-aliased bicubic interpolation which
  #   blends adjacent IDs producing fractional values the LUT cannot look up.
  scale_factors <- 2L ^ seq_len(max_pyramid_levels)

  ome_zarr$writer$write_image(
    image = arr,
    group = z_root,
    axes = "cyx",
    scale_factors = scale_factors,
    method = "nearest",
    storage_options = list(
      chunks = reticulate::tuple(1L, 512L, 512L),
      dimension_separator = "/"
    )
  )

  invisible()
}

pad_image_height <- function(image_arr, target_height) {
  current_height <- dim(image_arr)[1]
  if (current_height == target_height) return(image_arr)

  # grab bottom row and repeat to pad bottom of image
  bottom_row <- image_arr[nrow(image_arr), , , drop = FALSE]

  # replace pixels where where green channel is below cutoff
  # with average color of other pixels in row
  # cutoff determined empirically looking at RGB values of slides
  # goal is to not padd with pixels representing tissue (less green)
  green_cutoff <- 210 / 255
  maybe_tissue <- bottom_row[, , 2] < green_cutoff

  if (all(maybe_tissue)) {
    # use pale lavender for padding
    pad_color <- c(241, 234, 241) / 255

  } else {
    # reshape to a 2D matrix (pixels as rows, 3 channels as columns)
    pixel_matrix <- matrix(bottom_row[, !maybe_tissue, ], ncol = 3)

    # use average color of non-tissue pixels for padding
    pad_color <- colMeans(pixel_matrix)

  }

  # set RGB channels
  bottom_row[1, , 1] <- pad_color[1]
  bottom_row[1, , 2] <- pad_color[2]
  bottom_row[1, , 3] <- pad_color[3]


  padding <- bottom_row[rep(1, target_height - current_height), , ]
  padded_image_arr <- abind::abind(image_arr, padding, along = 1)

  message(
    "Padded image from height ", current_height,
    " to ", dim(padded_image_arr)[1]
  )

  return(padded_image_arr)
}

create_sample_file <- function(
  api_url, experiment_id, sample_id, file_type,
  file_size, sample_file_id, overwrite_existing, auth_jwt
) {
  url <- paste0(
    api_url,
    "/v2/experiments/", experiment_id,
    "/samples/", sample_id,
    "/sampleFiles/", file_type
  )

  body <- list(
    sampleFileId = sample_file_id,
    size = file_size,
    uploadStatus = "uploaded",
    overwriteExisting = overwrite_existing
  )

  response <- httr::POST(
    url,
    body = body,
    encode = "json",
    httr::add_headers(
      "Content-Type" = "application/json",
      "Authorization" = auth_jwt
    )
  )

  if (httr::status_code(response) >= 400) {
    stop("API post to create sample file failed with status code: ", httr::status_code(response))
  }
}

upload_obj2s_images_to_s3 <- function(pipeline_config, input, experiment_id, scdata) {
  # use the
  image_ids <- input$sampleIds
  names(image_ids) <- input$sampleNames

  # associate obj2s images with single sample id for experiment
  sample_id <- input$obj2sSampleId
  image_names <- Seurat::Images(scdata)

  for (image_name in image_names) {
    image_arr <- scdata[[image_name]]@image
    image_id <- image_ids[image_name]

    upload_image_to_s3(
      pipeline_config,
      input,
      experiment_id,
      image_arr,
      image_name,
      image_id,
      sample_id,
      overwrite_existing = FALSE
    )
  }
}

put_object_in_s3 <- function(
  pipeline_config,
  bucket,
  object,
  key,
  tagging = NULL
) {
  message(
    "- bucket: ", bucket, "\n",
    "- key: ", key
  )

  s3 <- paws::s3(config = pipeline_config$aws_config)

  retry_count <- 0
  max_retries <- 2

  while (retry_count < max_retries) {
    tryCatch(
      {
        response <- s3$put_object(
          Bucket = bucket,
          Key = key,
          Body = object,
          Tagging = switch(!is.null(tagging),
            tagging
          )
        )
        message("... object successfully uploaded to S3.")
        return(response)
      },
      error = function(e) {
        retry_count <<- retry_count + 1
        message(sprintf(
          "... upload failed. Retrying (%d/%d)", retry_count, max_retries
        ))
        # exponential back off, preventing sending a new request too fast
        Sys.sleep(2^retry_count)
      }
    )
  }
  stop("Failed to upload object to S3 after maximum number of retries.")
}

#' Upload a file to S3 using multipart upload
#'
#' @param pipeline_config A Paws S3 config object, e.g. from `paws::s3()`.
#' @param object The path to the file to be uploaded.
#' @param bucket The name of the S3 bucket to be uploaded to, e.g. `my-bucket`.
#' @param key The name to assign to the file in the S3 bucket e.g. `path/file`.
put_object_in_s3_multipart <- function(pipeline_config, bucket, object, key) {
  message(
    "- object: ", object, "\n",
    "- bucket: ", bucket, "\n",
    "- key: ", key
  )

  s3 <- paws::s3(config = pipeline_config$aws_config)
  multipart <- s3$create_multipart_upload(
    Bucket = bucket,
    Key = key
  )
  resp <- NULL
  on.exit({
    if (is.null(resp) || inherits(resp, "try-error")) {
      s3$abort_multipart_upload(
        Bucket = bucket,
        Key = key,
        UploadId = multipart$UploadId
      )
    }
  })
  resp <- try({
    parts <- upload_multipart_parts(s3, bucket, object, key, multipart$UploadId)
    s3$complete_multipart_upload(
      Bucket = bucket,
      Key = key,
      MultipartUpload = list(Parts = parts),
      UploadId = multipart$UploadId
    )
  })
  return(resp)
}

# The limit for part numbers is 10000
# aws cli uses 8mb default, and recommends up to 64mb for files >1G
upload_multipart_parts <- function(s3, bucket, object, key, upload_id) {
  file_size <- file.size(object)
  megabyte <- 2^20
  gigabyte <- 2^30
  part_size <- ifelse(file_size > gigabyte, 64 * megabyte, 8 * megabyte)
  num_parts <- ceiling(file_size / part_size)

  con <- base::file(object, open = "rb")
  on.exit({
    close(con)
  })
  parts <- list()
  tstart <- Sys.time()
  for (i in 1:num_parts) {
    part <- readBin(con, what = "raw", n = part_size)

    if (i %% 10 == 0 || i == num_parts) {
      message("... uploading part ", i, "/", num_parts)
    }

    # Upload part with retries and exponential backoff
    retry_count <- 0
    max_retries <- 3
    part_resp <- NULL

    while (retry_count < max_retries) {
      tryCatch(
        {
          part_resp <- s3$upload_part(
            Body = part,
            Bucket = bucket,
            Key = key,
            PartNumber = i,
            UploadId = upload_id
          )
          break # Exit loop on successful upload
        },
        error = function(e) {
          retry_count <<- retry_count + 1
          message("... error uploading part ", i, ": ", e$message)

          if (retry_count < max_retries) {
            wait_time <- 2^retry_count
            message(
              "... retrying part ", i,
              " (attempt ", retry_count, "/", max_retries, ")"
            )
            Sys.sleep(wait_time)
          }
        }
      )
    }

    if (is.null(part_resp)) {
      stop("Failed to upload part ", i, " after ", max_retries, " attempts")
    }

    parts <- c(parts, list(list(ETag = part_resp$ETag, PartNumber = i)))
  }
  ttask <- format(Sys.time() - tstart, digits = 2)
  message("... upload time: ", ttask)
  return(parts)
}


#' Load cellsets object from s3 and flatten it for easier manipulation in R
#'
#' @param s3 paws::s3 object
#' @param pipeline_config list
#' @param experiment_id character
#'
#' @return cellsets list
#' @export
#'
load_cellsets <- function(s3, pipeline_config, experiment_id) {
  message("loading cellsets file")

  bucket <- pipeline_config$cell_sets_bucket

  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = experiment_id
  )

  obj <- jsonlite::fromJSON(rawConnection(body), flatten = T)
  return(obj)
}

#' Load cellsets object flattened (as is done by load_cellsets)
#' and change it back to the original shape
#' This is useful to prepare the cell sets to upload back into s3
#'
#' @param s3 paws::s3 object
#' @param pipeline_config list
#' @param experiment_id character
#'
#' @return cellsets list
#' @export
#'
unflatten_cell_sets <- function(cell_sets) {
  formatted <- list()

  for (i in seq_along(cell_sets$key)) {
    cell_class <- list(
      key = cell_sets$key[[i]],
      name = cell_sets$name[[i]],
      rootNode = cell_sets$rootNode[[i]],
      type = cell_sets$type[[i]],
      children = list()
    )

    for (j in seq_along(cell_sets$children[[i]]$key)) {
      children <- cell_sets$children[[i]]

      cell_set <- list(
        key = children$key[[j]],
        name = children$name[[j]],
        rootNode = FALSE,
        color = children$color[[j]],
        # Take the type from the parent, some of the child cell sets don't have "type"
        # Some others do, but all parents do
        type = cell_sets$type[[i]],
        cellIds = ensure_is_list_in_json(children$cellIds[[j]])
      )

      cell_class$children <- append(cell_class$children, list(cell_set))
    }

    formatted <- append(formatted, list(cell_class))
  }

  return(formatted)
}


#' Bind columns not creating rows if there's an empty data.table
#'
#' `cbind` on `data.table` adds a row if binding an empty data.table to a non-empty
#' one. We do not want that behavior when parsing cellsets, because it implies
#' the "creation" of a cell that does not exists (i.e. when binding scratchpad
#' cellsets slots of an experiment without custom cellsets)
#'
#' @param dt data.table
#' @param ... columns to add
#'
#' @return data.table with new columns
#' @export
#'
safe_cbind <- function(dt, ...) {
  if (nrow(dt) > 0) {
    dt <- cbind(dt, ...)
  }
  return(dt)
}


#' add cellset type column to cellsets data.table
#'
#' helper to correctly name the cellset type column. some cellsets already
#' contain a "type" slot, which complicates matters, so we chose `cellset_type`,
#'
#' @param dt data.table
#' @param col string of corresponding cellset type
#'
#' @return data.table with cellset_type column
#' @export
#'
cbind_cellset_type <- function(dt, col) {
  dt <- safe_cbind(dt, cellset_type = col)
}


is_uuid <- function(x) {
  uuid_regex <- "^\\b[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}\\b$"
  return(grepl(uuid_regex, x))
}


get_cellset_type <- function(key, type) {
  cellset_type <- switch(key,
    louvain = "cluster",
    scratchpad = "scratchpad",
    sample = "sample",
    ifelse(is_uuid(key), "sctype", "metadata")
  )

  return(cellset_type)
}


#' Parse cellsets object to data.table
#'
#' Gets the cellsets list and converts it to a tidy data.table
#'
#' @param cellsets list
#'
#' @return data.table of cellset keys, names and corresponding cell_ids
#' @export
#'
parse_cellsets <- function(cellsets) {
  dt_list <- cellsets$cellSets$children
  cellset_types <- purrr::map2_chr(
    cellsets$cellSets$key,
    cellsets$cellSets$type,
    get_cellset_type
  )

  lapply(dt_list, data.table::setDT)

  dt_list <- purrr::map2(dt_list, cellset_types, cbind_cellset_type)

  # fill columns in case there are empty cellset classes
  dt <- data.table::rbindlist(dt_list, fill = TRUE)

  # unnest, and change column name
  dt <- dt[, setNames(.(unlist(cellIds)), "cell_id"), by = .(key, name, cellset_type)]
  data.table::setnames(dt, "cellset_type", "type")
  return(dt)
}

get_s3_rds <- function(bucket, key, aws_config) {
  s3 <- paws::s3(config = aws_config)

  c(body, ...rest) %<-% s3$get_object(
    Bucket = bucket,
    Key = key
  )

  conn <- gzcon(rawConnection(body))
  object <- readRDS(conn)

  close(conn)

  return(object)
}

put_s3_rds <- function(bucket, key, aws_config, rds) {
  s3 <- paws::s3(config = aws_config)

  file <- tempfile()
  saveRDS(rds, file = file)

  s3$put_object(
    Bucket = bucket,
    Key = key,
    Body = file
  )
}
