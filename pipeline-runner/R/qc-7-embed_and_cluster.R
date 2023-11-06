#' run clustering
#'
#'
#' @param scdata seurat object
#' @param config list with clustering parameters
#' @param sample_id character
#' @param cells_id list of cell ids that passed all filters so far
#' @param task_name character
#'
#' @return list
#' @export
#'
embed_and_cluster <-
 function(scdata,
           config,
           sample_id,
           cells_id,
           task_name = "configureEmbedding",
           ignore_ssl_cert = FALSE) {

    clustering_method <- config$clusteringSettings$method
    methodSettings <-
      config$clusteringSettings$methodSettings[[clustering_method]]
    cellSets <-
      runClusters(clustering_method, methodSettings$resolution, scdata)

    formated_cell_sets <-
      format_cluster_cellsets(cellSets, clustering_method, scdata@misc$color_pool)

    replace_cell_class_through_api(
      formated_cell_sets,
      config$api_url,
      scdata@misc$experimentId,
      clustering_method,
      config$auth_JWT,
      ignore_ssl_cert
    )

    # add cl metadata if any
    if (!is.null(config$metadataS3Path)) {
      cl_metadata_cellsets <- make_cl_metadata_cellsets(scdata, config)

      # remove previously uploaded cell level metadata (to allow user to replace the file)
      replace_cl_metadata_through_api(
        cl_metadata_cellsets,
        config$api_url,
        scdata@misc$experimentId,
        config$auth_JWT,
        ignore_ssl_cert
      )
    }


    result <- list(
      data = scdata,
      new_ids = cells_id,
      config = config,
      plotData = list()
    )

    return(result)
  }


format_cluster_cellsets <- function(cell_sets,
                                    clustering_method,
                                    color_pool,
                                    name = paste0(clustering_method, " clusters")) {
  message("Formatting cluster cellsets.")

  # careful with capital l on type for the key.
  cell_sets_object <-
    list(
      key = clustering_method,
      name = name,
      rootNode = TRUE,
      type = "cellSets",
      children = list()
    )
  for (cluster in sort_cluster_names(unique(cell_sets$cluster))) {
    cells <- cell_sets[cell_sets$cluster == cluster, "cell_ids"]
    is.num <- !is.na(as.numeric(cluster))
    set_name <- ifelse(is.num, paste("Cluster", cluster), cluster)

    new_set <- list(
      key = paste0(clustering_method, "-", cluster),
      name = set_name,
      rootNode = FALSE,
      type = "cellSets",
      color = color_pool[1],
      cellIds = ensure_is_list_in_json(unname(cells))
    )
    color_pool <- color_pool[-1]
    cell_sets_object$children <-
      append(cell_sets_object$children, list(new_set))
  }
  return(cell_sets_object)
}



#' replace cell class by key
#'
#' Replaces an entire cell class (group of cellsets) by key. Used to replace
#' clustering cellsets after re-clustering.
#'
#' @param cell_class_object list
#' @param api_url character
#' @param experiment_id character
#' @param cell_class_key character
#' @param auth_JWT character
#' @param ignore_ssl_cert boolean
#'
#' @return NULL
#' @export
#'
replace_cell_class_through_api <-
  function(cell_class_object,
           api_url,
           experiment_id,
           cell_class_key,
           auth_JWT,
           ignore_ssl_cert) {
    message("updating cellsets through API")
    httr_query <- paste0("$[?(@.key == \"", cell_class_key, "\")]")

    if (ignore_ssl_cert) {
      httr::set_config(httr::config(ssl_verifypeer = 0L))
    }

    body <- list(list(
        "$match" = list(query = httr_query, value = list("$remove" = TRUE))
      ),
      list("$prepend" = cell_class_object))

    httr::PATCH(
      paste0(api_url, "/v2/experiments/", experiment_id, "/cellSets"),
      body = body,
      encode = "json",
      httr::add_headers("Content-Type" = "application/boschni-json-merger+json",
                        "Authorization" = auth_JWT)
    )
  }


#' Replace and Append Cell Metadata Through API
#'
#' This function deletes cell-level metadata cellsets in an existing cellsets.json
#' file based on specified types. It then appends new cell-level metadata cellsets.
#'
#' @param cl_metadata_cellsets A list that will be appended to the existing cell sets.
#' @param api_url The base URL of the API.
#' @param experiment_id The unique identifier for the experiment to be modified.
#' @param auth_JWT The JSON Web Token used for authentication.
#' @param ignore_ssl_cert Boolean flag to indicate whether to ignore SSL certificate verification.
#'
#' @return The response from the API after performing the PATCH operation.
#'
#' @export
replace_cl_metadata_through_api <-
  function(cl_metadata_cellsets,
           api_url,
           experiment_id,
           auth_JWT,
           ignore_ssl_cert) {

    if (ignore_ssl_cert) {
      httr::set_config(httr::config(ssl_verifypeer = 0L))
    }

    appends <- list()
    for (i in seq_along(cl_metadata_cellsets)) {
      appends <- append(appends, list(list("$append" = cl_metadata_cellsets[[i]])))
    }

    httr::PATCH(
      paste0(api_url, "/v2/experiments/", experiment_id, "/cellSets"),
      body = appends,
      encode = "json",
      httr::add_headers("Content-Type" = "application/boschni-json-merger+json",
                        "Authorization" = auth_JWT)
    )
  }

#' Download and load cell-level metadata tsv file
#'
#' @param config list with AWS configuration and metadataS3Path
#'
#' @return data.table of cell-level metadata
#' @export
#'
download_cl_metadata_file <- function(config) {

  s3 <- paws::s3(config = config$aws_config)
  s3_path <- basename(config$metadataS3Path)

  # TODO figure out best way to temporarily store file, or read from S3 directly
  file_path <- paste0(basename(s3_path), ".tsv.gz")
  message("downloading cell-level metadata file to ", file_path)
  download_and_store(config$cl_metadata_bucket, s3_path, file_path, s3)
  cl_metadata <- data.table::fread(file_path)

  return(cl_metadata)
}


#' Makes cell-level metadata cellsets
#'
#' Main cell-level metadata cellsets function. Downloads, parses, de-duplicates
#' barcodes, detects variable types and formats cell-level metadata cellsets.
#' Returning a list ready to add to `cellsets.json` file.
#'
#' @param scdata SeuratObject
#' @param config list with AWS config and metadataS3Paths
#'
#' @return list of cell-level metadata cellsets
#' @export
#'
make_cl_metadata_cellsets <- function(scdata, config) {

  cl_metadata <- download_cl_metadata_file(config)

  # TODO deduplicate using sample column (if present) in the cell-level metadata file
  # TODO add "duplicated" variable to cl_metadata table to create cellset for duplicated barcodes
  #cl_metadata <- deduplicate_cl_metadata(scdata, cl_metadata)

  # extract barcode - cell_id (keep sample column for variable type detection)
  barcode_cell_id_map <- get_cell_id_barcode_map(scdata)

  # makes cell_id - metadata table
  cl_metadata <-
    make_cl_metadata_table(cl_metadata, barcode_cell_id_map)

  var_types <- detect_variable_types(cl_metadata)

  # creates cell-level metadata cellsets, setting the correct type
  cl_metadata_cellsets <-
    format_cl_metadata_cellsets(cl_metadata, var_types, scdata@misc$color_pool)

  return(cl_metadata_cellsets)
}


#' extract cell_id - barcode map
#'
#' Helper to get required variables from the SeuratObject
#'
#' @param scdata SeuratObject
#'
#' @return data.table containing cells_id, barcode, samples
#' @export
#'
get_cell_id_barcode_map <- function(scdata) {
  message("extracting cell_id - barcode map")
  # TODO add samples when able to get sample names from database, scdata contains sample ids only
  cols_to_keep <- c("cells_id")
  data.table::as.data.table(scdata@meta.data[cols_to_keep], keep.rownames = "barcode")
}


#' make cell-level metadata - cell_ids table
#'
#' Joins user supplied tsv file with the barcode-cell_id map, so that each cell_id
#' is assigned to the corresponding values of the cell-level metadata variables
#' provided.
#'
#' If a samples column is available, join is performed using the primary key
#' samples + barcode, to avoid duplicated barcodes.
#'
#' @param cl_metadata data.table with cell-level metadata
#' @param barcode_cell_ids data.table of barcodes, cell_ids and samples from scdata
#'
#' @return data.table of cell_ids and cell-level metadata
#' @export
make_cl_metadata_table <- function(cl_metadata, barcode_cell_ids) {

  join_cols <- c("barcode")
  if (!"samples" %in% names(cl_metadata)) {
    join_cols <- setdiff(join_cols, "samples")
    # remove samples from barcode-cell_id map if not used for join
    barcode_cell_ids <- barcode_cell_ids[, -"samples"]
  }

  cl_metadata[barcode_cell_ids, , on = join_cols]
}


#' check if variable is acceptable as cell-level metadata
#'
#' Tests values using several heurisitics to determine if a variable has a reasonable
#' number of values to be made into cellsets. Cellsets are categorical, therefore
#' continuous or high-cardinality variables are filtered out.
#'
#' @param check_vals vector of arbitrary type
#'
#' @return boolean specifying if variable is desirable or not
#' @export
#'
find_clm_columns <- function(check_vals) {
  # TODO refactor along pipeline::find_cluster_columns

  # skip if too few or way too many values
  value_counts <- table(check_vals)
  n.vals <- length(value_counts)
  if (n.vals > 500) {
    return(FALSE)
  }

  # skip if more than 1/3 of values are repeated fewer than 4 times
  nreps_lt4 <- sum(value_counts < 4)
  if (nreps_lt4 > n.vals / 3) {
    return(FALSE)
  }
  return(TRUE)
}


#' Detect cell-level metadata variable types
#'
#' detect cell level metadata types
#'   - group (CLMPerSample, like metadataCategorical)
#'   - "cellset type" (CLM, like cellSets)
#'   - excludes continuous/high cardinality (to avoid infinite cellsets)
#'
#' @param cl_metadata data.table
#'
#' @return list of cell-level metadata types
#' @export
#'
detect_variable_types <- function(cl_metadata) {

  # can only find group columns when samples are available
  if ("samples" %in% names(cl_metadata)) {
    # do not remove dups; if a user uploads some metadata it should be there
    clm_per_sample_cols <- find_group_columns(cl_metadata, remove.dups = F)
  } else {
    clm_per_sample_cols <- character()
  }

  undesirable_cols <- vapply(cl_metadata[, -..clm_per_sample_cols], find_clm_columns, logical(1))
  clm_cols <- names(undesirable_cols[unlist(undesirable_cols)])

  # remove samples var, useless from this point on
  clm_cols <- grep("^samples$", clm_cols, value = T, invert = T)

  return(list(CLM = clm_cols, CLMPerSample = clm_per_sample_cols))
}


#' Transform a variable into a cell-level metadata cell class
#'
#' The output is a list that contains one children cellset per value of the
#' variable to transform.
#'
#' @param variable character - name of variable to transform
#' @param type character - type of cell-level metadata variable
#' @param cl_metadata data.table
#' @param color_pool list
#'
#' @return list - correctly formatted cell-level metadata cellclass
#' @export
#'
make_cl_metadata_cellclass <- function(variable, type, cl_metadata, color_pool) {

  cl_metadata_cellset <- list(
    key = uuid::UUIDgenerate(),
    name = as.character(variable),
    rootNode = TRUE,
    type = type,
    children = list()
  )

  values <- unique(cl_metadata[[variable]])

  for (i in seq_along(values)) {
    cl_metadata_cellset$children[[i]] <- list(
      key = uuid::UUIDgenerate(),
      name = as.character(values[i]),
      rootNode = FALSE,
      type = type,
      color = color_pool[1],
      cellIds = ensure_is_list_in_json(cl_metadata[get(variable) == values[i], cells_id])
    )
    color_pool <- color_pool[-1]
  }

  return(cl_metadata_cellset)

}


#' Formats cell-level metadata cellsets
#'
#' @param cl_metadata data.table
#' @param var_types list of variable types
#' @param color_pool list
#'
#' @return list of cell-level metadata cell classes
#' @export
#'
format_cl_metadata_cellsets <-
  function(cl_metadata,
           var_types,
           color_pool) {
    # explicitly exclude cells_id var from cellset creation
    vars_to_cellset <-
      setdiff(names(cl_metadata), c("barcode", "cells_id"))
    # only create cellsets for vars that passed type detection
    vars_to_cellset <- lapply(var_types, intersect, vars_to_cellset)

    cl_metadata_cellsets <- list()
    for (i in seq_along(vars_to_cellset)) {
      cellsets <-
        lapply(
          vars_to_cellset[[i]],
          make_cl_metadata_cellclass,
          names(vars_to_cellset)[i],
          cl_metadata,
          color_pool
        )
      cl_metadata_cellsets <- append(cl_metadata_cellsets, cellsets)
    }

    return(cl_metadata_cellsets)
  }

#' Sort cluster names
#'
#' Sorts cluster names naturally, i.e. Cluster 1, Cluster 2, Cluster 10
#'
#' @param strings cluster names
#'
#' @return sorted vector
#' @export
#'
sort_cluster_names <- function(strings) {
  # extract letters and digits
  char <- gsub("\\d", "", strings)
  nums <- gsub("\\D", "", strings)

  sorted_indices <- order(char, as.integer(nums))

  return(strings[sorted_indices])
}
