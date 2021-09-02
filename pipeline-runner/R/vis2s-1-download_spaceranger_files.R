download_spaceranger_files <- function(input, pipeline_config, prev_out = list()) {
    project_id <- input$projectId
    sample_names <- input$sampleNames
    sample_uuids <- input$sampleIds

    s3 <- paws::s3(config = pipeline_config$aws_config)

    spaceranger_fnames <- c(
        "filtered_feature_bc_matrix.h5",
        "spatial/aligned_fiducials.jpg",
        "spatial/detected_tissue_image.jpg",
        "spatial/scalefactors_json.json",
        "spatial/tissue_hires_image.png",
        "spatial/tissue_lowres_image.png",
        "spatial/tissue_positions_list.csv"
        )

    unlink("/input", recursive = TRUE)

    for (sample in sample_uuids) {
        message("\nSample --> ", sample)

        for (spaceranger_fname in spaceranger_fnames) {
            vis_key <- file.path(project_id, sample, spaceranger_fname)

            message("Visium key: ", vis_key)

            sample_name <- sample_names[[match(sample, sample_uuids)]]
            # Preparing directories
            local_dir <- file.path("/input", sample)
            dir.create(file.path(local_dir, 'spatial'), recursive = TRUE)

            local_fpath <- file.path(local_dir, spaceranger_fname)

            # Download the file and store the output in a variable
            c(body, ...rest) %<-% s3$get_object(
                Bucket = pipeline_config$originals_bucket,
                Key = vis_key
            )

            # Write output to file
            writeBin(body, con = local_fpath)
        }
    }

    config <- list(
        name = input$experimentName,
        samples = input$sampleIds,
        organism = input$organism,
        input = list(type = input$input$type)
    )

    if ("metadata" %in% names(input)) {
        config$metadata <- input$metadata
    }

    prev_out$config <- config
    res <- list(
        data = list(),
        ouput = prev_out)

    message("\nDownloading of spaceranger files step complete.")
    return(res)
}
