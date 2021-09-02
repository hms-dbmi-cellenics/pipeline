process_spaceranger_files <- function(input, pipeline_config, prev_out) {
    scdata_list <- list()

    samples <- config$samples
    message("Samples to include in the analysis:\n- ", paste(samples, collapse = "\n- "))
    message("Loading 10x visium data set from input folder.")

    for (sample in samples) {
        sample_dir <- file.path("/input", sample)

        # load and normalize
        scdata <- Seurat::Load10X_Spatial(sample_dir)
        scdata <- Seurat::SCTransform(scdata, assay = "Spatial", verbose = TRUE)

        # cluster and run umap
        scdata <- Seurat::RunPCA(scdata, assay = "SCT", verbose = FALSE)
        scdata <- Seurat::FindNeighbors(scdata, reduction = "pca", dims = 1:30)
        scdata <- Seurat::FindClusters(scdata, verbose = FALSE)
        scdata <- Seurat::RunUMAP(scdata, reduction = "pca", dims = 1:30, umap.method = "umap-learn")
        scdata_list[[sample]] <- scdata
    }

    prev_out$scdata_list <- scdata_list

    res <- list(
        data = list(),
        output = prev_out)

    message("\nProcessing of spaceranger files step complete.")
    return(res)
}
