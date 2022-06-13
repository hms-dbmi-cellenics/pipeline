mock_cellranger_files <- function(counts, features, sample_dir) {

  # save features
  features_path <- file.path(sample_dir, "features.tsv")
  write.table(features, features_path, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  R.utils::gzip(features_path)

  # save barcodes
  barcodes <- colnames(counts)
  barcodes_path <- file.path(sample_dir, "barcodes.tsv")
  writeLines(barcodes, barcodes_path)
  R.utils::gzip(barcodes_path)

  # save Matrix
  if (is(counts, "data.frame")) counts <- as.matrix(counts)
  sparse.mat <- Matrix::Matrix(counts, sparse = TRUE)
  matrix_path <- file.path(sample_dir, "matrix.mtx")
  Matrix::writeMM(sparse.mat, matrix_path)
  R.utils::gzip(matrix_path)
}


mock_counts <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc_raw.txt", package = "Seurat"),
    as.is = TRUE
  )
}

local_cellranger_experiment <- function(counts, features, experiment_dir, sample_dir, env = parent.frame()) {

  sample_path <- file.path(experiment_dir, sample_dir)
  dir.create(sample_path, recursive = T)
  mock_cellranger_files(counts, features, sample_path)
  withr::defer(unlink(experiment_dir, recursive = T), envir = env)

}

mock_lists <- function(){
  counts <- mock_counts()
  symbols <- row.names(counts)
  ensids <- paste0("ENSFAKE", seq_len(nrow(counts)))
  features <- data.frame(input = ensids, name = symbols)
  rownames(counts) <- ensids

  counts_list <- list(sample1=counts,sample2=counts)
  annot_list <- list(sample1=features,sample2=features)
  return(list(counts_list=counts_list,annot_list=annot_list))
}

test_that("format_annot keeps unique rows", {
  annot_list <- list(
    sample1 = data.frame(input = 1:5, name = paste0("gene", 1:5)),
    sample2 = data.frame(input = 1:5, name = paste0("gene", 1:5))
  )

  annot <- format_annot(annot_list)

  expect_s3_class(annot, "data.frame")
  expect_true(nrow(annot) == nrow(annot_list$sample1))
})


test_that("format_annot deduplicates name column", {
  annot_list <- list(
    sample1 = data.frame(input = 1:6, name = paste0("gene", c(1, 1:5)))
  )

  annot <- format_annot(annot_list)

  expect_true(length(annot$name) == length(unique(annot$name)))
})


test_that("format_annot removes duplicated input (Ensembl IDs) column", {
  ensids <- c(1, 1:4)
  annot_list <- list(
    sample1 = data.frame(input = ensids, name = paste0("gene", 1:5))
  )

  annot <- format_annot(annot_list)

  expect_equal(
    length(unique(ensids)),
    length(annot$input)
  )
})


test_that("load_user_files loads a 10x count matrix", {
  counts <- mock_counts()
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)


  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  expect_true("counts_list" %in% names(out))
  expect_true(sample %in% names(out$counts_list))

  expect_s4_class(out$counts_list[[1]], "dgCMatrix")
})


test_that("load_user_files generates feature annotation for 10x data", {
  counts <- mock_counts()
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  expect_true("annot" %in% names(out))
  expect_true(
    all(c("input", "name", "original_name") %in% colnames(out$annot))
  )

})


test_that("load_user_files deduplicates gene symbols for 10x data", {
  counts <- mock_counts()

  symbols <- row.names(counts)
  symbols[1:5] <- "DUPLICATED"

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = symbols
  )

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  annot <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output$annot

  # unique gene names is same as number of gene names
  expect_length(unique(annot$name), length(symbols))

  # unique original names is same as unique gene names
  expect_length(unique(annot$original_name), length(unique(symbols)))
})


test_that("load_user_files uses appropriate feature columns for 10x data", {
  counts <- mock_counts()

  symbols <- row.names(counts)
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = symbols
  )

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  # ensembl ids are counts row names
  expect_equal(
    features$ensid,
    row.names(out$counts_list[[1]])
  )

  # ensembl ids are in column 'input' of annot
  expect_equal(out$annot$input, features$ensid)

  # symbols are in column 'name' of annot
  expect_equal(out$annot$name, symbols)

})


test_that("load_user_files loads 10x multisample experiments", {
  counts <- mock_counts()

  symbols <- row.names(counts)
  ensids <- paste0("ENSFAKE", seq_len(nrow(counts)))

  features <- data.frame(ensid = ensids, symbol = symbols)

  # most overlapping, some unique
  ensids2 <- ensids
  symbols2 <- symbols
  ensids2[1:20] <- paste0(ensids2[1:20], "_2")
  symbols2[1:20] <- paste0(symbols2[1:20], "_2")

  features2 <- data.frame(ensid = ensids2, symbols = symbols2)

  experiment_dir <- "./experiment_1"
  samples <- c("sample_a", "sample_b")

  local_cellranger_experiment(counts, features, experiment_dir, samples[1])
  local_cellranger_experiment(counts, features2, experiment_dir, samples[2])


  prev_out <- list(config = list(samples = samples, input = list(type = "10x")))
  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  # loaded both
  expect_equal(names(out$counts_list), samples)

  # kept only unique ensembl ids
  unique_ensids <- unique(c(ensids, ensids2))
  expect_length(out$annot$input, length(unique_ensids))

  # removed duplicated ensembl ids
  expect_lt(
    length(unique_ensids),
    length(c(ensids, ensids2))
  )

})

test_that("read_10x_files returns error if files missing", {
  counts <- mock_counts()
  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"
  sample_dir <- file.path(experiment_dir, sample)

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  files <- c("barcodes.tsv.gz", "matrix.mtx.gz")

  # remove files one by one renaming
  for (file in files) {
    file.rename(file.path(sample_dir, file), file.path(sample_dir, "blah"))
    expect_error(load_user_files(NULL, NULL, prev_out, experiment_dir), "file missing")
    file.rename(file.path(sample_dir, "blah"), file.path(sample_dir, file))
  }

  file <- "features.tsv.gz"

  file.rename(file.path(sample_dir, file), file.path(sample_dir, "blah"))
  expect_error(supressWarnings(load_user_files(NULL, NULL, prev_out, experiment_dir), "cannot open the connection"))
  file.rename(file.path(sample_dir, "blah"), file.path(sample_dir, file))

})

test_that("get_feature_types properly determines types", {
  annot_list <- list(
    sample1 = data.frame(ENSID = paste0("ENS", 1:5), SYMBOL = paste0("gene", 1:5)),
    sample2 = data.frame(ENSID = paste0("gene", 1:5), SYMBOL = paste0("gene", 1:5)),
    sample3 = data.frame(ENSID = paste0("ENS", 1:5), SYMBOL = paste0("ENS", 1:5)),
    sample4 = data.frame(ENSID = paste0("gene", 1:5), SYMBOL = paste0("ENS", 1:5))
  )

  expect_equal(pipeline:::get_feature_types(annot_list[["sample1"]]),1)
  expect_equal(pipeline:::get_feature_types(annot_list[["sample2"]]),0)
  expect_equal(pipeline:::get_feature_types(annot_list[["sample3"]]),2)
  expect_equal(pipeline:::get_feature_types(annot_list[["sample4"]]),-1)
})

test_that("get_feature_types identifies mixed columns", {
  annot_list <- list(
    sample1 = data.frame(ENSID = c(paste0("ENS", 1:5),paste0("gene",1:5)), SYMBOL = c(paste0("ENS", 1:4),paste0("gene",1:6))),
    sample2 = data.frame(ENSID = c(paste0("ENS", 1:4),paste0("gene",1:6)), SYMBOL = paste0("gene", 1:10))
  )

  expect_equal(pipeline:::get_feature_types(annot_list[["sample1"]]),1)
  expect_true(pipeline:::get_feature_types(annot_list[["sample2"]])==0)
})


test_that("equalize_annotation_types does nothing if all annotation types are the same",{
  input <- mock_lists()

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  counts_list <- input$counts_list
  annot_list <- input$annot_list

  res <- pipeline:::equalize_annotation_types(annot_list,counts_list,features_types_list, samples)

  expect_equal(res[[1]],counts_list)
  expect_equal(res[[2]],annot_list)
})

test_that("equalize_annotations does nothing if there are no samples with annotations",{
  input <- mock_lists()

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  samples <- c("sample1","sample2")
  for(sample in samples){
    input$annot_list[[sample]]$input <- input$annot_list[[sample]]$name
    rownames(input$counts_list[[sample]])  <- input$annot_list[[sample]]$input
  }

  counts_list <- input$counts_list
  annot_list <- input$annot_list

  res <- pipeline:::equalize_annotation_types(annot_list, counts_list, features_types_list, samples)

  expect_equal(res[[1]],counts_list)
  expect_equal(res[[2]],annot_list)
})


test_that("equalize_annotation_types infers gene ids from symbols and corrects counts rownames",{
  input <- mock_lists()

  sample2_annot <- input$annot_list$sample2
  sample2_annot$input <- sample2_annot$name
  input$annot_list$sample2 <- sample2_annot
  rownames(input$counts_list$sample2) <- sample2_annot$input

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expect_equal(res$annot_list$sample2,input$annot_list$sample1)
  expect_equal(rownames(res$counts_list$sample2),res$annot_list$sample2$input)
})

test_that("equalize_annotation_types infers gene symbols from ids",{
  input <- mock_lists()

  sample2_annot <- input$annot_list$sample2
  sample2_annot$name <- sample2_annot$input
  input$annot_list$sample2 <- sample2_annot
  rownames(input$counts_list$sample2) <- sample2_annot$input

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expect_equal(res$annot_list$sample2,input$annot_list$sample1)
  expect_equal(rownames(res$counts_list$sample2),input$annot_list$sample2$input)
})

test_that("equalize_annotation_types infers ids with incomplete match",{
  input <- mock_lists()

  sample2_annot <- input$annot_list$sample2
  sample2_annot$name[1:nrow(sample2_annot)%%2==1] <- paste0("gene",(1:(nrow(sample2_annot)/2)))
  sample2_annot$input <- sample2_annot$name
  rownames(input$counts_list$sample2) <- sample2_annot$input
  input$annot_list$sample2 <- sample2_annot

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expected_annot <- input$annot_list$sample1
  expected_annot$input[1:nrow(expected_annot)%%2==1] <- paste0("gene",(1:(nrow(expected_annot)/2)))
  expected_annot$name[1:nrow(expected_annot)%%2==1] <- paste0("gene",(1:(nrow(expected_annot)/2)))

  expect_equal(res$annot_list$sample2,expected_annot)
  expect_equal(rownames(res$counts_list$sample2),res$annot_list$sample2$input)
})

test_that("equalize_annotation_types infers symbols with incomplete match and doesnt modify ids",{
  input <- mock_lists()

  sample2_annot <- input$annot_list$sample2
  sample2_annot$name <- sample2_annot$input
  sample2_annot$input[1:nrow(sample2_annot)%%2==1] <- paste0("gene",(1:(nrow(sample2_annot)/2)))
  rownames(input$counts_list$sample2) <- sample2_annot$input
  input$annot_list$sample2 <- sample2_annot

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expected_annot <- input$annot_list$sample1
  expected_annot$input[1:nrow(expected_annot)%%2==1] <- paste0("gene",(1:(nrow(expected_annot)/2)))
  expected_annot$name[1:nrow(expected_annot)%%2==1] <- input$annot_list$sample1$input[1:nrow(expected_annot)%%2==1]

  expect_equal(res$annot_list$sample2,expected_annot)
  expect_equal(rownames(res$counts_list$sample2),input$annot_list$sample2$input)
  expect_equal(res$annot_list$sample2$input,input$annot_list$sample2$input)
})

test_that("equalize_annot properly infers ids with more than 2 samples",{
  input <- mock_lists()

  annot <- input$annot_list$sample1
  counts <- input$counts_list$sample1

  input$counts_list$sample3 <- counts[1:(nrow(counts)/2),]
  input$annot_list$sample3 <- annot[1:(nrow(annot)/2),]

  input$counts_list$sample1 <- counts[(nrow(counts)/2):nrow(counts),]
  input$annot_list$sample1 <- annot[(nrow(annot)/2):nrow(annot),]

  sample2_annot <- input$annot_list$sample2
  sample2_annot$name[1:nrow(sample2_annot)%%2==1] <- paste0("gene",(1:(nrow(sample2_annot)/2)))
  sample2_annot$input <- sample2_annot$name
  rownames(input$counts_list$sample2) <- sample2_annot$input
  input$annot_list$sample2 <- sample2_annot

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expected_annot <- annot
  expected_annot$input[1:nrow(expected_annot)%%2==1] <- paste0("gene",(1:(nrow(expected_annot)/2)))
  expected_annot$name[1:nrow(expected_annot)%%2==1] <- paste0("gene",(1:(nrow(expected_annot)/2)))

  expect_equal(res$annot_list$sample2,expected_annot)
  expect_equal(rownames(res$counts_list$sample2),res$annot_list$sample2$input)
})

test_that("duplicated genes dont lead to any rowname duplication",{
  input <- mock_lists()

  annot <- input$annot_list$sample1
  counts <- input$counts_list$sample1

  sample2_annot <- input$annot_list$sample2
  sample2_annot$name[1:nrow(sample2_annot)%%2==1] <- paste0("gene",(1:(nrow(sample2_annot)/2)))
  sample2_annot$input <- sample2_annot$name

  sample2_annot[1,1] <- sample2_annot[2,1]
  sample2_annot[1,1] <- annot[2,1]


  rownames(input$counts_list$sample2) <- sample2_annot$input
  input$annot_list$sample2 <- sample2_annot

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expected_annot <- annot
  expected_annot$input[1:nrow(expected_annot)%%2==1] <- paste0("gene",(1:(nrow(expected_annot)/2)))
  expected_annot$name[1:nrow(expected_annot)%%2==1] <- paste0("gene",(1:(nrow(expected_annot)/2)))

  expected_annot$input[[1]] <- "ENSFAKE2"
  expected_annot$input[[2]] <- "ENSFAKE2.1"

  expect_equal(res$annot_list$sample2,expected_annot)
  expect_equal(rownames(res$counts_list$sample2),expected_annot$input)
})

test_that("Mislabeling of features types to symbols results in no changes",{
  input <- mock_lists()

  sample2_annot <- input$annot_list$sample2
  sample2_annot$input[1:(nrow(sample2_annot)/2+1)] <- paste0("ab",1:(nrow(sample2_annot)/2+1))
  sample2_annot$name[1:(nrow(sample2_annot)/2+1)] <- paste0("ab",1:(nrow(sample2_annot)/2+1))
  input$annot_list$sample2 <- sample2_annot
  rownames(input$counts_list$sample2) <- sample2_annot$input

  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expect_equal(res$annot_list$sample2,input$annot_list$sample2)
  expect_equal(rownames(res$counts_list$sample2),input$annot_list$sample2$input)
})

test_that("Mislabeling of features types to ids results in no changes",{
  input <- mock_lists()

  sample2_annot <- input$annot_list$sample2
  sample2_annot$input[1:(nrow(sample2_annot)/2+1)] <- paste0("ENS",1:(nrow(sample2_annot)/2+1))
  sample2_annot$name[1:(nrow(sample2_annot)/2+1)] <- paste0("ENS",1:(nrow(sample2_annot)/2+1))
  input$annot_list$sample2 <- sample2_annot
  rownames(input$counts_list$sample2) <- sample2_annot$input


  features_types_list <- list(sample1=get_feature_types(input$annot_list$sample1),sample2=get_feature_types(input$annot_list$sample2))

  res <- pipeline:::equalize_annotation_types(input$annot_list,input$counts_list,features_types_list,samples=list("sample1","sample2"))

  expect_equal(res$annot_list$sample2,input$annot_list$sample2)
  expect_equal(rownames(res$counts_list$sample2),input$annot_list$sample2$input)
})


test_that("read_10x_files removes rows with empty feature names both in count matrix and annotation if present and < 0.1%", {
  # mock count matrix replicating it 10 times to mock a matrix with < 0.1% of empty features
  counts <- mock_counts()[rep(seq_len(nrow(mock_counts())), each = 10), ]
  rownames(counts)[2] <- ""
  rownames(counts)[3] <- ".1"

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )
  features[2:3, 1:2] <- ""

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(length(which(rownames(counts_list[[1]]) == "")), 0)
  expect_equal(length(which(annot[, 1] == "")), 0)
})


test_that("read_10x_files removes single row with empty feature names both in count matrix and annotation if present and < 0.1%", {
  # mock count matrix replicating it 10 times to mock a matrix with < 0.1% of empty features
  counts <- mock_counts()[rep(seq_len(nrow(mock_counts())), each = 10), ]
  rownames(counts)[2] <- ""

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )
  features[2, 1:2] <- ""

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output
  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(length(which(rownames(counts_list[[1]]) == "")), 0)
  expect_equal(length(which(annot[, 1] == "")), 0)
})


test_that("read_10x_files doesn't remove any rows with empty feature names both in count matrix and annotation if present and >= 0.1%", {
  counts <- mock_counts()
  rownames(counts)[2] <- ""
  rownames(counts)[3] <- ".1"

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )
  features[2:3, 1:2] <- ""

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(nrow(counts), nrow(counts_list[[1]]))
  # expect_equal(nrow(features), nrow(annot))  # decomment this line and delete the following line when make unique will be added in format_annot [BIOMAGE-1817]
  expect_equal(length(which(annot[, 1] == "")), 1)
})


test_that("read_10x_files doesn't remove any rows if no rows with empty rownames are present", {
  counts <- mock_counts()

  features <- data.frame(
    ensid = paste0("ENSFAKE", seq_len(nrow(counts))),
    symbol = row.names(counts)
  )

  experiment_dir <- "./experiment_1"
  sample <- "sample_a"

  local_cellranger_experiment(counts, features, experiment_dir, sample)

  prev_out <- list(config = list(samples = sample, input = list(type = "10x")))

  out <- load_user_files(NULL, NULL, prev_out, experiment_dir)$output

  counts_list <- out$counts_list
  annot <- out$annot

  expect_equal(length(which(rownames(counts_list[[1]]) == "")), 0)
  expect_equal(length(which(annot[, 1] == "")), 0)
  expect_equal(nrow(counts), nrow(counts_list[[1]]))
  expect_equal(nrow(features), nrow(annot))
})
